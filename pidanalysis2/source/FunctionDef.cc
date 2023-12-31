#include "FunctionDef.hh"

double distance_func(double* xx, double* par) {
    double x = xx[0];
    double y = xx[1];
    double z = par[0];

    // Emission Point
    double xE = par[1];
    double yE = par[2];
    double zE = par[3];

    // Perpendicular of back of Radiator
    double nx = par[4];
    double ny = par[5];
    double nz = par[6];

    // Center & Radius of Spherical Mirror
    double x0Mirror = par[7];
    double y0Mirror = par[8];
    double z0Mirror = par[9];
    double Radius = par[10];

    // Refractive Index
    double n_Aero = par[11];
    double n_Air = par[12];

    // Detection Point
    double xD = par[13];
    double yD = par[14];
    double zD = par[15];

    Vector3D rP(x, y, z);
    Vector3D rE(xE, yE, zE);
    Vector3D n_perp(nx, ny, nz);
    Vector3D r0Mirror(x0Mirror, y0Mirror, z0Mirror);
    Vector3D rD(xD, yD, zD);

    // Refraction
    Vector3D d_refraction = Vector3D::refraction(rP - rE, n_perp, n_Aero, n_Air);

    // Cross-Point of refraction light and Spherical Mirror
    Vector3D rPM0 = rP - r0Mirror;
    double t_c = -1.0 * Vector3D::dot(rPM0, d_refraction) +
                 std::sqrt(Radius * Radius - rPM0.norm2() +
                           std::pow(Vector3D::dot(rPM0, d_refraction), 2.0));
    Vector3D rM = rP + t_c * d_refraction;

    // Reflect on Spherical Mirror
    Vector3D d_reflect = Vector3D::reflect(d_refraction, rM - r0Mirror);

    // Distance between detection point and reflect light
    Vector3D cross_tmp = Vector3D::cross(rD - rM, d_reflect);

    return cross_tmp.norm() / d_reflect.norm();
}

double quantumEff(double* wlen, double* par) {
    double wlength = wlen[0];  // nm

    double p0 = -8.18311;
    double p1 = 0.0558217;
    double p2 = -0.000132938;
    double p3 = 1.40711e-07;
    double p4 = -6.2958e-11;
    double p5 = 7.41335e-15;

    if (wlength < 315 || wlength > 900) {
        return 0.0;
    } else {
        return p0 + p1 * wlength + p2 * pow(wlength, 2.0) + p3 * pow(wlength, 3.0) +
               p4 * pow(wlength, 4.0) + p5 * pow(wlength, 5.0);
    }
}

double refractive_index(double* wlen, double* par) {
    double wlength = wlen[0];  // nm
    double a0 = par[0];
    double wlength0 = par[1];  // nm

    return std::sqrt(1 + a0 * std::pow(wlength, 2.0) /
                             (std::pow(wlength, 2.0) - std::pow(wlength0, 2.0)));
}

double transmittance(double* wlen, double* par) {
    double wlength = wlen[0];  // nm
    double A = par[0];
    double C = par[1];      // Clarity factor [um^4/cm]
    double thick = par[2];  // cm
    return A * std::exp(-C * thick / std::pow(wlength * 1e-3, 4.0));
}

double transmittance2(double* wlen, double* par) {
    double wlength = wlen[0];  // nm
    double A = par[0];
    double d = par[1];      // diameter of scattering center [um]
    double thick = par[2];  // cm
    double a0 = par[3];
    double wlength0 = par[4];
    double R_index2 =
        1 + a0 * std::pow(wlength, 2.0) / (std::pow(wlength, 2.0) - std::pow(wlength0, 2.0));

    double rho_ap = 0.01;  // g/cm^3
    double rho_si = 2.2;   // g/cm^3
    return A * std::exp(-1 * 1e4 * 4 * pow(M_PI, 4.0) * rho_ap * pow(d, 3.0) * thick *
                        pow(wlength * 1e-3, -4.0) * (R_index2 - 1) / (R_index2 + 2) / rho_si);
}

double allDetectedPhoton(double mass, double momentum, double thickness,
                         TF1* qeff, TF1* transmit, TF1* r_index, double geo_eff, double wlenmin,
                         double wlenmax) {
    double alpha = 1.0 / 137;  // fine structure constant
    //double mass_list[4] = {139.57018, 493.677, 938.272081, 0.51099895};  // Pi, K , proton, electron [MeV/c^2]

    double beta = momentum / std::sqrt(std::pow(momentum, 2.0) +
                                       std::pow(mass * 1e-3, 2.0));  // v/c
    double reflectivity_mirror = 0.89;

    // search wavelength of 1/βn = 1
    double wlenmin2, wlenmax2;
    if(r_index->Eval(wlenmax)>1/beta){
        wlenmin2 = wlenmin;
        wlenmax2 = wlenmax;
    }else if(r_index->Eval(wlenmin)>1/beta){
        wlenmin2 = wlenmin;
        wlenmax2 = r_index->GetX(1/beta);
    }else{
        return 0.0;
    }

    double deltawlen = 1.0;  // nm
    int nstep = (int)((wlenmax2 - wlenmin2) / deltawlen);

    //double thickness = 10.0;  // cm
    double deltathick = 1.0;  // cm
    int nstep_thick = (int)(thickness / deltathick);

    double sum_temp = 0.0;
    double wlenlow;  //, wlenhigh;
    double NumOfPhoton = 0.0;
    for (int j = 0; j < nstep_thick; j++) {
        transmit->SetParameter(2, deltathick / 2.0 + j * deltathick);
        sum_temp = 0.0;
        for (int i = 0; i < nstep; i++) {
            wlenlow = wlenmin2 + i * deltawlen;
            if (i == 0 || i == nstep) {
                sum_temp +=
                    geo_eff * qeff->Eval(wlenlow) * reflectivity_mirror * transmit->Eval(wlenlow) *
                    (1.0 - std::pow(beta * r_index->Eval(wlenlow), -2.0)) / std::pow(wlenlow, 2.0);
            } else {
                sum_temp += geo_eff * 2.0 * qeff->Eval(wlenlow) * reflectivity_mirror *
                            transmit->Eval(wlenlow) *
                            (1.0 - std::pow(beta * r_index->Eval(wlenlow), -2.0)) /
                            std::pow(wlenlow, 2.0);
            }
        }
        NumOfPhoton += 2.0 * M_PI * alpha * deltathick * 1e7 * 0.5 * deltawlen * sum_temp;
    }

    return NumOfPhoton;
}

double angle_probDis(double theta, double theta_theory, double sigma){
    return std::exp(-0.5*std::pow((theta-theta_theory)/sigma,2.0))/(std::pow(2.0*M_PI,3.0/2)*sigma);
}