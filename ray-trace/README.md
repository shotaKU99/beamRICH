# Ray traceing


## 実行方法

```bash
$ mkdir build
$ cd build
/build$ cmake ../
/build$ make -j4
/build$ cd ..
$ ./raytrace
```

ROOT のパスを通しておく必要がある。`thisroot.sh` を実行しておく。

## Yaml-cpp の使い方

### インストール方法

```bash
$ git clone https://github.com/jbeder/yaml-cpp.git
$ cd yaml-cpp
$ mkdir build
$ cd build
$ cmake -DYAML_BUILD_SHARED_LIBS=ON ..
$ make -j4
$ sudo make install
```

共有ライブラリへのパスを通す
`libyaml-cpp.so.0.X` が `/usr/local/bin` にある。

`.bashrc` にパスを追加する
```
export LD_LIBRARY_PATH=/usr/local/bin:$LD_LIBRARY_PATH
```

### 実行


```bash
g++ -o YamlCpp YamlCpp.cc -lyaml-cpp
```
のように実行する
