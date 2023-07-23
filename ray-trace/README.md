# Ray traceing


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