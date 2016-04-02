# GLSL compiler

## Revision

* 2013.11.13  Changed to use Mesa package from LSGL repo.
* 2013.10.22  Added note on python version.
* 2013.10.17  Added sparc64 build of glsl compiler.

## Requirements

* Python 2.4 or later(2.6 preferred)

## Setup(for developers)

### Patching

Apply a patch to MESA's GLSL Compiler to add LSGL extension and also fix a bug of S-exp print for structure.
(Tested Mesa 9.0.1)

    $ cd $lsgl/glsl/deps
    $ tar -jxvf Mesa-9.0.1.tar.bz2
    $ cd Mesa-9.0.1
    $ patch -p1 < ../../patch/mesa_glsl_fix.patch

### Mesa build on K native

    1. Install scons 2.3.0
    2. cd $lsgl/glsl/deps/Mesa-9.0.1/src/glsl
    3. scons -u   # Use K(Sparc) gcc
    4. rename $lsgl/glsl/deps/Mesa-9.0.1/build/linux-debug/glsl/glsl2 to glsl_compiler
   

### Mesa build on x86
(Tested Mesa 9.0.1 with MacOSX 10.7.5)

    1.brewなどでsconsをいれる
     brew install scons
    2.libxml2をダウンロード/解凍
    3.libxml2をビルド
     ./configure --prefix=/usr/local --with-python
    4.libxml2のpythonのプラグインをビルド
     cd python
     python setup.py build
     sudo python setup.py install
    5.上記patchを当てる
      Mesa-9.0.1/src/glsl/以下など
    6.Mesaをビルド
     cd Mesa-9.0.1
     scons build=release gles=yes
    7.glslコンパイラを取り出す
     cd build/darwin-x86_64/glsl
     cp glsl2 ../path/to/LSGL/glsl/bin/macosx64/
    8.LSGLのglsl_config.pyを設定する

### Mesa build on Windows

1. Install python 2.7.6
2. Install scons 2.2.0
3. Install flex and bison: win_flex_bison 2.4.1 http://sourceforge.net/projects/winflexbison/
  (rename win_bison -> bison and win_flex -> flex)
4. Open Visual Studio command line shell(VS2012-x86_64 or VS2013-x86_64 preferred)
5. Build patched Mesa GLSL compiler.

    > cd Mesa-9.0.1/src/glsl
    > scons build=release machine=x86_64 -u
    > cp ..\\..\\build\\windows-x86_64\\glsl\\glsl2.exe /path/to/lsgl/glsl/bin/windows_x64/glsl_compiler.exe

### Configuration

edit glsl_config.py to match your evironment.

## Usage

    ./glslc [options] input.frag

    options: -S          output C/C++ source code
             -o          specify output shader dll(default "shader.so")
             -v          verbose mode

## How it works

Mesa の GLSL コンパイラにより, GLSL ソースのパースと抽象木レベルでの最適化, 中間言語表現(S 式)への変換を行う.

その後, python スクリプトによるトランスレータ(glsl_translate.py) により S 式を読み込み C/C++ コードへと変換する.

C/C++ コードをターゲット環境用の C/C++ コンパイラでコンパイルし, シェーダモジュール(DLL)を作成する.

glslc(python で記述) では, 上記の処理を一式行う. 従って, ユーザから見れば glslc のみを利用すればよい.
