# SURFACE 利用手順書

## Revisions

* Jan 13, 2015 : SURFACE に名称変更
* Sep 02, 2014 : Linux 環境での設定について追記
* Mar 04  2014 : 日本語版
* Feb 06, 2014 : Introduce --platform option to replace --cross-k option.
* Jan 06, 2014 : GLES unit test is now built with --with-unittest option.
* Nov 21, 2013 : Refresh GLSL build instruction. Add directory layout. Add build option --with-glsl-opt
* Nov 11, 2013 : Add note on library names.
* Oct 29, 2013 : Update build environment. Add note on Intel C/C++ compiler.
* Mar 07, 2013 : Add note on Windows.
* Feb 13, 2012 : Initial.

## Abstract

SURFACE(Scalable and Ubiquitous Rendering Framework for Advanced Computing Environments) は 2 つのモジュールで構成されます.

* `libLSGLES`, GLES v2.0 API 互換インターフェイスを持つレンダラライブラリ
* `glslc`, オフライン GLSL コンパイラ. 一部 Python で記述

LSGL を利用したアプリケーションを動作させるには, 両方とも適切にセットアップする必要があります.

## Requirements

* premake4(`tools` ディレクトリにプレビルドバイナリがあります)
* C/C++ コンパイラ
  * gcc 4.4 or later(For Linux or Mac)
  * Intel icc 13.0(For Linux)
  * clang 3.2 or later(For Linux or Mac)
  * Fujitsu's C/C++ compiler(for K/FX10)
  * Visual Studio 2010 or later(for Windows)
  * Intel C++ 13.0 or later(for Linux)
  * MinGW x64 http://www.drangon.org/mingw/ (Windows でシェーダコードをコンパイルするときに利用)
* Python 2.6 or 2.7(Linux, K/FX10, MacOSX, Windows 共通)

## Working environments

LSGL は, 以下の環境での動作を確認しています.

* Linux 64bit 
  * RHEL 6.1, CentoS 6.5.
  * OpenMPI および Intel MPI
* K/FX10 : K-1.2.0-15.
* Windows7 64bit : Visual Studio 2010, Visual Studio 2013.
* MacOSX : 10.8 以降.


## Directory layout

    gles/         GLES v2.0 compatible API implmenetation
    glsl/         GLSL compiler implementation
    render/       Raytracing render core
    include/      GLES header files
    examples/     Examples using LSGL library
    test/         Test directory
    doc/          Document 
    tools/        Tools directory 

## libLSGLES のビルド.

libLSGLES のビルドには premake4 を利用しています.
premake4 を用いた, 各環境での Makefile を生成するブートストラップスクリプトを用意しています.

### K/FX10

フロントエンドノードでクロスコンパイル(MPI サポートなし)

    $ cd gles
    $ ../scripts/setup_k.sh
    $ make

フロントエンドノードでクロスコンパイル(MPI サポートあり)

    $ cd gles
    $ ../scripts/setup_k_mpi.sh
    $ make

### Linux x86

    $ cd gles
    $ ../scripts/setup_linux_intel.sh
    $ make

### Linux x86(MPI サポートあり)

    $ cd gles
    $ ../scripts/setup_linux_intel_mpi.sh
    $ make

### Linux x86(Intel C/C++ compiler 利用)

    $ cd gles
    $ ../scripts/setup_linux_intel_icc.sh
    $ make

### MacOSX

    $ cd gles
    $ ../scripts/setup_macosx.sh
    $ make

### Windows

vcbuild.bat でターゲットとする Visual Studio のバージョンを指定したのち,

    > cd gles
    > vcbuild.bat

ソリューションファイルが生成されるので Visual Stuio で開いてビルドを行います.

### Bulid オプション

より明示的なカスタマイズのために, premake4 に以下のオプションを用意しています.

    --platform="k-cross"

K/FX10 クロスコンパイル用の設定を行います.

    --platform="k-cross-mpi"

K/FX10 で MPI 利用のクロスコンパイル用の設定を行います.

    --with-debugtrace

GLES API のデバッグトレース機能を有効にします. 各 GLES API コールで, API がどの値で呼び出されたかの情報が stdout に出力されるようになります.

    --with-unittest

ユニットテストをビルドするようにします.

    --with-openmp

OpenMP サポートを有効にします(Mac, Linux, K/FX10).

    --with-mpi

MPI サポートを有効にします(Mac, Linux, K/FX10). ライブラリ名に `_mpi` プレフィックスが付くようになります.

    --with-screen-parallel

MPI 有効オプション(--with-mpi)と組み合わせて, スクリーン並列による並列化を行うようにします. ライブラリ名に `_screen_parallel` プレフィックスが付くようになります.

    --with-glsl-opt

GLSL コード生成を最適化します(実験的機能)

### ライブラリ名のまとめ

* libLSGL : MPI と screen parallel フラグが無効のリリースビルド(デフォルト構成)
* libLSGL_mpi : MPI 有効のリリースビルド
* libLSGL_mpi_screen_parallel: MPI + screen parallel が有効のリリースビルド.
* libLSGLd : MPI と screen parallel が無効のデバッグビルド.
* libLSGLd_mpi : MPI 有効のデバッグビルド.
* libLSGLd_mpi_screen_parallel: MPI + screen parallel flag が有効のデバッグビルド.

## GLSL コンパイラのビルドと設定.

すでにプレビルド済みのバイナリを `glsl/bin` に配置してあります. ビルドについては `glsl/README.md` を参照してください.

### GLSL コンパイラ設定

`glShaderBinary`, `glShaderSource` / `glCompileShader` が適切に動作するように, 以下の手順で GLSL コンパイラを適切に設定するのを忘れずに行ってください.

`glsl/glslc` が GLSL コンパイラのフロントエンドになります(Windows の場合は glslc.bat).

`glsl/glsl_config.py` を編集し, ユーザの環境ごとに設定してください(glslc へのパスや, C++ コンパイラのパスなど). Windows の場合は C++ コンパイラには Visual Studio のコンパイラ(cl.exe)ではなく, MinGW x64(gcc) コンパイラを指定してください.

環境変数 `GLSL_COMPILER` を適切にセットします. この環境変数は, `glShaderSource` が内部で GLSL コンパイラを呼び出すときに利用します.

    $ export GLSL_COMPILER=/path/to/lsgl/glsl/glslc

(Windows の場合は glslc.bat へのパスとなります).

設定が適切に行われているか, 以下のテストコードをコンパイルして確認します.

    $ cat input.frag
    #ifdef GL_ES
    precision mediump float;
    #endif
    
    void main()
    {
      gl_FragColor = vec4(1.0, 1.0, 1.0, 1.0);
    }
    
    $ $GLSL_COMPILER input.frag 

設定が適切に行われていれば, `shader.so` ファイルがカレントディレクトリに生成されます.

## サンプルアプリケーション

サンプルプログラムが `examples` ディレクトリに配置してあります.

たとえば, particle_render サンプルをビルドするには以下のようにします.

    $ cd examples/particle_render
    $ make
    # symlink DLL path(MacOSX)
    $ ln -s ../../gles/libLSGLES.dylib libLSGLES.dylib
    # symlink DLL path(Linux)
    $ ln -s ../../gles/libLSGLES.so libLSGLES.so
    # Add current directory to DLL path
    $ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./
    # compile shader
    $ $GLSL_COMPILER input.frag
    $ ./a.out

インストールと設定が適切に行われていれば, レンダリング結果として `colorbuf.tga` ファイルが生成されます.

## Unit testing

ユニットテストを有効で libLSGLES をビルドすると, `gles/` ディレクトリにユニットテストプログラムが生成されます. このユニットテストプログラムの実行は単に

    $ ./test_LSGL
    ...

として実行します.

以上.
