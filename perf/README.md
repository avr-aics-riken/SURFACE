# SURFACE benchmark suite

Do renderer benchmark with hayai benchmark library. https://bruun.co/2012/02/07/easy-cpp-benchmarking

## Setup

If it is your first time of git clone, checkout dependenccy libraries

    $ cd ${SURFACE.git}/perf/
    $ mkdir deps
    $ cd deps
    $ git clone https://github.com/nickbruun/hayai.git

Or run `setup.sh` which does same thing in the above commands.

    $ ./setup.sh

## Build

### MacOSX

Edit compiler line in bootstrap/macosx.sh , then

    $ ./bootstrap/macosx.sh
    $ make

### Linux

T.B.W.

## Run benchmark

$ ./render-benchmark
