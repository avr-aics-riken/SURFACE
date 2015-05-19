#!/bin/sh

# Use jemalloc
LD_PRELOAD=$HOME/local/lib/libjemalloc.so ./render-benchmark
