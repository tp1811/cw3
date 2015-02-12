#!/bin/sh
for K in 1 2 4 8 16 32 64 128 256 512 1024 2048 4096 8192 16384
do
    echo $K
    export HPCE_FFT_RECURSION_K=$K
    ./time_fourier_transform hpce.tp1811.fast_fourier_transform_taskgroup
	time_fourier_transform hpce.tp1811.fast_fourier_transform_taskgroup >> taskgroup_opt.txt
done
