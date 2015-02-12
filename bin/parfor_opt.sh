#!/bin/sh
for K in 1 2 4 8 16 32
do
    echo $K
    export HPCE_FFT_LOOP_K=$K
    ./time_fourier_transform hpce.tp1811.fast_fourier_transform_parfor
	time_fourier_transform hpce.tp1811.fast_fourier_transform_parfor >> parfor_opt.txt
done