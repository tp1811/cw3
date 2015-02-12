#!/bin/sh

./time_fourier_transform hpce.tp1811.direct_fourier_transform_parfor
time_fourier_transform hpce.tp1811.direct_fourier_transform_parfor >> comparison.txt

./time_fourier_transform hpce.tp1811.fast_fourier_transform_combined
time_fourier_transform hpce.tp1811.fast_fourier_transform_combined >> comparison.txt

./time_fourier_transform hpce.tp1811.fast_fourier_transform_opt
time_fourier_transform hpce.tp1811.fast_fourier_transform_opt >> comparison.txt
