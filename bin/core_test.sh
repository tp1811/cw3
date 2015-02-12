#!/bin/sh

./time_fourier_transform hpce.tp1811.fast_fourier_transform_opt 1 3
time_fourier_transform hpce.tp1811.fast_fourier_transform_opt >> core_test.txt

./time_fourier_transform hpce.tp1811.fast_fourier_transform_opt 16 3
time_fourier_transform hpce.tp1811.fast_fourier_transform_opt >> core_test.txt

./time_fourier_transform hpce.tp1811.fast_fourier_transform_opt 64 3
time_fourier_transform hpce.tp1811.fast_fourier_transform_opt >> core_test.txt

./time_fourier_transform hpce.tp1811.fast_fourier_transform_opt 256 3
time_fourier_transform hpce.tp1811.fast_fourier_transform_opt >> core_test.txt