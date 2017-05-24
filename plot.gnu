#!/usr/bin/gnuplot

set terminal png

set output 'Test_Case_1.png'
set xlabel 'Incident Angle'
set ylabel 'Reflectance'
plot 'First_Case_Fresnel.txt' u 1:3 w l title 'Fresnel', \
'FDTD_TEST/Ref_Trans_Air_n_1.5_0.1i_Air_d_400nm_lambda_500nm.txt' u 1:2 w l title 'FDTD'
