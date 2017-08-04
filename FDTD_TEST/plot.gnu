#!/usr/bin/gnuplot

set terminal postscript enhanced color 'Helvetica' 25
set output 'Reflectance.eps'

set xlabel 'Angle of Incidence (Deg)'
set ylabel 'Reflectance'
plot 'File.txt' u 1:2 w l lw 3 title 'TMM Reflectance', \
'FDTD_File_Ref.txt' u 1:2 w l lw 3 title 'Lumerical Reflectance'
