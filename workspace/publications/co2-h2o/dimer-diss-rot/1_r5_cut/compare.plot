set terminal postscript eps enhanced color font "Helvetica,22"
set out 'cut_compare.eps'
set encoding iso_8859_1


set xlabel 'r_{C--O} ({\305})' font "Helvetica,22"
set ylabel 'Energy (cm^{-1})' font "Helvetica,22"
set key at 15,800

plot 'rigid'  with lines lt rgb "black" lw 5 title 'Rigid cut' ,'rlx' with lines lt rgb "violet" lw 5 title 'Relaxed cut'
