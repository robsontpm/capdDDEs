set terminal png size 1600,1200 lw 3 font ",20"
set output 'numeric.png'
set xrange [-1:20]
plot \
    'numeric-solutionddes-plot.dat' using 1:3 with lines title "solution",\
    'numeric-derivativeddes-plot.dat' using 1:3 with lines title "derivative"
