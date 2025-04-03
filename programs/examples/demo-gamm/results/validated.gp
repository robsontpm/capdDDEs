set terminal png size 1600,1200 lw 3 font ",20"
set output 'validated.png'
set xrange [-1:20]
plot \
    'validated-solutionddes-plot.dat' using 1:($3-$4) with lines title "solution, lower bound", \
    'validated-solutionddes-plot.dat' using 1:($3+$4) with lines title "solution, upper bound"
