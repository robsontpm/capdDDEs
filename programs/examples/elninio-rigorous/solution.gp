set terminal png size 800,600
set output 'solution.png'
plot 'solution.dat' using 1:3 with lines, \
     'solution.dat' using 1:($3-$4):($3+$4) with filledcu notitle