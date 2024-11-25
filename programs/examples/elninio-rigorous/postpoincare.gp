set terminal png size 800,600
set output 'postpoincare.png'
plot 'postpoincare.dat' using 1:3 with lines, \
     'postpoincare.dat' using 1:($3-$4):($3+$4) with filledcu notitle