set terminal png size 800,600
set output 'prepoincare.png'
plot 'prepoincare.dat' using 1:3 with lines, \
     'prepoincare.dat' using 1:($3-$4):($3+$4) with filledcu notitle