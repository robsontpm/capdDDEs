set terminal png size 1600,1200
set output 'phase-0-initial--ddes-plot.png'
plot 'phase-0-initial--ddes-plot.dat' using 1:($3-$4) with lines, 'phase-0-initial--ddes-plot.dat' using 1:($3+$4) with lines
