#!/bin/bash

cd bin
./$BINDIR/draw-nonrig --init-steps=-100 --draw-steps=-10 --kappa=11 --tau=1.2 --gamma=3 --alpha=1 '--ivp={1.1}' --dirpath=circle
./$BINDIR/draw-nonrig --init-steps=-100 --draw-steps=-10 --kappa=11 --tau=1.2 --gamma=3 --alpha=1 '--ivp={1.5}' --dirpath=clover