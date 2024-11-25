let SessionLoad = 1
let s:so_save = &g:so | let s:siso_save = &g:siso | setg so=0 siso=0 | setl so=-1 siso=-1
let v:this_session=expand("<sfile>:p")
silent only
silent tabonly
cd ~/Repositories/capdDDEs
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
let s:shortmess_save = &shortmess
if &shortmess =~ 'A'
  set shortmess=aoOA
else
  set shortmess=aoO
endif
badd +2 ~/Repositories/capdDDEs/.gitignore
badd +10048 term://~/Repositories/capdDDEs//1856:/usr/bin/zsh
badd +6 ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/Makefile
badd +28 ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/main.cpp
badd +79 ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/equation.h
badd +1 ~/Repositories/capdDDEs/programs/examples/rossler-ode-vs-dde-code/dde-vs-ode-code.cpp
badd +0 ~/Repositories/capdDDEs/programs/examples/rossler-ode-vs-dde-code/setup.h
badd +101 ~/Repositories/capdDDEs/include/capd/ddeshelper/DDEHelperNonrigorous.h
badd +612 ~/Repositories/capdDDEs/include/capd/ddeshelper/DDEHelperRigorous.h
badd +686 ~/Repositories/capdDDEs/include/capd/ddes/DDESolutionCurve.h
badd +67 ~/Repositories/capdDDEs/include/capd/ddes/DDESolutionCurve_PiecesManip.hpp
badd +236 ~/Repositories/capdDDEs/include/capd/ddes/DDEPoincareMap.h
badd +27 ~/Repositories/capdDDEs/bin/capd/include/capd/map/lib.h
badd +143 ~/Repositories/capdDDEs/bin/capd/include/capd/map/Map.h
badd +302 ~/Repositories/capdDDEs/bin/capd/include/capd/autodiff/NodeType.h
badd +682 ~/Repositories/capdDDEs/external/capd/capdAlg/include/capd/filib/Interval.h
badd +70 ~/Repositories/capdDDEs/bin/capd/include/capd/fadbad/fadbad.h
badd +36 ~/Repositories/capdDDEs/bin/capd/include/capd/intervals/Interval.h
badd +51 ~/Repositories/capdDDEs/bin/capd/include/capd/intervals/lib.h
badd +1 ~/Repositories/capdDDEs/bin/capd/include/capd/filib/Interval.h
badd +118 ~/Repositories/capdDDEs/include/capd/ddes/DDETaylorSolver.h
badd +127 ~/Repositories/capdDDEs/include/capd/ddes/DDESolutionCurve_Move.hpp
argglobal
%argdel
set stal=2
tabnew +setlocal\ bufhidden=wipe
tabnew +setlocal\ bufhidden=wipe
tabnew +setlocal\ bufhidden=wipe
tabnew +setlocal\ bufhidden=wipe
tabnew +setlocal\ bufhidden=wipe
tabnew +setlocal\ bufhidden=wipe
tabnew +setlocal\ bufhidden=wipe
tabrewind
edit ~/Repositories/capdDDEs/include/capd/ddeshelper/DDEHelperRigorous.h
argglobal
balt ~/Repositories/capdDDEs/include/capd/ddes/DDEPoincareMap.h
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 444 - ((26 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 444
normal! 018|
tabnext
edit ~/Repositories/capdDDEs/include/capd/ddes/DDESolutionCurve.h
let s:save_splitbelow = &splitbelow
let s:save_splitright = &splitright
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
let &splitbelow = s:save_splitbelow
let &splitright = s:save_splitright
wincmd t
let s:save_winminheight = &winminheight
let s:save_winminwidth = &winminwidth
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 30 + 105) / 210)
exe 'vert 2resize ' . ((&columns * 179 + 105) / 210)
argglobal
enew
file NvimTree_8
balt ~/Repositories/capdDDEs/include/capd/ddes/DDESolutionCurve.h
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal nofen
wincmd w
argglobal
balt ~/Repositories/capdDDEs/include/capd/ddes/DDESolutionCurve_Move.hpp
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 191 - ((46 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 191
normal! 011|
wincmd w
exe 'vert 1resize ' . ((&columns * 30 + 105) / 210)
exe 'vert 2resize ' . ((&columns * 179 + 105) / 210)
tabnext
edit ~/Repositories/capdDDEs/include/capd/ddes/DDETaylorSolver.h
let s:save_splitbelow = &splitbelow
let s:save_splitright = &splitright
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
let &splitbelow = s:save_splitbelow
let &splitright = s:save_splitright
wincmd t
let s:save_winminheight = &winminheight
let s:save_winminwidth = &winminwidth
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 30 + 105) / 210)
exe 'vert 2resize ' . ((&columns * 179 + 105) / 210)
argglobal
enew
file NvimTree_4
balt ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/main.cpp
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal nofen
wincmd w
argglobal
balt ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/main.cpp
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 118 - ((5 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 118
normal! 064|
wincmd w
exe 'vert 1resize ' . ((&columns * 30 + 105) / 210)
exe 'vert 2resize ' . ((&columns * 179 + 105) / 210)
tabnext
edit ~/Repositories/capdDDEs/include/capd/ddes/DDEPoincareMap.h
argglobal
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 303 - ((35 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 303
normal! 053|
tabnext
edit ~/Repositories/capdDDEs/programs/examples/rossler-ode-vs-dde-code/setup.h
let s:save_splitbelow = &splitbelow
let s:save_splitright = &splitright
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
let &splitbelow = s:save_splitbelow
let &splitright = s:save_splitright
wincmd t
let s:save_winminheight = &winminheight
let s:save_winminwidth = &winminwidth
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 104 + 105) / 210)
exe 'vert 2resize ' . ((&columns * 105 + 105) / 210)
argglobal
balt ~/Repositories/capdDDEs/programs/examples/rossler-ode-vs-dde-code/dde-vs-ode-code.cpp
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 39 - ((38 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 39
normal! 019|
wincmd w
argglobal
if bufexists(fnamemodify("~/Repositories/capdDDEs/programs/examples/rossler-ode-vs-dde-code/dde-vs-ode-code.cpp", ":p")) | buffer ~/Repositories/capdDDEs/programs/examples/rossler-ode-vs-dde-code/dde-vs-ode-code.cpp | else | edit ~/Repositories/capdDDEs/programs/examples/rossler-ode-vs-dde-code/dde-vs-ode-code.cpp | endif
if &buftype ==# 'terminal'
  silent file ~/Repositories/capdDDEs/programs/examples/rossler-ode-vs-dde-code/dde-vs-ode-code.cpp
endif
balt ~/Repositories/capdDDEs/programs/examples/rossler-ode-vs-dde-code/setup.h
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 155 - ((47 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 155
normal! 0
wincmd w
exe 'vert 1resize ' . ((&columns * 104 + 105) / 210)
exe 'vert 2resize ' . ((&columns * 105 + 105) / 210)
tabnext
edit ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/equation.h
let s:save_splitbelow = &splitbelow
let s:save_splitright = &splitright
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
let &splitbelow = s:save_splitbelow
let &splitright = s:save_splitright
wincmd t
let s:save_winminheight = &winminheight
let s:save_winminwidth = &winminwidth
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 104 + 105) / 210)
exe 'vert 2resize ' . ((&columns * 105 + 105) / 210)
argglobal
balt ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/main.cpp
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 67 - ((19 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 67
normal! 057|
wincmd w
argglobal
if bufexists(fnamemodify("~/Repositories/capdDDEs/programs/examples/elninio-rigorous/main.cpp", ":p")) | buffer ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/main.cpp | else | edit ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/main.cpp | endif
if &buftype ==# 'terminal'
  silent file ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/main.cpp
endif
balt ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/equation.h
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 28 - ((23 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 28
normal! 081|
wincmd w
exe 'vert 1resize ' . ((&columns * 104 + 105) / 210)
exe 'vert 2resize ' . ((&columns * 105 + 105) / 210)
tabnext
argglobal
if bufexists(fnamemodify("term://~/Repositories/capdDDEs//1856:/usr/bin/zsh", ":p")) | buffer term://~/Repositories/capdDDEs//1856:/usr/bin/zsh | else | edit term://~/Repositories/capdDDEs//1856:/usr/bin/zsh | endif
if &buftype ==# 'terminal'
  silent file term://~/Repositories/capdDDEs//1856:/usr/bin/zsh
endif
balt ~/Repositories/capdDDEs/programs/examples/elninio-rigorous/Makefile
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 10048 - ((47 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 10048
normal! 055|
tabnext
edit ~/Repositories/capdDDEs/bin/capd/include/capd/filib/Interval.h
let s:save_splitbelow = &splitbelow
let s:save_splitright = &splitright
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
let &splitbelow = s:save_splitbelow
let &splitright = s:save_splitright
wincmd t
let s:save_winminheight = &winminheight
let s:save_winminwidth = &winminwidth
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
exe 'vert 1resize ' . ((&columns * 30 + 105) / 210)
exe 'vert 2resize ' . ((&columns * 179 + 105) / 210)
argglobal
enew
file NvimTree_7
balt ~/Repositories/capdDDEs/bin/capd/include/capd/fadbad/fadbad.h
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal nofen
wincmd w
argglobal
balt ~/Repositories/capdDDEs/bin/capd/include/capd/fadbad/fadbad.h
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
let &fdl = &fdl
let s:l = 588 - ((23 * winheight(0) + 24) / 48)
if s:l < 1 | let s:l = 1 | endif
keepjumps exe s:l
normal! zt
keepjumps 588
normal! 025|
wincmd w
exe 'vert 1resize ' . ((&columns * 30 + 105) / 210)
exe 'vert 2resize ' . ((&columns * 179 + 105) / 210)
tabnext 7
set stal=1
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0 && getbufvar(s:wipebuf, '&buftype') isnot# 'terminal'
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20
let &shortmess = s:shortmess_save
let &winminheight = s:save_winminheight
let &winminwidth = s:save_winminwidth
let s:sx = expand("<sfile>:p:r")."x.vim"
if filereadable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &g:so = s:so_save | let &g:siso = s:siso_save
set hlsearch
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
