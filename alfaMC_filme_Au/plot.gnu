#set terminal x11 size 800,800

set title 'Energy in detector'
set xlabel "E (MeV)"
set ylabel "events"
#set logscale y
plot 'Edet.csv' u 1:2 w histeps lc 2
pause mouse 'click mouse'


set title 'Energy in target'
set xlabel "E (MeV)"
set ylabel "events"
#set logscale y
plot 'Etar.csv' u 1:2 w histeps lc 2
pause mouse 'click mouse'

