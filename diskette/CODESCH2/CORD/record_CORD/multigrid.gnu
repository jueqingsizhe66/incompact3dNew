#set ter jpeg
set term postscript eps enhanced color
set out 'multifuck.eps'
#set xrange [-pi:pi]
 
# Uncomment the following to line up the axes
# set lmargin 6
 
# Gnuplot recommends setting the size and origin before going to
# multiplot mode
# This sets up bounding boxes and may be required on some terminals
#set size 1.5,1.5
#set origin 0,0
 
# Done interactively, this takes gnuplot into multiplot mode
# and brings up a new prompt ("multiplot >" instead of "gnuplot >")
set style data lines  # lines    points    linespoints
set key left top
set xlabel 'x'
set ylabel 'y'

#set xrange [0:1]
#set yrange [1:50]
set multiplot
 
# plot the first graph so that it takes a quarter of the screen
set size 0.33,0.4
set origin 0,0.5
#set label "Cord_a"
plot './coord1a.plo'
 
# plot the second graph so that it takes a quarter of the screen
set size 0.33,0.4
set origin 0,0
#set label "Cord_b"


plot './coord1b.plo'
 
# plot the third graph so that it takes a quarter of the screen
set size 0.33,0.4
set origin 0.33,0.5
#set label "Cord_c"
plot './coord1c.plo'
 
# plot the fourth graph so that it takes a quarter of the screen
set size 0.33,0.4
set origin 0.33,0
#set label "Cord_d"
plot './coord1d.plo'

 # plot the fifth graph so that it takes a quarter of the screen
set size 0.33,0.4
set origin 0.66,0.5
#set label "Cord_e"
plot './coord1e.plo'
 
# plot the sixth graph so that it takes a quarter of the screen
set size 0.33,0.4
set origin 0.66,0
#set label "Cord_f"
plot './coord1f.plo'
 
# On some terminals, nothing gets plotted until this command is issued
unset multiplot
 
# remove all customization
reset


#  多曲线打印
#plot y1(x) t "{/Symbol b} in low-{/Symbol w} regime" lt 2, \
#     y2(x) t "{/Symbol b} in high-{/Symbol w} regime" lt 3 , \
#     "data.txt" t "Exponential" lw 2 lt 1


