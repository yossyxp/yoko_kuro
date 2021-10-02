set term pngcairo size 800,800
#set term pdfcairo enhanced size 4in, 4in
#set object 1 rect behind from screen 0,0 to screen 1,1 fc rgb "#000000" fillstyle solid 1.0
#set border lc rgb "white"
#set key tc rgb "white"
unset key
#set xlabel "x" tc rgb "white"
#set ylabel "y" tc rgb "white" offset 1, 0
#set palette rgbformulae 22, 13, -31
#set style fill border 6
set xrange[-0.05:0.05]
set yrange[-0.05:0.05]
set isosamples 1000
do for [i=0:3]{
outfile = sprintf("yoko_kuro_mfs2%06d.png", i)
set output outfile
plot for[j=0:i] sprintf("yoko_kuro_mfs%06d.dat", 30 * j) using 1:2 w l lw 1.5 lc "blue"
#plot sprintf("yoko_kuro%06d.dat", i) using 1:2 lw 2 lc red
#plot sprintf("yoko_kuro%06d.dat", i) using 1:2 with filledcurves lc "white" linewidth 3
}