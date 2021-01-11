#set term pngcairo size 1000,1000
set term pdfcairo enhanced size 4in, 4in
set object 1 rect behind from screen 0,0 to screen 1,1 fc rgb "#000000" fillstyle solid 1.0
set border lc rgb "white"
set key tc rgb "white"
unset key
set xlabel "x" tc rgb "white"
set ylabel "y" tc rgb "white" offset 1, 0
set palette rgbformulae 22, 13, -31
set style fill border 6
set xrange[-20:20]
set yrange[-20:20]
set isosamples 1000
outfile = "snow.pdf"
#outfile = "snow.png"
set output outfile
plot for[i=0:200] sprintf("snow%06d.dat", i) using 1:2 with lines lc "skyblue" lw 0.5	
#plot for[i=0:123] sprintf("snow%06d.dat", i) using 1:2 with filledcurves lc "white" lw 3