# se corre gnuplot -persist plot.gp  
#
#
set key left top Left
set title "Tiempos de ordenamiento"
set xlabel "Numero de datos"
set ylabel "Tiempo (ms)"
set yrange [0.01:1e7]
set logscale xy
plot "tiempos.data" u 1:2 title "Selection" with lines, \
       "" u 1:3 title "Insertion" with lines, \
       "" u 1:4 title "Bubble" with lines, \
       "" u 1:5 title "Shell Knuth 1959" with lines, \
       "" u 1:6 title "Shell Wiki 1959" with lines, \
       "" u 1:7 title "Shell Frank y Lazarus, 1960" with lines, \
       "" u 1:8 title "Quick Básico" with lines, \
       "" u 1:9 title "Quick Insertion" with lines, \
       "" u 1:10 title "Quick  Mediana" with lines, \
       "" u 1:11 title "Merge Basico" with lines, \
       "" u 1:12 title "Merge Basico Mejorado" with lines, \
       "" u 1:13 title "Merge Bottom Up" with lines, \
       "" u 1:14 title "Heap Sort" with lines
set terminal jpeg enhanced
set output "salida.jpg"
replot
