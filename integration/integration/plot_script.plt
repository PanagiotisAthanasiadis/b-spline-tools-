set key autotitle columnhead
set title "Aproximation of sin(2*pi*t)"
plot for [i=2:*] "Aproximations.txt" using 1:i with lines
