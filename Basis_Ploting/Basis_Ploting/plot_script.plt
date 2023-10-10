set key autotitle columnhead
set title "Quadratic basis functions"
plot for [i=2:*] "data.txt" using 1:i with lines
