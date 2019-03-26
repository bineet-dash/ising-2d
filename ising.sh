g++ ising2dim.cpp -o ising2d

read -p "Enter lattice size: " size
read -p "Enter number of sweeps: " no_sweeps

./ising2d  $size $no_sweeps
# gnuplot -p -e "set xlabel \"Temperature (T) \"; set ylabel \"average magnetization(M) \"; set autoscale xy; plot 'data.txt' using 1:3 with lines notitle;  quit"
# gnuplot -p -e "set xlabel \"Temperature (T) \"; set ylabel \"average Energy(E) \"; set autoscale xy; plot 'data.txt' using 1:2 with lines notitle;  quit"

# read -p "Enter lattice size again: " size
gnuplot -e "latticeSize=${size}" plot_results.gnuplot

rm ising2d
