#!/bin/gnuplot
# Usage: gnuplot plot.p
# Requirements: That there is a folder named output_data with the used files in

cd "output_data"
set terminal png size 1200,1000

set yrange [-200:-20]
set output '../carrier_sensing_range.png'
set title "Received Signal Strength at Range in Corresponding Environment"
set label 1 gprintf("%gm",28) at 25, -84
set label 2 gprintf("%gm",50) at 46.4, -84 
set xlabel "Distance (m)"
set ylabel "Received power (dBm)"
set arrow 1 from 0,-82 to 100,-82 nohead lw 2 lc 5
set arrow 2 from 28.3,-20.2 to 28.3,-200 nohead lw 2 lc 1 dt 0
set arrow 3 from 50,-20.2 to 50,-200 nohead front lw 2 lc 3 dt 0
plot "office.csv" title "Office/Residential" with lines lw 2, "LoS.csv" title "Line of Sight" with lines lw 2, "shoppingmall.csv" title "Shoppingmall" with lines lw 2

reset
set output '../los_area_data_rate.png'
set title "Area Data Rate Line of Sight"
set xlabel "Deployment Density (AP/m^{2})"
set ylabel "Area Data Rate (Mbps/m^{2})"
set yrange [0:0.3]
plot "los1" title "1 channel" with linespoints, "los2" title "2 channels" with linespoints, "los3" title "3 channels" with linespoints


reset
set output '../home_area_data_rate.png'
set title "Area Data Rate Residential Home"
set xlabel "Deployment Density (AP/m^{2})"
set ylabel "Area Data Rate (Mbps/m^{2})"
set yrange [0:0.3]
plot "home1" title "1 channel" with linespoints, "home2" title "2 channels" with linespoints, "home3" title "3 channels" with linespoints

reset
set output '../mall_area_data_rate.png'
set title "Area Data Rate Shopping Mall"
set xlabel "Deployment Density (AP/m^{2})"
set ylabel "Area Data Rate (Mbps/m^{2})"
set yrange [0:0.3]
plot "mall1" title "1 channel" with linespoints, "mall2" title "2 channels" with linespoints, "mall3" title "3 channels" with linespoints

reset
set table '/tmp/stats.dat'
set yrange [0:0.12]
plot "random_mall[20]"
unset table
plot '/tmp/stats.dat'
max1_y = GPVAL_DATA_Y_MAX
max1_x = GPVAL_DATA_X_MAX
set table '/tmp/stats.dat'
plot "random_mall[20, 20]"
unset table
plot '/tmp/stats.dat'
max2_y = GPVAL_DATA_Y_MAX
max2_x = GPVAL_DATA_X_MAX
set table '/tmp/stats.dat'
plot "random_mall[20, 20, 20]"
unset table
plot '/tmp/stats.dat'
max3_y = GPVAL_DATA_Y_MAX
max3_x = GPVAL_DATA_X_MAX
set label 1 gprintf("%.4fMbps/m^2", max1_y) at max1_x+0.001, max1_y+0.001
set label 2 gprintf("%.4fMbps/m^2", max2_y) at max2_x+0.001, max2_y+0.001
set label 3 gprintf("%.4fMbps/m^2", max3_y) at max3_x+0.001, max3_y+0.001
set output '../mall_random_area_data_rate.png'
set title "Area Data Rate Shopping Mall"
set xlabel "Deployment Density (AP/m^{2})"
set ylabel "Area Data Rate (Mbps/m^{2})"
plot "random_mall[20]" title "1 Channel" with linespoints, "random_mall[20, 20]" title "2 Channels" with linespoints, "random_mall[20, 20, 20]" title "3 Channels" with linespoints

reset
set table '/tmp/stats.dat'
set yrange [0:0.3]
plot "random_home[20]"
unset table
plot '/tmp/stats.dat'
max1_y = GPVAL_DATA_Y_MAX
max1_x = GPVAL_DATA_X_MAX
set table '/tmp/stats.dat'
plot "random_home[20, 20]"
unset table
plot '/tmp/stats.dat'
max2_y = GPVAL_DATA_Y_MAX
max2_x = GPVAL_DATA_X_MAX
set table '/tmp/stats.dat'
plot "random_home[20, 20, 20]"
unset table
plot '/tmp/stats.dat'
max3_y = GPVAL_DATA_Y_MAX
set output '../home_random_area_data_rate.png'
set title "Area Data Rate Residential Home 20MHz Channel Bandwidth"
set label 1 gprintf("%.4fMbps/m^2", max1_y) at max1_x+0.001, max1_y+0.001
set label 2 gprintf("%.4fMbps/m^2", max2_y) at max2_x+0.001, max2_y+0.001
set label 3 gprintf("%.4fMbps/m^2", max3_y) at max3_x+0.001, max3_y+0.001
set xlabel "Deployment Density (AP/m^{2})"
set ylabel "Area Data Rate (Mbps/m^{2})"
plot "random_home[20]" title "1 Channel" with linespoints, "random_home[20, 20]" title "2 Channels" with linespoints, "random_home[20, 20, 20]" title "3 Channels" with linespoints

reset
set yrange [0:0.3]
set table '/tmp/stats.dat'
plot "random_home[20, 20, 20]"
unset table
plot '/tmp/stats.dat'
max1_y = GPVAL_DATA_Y_MAX
max1_x = GPVAL_DATA_X_MAX
set table '/tmp/stats.dat'
plot "random_home[20, 40]"
unset table
plot '/tmp/stats.dat'
max2_y = GPVAL_DATA_Y_MAX
max2_x = GPVAL_DATA_X_MAX
set output '../home_random_area_data_rate_20vs40.png'
set title "Area Data Rate Comparison 20+20+20MHz vs 40+20MHz Bandwidth Channels"
set label 1 gprintf("%.4fMbps/m^2", max1_y) at max1_x+0.001, max1_y+0.001
set label 2 gprintf("%.4fMbps/m^2", max2_y) at max2_x+0.001, max2_y+0.001
set xlabel "Deployment Density (AP/m^{2})"
set ylabel "Area Data Rate (Mbps/m^{2})"
plot "random_home[20, 20, 20]" title "20, 20, 20 Bands" with linespoints, "random_home[20, 40]" title "20, 40 Bands" with linespoints

reset
set table '/tmp/stats.dat'
set yrange [0:0.03]
plot "random_los[20]"
unset table
plot '/tmp/stats.dat'
max1_y = GPVAL_DATA_Y_MAX
max1_x = GPVAL_DATA_X_MAX
set table '/tmp/stats.dat'
plot "random_los[20, 20]"
unset table
plot '/tmp/stats.dat'
max2_y = GPVAL_DATA_Y_MAX
max2_x = GPVAL_DATA_X_MAX
set table '/tmp/stats.dat'
plot "random_los[20, 20, 20]"
unset table
plot '/tmp/stats.dat'
max3_y = GPVAL_DATA_Y_MAX
max3_x = GPVAL_DATA_X_MAX
set output '../los_random_area_data_rate.png'
set title "Area Data Rate Line of Sight"
set label 1 gprintf("%.4fMbps/m^2", max1_y) at max1_x+0.001, max1_y+0.001
set label 2 gprintf("%.4fMbps/m^2", max2_y) at max2_x+0.001, max2_y+0.001
set label 3 gprintf("%.4fMbps/m^2", max3_y) at max3_x+0.001, max3_y+0.001
set xlabel "Deployment Density (AP/m^{2})"
set ylabel "Area Data Rate (Mbps/m^{2})"
plot "random_los[20]" title "1 Channel" with linespoints, "random_los[20, 20]" title "2 Channels" with linespoints, "random_los[20, 20, 20]" title "3 Channels" with linespoints

reset
set output '../data_rate_user.png'
set title 'Data Rate Attainable for Every User of an AP'
set xlabel 'Users [N]'
set ylabel 'Average Datarate Mbps/User'
plot "avg_user_datarate_per_AP" smooth bezier title "Average Datarate" with lines
