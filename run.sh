#!/bin/bash

# Dataset:
dataset_title="long_run_2000d0p1"

# For this simulation

# 4-12 MHz is used as the detuning (~ unity detuning in Gamma)
# 75um diameter mot radius (~48.07 wavelengths)
# velocity is manually chosen to get the energy close to escaping, but not quite escaping (was done by checking simulations)
# threshhold was chosen randomly by what seems to work OK... not really certain how to pick it
# motradius was also chosen somewhat randomly, probably could figure out a value based on the beam diameters

reps=2000
for ((i=1;i<=reps;i++)); do
echo Starting rep $i...

# Parameters
rabi_frequency=0.707
detuning=-0.1
wavenumber=1.0
decay_rate=1.0
vx=0.0
vy=0.0
vz=0.0
tstep=0.1
tfinal=30000000
save_interval=1000
seed=$i
focalradius=48.07
focalthreshhold=54
motradius=250

# File Settings
daystring=$(date +%Y%m%d)
timestring=$(date +%H%M%S)
datafolder="./$daystring/$dataset_title/$timestring"
simulation_command="./mot $rabi_frequency $detuning $wavenumber $decay_rate $vx $vy $vz $tstep $tfinal $save_interval $seed"
main_data="main_data"
emission_data="emission_data"
png_plot="plot.png"

echo "Running MOT simulation..."

time ./mot $rabi_frequency $detuning $wavenumber $decay_rate $vx $vy $vz $focalradius $focalthreshhold $motradius $tstep $tfinal $save_interval $seed > $main_data

echo "Simulation complete."


echo "Plotting..."

# Create a png:
gnuplot -e "set term pngcairo size 1024,768; set output '$png_plot'; datafile='$main_data'" motplot.gnuplot

# Plot to screen:
#gnuplot --persist -e "set term 'x11'; datafile='$main_data'" motplot.gnuplot

echo "Plotting Complete."



echo "Saving data..."

if [ ! -d "$datafolder" ]; then
  mkdir -p $datafolder
fi

echo $simulation_command > $datafolder/simulation_command.sh
cp run.sh $datafolder/
cp log $datafolder/
cp motplot.gnuplot $datafolder/
cp $main_data $datafolder/
cp $png_plot $datafolder/
cp $png_plot $daystring/$daystring-$timestring-$dataset_title.png

echo "Saving Data Complete."

done