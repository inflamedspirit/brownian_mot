N = 1


#COLLECT PLOT DATA For Multiple Axes Labels:

output_term = GPVAL_TERM
set term unknown

plot datafile every N u 1:2, datafile every N u 1:3, datafile every N u 1:4
position_ymin = GPVAL_Y_MIN
position_ymax = GPVAL_Y_MAX
position_xmin = GPVAL_X_MIN
position_xmax = GPVAL_X_MAX

plot datafile every N u 1:(sqrt(($2)*($2)+($3)*($3)+($4)*($4)))
radius_ymin = GPVAL_Y_MIN
radius_ymax = GPVAL_Y_MAX

set term x11

#FINISH COLLECTING PLOT DATA


set multiplot layout 2, 2

# title "                                                                                                                                                                      MOT Simulation command: ".titletext
set tmargin 5


#START = 0.0
#STOP = 0.4
#set xrange [START:STOP]

#Conversions:
GammaLambda = 29.7405          # Multiply to convert velocity to velocity in m/s
GammaInverse = 2.6235*10**(-8) # Multiply to convert time to time in s
Lambda = 7.80241*10**(-7)      # Multiply to convert length to length in m
EnergyConv = 5.10858*10**7     # Multiply to convert energy to energy in h^2k^2/2m for k=2pi/lambda


set x2range [GammaInverse*position_xmin : GammaInverse*position_xmax]
set y2range [Lambda*position_ymin : Lambda*position_ymax]
set y2tics nomirror
set ytics nomirror
set x2tics nomirror
set xtics nomirror

set title "Position (multiples of 780nm / m)"
set size 0.5, 0.5
set origin 0.0, 0.5
plot datafile every N u 1:2 w l t "x", datafile every N u 1:3 w l t "y", datafile every N u 1:4 w l t "z"
#plot datafile every N u ($1*GammaInverse):2 w l t "x", datafile every N u ($1*GammaInverse):3 w l t "y", datafile every N u ($1*GammaInverse):4 w l t "z"

unset y2range
unset y2tics

#set title "Velocity (multiples of 780nm*Gamma / m/s)"
#set size 0.5, 0.5
#set origin 0.0, 0.0
#plot datafile every N u 1:5 w l t "vx", datafile every N u 1:6 w l t "vy", datafile every N u 1:7 w l t "vz"

set title "State"
set ytics ("foc vol." 0, "diffusing." 1, "returned." 2, "ejected." 3)
set size 0.5, 0.5
set origin 0.0, 0.0
plot datafile every N u 1:9 w l t "Escape State"

unset ytics
set ytics scale default
set y2range [Lambda*radius_ymin : Lambda*radius_ymax]
set y2tics nomirror
set ytics nomirror
set title "Radius (multiples of 780nm / m)"
set size 0.5, 0.5
set origin 0.5, 0.5
set yrange [*:*]
plot datafile every N u 1:(sqrt(($2)*($2)+($3)*($3)+($4)*($4))) w l t "|r|"



unset y2range
unset y2tics
set title "Energy (recoil energies, hk^2/2m)"
set size 0.5, 0.5
set origin 0.5, 0.0
set yrange [*:*]
plot datafile every N u 1:($8*EnergyConv) w l t "H", datafile every N u 1:(EnergyConv*0.5*(($5)**2+($6)**2+($7)**2)) w l t "KE", datafile every N u 1:($8*EnergyConv-EnergyConv*0.5*(($5)**2+($6)**2+($7)**2)) w l t "V"