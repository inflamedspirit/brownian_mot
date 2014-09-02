N = 1


# title "                                                                                                                                                                      MOT Simulation command: ".titletext
#set tmargin 5


#Conversions:
GammaLambda = 29.7405          # Multiply to convert velocity to velocity in m/s
GammaInverse = 2.6235*10**(-8) # Multiply to convert time to time in s
Lambda = 7.80241*10**(-7)      # Multiply to convert length to length in m
EnergyConv = 5.10858*10**7     # Multiply to convert energy to energy in h^2k^2/2m for k=2pi/lambda


#set x2range [GammaInverse*position_xmin : GammaInverse*position_xmax]
set x2tics nomirror
set xtics nomirror

set title "Fluorescence Rate (photons/sec)"

plot datafile every N u 1:10 w l t "photons"
