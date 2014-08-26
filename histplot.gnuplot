scale=2.62350132*10**(-8)

#lorentzian(x) = (g/2.0)/((x-xo)**2.0 + g*g/4.0)*h
#xo = 0.1
#g = 1
#h = 1
#fit lorentzian(x) 'temp2.dat' via g, h, xo
#
#lorentzian2(x) = (g2/2.0)/((x-xo2)**2.0 + g2*g2/4.0)*h2
#xo2 = 0.1
#g2 = 1
#h2 = 1
#fit lorentzian2(x) 'temp1.dat' via g2, h2, xo2

#exponential(x) = a*exp(-b*x)
#a = 1
#b = 1
#fit exponential(x) 'temp1.dat' via a, b
#
#exponential2(x) = a2*exp(-b2*x)
#a2 = 1
#b2 = 1
#fit exponential2(x) 'temp2.dat' via a2, b2
#
#levy(x) = s*sqrt(c/6.282)*exp(-c/(2.0*(x-mu)))/(x-mu)**(3.0/2.0)
#s=10
#c=5
#mu=-0.0342234
#fit levy(x) 'temp1.dat' via c, mu, s
#
#levy2(x) = s2*sqrt(c2/6.282)*exp(-c2/(2.0*(x-mu2)))/(x-mu2)**(3.0/2.0)
#s2=0.87
#c2=2.32223
#mu=-0.132123
#fit levy2(x) 'temp2.dat' via c2, mu2, s2
#

pow32(x) = q*x**(-o)
fit [0.2:*] pow32(x) 'temp1.dat' via q, o

pow32two(x) = q2*x**(-o2)
fit [0.001:*] pow32two(x) 'temp2.dat' via q2, o2

#set xrange [0:0.02]
#set yrange [0:100]
set logscale xy

set term x11
set output 'hist.png'
plot 'temp1.dat' u 1:2 w points t "Escape Times", 'temp2.dat' w points t "Return Times", pow32(x), pow32two(x)
# levy(x), levy2(x), exponential(x), lorentzian(x),  lorentzian2(x),



set term png
set output 'histlog.png'
set logscale xy
plot 'temp1.dat', 'temp2.dat'




