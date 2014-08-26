ejectiondata = "times_ejection"
returndata = "times_return"

scale=2.62350132*10**(-8)
bin_width=160000
bin_number(x)=floor(x/bin_width)
rounded(x) = bin_width*(bin_number(x))

set table 'temp1.dat'
plot ejectiondata u (scale*rounded($1)):(1) smooth frequency
unset table

set table 'temp2.dat'
plot returndata u (scale*rounded($1)):(1) smooth frequency
unset table

