ejectiondata = "times_ejection"
returndata = "times_return"

scale=2.62350132*10**(-8)
bin_width=80000
bin_number(x)=floor(x/bin_width)
rounded(x) = bin_width*(bin_number(x))

set output "| head -n -2 > temp1.dat"
set table
plot ejectiondata u (scale*rounded($1)):(1) smooth frequency
unset table

set output "| head -n -2 > temp2.dat"
set table
plot returndata u (scale*rounded($1)):(1) smooth frequency
unset table

