#!/usr/bin/perl


use strict;

# Output files
our $stdoutfile;
our $stderrfile;

# Declare Parameters to be loaded:
our $rabi_frequency;
our $detuning;
our $wavenumber;
our $step_size;
our $vx;
our $vy;
our $vz;
our $tstep;
our $tfinal;
our $save_interval;
our $seed;
our $focalradius;
our $focalthreshhold;
our $motradius;

# Load Parameters
if(@ARGV != 1 ){
    die "Needs parameter file as argument.";
}
my $paramfile = $ARGV[0];
require $paramfile;

# Launch Mot
system("time ./mot $rabi_frequency $detuning $wavenumber $step_size $vx $vy $vz $focalradius $focalthreshhold $motradius $tstep $tfinal $save_interval $seed 1> $stdoutfile 2> $stderrfile ");
