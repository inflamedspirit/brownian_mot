#!/usr/bin/perl


use strict;

if( @ARGV != 1 ) {
    die("Needs number of paramfiles to generate as input.");
}

# File and job variables
my $num_jobs = $ARGV[0];
my $stdoutbase = "stdout";
my $stderrbase = "stderr";

# Unused Parameters
my $rabi_frequency=0.707;
my $detuning=-0.12;
my $wavenumber=1.0;
my $vx=0.0;
my $vy=0.0;
my $vz=0.0;

# Base Parameters
my $tstep=0.1;
my $tfinal=30000000;
my $save_interval=100000;
my $step_size=0.1;
my $seed=0;
my $focalradius=48.07;
my $focalthreshhold=1.1*$focalradius;
my $motradius=20.0*$focalradius;

my $j;
for ($j=1; $j<=$num_jobs; $j++) {

# Adjusted parameters
    $seed=$seed+1;

    open(my $paramfh, ">PARAMETER.$j");
    print $paramfh "# file automatically generated by makeparameters.pl\n";
    print $paramfh "\$stdoutfile='$stdoutbase.$j';\n";
    print $paramfh "\$stderrfile='$stderrbase.$j';\n";
    print $paramfh "\n";
    print $paramfh "# PARAMETERS.$j\n";
    print $paramfh "\$rabi_frequency=$rabi_frequency;\n";
    print $paramfh "\$detuning=$detuning;\n";
    print $paramfh "\$wavenumber=$wavenumber;\n";
    print $paramfh "\$step_size=$step_size;\n";
    print $paramfh "\$vx=$vx;\n";
    print $paramfh "\$vy=$vy;\n";
    print $paramfh "\$vz=$vz;\n";
    print $paramfh "\$tstep=$tstep;\n";
    print $paramfh "\$tfinal=$tfinal;\n";
    print $paramfh "\$save_interval=$save_interval;\n";
    print $paramfh "\$seed=$seed;\n";
    print $paramfh "\$focalradius=$focalradius;\n";
    print $paramfh "\$focalthreshhold=$focalthreshhold;\n";
    print $paramfh "\$motradius=$motradius;\n";
    print $paramfh "\n";
    close($paramfh);

}
