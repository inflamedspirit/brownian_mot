#!/usr/bin/perl -w
#use File::Glob;

my $ejection_time_file = "times_ejection";
my $return_time_file = "times_return";

my @logfiles = glob 'stderr.*';

my @ejection_time;
my @return_time;


print "Scanning through @logfiles stderr logfiles...\n";

for(@logfiles){

    open FILE, "<$_" or die "Error, no data file.\n";


    while(<FILE>){
	if(/return_time/){
	    $_=<FILE>;
	    if(/diffusion_time\s*=\s*(\d+.\d+)/){
		push(@return_time,$1);
	    }
	}

	if(/ejection_time/){
	    $_=<FILE>;
	    if(/diffusion_time\s*=\s*(\d+.\d+)/){
		push(@ejection_time,$1);
	    }
	}


    }
    close FILE;
}


open FILE, ">$ejection_time_file" or die "Error, output ejection times.";
for(@ejection_time){
    print FILE "$_\n";
}
close FILE;

open FILE, ">$return_time_file" or die "Error, output return times.";
for(@return_time){
    print FILE "$_\n";
}
close FILE;



