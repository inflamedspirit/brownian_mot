#!/usr/bin/perl
#PBS -N motsim
#PBS -e motsim.err
#PBS -o motsim.out
#PBS -m abe
#PBS -q normal
#PBS -l nodes=4:ppn=8
#PBS -l walltime=24:00:00

# -N Name
# -e stderr file for THIS script
# -o stdout file for THIS script
# -m send emails (abe = all the emails)
# -q not sure
# -l nodes and processers per node (always choose 8 ppn)
# -l time for job walltime=hh:mm:ss

### above PBS settings are: name = sample, std error -> sample.err,
###   std out -> sample.log, mail notification to the job owner,
###   queue = normal, request for 16 processors on 2 nodes (i.e.
###   8 jobs per node)

use POSIX ":sys_wait_h";
use strict;


# data output for universal stuff
my $ejection_time_file = "times_ejection";
my $return_time_file = "times_return";
my @ejection_time;
my @return_time;

### set locations
#     assume the executable names are of the form
#     a.out.1, a.out.2, ..., a.out.n
#     or
#     executable file.1, executable file.2, ...
my $workdir_base = $ENV{'PBS_O_WORKDIR'};
# on the next line: put everything in the command,
#   the dot and last number will be appended to this string
#   (this could also be a.out for the other example above)
my $num_jobs = 8000; # number of jobs to run
system("cd $workdir_base; ./makeparameters.pl $num_jobs");
my $job_base = "$workdir_base/runparam.pl PARAMETER";

open FILEA, ">$workdir_base/$ejection_time_file" or die "Error, output ejection times.";
open FILEB, ">$workdir_base/$return_time_file" or die "Error, output return times.";



### get node file name from shell
my $pbs_nodefile = $ENV{'PBS_NODEFILE'};

### log stuff to standard error
my ($hostname, $number_of_procs);
chomp($hostname = `hostname`);
chomp($number_of_procs = `wc -l <$pbs_nodefile`); $number_of_procs =~ s/\s//g;
print STDERR "Working base directory is $workdir_base\n";
print STDERR "Running on host $hostname\n";
print STDERR "Time is " . localtime(time) . "\n";
print STDERR "Allocated $number_of_procs nodes\n";

### make list of jobs to do, and set up other lists
my $j; my @job_array; my @child_array; my @running_job_array;
for ($j=1; $j<=$num_jobs; $j++) {
  push @job_array, $j
}

### make a list of available processors
my  @proc_array = ();
open(procfile, "<$pbs_nodefile");
while(<procfile>) { chomp; push @proc_array, $_; }
close(procfile);

### make a job table to keep track of which processors are working;
#     a value of 1 means busy, a value of 0 means idle
my @busy_table = ();
for ($j=0; $j<$number_of_procs; $j++) { push @busy_table, 0; }

### run the jobs!
my $running_jobs = 0;
while ( $#job_array >= 0 or $running_jobs > 0 ) {

  # spawn new jobs, if needed
    while ( $running_jobs < $number_of_procs and $#job_array >= 0  ) {
	my $job = shift @job_array;
	my $job_command = "$job_base.$job"; #OKAY... is is clever, the $job appends to the input parameter file

    # find first free job slot and launch
	$j = 0; while($busy_table[$j] == 1) {$j++;}
	if ( $j >= $number_of_procs ) {
	    die "Internal error with allocation table";
	}
	$busy_table[$j] = 1;
	print STDERR "on $proc_array[$j]:\n";
	$job_command = "ssh -n $proc_array[$j] 'cd $workdir_base; $job_command'";
    print STDERR "Launching job #$job: \n$job_command\n";
	print STDERR "\n";

	my $child = fork();
	if ($child == 0) {
	    exec($job_command);
	} else {
	    $child_array[$j] = $child;
	    $running_job_array[$j] = $job;
	}
	$running_jobs++;
    }

  # check for finished jobs
    for ($j=0; $j<=$#child_array; $j++) {
	if ( $busy_table[$j] == 1 ) {
	    my $child = waitpid($child_array[$j], WNOHANG);
	    if ( $child == -1 ) {
		my $exitstatus = int($?+0.49)/256;
        print STDERR "Finished job #$running_job_array[$j]";
		my $finished_job = $running_job_array[$j];

		

		open FILE, "<$workdir_base/stderr.$finished_job" or die "no error file for $finished_job\n";
		while(<FILE>){
		    if(/return_time/){
			$_=<FILE>;
			if(/diffusion_time\s*=\s*(\d+.\d+)/){
			    push(@return_time,$1);
			    print FILEA "$1\n";
			}
		    }

		if(/ejection_time/){
		    $_=<FILE>;
		    if(/diffusion_time\s*=\s*(\d+.\d+)/){
			push(@ejection_time,$1);
			print FILEB "$1\n";
		    }
		}

		}
		close FILE;

		if ( $exitstatus != 0 ) {print STDERR ", exit status = $exitstatus";}
		print STDERR "\n\n";
		$child_array[$j] = -1;
		$running_job_array[$j] = -1;
		$running_jobs--;
		$busy_table[$j] = 0;
	    }
	}
    }



    



  # wait a bit so we don't overload the processor with admin stuff
    sleep 1;

}

close FILEA;
close FILEB;
