#!/usr/bin/perl

######################################################################################
#
# correlation.pl - parallel processing of correlation on binary genotype data
#
#
#    Copyright (C) 2012  Katie Elizabeth Hyma (keh233@cornell.edu)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
######################################################################################


use strict;
use warnings;
use Statistics::R;
use Statistics::Descriptive;
use POSIX qw(WNOHANG);
use Getopt::Long;


my $usage = "perl correlation.pl [options]\n\n";
$usage = $usage."\t-b binaryfile (required) the progeny binary file, the result of get_pseudotestcross_markers.pl\n";
$usage = $usage."\t-t number of threads to use (defaults to 1)\n";
$usage = $usage."\t-o output file base (defaults to binaryfile name)\n";
$usage = $usage."\t-d temporary directory to hold output (defaults to \"temp\")\n";
$usage = $usage."\t-s chunk size for processing markers (defaults to 5,000) change depending on amount of RAM and threads available\n";

my ($infile, $maxthreads, $tempdir, $outfile, $chunk);

if (scalar(@ARGV==0)) {
	die "\n\n$usage\n\n";
}

GetOptions(
	'b:s' => \$infile,
	't:i' => \$maxthreads,
	'o:s' => \$outfile,
	'd:s' => \$tempdir,
	's:i' => \$chunk,
);

unless (defined $infile) {
	die "\n\nmust specify binary input file with option -i\n\n$usage\n\n";
}
unless (defined $maxthreads) {
	$maxthreads = 1;
}
unless (defined $outfile) {
	$outfile = $infile.".r2";
}
unless (defined $tempdir) {
	$tempdir = "temp";
}
if (-e $outfile) {
	die "\n\noutput file $outfile already exists\n\n";
}
if (-d $tempdir) {
	die "\n\ntemporary directory $tempdir already exists\n\n";
}
unless (defined $chunk) {
	$chunk = 5000;
}

if(scalar(@ARGV) >0) { die "\n\nunknown options @ARGV\n\n$usage\n\n"; }

my ($logfile, $indexfile) = ("$outfile.log", "$outfile.markerlist");
open LOG, ">$logfile";
open INDEX, ">$indexfile";

print LOG "$0\n\noptions in effect:\n-b $infile\n-t $maxthreads\n-o $outfile\n-d $tempdir\n-s $chunk\n";
my $time = localtime();
my $inmd5 = `md5sum $infile `;
print LOG "\ninput file $infile has md5sum:\n $inmd5\n";
print LOG "script executed at localtime: $time\n\n";

system "mkdir $tempdir";

my $status = Statistics::R->new();
$status->run(qq`read.table(\"$infile\", header=T)->s3`);
$status->run(qq`s3[,6:dim(s3)[2]]->datas3`);
$status->run(qq`rownames(datas3)<-s3\$ID`);
$status->run(qq`datas3t<-t(datas3)`);
#loop over blocks to do the correlation!!  
my $nmarkers = $status->get('dim(datas3t)[2]');
my $names = $status->get('colnames(datas3t)');
my $chrs = $status->get('s3[,1]');
my $pos = $status->get('s3[,2]');
$status->stop();

print LOG "$nmarkers total markers\n\n";

for (my $i = 0; $i < $nmarkers; $i++) {
	print INDEX "@{$names}[$i]\t@{$chrs}[$i]\t@{$pos}[$i]\n";
}
close INDEX;

my @subs = ();
my %paste;
my %paste2;

my ($last, $lasti)  = (1,1);
for (my $n=1; $n<=$nmarkers - $chunk; $n=$n+$chunk) {
	for (my $i=1;$i<=$nmarkers-$chunk; $i=$i+$chunk) {
		my $t = 0;
		unless ($n == $i) {
			$t = 1;
		}
		my @coords = ($n, $n+$chunk-1, $i, $i+$chunk-1, $t); #0 or 1 value indicates whether transposition is necessary
		push (@subs, [@coords]);
		$last = $n + $chunk;
		$lasti = $i + $chunk;
		my $filename = $coords[0]."_".$coords[1]."_".$coords[2]."_".$coords[3]; 
		$paste{$coords[0]}{$coords[2]} = "$tempdir/$filename.r";
		$paste2{$coords[0]}{$coords[2]} = "$tempdir/$filename.r2";
		if ($t ==1) {
			my $tfilename = $coords[2]."_".$coords[3]."_".$coords[0]."_".$coords[1]; 
			$paste{$coords[2]}{$coords[0]} = "$tempdir/$tfilename.r";
			$paste2{$coords[2]}{$coords[0]} = "$tempdir/$tfilename.r2";
		}
	}
	#add in the last bit
	my @lastcoords = ($n, $n+$chunk-1, $lasti, $nmarkers, 1);
	push (@subs, [@lastcoords]);
	my $lastfilename = $lastcoords[0]."_".$lastcoords[1]."_".$lastcoords[2]."_".$lastcoords[3]; 
	$paste{$lastcoords[0]}{$lastcoords[2]} = "$tempdir/$lastfilename.r";
	$paste2{$lastcoords[0]}{$lastcoords[2]} = "$tempdir/$lastfilename.r2";
	my $tfilename = $lastcoords[2]."_".$lastcoords[3]."_".$lastcoords[0]."_".$lastcoords[1]; 
	$paste{$lastcoords[2]}{$lastcoords[0]} = "$tempdir/$tfilename.r";
	$paste2{$lastcoords[2]}{$lastcoords[0]} = "$tempdir/$tfilename.r2";
}

my @lastcoords = ($last, $nmarkers, $last, $nmarkers, 0);
push (@subs, [@lastcoords]);
my $lastfilename = $lastcoords[0]."_".$lastcoords[1]."_".$lastcoords[2]."_".$lastcoords[3]; 
$paste{$lastcoords[0]}{$lastcoords[2]} = "$tempdir/$lastfilename.r";
$paste2{$lastcoords[0]}{$lastcoords[2]} = "$tempdir/$lastfilename.r2";

my $ntasks = scalar(@subs);

print LOG "$ntasks total tasks\n";
print LOG "task\tstart_x\tstop_x\tstart_y\tstop_y\ttranspose\n";

for (my $n=0; $n<scalar(@subs); $n++) {
	print LOG "$n\t@{$subs[$n]}\n";
}

#fork children (this takes a long time)
fork_child_processes(\&corr, $maxthreads, $ntasks, $infile);

print LOG "\ncalculating correlations:\n";
$time = localtime();
print LOG "localtime: $time\n\n";


sub corr {
	my ($infile, $i) = @_;
	$i = $i-1;
	my $status = Statistics::R->new();
	$status->run(qq`library(WGCNA)`);
	$status->run(qq`allowWGCNAThreads(nThreads=$maxthreads)`);
	$status->run(qq`read.table(\"$infile\", header=T)->s3`);
	$status->run(qq`s3[,6:dim(s3)[2]]->datas3`);
	$status->run(qq`rownames(datas3)<-s3\$ID`);
	$status->run(qq`datas3t<-t(datas3)`);
	my ($start1, $stop1, $start2, $stop2, $t) = @{$subs[$i]}; 
	print LOG "$i\t$start1\t$stop1\t$start2\t$stop2\t$t\n";
	$status->run(qq`cor(datas3t[,$start1:$stop1],datas3t[,$start2:$stop2], use=\"pairwise.complete.obs\") -> subr`);
	$status->run(qq`subr^2 -> subr2`);
	my $filename = $start1."_".$stop1."_".$start2."_".$stop2; 
	$status->run(qq`write.table(subr2, \"$tempdir/$filename.r2\", sep=\"\t\", quote=FALSE, col.names=FALSE, row.names=FALSE)`);
	if ($t == 1) {
		my $tfile = $start2."_".$stop2."_".$start1."_".$stop1;
		$status->run(qq`write.table(t(subr2), \"$tempdir/$tfile.r2\", sep=\"\t\", quote=FALSE, col.names=FALSE, row.names=FALSE)`);
	}
	$status->stop();
	exit;
}

	
#Now put the files back together...
my %cat;
my %cat2;

print LOG "pasting:\n";
$time = localtime();
print LOG "localtime: $time\n\n";

foreach my $n (sort {$a <=> $b} keys %paste2) {
	$cat2{$n} = "$tempdir/$n.r2";
	my $count = 0;
	foreach my $i (sort {$a <=> $b} keys %{$paste2{$n}}) {
		print LOG "pasting: $paste2{$n}{$i}\n";
		if ($count ==0) { 
			$status = system "paste $paste2{$n}{$i} > $tempdir/$n.r2";
		}
		elsif ($count > 0) {
			$status = system "paste $tempdir/$n.r2 $paste2{$n}{$i} > $tempdir/temp";
			$status = system "mv $tempdir/temp $tempdir/$n.r2";
		}
		$count++;
	}
}


#now the column together	

print LOG "\ncat-ing:\n";
$time = localtime();
print LOG "localtime: $time\n\n";

foreach my $n (sort {$a <=> $b} keys %cat2) {
	print LOG "cating $cat2{$n}\n";
	$status = system "cat $cat2{$n} >> $outfile";
}

$status = system "rm -r $tempdir";

print LOG "\ndone\n";

$time = localtime();
my $outmd5 = `md5sum $outfile `;
my $outlistmd5 = `md5sum $outfile.markerlist `;
print LOG "\nproduced the following output files:\n";
print LOG "$outmd5$outlistmd5\n\n";
print LOG "finished at localtime: $time\n\n";

close LOG;

sub fork_child_processes {

	#fork child processes
	my $child = shift;
	my $maxthreads = shift;
	my $ntasks = shift;
	my @params = @_;
	
	my @procs;
	my $task = 0;
	for(my $i=1; $i<=$maxthreads; $i++)
	{
		$task++;
		unless ($task > $ntasks) {
		
			print STDERR "starting child $i task $task ";
			my $pid = fork();
			if($pid < 0)
			{
			#error
				print STDERR "\n\nERROR: Cannot fork child $i\n";
				for(my $j=0; $j<=$#procs; $j++)
				{
					system("kill -9 " . $procs[$j]);
				}
				exit;
			}
			if($pid == 0)
			{
			#child code
				$child->(@params, $task);
				#child_exec($i);
				exit;
			}
			#master - continue, $pid contains child pid
			$procs[$i-1] = $pid;
			print STDERR " pid $pid\n";
		}
	}
	
	#waiting for child processes to finish and execute remaining tasks
	while(1)
	{
		sleep(1);
		my $n=0;
		for(my $i=0; $i<=$#procs; $i++)
		{
			if($procs[$i] != 0)
			{
	        		my $kid = waitpid($procs[$i], WNOHANG);
	        		#print "STATUS $kid " . $procs[$i] . " ($n)\n";
				if($kid <= 0)
				{
					#process exists
					$n++;
				}
				else
					{
					print STDERR "Child " . ($i+1) . " finished (pid=" . $procs[$i] . ")\n";
					$procs[$i] = 0 ;
					if($task < $ntasks)
					{
						$task++;
						my $pid = fork();
						if($pid < 0)
						{
						#error
							print STDERR "\n\nERROR: Cannot fork child $i\n";
							for(my $j=0; $j<=$#procs; $j++)
							{
								system("kill -9 " . $procs[$j]);
							}
							exit;
						}
						if($pid == 0)
						{
						#child code
							#child_exec($i);
							$child->(@params, $task);
							exit;
						}
						#master - continue, $pid contains child pid
						$procs[$i-1] = $pid;
						
						print STDERR " child " . ($i+1) . " restarted for task $task with pid $pid\n";
					}
				}
			}	
		}
		if($n==0){last;}
	}
	
	#print "ALL DONE\n";
}
	
sub child_exec {
	(my $num) = @_;	

	my $seed = time * $num * $num;
	srand($seed);

	my $nsec = int(rand(20)) + 5;
	sleep($nsec);	
}
