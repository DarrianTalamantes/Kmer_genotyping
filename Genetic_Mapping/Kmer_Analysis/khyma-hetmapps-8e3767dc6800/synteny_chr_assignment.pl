#!/usr/bin/perl

######################################################################################
#
# synteny_chr_assignment.pl tests accuracy of chromosome alignment
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
######################################################################################


use strict;
use warnings;
use Statistics::Descriptive;
use Data::Dumper;
use Getopt::Long;


my $usage = "perl synteny_chromosome_assignment.pl [options]\n\n";
$usage = $usage."\t-c correlation file (required), output from correlation_forked.pl\n";
$usage = $usage."\t-m marker list (required), output from correlation_forked.pl\n";
$usage = $usage."\t-b binary file (required), input to correlation_forked.pl\n";
$usage = $usage."\t-chr chromosome info file, a file with the columns chromosome# (integer values only)\n";
$usage = $usage."\t\t and length (only chromosomes you want to create linkage groups for)\n";
$usage = $usage."\t-diff ratio cutoff for correlation (default = 2), the number of times greater the mean correlation\n";
$usage = $usage."\t\tto the top ranking chr must be than the second ranking chr\n";
$usage = $usage."\t-o outfilebase (default = correlationfilename.diff[diff_cutoff])\n";
$usage = $usage."\t-v verbose (print info about each marker to logfile in addition to summary\n";

my ($corrfile, $markerlist, $binaryfile, $chrfile, $diff_cutoff, $outprefix, $verbose);

if (scalar(@ARGV==0)) {
	die "\n\n$usage\n\n";
}

GetOptions(
	'c:s' => \$corrfile,
	'm:s' => \$markerlist,
	'b:s' => \$binaryfile,
	'chr:s' => \$chrfile,
	'diff:f' => \$diff_cutoff,
	'o:s' => \$outprefix,
	'v' => \$verbose,
);


unless (defined $corrfile) {
	die "\n\nmust specify input corrleation file with option -c\n\n$usage\n\n";
}
unless (defined $chrfile) {
	die "\n\nmust specify chromosome info file with option -i\n\n$usage\n\n";
}
unless (defined $binaryfile) {
	die "\n\nmust specify binary file with option -b\n\n$usage\n\n";
}
unless (defined $diff_cutoff) {
	$diff_cutoff = 2;
}
unless (defined $outprefix) {
	$outprefix = $corrfile.".diff".$diff_cutoff;
}
if (-e "$outprefix.binary") {
	die "\n\noutput file $outprefix.binary already exists\n\n$usage\n\n";
}
unless (-e $corrfile) {
	die "\n\ncorrelation file $corrfile does not exist\n\n";
}
unless (-e $markerlist) {
	die "\n\nmarker list file $markerlist does not exits\n\n";
}
unless (-e $binaryfile) {
	die "\n\nbinary file $binaryfile does not exist\n\n";
}
unless (-e $chrfile) {
	die "\n\nchromosome info file $chrfile does not exist\n\n";
}


open CHRFILE, $chrfile || die "\n\ncould not open chromosome info file $chrfile\n\n";;
my %chrs;
while (my $line = <CHRFILE>) {
	chomp $line;
	my ($chr, $length) = split('\s+', $line);
	$chrs{$chr} = $length;
}
close CHRFILE;

open IN, "$corrfile" || die "\n\ncould not open correlation file $corrfile\n\n";
open MARKERS, "$markerlist" || die "\n\ncould not open marker list file $markerlist\n\n";

open LOG, ">$outprefix.log";

print LOG "$0\n\noptions in effect:\n-c $corrfile\n-m $markerlist\n-b $binaryfile\n-chr $chrfile\n-diff $diff_cutoff\n-o $outprefix\n";
if ($verbose) {
	print LOG "-v\n";
}

my $inmd5 = `md5sum $corrfile `;
my $inlist = `md5sum $markerlist `;
my $inchrinfomd5 = `md5sum $chrfile `;
print LOG "\ninput files:\n$inmd5$inlist$inchrinfomd5";
my $time = localtime();
print LOG "\nscript executed at localtime: $time\n\n";

if ($verbose) {
	print LOG "SNP\tCHR\tPOS\tSTATUS\tALIGNED_CHR\tMEAN_CORR_2_ALIGNED_CHR\tFIRST_CLUSTERED_CHR\tMEAN_CORR_2_FIRST_CLUSTERED_CHR\tSECOND_CLUSTERED_CHR\tMEAN_CORR_2_SECOND_CLUSTERED_CHR\tDIFFERENCE\tKEEP\n";
}

my (%row2chr, %row2id, %row2pos, $row, @cols, @coldata, @rank, @keepids);
	
#read the .corr file, for each remaining line and get the chromosome 
my $i=0; #note we are starting at line 0...

while (<MARKERS>) {
    chomp;
	my @line = split("\t", $_);
	my $chr = $line[1];
    $row2chr{$i}=$chr;
	$row2id{$i} = $line[0];
	$row2pos{$i} = $line[2];
	$i++;
}
close MARKERS;

open IN, "$corrfile";
	
my ($j, $chrom, $chrom2, $k, $truemembers, $falsemembers, $undetect, $chromosome) = (0,0,0,0,0,0,0,0);
my ($pos, $r) = ("","");
my %pos2value;
my $m=0; # m is the line number
my ($agree, $disagree, $unresolved) = (0,0,0);
while (<IN>) { 
	my @tmp = ();
	chomp;
	@cols=split ("\t", $_);
	my $chr = $row2chr{$m};
	my $id = $row2id{$m};
	my $pos = $row2pos{$m};
	my @coldata = @cols;
	#get the correlation value lookup hash
	my %bychr = ();
	for ($k=0;$k<=$#coldata;$k++) {
		#find out which chr it is on and add it to the vector of correlations
		my $chr2 = $row2chr{$k};
		unless ($coldata[$k] eq "NA") {push (@{$bychr{$chr2}}, $coldata[$k]);}
	}
	#make a table with the mean correlation to other markers on all chromosomes
	my @tr; 
	for my $c2 (sort { $a <=> $b } keys %chrs) {
		push (@tr, $c2);
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data(@{$bychr{$c2}});
		my $meancorr = $stat->mean();
		if (defined $meancorr) {
			if ($meancorr eq "n/a") {
				$meancorr = 0.00001; #avoid sorting and divide by zero errors
			}
			else {
				$meancorr = $meancorr + 0.00001; #avoid sorting and divide by zero errors
			}
		}
		else { 
			$meancorr = 0.00001; #avoid sorting and divide by zero errors
		}
		push(@tmp, $meancorr);
	}
	$m++;

	#now that the line is processed, let's print the chrassign
	my @sorted = sort {$b <=> $a} @tmp;
	my @index = grep {$tmp[$_] == $sorted[0] } 0 .. $#tmp;
	my @second = grep {$tmp[$_] == $sorted[1] } 0 .. $#tmp;
	my $diff = sprintf "%.2f", $sorted[0]/$sorted[1];
	$sorted[0] = sprintf "%.2f", $sorted[0];
	$sorted[1] = sprintf "%.2f", $sorted[1];
	my $alignedcorr;
	if (defined $tmp[$chr-1]) {
		$alignedcorr = sprintf "%.4f", $tmp[$chr-1];
	}
	else { 
		$alignedcorr = "NA";
	}
	#if ($alignedcorr eq "") {$alignedcorr = "NA";}
	if (scalar(@index) > 1) { ## there is more than one first best match, so it is unresolved
		$unresolved++;
		if ($verbose) {
			print LOG "$id\t$chr\t$pos\t0\t0\t1\t$chr\t$alignedcorr\tNA\t$sorted[0]\tNA\t$sorted[1]\t$diff\t0\n";
		}
	}
	if (scalar(@index) == 1) { #there is just one best match, check it
		if ($diff >= $diff_cutoff) { #best match and second best are more than the diff cutoff, it is resolved
			if ($chr == $index[0]+1) { #agrees
				$agree++;
				if ($verbose) {
					print LOG "$id\t$chr\t$pos\tagree\t$chr\t$alignedcorr\t$tr[$index[0]]\t$sorted[0]\t$tr[$second[0]]\t$sorted[1]\t$diff\t1\n";
				}
				#print the binary file
				push @keepids, $id;
			}
			if ($chr != $index[0]+1) { # disagrees
				$disagree++;
				if ($verbose) {
					print LOG "$id\t$chr\t$pos\tdisagree\t$chr\t$alignedcorr\t$tr[$index[0]]\t$sorted[0]\t$tr[$second[0]]\t$sorted[1]\t$diff\t0\n";
				}
			}
		}
		if ($diff < $diff_cutoff) { #didn't meet the cutoff, it is unresolved
			$unresolved++;
			if ($verbose) {
				print LOG "$id\t$chr\t$pos\tunresolved\t$chr\t$alignedcorr\t$tr[$index[0]]\t$sorted[0]\t$tr[$second[0]]\t$sorted[1]\t$diff\t0\n";
			}
		}			
	}	
}
	
close IN;

print LOG "\nat diff level $diff_cutoff:\n$agree agree\n$disagree disagree\n$unresolved unresolved\n\n";

open BINARYFILE, "$binaryfile" || die "\n\ncould not open binary file $binaryfile\n\n";
my $header = <BINARYFILE>;

open OUT, ">$outprefix.binary";

print OUT $header;
while (my $line = <BINARYFILE>) {
	my @line = split("\t", $line);
	my $id = $line[2];
	if ($id ~~ @keepids) {
		print OUT $line;
	}
}
close OUT;

my $outmd5 = `md5sum $outprefix.binary `;
print LOG "\nproduced the following output file:\n";
print LOG "$outmd5\n\n";
$time = localtime();
print LOG "finished at localtime: $time\n\n";

close LOG;
	
