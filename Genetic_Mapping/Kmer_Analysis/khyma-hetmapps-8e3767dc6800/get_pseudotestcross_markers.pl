#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

######################################################################################
#
#    get_pseudotestcross_markers.pl 
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

######################################################################################
## set it up #########################################################################
######################################################################################

my $usage = "perl get_pseudotestcross_markers.pl [options]\n\n";
$usage = $usage."\t-vcf vcf (required) the vcf file, with only true progeny and progenitors, and in major/minor allele format\n";
$usage = $usage."\t-i indfile (optional) the file indicating parents and grandparents.\n";
$usage = $usage."\t\ttwo tab delimited columns, first column has the id (must match vcf file exactly)\n";
$usage = $usage."\t\tsecond column has the designation (either parent or grandparent)\n";
$usage = $usage."\t\tall others are assumed to be progeny.\n";
$usage = $usage."\t\tNote: if you create this file in windows and transfer it to a linux platform\n"; 
$usage = $usage."\t\tbe sure to save in tab delimited format and use the command dos2unix <filename> after transferring\n";
$usage = $usage."\t-o outfilebase. The base name of output files\n";
$usage = $usage."\t-error_maf sequencing error threshold (optional, default = 0.05).\n";
$usage = $usage."\t\t The maximum allele frequency an allele can be at to be considered a sequencing error\n";
$usage = $usage."\t-MAF (optional, default = 0.25) the expected minor allele frequency. Should be .25 for pseudo testcross markers\n";
$usage = $usage."\t-tol tolerance (optional, default = 0.125). The tolerance away from the expected allele or genotype frequencies (one-sided)\n";
$usage = $usage."\t\ti.e. tolerance of .1 for MAF of 0.5 wouild give a possible range of 0.4 to 0.6. Can't be set higher than .125\n";
$usage = $usage."\t-geno (optional, default = 0.5) minimum genotyping rate (after errror correction) to retain SNP\n";
$usage = $usage."\t-minGQ (optional, default = 98) minimum genotype quality score to retain a homozygote,\n";
$usage = $usage."\t\t also used for filtering progenitor genotypes\n";
$usage = $usage."\t-maxerr (optional, default = 0.05) maximum proportion of errors\n";
$usage = $usage."\t\t(not AA, AB, or BB genotypes and BB with GQ > minGQ) to total genotypes\n";
$usage = $usage."\t-v (optional, default is non-verbose) print out the result for each marker to the log file in addition to the summary\n";


my ($vcf, $indfile, $outfile, $error_threshold, $tolerance, $minGeno, $minGQ, $maxerr, $verbose, $MAF);

if (scalar(@ARGV)==0) {
	die "\n\n$usage\n\n";
}

GetOptions(
	'vcf:s' => \$vcf,
	'i:s' => \$indfile,
	'o:s' => \$outfile,
	'error_maf:f' => \$error_threshold,
	'tol:f' => \$tolerance,
	'geno:f' => \$minGeno,
	'minGQ:i' => \$minGQ,
	'maxerr:f' => \$maxerr,
	'v' => \$verbose,
	'MAF:f' => \$MAF
);

if (scalar(@ARGV) >0) {
	die "\n\nunknown options @ARGV\n\n$usage\n\n";
}
#check the params
unless (defined $vcf) {
	die "\n\nmust specify the vcf file with -v\n\n$usage\n\n";
}
unless (defined $outfile) {
	die "\n\nplease specify the outfile base name with -o\n\n$usage\n\n";
}
unless (defined $error_threshold) {
	$error_threshold = 0.05;
}
unless (defined $MAF) {
	$MAF = 0.25;
}
unless (defined $tolerance) {
	$tolerance = 0.125;
}
unless (defined $minGeno) {
	$minGeno = 0.5;
}
unless (defined $minGQ) {
	$minGQ = 98;
}
unless (defined $maxerr) {
	$maxerr = 0.05;
}

#check for overlap in ranges
unless (($MAF + $tolerance) <= (1 - $MAF - $tolerance)) {
	die "\n\nmajor and minor allele frequency ranges overlap for specified MAF and tolerance with MAF = $MAF and tolerance = $tolerance\n\n";
}


#check if input file exists
unless (-e $vcf) {
	die "\n\nvcf file $vcf does not exist\n\n$usage\n\n";
}

#open the files to write to
if (-e "$outfile.vcf") {
	die "\n\n$outfile.vcf already exists\n\n$usage\n\n";
}
if (-e "$outfile.binary") {
	die "\n\n$outfile.binary already exists\n\n$usage\n\n";
}
if (-e "$outfile.log") {
	die "\n\n$outfile.log already exists\n\n$usage\n\n";
}
if (-e "$outfile.markerinfo") {
	die "\n\n$outfile.markerinfo already exists\n\n$usage\n\n";
}

open LOG, ">>$outfile.log";
open OUTVCF, ">>$outfile.vcf";
open OUTBINARY, ">>$outfile.binary";
open OUTMARKERS, ">>$outfile.markerinfo";

print LOG "$0\n\noptions in effect:\n-vcf: $vcf\n-o: $outfile\n-error_maf: $error_threshold\n-MAF: $MAF\n-tol: $tolerance\n-geno: $minGeno\n-minGQ: $minGQ\n-maxerr: $maxerr\n";
if (defined $indfile) {
	print LOG "-i: $indfile\n";
}
print LOG "\n";
my $time = localtime();
my $inmd5 = `md5sum $vcf `;
print LOG "\ninput file $vcf has md5sum:\n$inmd5\n";
print LOG "script executed at localtime: $time\n\n";

######################################################################################
## check for parents and grandparents ################################################
######################################################################################

my %progenitors; # progenitors{name} = type
my %progenitors_index; #progenitors_index{num} = name

my ($parentfile, $grandparentfile);

if (defined $indfile) {
	print LOG "from $indfile:\n";
	open INDFILE, $indfile;
	while (my $line = <INDFILE>) {
		chomp $line; $line=~s/\s+$//; next LLL unless ($line=~/\w/);
		my @line = split("\t", $line);
		unless (scalar(@line) == 2) {
			die "the individual file $indfile must be a tab delimited file with two columns (id and type)\nwhere type can be 'parent', 'grandparent' or 'progeny'\n\nsaw: $line\n\n";
		}
		if ($line[1] eq 'progeny') {
			#skip it
		}
		elsif ($line[1] eq 'parent') { 
			$progenitors{$line[0]} = $line[1];
			$parentfile = "$outfile.parents.GQ$minGQ";
			print LOG "found parent $line[0]\n";
		}
		elsif ($line[1] eq 'grandparent') { 
			$progenitors{$line[0]} = $line[1];
			$grandparentfile = "$outfile.grandparents.GQ$minGQ";
			print LOG "found grandparent $line[0]\n";
		}
		else {
			die "invalid individual type $line[1] in indfile $indfile\n$line\nmust be 'parent', 'grandparent', or 'progeny'\n\n";
		}
	}
	close INDFILE;
}
print LOG "\n";

open VCF, $vcf || die "\n\ncould not open vcf file $vcf\n\n";

if (defined $parentfile) {
	open OUTPARENTVCF, ">>$parentfile.vcf";
	open OUTPARENTBINARY, ">>$parentfile.binary";
}
if (defined $grandparentfile) {
	open OUTGRANDPARENTVCF, ">>$grandparentfile.vcf";
	open OUTGRANDPARENTBINARY, ">>$grandparentfile.binary";
}

######################################################################################
## find the marker types #############################################################
######################################################################################

my $totalmarkers = 0;
my $pms = 0;
my @inds;

if ($verbose) { print LOG "MARKER_ID\tGenotyping_rate\tNalleles_sequenced\tNalleles_kept\tMajor_Minor_allele_freq\tsegregating_genotypes\tgenotype_frequencies\tgenotyping_error_rate\tgenotype_correction_rate\tgenotype_masking_rate\tRESULT\n";}


my %outcome_summary;
$outcome_summary{PtM} = 0;
$outcome_summary{nonPtM} = 0;
$outcome_summary{PtM_filtered_err_rate} = 0;
$outcome_summary{PtM_filtered_genorate} = 0;
$outcome_summary{invariant} = 0;
$outcome_summary{filtered_genorate} = 0;

LINE: while (my $line = <VCF>) {
	my %outcome;
	$outcome{id} = ""; 
	$outcome{genorate} = ""; 
	$outcome{alleles_seqd} = ""; 
	$outcome{alleles_kept} = ""; 
	$outcome{allelefreqs} = ""; 
	$outcome{genolist} = ""; 
	$outcome{genofreqs} = ""; 
	$outcome{hetfreq} = ""; 
	$outcome{err} = ""; 
	$outcome{corr} = ""; 
	$outcome{masked} = "";
	$outcome{outcome} = "";

	# first deal with the headers
	if ($line =~ /^\#CHROM/) {
		chomp $line;
		my @line = split ("\t", $line);
		print OUTMARKERS "$line[0]\t$line[1]\t$line[2]\n";
		@inds = @line[9..(scalar(@line)-1)];
		#check that the defined progenitors are in the file
		foreach my $k (keys %progenitors) {
			unless ($k ~~ @inds) {
				die "\n\ncould not find progenitor $k in VCF file $vcf\n\n";
			}
		}
		#print the proper info to the output files
		#VCF
		for (my $i = 0; $i<=8;$i++) {
			unless ($i==0) {
				print OUTVCF "\t";
			}
			print OUTVCF $line[$i];
			if (defined $parentfile) {
				unless ($i==0) {
					print OUTPARENTVCF "\t";
				}
				print OUTPARENTVCF $line[$i];
			}
			if (defined $grandparentfile) {
				unless ($i==0) {
					print OUTGRANDPARENTVCF"\t";
				}
				print OUTGRANDPARENTVCF $line[$i];
			}
			#binary
			if ($i <=4) {
				my $binaryout = $line[$i];
				$binaryout =~ s/\#//;
				unless ($i==0) {
					print OUTBINARY "\t";
				}
				print OUTBINARY $binaryout;
				if (defined $parentfile) {
					unless ($i==0) {
						print OUTPARENTBINARY "\t";
					}
					print OUTPARENTBINARY $binaryout;
				}
				if (defined $grandparentfile) {
					unless ($i==0) {
						print OUTGRANDPARENTBINARY "\t";
					}
					print OUTGRANDPARENTBINARY $binaryout;
				}
			}
		}
	
		for (my $i = 0; $i<=$#inds; $i++) {
			my $ind = $inds[$i];
			if (exists $progenitors{$ind}) {
				$progenitors_index{$i} = $ind;
				if ($progenitors{$ind} eq 'parent') {
					print OUTPARENTVCF "\t", $ind;
					print OUTPARENTBINARY "\t", $ind;
				}
				elsif ($progenitors{$ind} eq 'grandparent') {
					print OUTGRANDPARENTVCF "\t", $ind;
					print OUTGRANDPARENTBINARY "\t", $ind;
				}
				else {
					die "progenitor discrepancy for $ind type $progenitors{$ind}\n";
				}
			}
			else { #assume it is progeny
				print OUTVCF "\t", $ind;
				print OUTBINARY "\t", $ind;
			}
		}
		print OUTVCF "\n";
		print OUTBINARY "\n";
		if (defined $parentfile) {
			print OUTPARENTVCF "\n";
		}
		if (defined $grandparentfile) {
			print OUTGRANDPARENTVCF "\n";
		}
		if (defined $parentfile) {
			print OUTPARENTBINARY "\n";
		}
		if (defined $grandparentfile) {
			print OUTGRANDPARENTBINARY "\n";
		}
	}
	elsif ($line =~ /^\#/) {
		print OUTVCF $line;
		if (defined $parentfile) {
			print OUTPARENTVCF $line;
		}
		if (defined $grandparentfile) {
			print OUTGRANDPARENTVCF $line;
		}
	}

	#then the actual markers - 
	else {
		chomp $line;
		my @line = split("\t", $line);
		$line=~s/\t1\/0:/\t0\/1:/g;
		my ($chr, $pos, $id, $major, $minor, $info, $extra) = ($line[0], $line[1], $line[2], $line[3], $line[4], $line[7], (join "\t", @line[5..8]));
		$outcome{id} = $id;
		my $len = scalar(@line);
		my @genotypes=@line[9..$len-1];
		my $ntaxa = scalar(@genotypes) - scalar(keys %progenitors);
		
		#find the number of alleles
		my @alt = split(",", $minor);
		my $nalt = scalar(@alt);
		my $nalleles = $nalt + 1;
		$outcome{alleles_seqd} = $nalleles;
		
		#initialize allele counts and geno counts
		my @allele_counts;
		for (my $i = 0; $i < $nalleles; $i++) {
			$allele_counts[$i] = 0;
		}

		my ($all, $total, $nhets) = (0,0,0);
		my %genos;
		my %genos_new;
		my @seqdgenos = ();

		#find the allele counts, genotype counts and frequencies
		for (my $i = 0; $i<=$#genotypes; $i++) {
			if (exists $progenitors_index{$i}) {
				#skip it for calculating
				next;
			}
			else {
				my $genotype = $genotypes[$i];
				my @g = split(":", $genotype);
				my $geno = $g[0];
				my $dp = $g[2];
				unless($geno eq "./.") {
					if (exists $genos{$geno}) {
						$genos{$geno} = $genos{$geno} +1;
					}
					else {
						$genos{$geno} = 1;
						push @seqdgenos, $geno;
					}
					my ($a1, $a2) = split("/", $geno);
					$allele_counts[$a1]++;
					$allele_counts[$a2]++;
					$all++;
				}
			}
		}
		my $seqdgenos = join(":", @seqdgenos);
		my $nseqdgenos = scalar(@seqdgenos);

		#convert counts to frequencies, removing putative sequencing errors (low frequency alleles) from nalleles
		my $nalleles_new=0;
		for (my $i = 0; $i < $nalleles; $i++) {
			my $freq = $allele_counts[$i]/($all*2+0.000000001);
			if ($freq > $error_threshold) {
				$nalleles_new++;
			}
		}

		#recalculate the number of alleles, genotypes, frequencies, hets and depth after removing putative sequencing errors

		@allele_counts = ();
		for (my $i = 0; $i < $nalleles_new; $i++) {
			$allele_counts[$i] = 0;
		}
		
		for (my $i = 0; $i<=$#genotypes; $i++) {
			if (exists $progenitors_index{$i}) {
				#skip it for calculating
				next;
			}
			else {
				my $genotype = $genotypes[$i];
				my @g = split(":", $genotype);
				my $geno = $g[0];
				my $dp = $g[2];
				unless($geno eq "./.") {
					my ($a1, $a2) = split("/", $geno);
					if ($a1 < $nalleles_new && $a2 < $nalleles_new) { 
						if (exists $genos_new{$geno}) {
							$genos_new{$geno} = $genos_new{$geno} +1;
						}
						else {$genos_new{$geno} = 1;}
						$allele_counts[$a1]++;
						$allele_counts[$a2]++;
						$total++;
						if ($a1 ne $a2) {
							$nhets++;
						}
					}
					elsif ($a1 < $nalleles_new && $a2 >= $nalleles_new) {
						if (exists $genos_new{$geno}) {
							$genos_new{$geno} = $genos_new{$geno} + 1;
						}
						else {
							$genos_new{$geno} = 1;
						}
						$allele_counts[$a1]++;
						$total++;
					}
					elsif ($a1 >= $nalleles_new && $a2 < $nalleles_new) {
						if (exists $genos_new{$geno}) {
							$genos_new{$geno} = $genos_new{$geno} + 1;
						}
						else {
							$genos_new{$geno} = 1;
						}
						$allele_counts[$a2]++;
						$total++;
					}
					else {
					}
				}
			}
		}

		my $nallelesseqd = $nalleles;
		$nalleles = $nalleles_new;
		$outcome{alleles_kept} = $nalleles_new;

		#convert counts to frequencies, removing putative sequencing errors (low frequency alleles) from nalleles
		my @allele_freqs = ();
		for (my $i = 0; $i < $nalleles_new; $i++) {
			$allele_freqs[$i] = sprintf "%.3f", $allele_counts[$i]/($total*2 + 0.000000001);
		}
		$outcome{allelefreqs} = join(":",@allele_freqs);

		#mask sequencing errors in calculation of total
		my $het_freq = sprintf "%.3f", $nhets/($total+0.000000001);
		$outcome{hetfreq} = $het_freq;
		my $genotyping_rate = sprintf "%.3f", $total/$ntaxa;
		$outcome{genorate} = $genotyping_rate;
		if ($genotyping_rate <= $minGeno) {
			$outcome{outcome} = "filtered_genorate"; 
			$outcome_summary{filtered_genorate}++;
		}
		
		#find the number of segregating genotypes	
		my $seggenos = 0;
		my @genolist = ();
		my @genocount = ();
		foreach my $g (sort keys %genos_new) {
			if (($genos{$g} / $total+0.000000001) > $error_threshold) {
				$seggenos++;
				$seqdgenos++;
				push @genolist, $g;
				my $rounded = sprintf "%.2f", $genos{$g}/$total;
				push @genocount, $rounded;
			}
		}
		my $genolist = join(":", @genolist);
		my $genocount = join(":", @genocount);

		$outcome{genolist} = $genolist;
		$outcome{genofreqs} = $genocount;

		if ($seggenos == 0) {$genolist = "./.";}

		#get genotype frequencies
		my ($aa, $ab, $bb) = (0,0,0);
		if (exists $genos{"0/0"}) {$aa = sprintf "%.2f", $genos{"0/0"}/($total + 0.000000001);}
		if (exists $genos{"0/1"}) {$ab = sprintf "%.2f", $genos{"0/1"}/($total + 0.000000001);}
		if (exists $genos{"1/1"}) {$bb = sprintf "%.2f", $genos{"1/1"}/($total + 0.000000001);}

		if ($nalleles ==1 ) {
			#unless ($outcome{outcome} ne "" && $outcome{outcome} ne "filtered_genorate") {
			unless ($outcome{outcome} eq "filtered_genorate") {
				$outcome{outcome} = "invariant";
				$outcome_summary{invariant}++;
			}
		}
		elsif ($nalleles == 2 ) {
			if ($seggenos >= 2) {
				if ($allele_freqs[0] > (1-$MAF) - $tolerance && $allele_freqs[0] < (1- $MAF) + $tolerance && $allele_freqs[1] > $MAF - $tolerance && $allele_freqs[1] < $MAF + $tolerance) {
					#this is a PM
					my ($totalgeno, $maskedgeno, $keptgeno, $errorsgeno, $missinggeno, $correctedgeno) = (0,0,0,0,0,0);
					#print the progeny file, and parents and grandparents if they exist!!
					my %outbinary;
					my %outvcf;
					my $minoralleles = $line[4];
					$minoralleles =~ s/,.*//g;
					$outbinary{info} = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$minoralleles\t";
					$outvcf{info} = "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$minoralleles\t$extra\t";
					$outbinary{progeny} = "";
					$outvcf{progeny} = "";
					for (my $n=0;$n<=$#genotypes ;$n++){
						my $type;
						if (exists $progenitors_index{$n}) {
							$type = $progenitors{$progenitors_index{$n}};
							unless (exists $outbinary{$type}) {
								$outbinary{$type} = "";
								$outvcf{$type} = "";
							}
						}
						else {
							$type = "progeny";
						}
						if ($genotypes[$n]=~/^(\d+)\/(\d+)\:/) {
							my ($a1, $a2) = ($1, $2);
							my ($GT, $AD, $DP, $GQ, $PL) = split(":",$genotypes[$n]);
							$AD =~ s/(\d+,\d+).*/$1/;
							if ($a1 == 0 && $a2 == 0) { #homozygous major, could be het, filtering on GQ for all
								if ($type eq "progeny") {
									if ($GQ >= $minGQ) {
										$outbinary{$type} = $outbinary{$type}."0\t"; $keptgeno++; $totalgeno++;
										$outvcf{$type} = $outvcf{$type}.$genotypes[$n]."\t"; #kept it
									}
									else {
										$outbinary{$type} = $outbinary{$type}."NA\t"; $maskedgeno++; $totalgeno++;
										$outvcf{$type} = $outvcf{$type}."./.:".$AD.":".$DP.":".$GQ.":".$PL."\t";
									}
								}
								elsif ($type eq "parent" || $type eq "grandparent") {
									if ($GQ >= $minGQ) {
										$outbinary{$type} = $outbinary{$type}."0\t";
										$outvcf{$type} = $outvcf{$type}.$genotypes[$n]."\t"; #kept it
									}
									else {
										$outbinary{$type} = $outbinary{$type}."NA\t";
										$outvcf{$type} = $outvcf{$type}."./.:".$AD.":".$DP.":".$GQ.":".$PL."\t";
									}
								}
								else {
									die "invalid type of individual \"$type\"\n\n";
								}
							}
							elsif ($a1 == 0 && $a2 == 1) { # het, filter on minGQ for progenitors
								if ($type eq "progeny") {
									$outbinary{$type} = $outbinary{$type}."1\t"; $keptgeno++; $totalgeno++;
									$outvcf{$type} = $outvcf{$type}.$genotypes[$n]."\t";
								}
								elsif ($type eq "parent" || $type eq "grandparent") {
									if ($GQ >= $minGQ) {
										$outbinary{$type} = $outbinary{$type}."1\t";
										$outvcf{$type} = $outvcf{$type}.$genotypes[$n]."\t"; #kept it
									}
									else {
										$outbinary{$type} = $outbinary{$type}."NA\t";
										$outvcf{$type} = $outvcf{$type}."./.:".$AD.":".$DP.":".$GQ.":".$PL."\t";
									}
								}
								else {
									die "invalid type of individual \"$type\"\n\n";
								}
							}
							elsif ($a1 == 1 && $a2 == 1) { # the minor allele homozygote can't exist
								if ($type eq "progeny") {
									if ($GQ >= $minGQ) {
										$outbinary{$type} = $outbinary{$type}."NA\t"; $errorsgeno++; $totalgeno++;
										$outvcf{$type} = $outvcf{$type}."./.:".$AD.":".$DP.":".$GQ.":".$PL."\t";
									}
									else {
										$outbinary{$type} = $outbinary{$type}."1\t"; $correctedgeno++; $totalgeno++; #change it to het if low qual
										$outvcf{$type} = $outvcf{$type}."0/1:".$AD.":".$DP.":".$GQ.":".$PL."\t";
									}
								}
								elsif ($type eq "parent" || $type eq "grandparent") {
									if ($GQ >= $minGQ) {
										$outbinary{$type} = $outbinary{$type}."NA\t";
										$outvcf{$type} = $outvcf{$type}."./.:".$AD.":".$DP.":".$GQ.":".$PL."\t";
									}
									else { #correct it
										$outbinary{$type} = $outbinary{$type}."1\t"; #change it to het if low qual
										$outvcf{$type} = $outvcf{$type}."0/1:".$AD.":".$DP.":".$GQ.":".$PL."\t";
									}
								}
								else {
									die "invalid type of individual \"$type\"\n\n";
								}
							}
							else {
								if ($type eq "progeny") {
									$outbinary{$type} = $outbinary{$type}."NA\t"; $errorsgeno++; $totalgeno++;
									$outvcf{$type} = $outvcf{$type}."./.:".$AD.":".$DP.":".$GQ.":".$PL."\t";
								}
								elsif ($type eq "parent" || $type eq "grandparent") {
									$outbinary{$type} = $outbinary{$type}."NA\t";
									$outvcf{$type} = $outvcf{$type}."./.:".$AD.":".$DP.":".$GQ.":".$PL."\t";
								}
								else {
									die "invalid type of individual \"$type\"\n\n";
								}
							} 
						}
						else { # set it to missing
							if ($type eq "progeny") {
								$outbinary{$type} = $outbinary{$type}."NA\t"; $missinggeno++, $totalgeno++;
								$outvcf{$type} = $outvcf{$type}.$genotypes[$n]."\t";
							}
							elsif ($type eq "parent" || $type eq "grandparent") {
								$outbinary{$type} = $outbinary{$type}."NA\t";
								$outvcf{$type} = $outvcf{$type}.$genotypes[$n]."\t";
							}
							else {
								die "invalid type of individual \"$type\"\n\n";
							}
						}
					}
					#check that we still keep the locus
					my $keepit = 1;
					my $endrate = sprintf "%.2f", $keptgeno / $totalgeno;  
					my $enderr = sprintf "%.2f", $errorsgeno / ($totalgeno-$missinggeno);
					my $correctrate = sprintf "%.2f", $correctedgeno / ($totalgeno - $missinggeno);
					my $maskedrate = sprintf "%.2f", $maskedgeno / ($totalgeno - $missinggeno);
					$outcome{genorate} = $endrate;
					$outcome{err} = $enderr;
					$outcome{corr} = $correctrate;
					$outcome{masked} = $maskedrate;

					if ($outcome{outcome} eq "filtered_genorate" || $outcome{outcome} eq "invariant") {
						$keepit = 0;
					}
					if ($outcome{outcome} ne "filtered_genorate" && $outcome{outcome} ne "invariant" && $endrate <= $minGeno) {
						$outcome{outcome} = "PtM_filtered_genorate";
						$outcome_summary{PtM_filtered_genorate}++;
						$keepit = 0;
					}
					if ($outcome{outcome} ne "filtered_genorate" && $outcome{outcome} ne "invariant" && $outcome{outcome} ne "PtM_filtered_genorate" && $enderr > $maxerr) {
						$outcome{outcome} = "PtM_filtered_err_rate";
						$outcome_summary{PtM_filtered_err_rate}++;
						$keepit = 0;
					}

					if ($keepit==1) {
						$outcome{outcome} = "PtM";
						$outcome_summary{PtM}++;
						foreach my $t (keys %outbinary) {
							unless ($t eq "info") {
								$outbinary{$t} =~ s/\t&//;
								$outvcf{$t} =~ s/\t&//g;
								$outbinary{$t} = $outbinary{$t}."\n";
								$outvcf{$t} = $outvcf{$t}."\n";
							}
						}
						print OUTBINARY $outbinary{info}.$outbinary{progeny};
						print OUTVCF $outvcf{info}.$outvcf{progeny};
						if (defined $parentfile) {
							print OUTPARENTBINARY $outbinary{info}.$outbinary{parent};
							print OUTPARENTVCF $outvcf{info}.$outvcf{parent};
						}
						if (defined $grandparentfile) {
							print OUTGRANDPARENTBINARY $outbinary{info}.$outbinary{grandparent};
							print OUTGRANDPARENTVCF $outvcf{info}.$outvcf{grandparent};
						}
						print OUTMARKERS "$chr\t$pos\t$id\n";

					}
				}
				else {
					unless ($outcome{outcome} eq "filtered_genorate" || $outcome{outcome} eq "invariant") {
						$outcome{outcome} = "nonPtM";
						$outcome_summary{nonPtM}++;
					}
				}
			}
			else {
				unless ($outcome{outcome} eq "filtered_genorate" || $outcome{outcome} eq "invariant") {
					$outcome{outcome} = "nonPtM";
					$outcome_summary{nonPtM}++;
				}
			}
		}
		else {
			unless ($outcome{outcome} eq "filtered_genorate" || $outcome{outcome} eq "invariant") {
				$outcome{outcome} = "nonPtM";
				$outcome_summary{nonPtM}++;
			}
		}
		if ($verbose) {
			foreach my $p (keys %outcome) {
				if ($outcome{$p} eq "") {
					$outcome{$p} = "NA";
				}
			}
			print LOG $outcome{id}."\t";
			print LOG $outcome{genorate}."\t";
			print LOG $outcome{alleles_seqd}."\t";
			print LOG $outcome{alleles_kept}."\t";
			print LOG $outcome{allelefreqs}."\t";
			print LOG $outcome{genolist}."\t";
			print LOG $outcome{genofreqs}."\t";
			print LOG $outcome{err}."\t";
			print LOG $outcome{corr}."\t";
			print LOG $outcome{masked}."\t";
			print LOG $outcome{outcome}."\n";
		}
	}
}	
close VCF;

print LOG "\n\nsummary of marker types found:\n\n";
print LOG "$outcome_summary{filtered_genorate} markers filtered due to low genotyping rate prior to testing marker type\n";
print LOG "$outcome_summary{invariant} invariant sites removed from dataset\n";
print LOG "$outcome_summary{nonPtM} markers passed initial filters but do not fit segregation patterns for pseudotestcross markers\n";
my $PtMtot = $outcome_summary{PtM} + $outcome_summary{PtM_filtered_genorate} + $outcome_summary{PtM_filtered_err_rate};
print LOG "\n$PtMtot total potential pseudotestcross markers found\n";
print LOG "$outcome_summary{PtM_filtered_genorate} of those markers were filtered due to low genotyping rate after masking low quality genotypes\n";
print LOG "$outcome_summary{PtM_filtered_err_rate} of those markers were filtered due to high genotyping error rate\n";
print LOG "$outcome_summary{PtM} total pseudotestcross markers remain and were output\n\n";

my $outmd5binary = ` md5sum $outfile.binary `;
my $outmd5vcf = ` md5sum $outfile.vcf `;
print LOG "generated the following output files:\n$outmd5binary$outmd5vcf";
if (defined $parentfile) {
	my $outmd5binaryparents = ` md5sum $outfile.parents.GQ$minGQ.binary `;
	my $outmd5vcfparents = ` md5sum $outfile.parents.GQ$minGQ.vcf `;
	print LOG $outmd5binaryparents;
	print LOG $outmd5vcfparents;
}
if (defined $grandparentfile) {
	my $outmd5binarygrandparents = ` md5sum $outfile.grandparents.GQ$minGQ.binary `;
	my $outmd5vcfgrandparents = ` md5sum $outfile.grandparents.GQ$minGQ.vcf `;
	print LOG $outmd5binarygrandparents;
	print LOG $outmd5vcfgrandparents;
}

$time = localtime();
print LOG "\n\nfinished at localtime: $time\n\n";

close LOG;

