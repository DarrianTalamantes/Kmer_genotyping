#!/usr/bin/perl
use Getopt::Long;


######################################################################################
#
# major_minor_VCF.pl - output in major/minor format and recalculate NS,DP,AF
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

my $usage = "perl major_minor_VCF.pl [options]\n\n";
$usage = $usage."\t-vcf the input vcf file (required), with only true F1 progeny and progenitors\n";
$usage = $usage."\t-o The base name of the output file. Defaults to inputfileprefix.majorminor.vcf\n";

my ($vcf, $outfile);

GetOptions(
	'vcf:s' => \$vcf,
	'o:s' => \$outfile,
);

unless (defined $vcf) {
	die "\n\nmust specify vcf input file with option -vcf\n\n$usage\n\n";
}
unless (defined $outfile) {
	my $prefix = $vcf;
	$prefix =~ s/\.vcf//;
	$outfile = $prefix.".majorminor.vcf";
}
if (scalar(@ARGV)>0) {
	die "\n\nunknown options @ARGV\n\n$usage\n\n";
}
if (-e $outfile) {
	die "\n\noutput file $outfile already exists\n\n";
}
if (-e "$outfile.log") {
	die "\n\nlogfile $outfile.log already exists\n\n";
}

# we we output in major/minor format rather than ref/alt format and recalculate NS,DP and AF
#new line reads: CHR POS ID REF ALT QUAL FILTER INFO(NS=;DP=;AF=) FORMAT(GT:AD:DP:GQ:PL) genotypes
	
open VCF, $vcf;
open OUTFILE, ">$outfile";
open LOGFILE, ">$outfile.log";

print LOGFILE "$0\n@ARGV\n\n";
my $time = localtime();
my $inmd5 = ` md5sum $vcf `;
print LOGFILE "input file $vcf has md5sum:\n$inmd5\n\n";
print LOGFILE "script executed at localtime: $time\n\n";

my @inds;
my ($all_missing, $monomorphic, $biallelic, $multiallelic, $difforder) = (0,0,0,0,0);	

	
while (<VCF>) {
	#change the header lines that need it (ID=AD, ID=PL, ID=AF);
	if ($_ =~ /^\#+FORMAT=<ID=AD/) {
		print OUTFILE "\#\#FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the major and minor alleles in the order listed\">\n";
	}
	elsif ($_ =~ /^\#+FORMAT=<ID=PL/) {
		print OUTFILE "\#\#FORMAT=<ID=PL,Number=3,Type=Float,Description=\"Normalized, Phred-scaled likelihoods for AA,AB,BB genotypes where A=major and B=minor; not applicable if site is not biallelic\">\n";
	}
	elsif ($_ =~ /^\#+INFO=<ID=AF/) {
		print OUTFILE "\#\#INFO=<ID=AF,Number=.,Type=Float,Description=\"Minor Allele Frequency\">\n";
	}
	elsif ($_ =~ /^\#+INFO=<ID=DP/) {
		print OUTFILE "\#\#INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (including filtered genotypes)\">\n";
	}
	elsif ($_ =~ /^\#/) {
		print OUTFILE $_;
	}
	else {
		chomp $_;
		my @line = split("\t", $_);
		my ($chr, $pos, $name, $ref, $alt, $qual, $filter, $info, $format) = @line[0..8];
		my @genotype_data = @line[9..scalar(@line-1)];
		my ($NS, $DP, $AF) = (0,0,0);
		my @refalt = ($ref, split(",", $alt));
		my $nalleles = scalar(@refalt);
		#loop through the genotypes and count up the number of times each genotype is seen to figure out if the alleles need changing
		my @allelecounts = (); #initialize the counts
		my @alleledepths;
		for (my $i=0; $i<$nalleles; $i++) {
			$allelecounts[$i] = 0;
		}
		for (my $i=0; $i<scalar(@genotype_data); $i++) { 
			#calculate the counts
			my @data = split(":", $genotype_data[$i]);
			if ($data[0] eq "./.") {
				#do nothing
			}
			else {
				$NS++;
				my @genotypes = split("/", $data[0]);
				$allelecounts[$genotypes[0]]++;
				$allelecounts[$genotypes[1]]++;
				#calculate NS, DP and AFs
				my @allele_depths = split(",", $data[1]);
				for (my $n=0;$n <scalar(@allele_depths); $n++) {
					$DP++;
				}
			}
			
		}
		my $nalleles_new = 0;
		my %alleles;
		#now that we have the counts, find out if the number and order matches ref/alt
		for (my $n=0; $n<$nalleles; $n++) {
			if ($allelecounts[$n] > 0) {
				$nalleles_new++;
				$alleles{$n} = $allelecounts[$n];
			}
			else {
				#do nothing
			}
		}
		if ($nalleles_new == 0) {
			$all_missing++;
		}
		elsif ($nalleles_new == 1) { #if it is now monomorphic, don't print it?
			$monomorphic++;
		}
		elsif ($nalleles_new  > 1) {
			if ($nalleles_new ==2) {$biallelic++;}
			if ($nalleles_new >2) {$multiallelic++;}
		
			my %remap_alleles = ();
			my %remap_alleles_reverse = ();
			for (my $n = 0; $n<$nalleles_new; $n++) {
				$remap_alleles{$n}{value} = 0;
				$remap_alleles{$n}{orig} = 9;
			}
			
			
			GO: foreach (my $j=0; $j<$nalleles; $j++) {
				my $assigned = 0;
				for (my $n=0; $n<$nalleles_new; $n++) {
					if ($allelecounts[$j] > $remap_alleles{$n}{value}) {
						for (my $z=$nalleles_new-1; $z>$n;$z--) {
							$remap_alleles{$z}{value} = $remap_alleles{$z-1}{value};
							$remap_alleles{$z}{orig} = $remap_alleles{$z-1}{orig};
						}
						$remap_alleles{$n}{value} = $allelecounts[$j];
						$remap_alleles{$n}{orig} = $j;
						next GO;
					}
					else{}
				}	
			}
			
			foreach my $keys (sort {$a <=> $b} keys %remap_alleles) {
				if ($remap_alleles{$keys}{orig} == 9) {
					delete $remap_alleles{$keys};
				}
			}
			
			foreach my $keys (sort {$a <=> $b} keys %remap_alleles) {
				$remap_alleles_reverse{$remap_alleles{$keys}{orig}} = $keys;
			}
			 
			#calculate AF and check for order and number
				
			$AF = sprintf("%.2f", $remap_alleles{1}{value} / ($NS*2));
			my $major = $refalt[$remap_alleles{0}{orig}];
			my $minor = $refalt[$remap_alleles{1}{orig}];
			#print  "oldalleles:$nalleles\tnewalleles:$nalleles_new\t@refalt\n";
			if ($nalleles_new == 2) {
				print OUTFILE "$chr\t$pos\t$name\t$major\t$minor\t$qual\t$filter\tNS=$NS;DP=$DP;AF=$AF\t$format\t";
			}
			elsif ($nalleles_new > 2) {
				print OUTFILE "$chr\t$pos\t$name\t$major\t";
				for (my $n=1; $n<$nalleles_new-1; $n++) {
					print OUTFILE "$refalt[$remap_alleles{$n}{orig}],"
				}
				print OUTFILE "$refalt[$remap_alleles{$nalleles_new -1}{orig}]\t$qual\t$filter\tNS=$NS;DP=$DP;AF=$AF\t$format\t";
			}
			
			
			my $order = 1;
			my $number = 1;
			
			if ($nalleles != $nalleles_new) {$number = 0;}
			for (my $n=0; $n<$nalleles_new-1; $n++) {
				unless ($remap_alleles{$n}{orig} < $remap_alleles{$n+1}{orig}) {
					$order = 0;
				}
			}
			
			if ($order == 0) {$difforder++;}

			#print out the genotype data
			for (my $n=0; $n<scalar(@genotype_data); $n++) {
			
				if ($number == 1 && $order == 1) { #nothing needs changing
					print OUTFILE "$genotype_data[$n]";
				}
				else { #either the number of order is different
					my ($gt, $ad, $dp, $gq, $pl) = split(":", $genotype_data[$n]);
					my ($gt1, $gt2) = split("/", $gt);
					my @ad = split(",", $ad);
					my ($newgt1, $newgt2) = ($remap_alleles_reverse{$gt1}, $remap_alleles_reverse{$gt2});
					
					unless ($gt eq "./.") {}
					if ($gt eq "./.") { #it is missing
						#print "$genotype_data[$n]\t";
						print OUTFILE "$gt:.:.:.:.";
						if ($n < scalar(@genotype_data) -1) {print OUTFILE "\t";}
						#else {print "\n";}
						next;
					}
					
					elsif ($newgt1 < $newgt2 || $newgt1 == $newgt2) { #if the new g1 is less than new g2, het or homozyg
						print OUTFILE "$newgt1/$newgt2:";
					}
					elsif ($newgt1 > $newgt2) { #het reversed
						print OUTFILE "$newgt2/$newgt1:";
					}
					
					
					for (my $j =0; $j<$nalleles_new-1; $j++) { #print the reordered alleles depths
						print OUTFILE "$ad[$remap_alleles{$j}{orig}],";
					}
					print OUTFILE "$ad[$remap_alleles{$nalleles_new-1}{orig}]:$dp:$gq:";
					
					my ($pl1, $pl2, $pl3) = split (",", $pl);
					
					if ($newgt1 == $newgt2) { 
						if ($newgt1 == 0) { #homozygous major
							if ($pl1 < $pl3) {
								print OUTFILE "$pl";
							}
							elsif ($pl1 > $pl3) {
								print OUTFILE "$pl3,$pl2,$pl1";
							}
						}
						elsif ($newgt1 != 0) { #homozygous alternative
							if ($pl1 <$pl3) {
								print OUTFILE "$pl3,$pl2,$pl1";
							}
							elsif ($pl1 > $pl3) {
								print OUTFILE "$pl";
							}
						}
					}
					elsif ($newgt1 != $newgt2) {
						# it is a het, do we need to switch the pl order?
						if ($remap_alleles{$newgt1}{orig} < $remap_alleles{$newgt2}{orig}) {
							print OUTFILE "$pl";
						}
						else {
							print OUTFILE "$pl3,$pl2,$pl1";
						}
					}
					else {print OUTFILE "$name\tumm, error!\n";}
								
				}
				if ($n < scalar(@genotype_data) -1) {
					print OUTFILE "\t";
				}	
				
			}
			print OUTFILE "\n";
		}
		
	}
}
			
close VCF;
close OUTFILE;
my $outmd5 = ` md5sum $outfile `;
print LOGFILE "generated output file $outfile with md5sum:\n$outmd5\n\n";
print LOGFILE "$all_missing sites with all missing data were not output to $outfile\n";
print LOGFILE "$monomorphic sites are monomorphic and were not output to $outfile\n";
my $outmarkers = $biallelic + $multiallelic;
print LOGFILE "\na total of $outmarkers sites were output to $outfile, of which:\n";
print LOGFILE "$biallelic sites are biallelic\n";
print LOGFILE "$multiallelic site are multiallelic\n";
print LOGFILE "$difforder sites changed major/minor order\n";
$time = localtime();
print LOGFILE "\n\nfinished at localtime: $time\n\n";

close LOGFILE;
			
			
####################################################################################################			


			
			
			
			
			
			
			
			




