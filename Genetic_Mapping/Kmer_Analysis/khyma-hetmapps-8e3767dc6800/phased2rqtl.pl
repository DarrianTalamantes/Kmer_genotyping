#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

######################################################################################
#
# phased2rqtl.pl separate minor allele contribution by phase 
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


my $usage = "perl phased2rqtl.pl\n\n";
$usage = $usage."\t-b binaryfile, required\n";
$usage = $usage."\t-l linkage group / marker file (phased.txt), required\n";
$usage = $usage."\t-c cross type (\"BC\" or \"4way\"), required\n\t\tif using BC format, input to Rqtl with the option genotypes=c(\"A\",\"B\")\n";
$usage = $usage."\t-map LGmap (mapping of LG / phase onto some meaningful names, optional)\n";
$usage = $usage."\t\tLGmap has 5 tab delimited columns: LG, PHASE, LABEL, PARENT_NUM, GRANDPARENT_NUM\n";
$usage = $usage."\t\tLG and PHASE must match the linkage group file,\n";
$usage = $usage."\t\tLABEL is the label that you would like to propogate down the pipeline\n";
$usage = $usage."\t\tif you want to join linkage groups this is the place to do it, prior to ordering\n";
$usage = $usage."\t\tyou can use the information from the GENOvPARENT and GENOvGRANDPARENT files to create the LGmap\n";
$usage = $usage."\t\tPARENT_NUM must be 1 or 2 (up to 2 parents per LG - expect 1 per LG\n"; 
$usage = $usage."\t\tGRANDPARENT_NUM must be 1 or 2 (up to 2 grandparents per LG, expect one per LG/PHASE\n";
$usage = $usage."\t-order new order, tab delimited file with ID, LG, POS, can use all.n.mstmap from phased2mstmap.pl (optional)\n";
$usage = $usage."\t-drop_markers markers to drop, file with one marker name per line (optional)\n";
$usage = $usage."\t\tIf the iterative ordering option in phased2mstmap was used, if iteration N is chosen,\n\t\tuse the output all.[N-1].markers2drop to drop the same markers that were dropped for ordering in iteration N.\n";
$usage = $usage."\t-drop_genotypes genotypes to drop, tab delimited file with markerID, taxaID (optional)\n";
$usage = $usage."\t\tIf the iterative ordering option in phased2mstmap was used, if iteration N is chose,\n\t\tuse the output all.[N-1].genotypes2drop to drop the same genotypes that were dropped for ordering in iteration N.\n";
$usage = $usage."\t-o outfilebase (defaults to binaryfile.rqtl)\n\t\toutput files are in the format \"csvsr\" for Rqtl\n";
$usage = $usage."\t-v verbose (print detailed output)\n";

my ($binaryfile, $markerfile, $crosstype, $lgmapfile, $order, $markers2drop, $genotypes2drop, $outfilebase, $verbose);

GetOptions(
	'b=s' => \$binaryfile,
	'l:s' => \$markerfile,
	'c:s' => \$crosstype,
	'map:s' => \$lgmapfile,
	'order:s' => \$order,
	'drop_markers:s' => \$markers2drop,
	'drop_genotypes:s' => \$genotypes2drop,
	'o:s' => \$outfilebase,
	'v' => \$verbose,
);



unless (defined $binaryfile){
	die "\n\n-b option is required\n\nusage: $usage\n\n";
}
unless (defined $markerfile) {
	die "\n\nplease define linkage group / marker file with -l\n\n$usage\n\n";
}
unless ($crosstype eq "BC" or $crosstype eq "4way") {
	die "invalid cross type $crosstype, please specify either \"BC\" or \"4way\"\b\b: $usage\n\n";
}
if (scalar(@ARGV)>0) {
	die "\nunknown options: @ARGV\n\nusage: $usage\n\n";
}
unless (defined $outfilebase) {
	$outfilebase = $binaryfile.".rqtl";
}
open LOG, ">$outfilebase.log";
print LOG "$0\n\noptions in effect:\n-b $binaryfile\n-l $markerfile\n-c $crosstype\n-o $outfilebase\n";
if (defined $lgmapfile) {
	print LOG "-map $lgmapfile\n";
}
if (defined $order) {
	print LOG "-order $order\n";
}
if (defined $markers2drop) {
	print LOG "-drop_markers $markers2drop\n";
}
if (defined $genotypes2drop) {
	print LOG "-drop_genotypes $genotypes2drop\n";
}
if (defined $verbose) {
	print LOG "-v\n";
}

my $time = localtime();
print LOG "\nscript executed at localtime: $time\n\n";

print LOG "using the following input files:\n";
my $binarymd5 = ` md5sum $binaryfile `;
print LOG $binarymd5;
my $markermd5 = ` md5sum $markerfile `;
print LOG $markermd5;

if (defined $lgmapfile) {
	my $lgmapmd5 = ` md5sum $lgmapfile `;
	print LOG $lgmapmd5;
}
if (defined $order) {
	my $ordermd5 = ` md5sum $order `;
	print LOG $ordermd5;
}

if (defined $markers2drop) {
	my $markerdropmd5 = ` md5sum $markers2drop `;
	print LOG $markerdropmd5;
}
if (defined $genotypes2drop) {
	my $genodropmd5 = ` md5sum $genotypes2drop `;
	print LOG $genodropmd5;
}

print LOG "\n\n";

my %lgmap;
if (defined $lgmapfile) {
	open LGMAP, $lgmapfile || die "could not open $lgmapfile\n";
	while (<LGMAP>) {
		chomp;
		if ($_ =~ /^\#/) {next;}
		my @line = split("\t", $_);
		my ($lg, $phase, $label, $parent, $grandparent) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
		$lgmap{$lg}{$phase}{parent} = $parent;
		$lgmap{$lg}{$phase}{grandparent} = $grandparent;
		$lgmap{$lg}{$phase}{label} = $label;
	}
	close LGMAP;
}

my %markerconv;
my %nloci;

my %markers2drop;
if (defined $markers2drop) {
	open MARKERDROP, $markers2drop || die "could not open $markers2drop\n";
	while (<MARKERDROP>) {
		chomp;
		$markers2drop{$_} = 1;
	}
	close MARKERDROP;
}

my %genotypes2drop;
if (defined $genotypes2drop) {
	open GENOTYPEDROP, $genotypes2drop || die "could not open $genotypes2drop\n";
	while (<GENOTYPEDROP>) {
		chomp;
		my ($marker, $id) = split ('\s+', $_);
		$genotypes2drop{$marker}{$id} = 1;
	}
	close GENOTYPEDROP;
}


open MARKERS, $markerfile || die "could not open markerfile $markerfile\n";
my $header = <MARKERS>;

while (<MARKERS>) {
	chomp; 
	my @line = split("\t", $_);
	my ($id, $chr, $pos, $lg, $phase) = ($line[0], $line[1], $line[2], $line[3], $line[4]);
	if ($id eq "ID") { next; }
	if ($lg =~ m/-0/) { next; }
	if (defined $lgmapfile) {
		my $parent = $lgmap{$lg}{$phase}{parent}; 
		my $grandparent = $lgmap{$lg}{$phase}{grandparent};
		my $label = $lgmap{$lg}{$phase}{label};
		$markerconv{$id}{parent} = $parent;
		$markerconv{$id}{grandparent} = $grandparent;
		$markerconv{$id}{label} = $label;
		if (exists $nloci{$label}) {
			$nloci{$label}++;
		}
		else {
			$nloci{$label} = 1
		}
		if ($parent ne 1 && $parent ne 2) {
			die "\n\nparent $parent not a valid option for id$id lg$lg phase$phase\n\n";
		}
			if ($grandparent ne 1 && $grandparent ne 2) {
			die "\n\ngrandparent $grandparent not a valid option for id$id lg$lg phase$phase\n\n";
		}

		else {
			if ($parent eq 1 && $grandparent eq 1) {
				if ($crosstype eq "4way") {
					$markerconv{$id}{0} = 6;
					$markerconv{$id}{1} = 5;
				}
				elsif ($crosstype eq "BC") {
					$markerconv{$id}{0} = "B";
					$markerconv{$id}{1} = "A";
				}
				$markerconv{$id}{LG} = $lgmap{$lg}{$phase}{label};
			}
			elsif ($parent eq 1 && $grandparent eq 2) {
				if ($crosstype eq "4way") {
					$markerconv{$id}{0} = 5;
					$markerconv{$id}{1} = 6;
				}
				elsif ($crosstype eq "BC") {
					$markerconv{$id}{0} = "A";
					$markerconv{$id}{1} = "B";
				}
				$markerconv{$id}{LG} = $lgmap{$lg}{$phase}{label};
			}
			elsif ($parent eq 2 && $grandparent eq 1) {
				if ($crosstype eq "4way") {
					$markerconv{$id}{0} = 8;
					$markerconv{$id}{1} = 7;
				}
				elsif ($crosstype eq "BC") {
					$markerconv{$id}{0} = "B";
					$markerconv{$id}{1} = "A";
				}
				$markerconv{$id}{LG} = $lgmap{$lg}{$phase}{label};
			}
			elsif ($parent eq 2 && $grandparent eq 2) {
				if ($crosstype eq "4way") {
					$markerconv{$id}{0} = 7;
					$markerconv{$id}{1} = 8;
				}
				elsif ($crosstype eq "BC") {
					$markerconv{$id}{0} = "A";
					$markerconv{$id}{1} = "B";
				}
				$markerconv{$id}{LG} = $lgmap{$lg}{$phase}{label};
			}
			else {
				$markerconv{$id}{0} = "-"; 
				$markerconv{$id}{1} = "-";
				$markerconv{$id}{LG} = "NA";
			}

		}
	}
	else {
		if (exists $nloci{$lg}) {
			$nloci{$lg}++;
		}
		else {
			$nloci{$lg} = 1
		}
		if ($phase eq 1) {
			if ($crosstype eq "4way") {
				$markerconv{$id}{0} = 1;
				$markerconv{$id}{1} = 2;
			}
			elsif ($crosstype eq "BC") {
				$markerconv{$id}{0} = "A";
				$markerconv{$id}{1} = "B";
			}
			$markerconv{$id}{LG} = $lg;
		}
		elsif($phase eq 2) {
			if ($crosstype eq "4way") {
				$markerconv{$id}{0} = 2;
				$markerconv{$id}{1} = 1;
			}
			elsif ($crosstype eq "BC") {
				$markerconv{$id}{0} = "B";
				$markerconv{$id}{1} = "A";
			}
			$markerconv{$id}{LG} = $lg;
		}
		else {
			$markerconv{$id}{0} = "-"; 
			$markerconv{$id}{1} = "-";
			$markerconv{$id}{LG} = $lg;
		}
	}
}


if (defined $order) {
	open ORDER, $order || die "could not open $order\n";
	while (<ORDER>) {
		chomp;
		if ($_ =~ /^\#/) { next;}
		my @line = split("\t", $_);
		my ($id, $lg, $pos) = ($line[0], $line[1], $line[2]);
		if (exists $markerconv{$id}) {
			$markerconv{$id}{orderLG} = $lg;
			$markerconv{$id}{orderpos} = $pos;
		}
	}
	close ORDER;
}


open BINARY, $binaryfile || die "could not open binary file $binaryfile\n";
my @prefix = split("/", $binaryfile);
my $prefix = $prefix[-1];
chomp $prefix;
open OUTGENO, ">$outfilebase.geno";
open OUTPHENO, ">$outfilebase.pheno";

my $head = <BINARY>;
chomp $head;
my @header = split("\t", $head);
my @inds = @header[5..scalar(@header)-1];

# create a lookup table for filtering genotypes
my %lookup;
for (my $i=0;$i<scalar(@inds) ;$i++) {
	$lookup{$inds[$i]}=$i;
}


print OUTGENO "id,,,";
print OUTPHENO "id,";

for (my $i=0; $i<scalar(@inds)-1; $i++) {
	print OUTGENO "$inds[$i],";
	print OUTPHENO "$inds[$i],";
}
print OUTGENO "$inds[-1]\n";
print OUTPHENO "$inds[-1]\n";

while (<BINARY>) {
	chomp;
	my @line = split("\t", $_);
	my ($chr2, $pos, $id, $major, $minor) = ($line[0], $line[1], $line[2], $line[3], $line[4]); ##fill me in
	unless (exists $markerconv{$id}) {
		if (defined $verbose) {
			print LOG "marker $id has no information\n";
		}
		next;
	}
	my @genos = @line[5..scalar(@line)-1];
	
	my $lg;
	if (defined $order) {
		if (exists $markerconv{$id}{orderLG}) {
			$lg = $markerconv{$id}{orderLG};
			$pos = $markerconv{$id}{orderpos};
		}
		else {
			if (defined $verbose) {
				print LOG "marker $id has no order information\n";
			}
			next;
		}
	}
	else { 
		$lg = $markerconv{$id}{LG};
	}
	#check if it is filtered
	if (exists $markers2drop{$id}) {
		if (defined $verbose) {
			print LOG "dropping marker $id\n";
		}
		next;
	}
	print OUTGENO "$id,$lg,$pos,";
	#look up the number for that marker
	for (my $i=0; $i<scalar(@genos)-1; $i++) {
		my $ind = $inds[$i];
		if ($genos[$i] eq 1 || $genos[$i] eq 0 && ! defined $genotypes2drop{$id}{$ind}) {
			print OUTGENO "$markerconv{$id}{$genos[$i]},";
		}
			else { 
			print OUTGENO "-,"; 
		}
	}
	my $ind = $inds[-1];
	if ($genos[-1] eq 1 || $genos[-1] eq 0 && ! defined $genotypes2drop{$id}{$ind}) {
		print OUTGENO "$markerconv{$id}{$genos[-1]}\n";
	}
	else { 
		print OUTGENO "-\n"; 
	}
}

#sort the files
system "head -n 1 $outfilebase.geno > $outfilebase.temp.geno";
system "tail -n +2 $outfilebase.geno | sort -t\",\" -k2,3 -d  >> $outfilebase.temp.geno";
system "mv $outfilebase.temp.geno $outfilebase.geno";

my $outmd5 = `md5sum $outfilebase.geno `;
print LOG "\nproduced the following output files:\n";
print LOG $outmd5;
$outmd5 = ` md5sum $outfilebase.pheno `;
print LOG $outmd5;
$time = localtime();
print LOG "\n\nfinished at localtime: $time\n\n";

close LOG;


close MARKERS;
close BINARY;
close OUTGENO;
close OUTPHENO;
