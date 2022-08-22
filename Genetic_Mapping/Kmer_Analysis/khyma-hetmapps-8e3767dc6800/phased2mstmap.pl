#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

######################################################################################
#
# phased2mstmap.pl take phased info and convert to mstmap format 
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

my $usage = "perl phased2mstmap.pl\n\n";
$usage = $usage."\t-x mstmap executable location (required if not in path)\n";
$usage = $usage."\t-b binaryfile required\n";
$usage = $usage."\t-od output directory\n";
$usage = $usage."\t-l linkage group / marker file (phased.txt)\n";
$usage = $usage."\t-m LGmap (mapping of LG / phase onto some meaningful names, optional)\n";
$usage = $usage."\t\tLGmap has 5 tab delimited columns: LG, PHASE, LABEL, PARENT_NUM, GRANDPARENT_NUM\n";
$usage = $usage."\t\tLG and PHASE must match the linkage group file,\n";
$usage = $usage."\t\t LABEL is the label that you would like to propogate down the pipeline\n";
$usage = $usage."\t\tif you want to join linkage groups this is the place to do it, prior to ordering\n";
$usage = $usage."\t\tyou can use the information from the GENOvPARENT and GENOvGRANDPARENT files to create the LGmap\n";
$usage = $usage."\t\tPARENT_NUM must be 1 or 2 (up to 2 parents per LG - expect 1 per LG\n"; 
$usage = $usage."\t\tGRANDPARENT_NUM must be 1 or 2 (up to 2 grandparents per LG, expect one per LG/PHASE\n";
$usage = $usage."\t-r default=1 number of times to do iterative dropping and re-ordering,\n\t\tsuspicious genotypes reported by MSTMap and duplicate markers are dropped, and -xo filter is applied\n";
$usage = $usage."\t-df distance_functon (kosambi or haldane) default=kosambi\n";
$usage = $usage."\t-cp cut_off_p_value default=2\n";
$usage = $usage."\t-md no_map_dist default=15\n";
$usage = $usage."\t-ms no_map_size default=5\n";
$usage = $usage."\t-mt missing threshold default = 0.5\n";
$usage = $usage."\t-of objective function default = ML\n", 
$usage = $usage."\t-xo markers creating double crossovers for greater than this proportion of individuals will be filtered\n\t\t(ie 0 allows no double crossovers, 1 allows all)\n\t\t(implemented when the r parameter is > 0, default = 0.2)\n";

my ($mstmapexe, $binaryfile, $outdir, $phasedir, $markerfile, $lgmapfile);
my ($distance_function, $cutoffp, $nomapdist, $nomapsize, $missingthreshold, $objectivefunction) = ("kosambi", 2, 15, 5, 0.5, "ML");
my ($maxPropDXO, $reorder) = (0.2,1);

GetOptions(
	'x:s'=> \$mstmapexe,
	'b:s' => \$binaryfile,
	'od:s' => \$outdir,
	'l:s' => \$markerfile,
	'm:s' => \$lgmapfile,
	'df:s' => \$distance_function,
	'cp:i' => \$cutoffp,
	'md:i' => \$nomapdist,
	'ms:f' => \$missingthreshold,
	'of:s' => \$objectivefunction,
	'r:i' => \$reorder,
	'xo:f' => \$maxPropDXO,
);

unless (defined $mstmapexe) {
	$mstmapexe = "mstmap";
	unless ($mstmapexe) {
		die "\n\ncould not execute mstmap, please specify the path to the executable with -x\n\n";
	}
}
unless (defined $binaryfile){
	die "\n-b option is required\n\nusage:\n\n$usage\n\n";
}
unless (defined $markerfile) {
	die "\nmust provide marker / lg file with option -l:\n\n$usage\n\n";
}
unless (defined $outdir) {
	die "\n\nmust provide output directory name with option -od:\n\n$usage\n\n";
}
if (scalar(@ARGV)>0) {
	die "\nunknown options: @ARGV\n\nusage: $usage\n\n";
}
if (-d $outdir) {
	die "output directory $outdir already exists\n";
}

if (scalar(@ARGV) > 0 ) {
	die "options @ARGV are not valid\n\n$usage\n\n";
}
mkdir $outdir;
open LOG, ">$outdir/log.txt";

print LOG "$0\n\noptions in effect:\n-x $mstmapexe\n-b $binaryfile\n-od $outdir\n-l $markerfile\n-df $distance_function\n-cp $cutoffp\n-md $nomapdist\n-ms $missingthreshold\n-of $objectivefunction\n-r $reorder\n-xo $maxPropDXO\n";

if ($lgmapfile) {
	print LOG "-m $lgmapfile\n";
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

print LOG "LG\titeration\tstarting_markers\tnomapsize_filt(current_iteration)\tkilled_missing_genos_filt(current_iteration)\tdXO_filt(for_next_iteration)\tdup_filt(for_next_iteration)\tgenotypes_filt(for_next_iteration)\n";

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
				$markerconv{$id}{0} = "A";
				$markerconv{$id}{1} = "B";
				$markerconv{$id}{LG} = $lgmap{$lg}{$phase}{label};
			}
			elsif ($parent eq 1 && $grandparent eq 2) {
				$markerconv{$id}{0} = "B";
				$markerconv{$id}{1} = "A";
				$markerconv{$id}{LG} = $lgmap{$lg}{$phase}{label};
			}
			elsif ($parent eq 2 && $grandparent eq 1) {
				$markerconv{$id}{0} = "A";
				$markerconv{$id}{1} = "B";
				$markerconv{$id}{LG} = $lgmap{$lg}{$phase}{label};
			}
			elsif ($parent eq 2 && $grandparent eq 2) {
				$markerconv{$id}{0} = "B";
				$markerconv{$id}{1} = "A";
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
			$markerconv{$id}{0} = "A";
			$markerconv{$id}{1} = "B";
			$markerconv{$id}{LG} = $lg;
		}
		elsif($phase eq 2) {
			$markerconv{$id}{0} = "B";
			$markerconv{$id}{1} = "A";
			$markerconv{$id}{LG} = $lg;
		}
		else {
			$markerconv{$id}{0} = "-"; 
			$markerconv{$id}{1} = "-";
			$markerconv{$id}{LG} = $lg;
		}
	}
}

my %markerGeno;

foreach my $working (keys %nloci) {
	if ($working =~ m/-0/) {
		next;
	}

	open BINARY, $binaryfile || die "could not open binary file $binaryfile\n";
	open OUTGENO, ">$outdir/$working.0.unordered";

	my $head = <BINARY>;
	chomp $head;
	my @header = split("\t", $head);
	my @inds = @header[5..scalar(@header)-1];

	my $nind = scalar(@inds);
	my $nloci = $nloci{$working};

	print OUTGENO "population_type DH\n";
	print OUTGENO "population_name $working\n";
	print OUTGENO "distance_function $distance_function\n";
	print OUTGENO "cut_off_p_value $cutoffp\n";
	print OUTGENO "no_map_dist $nomapdist\n";
	print OUTGENO "no_map_size $nomapsize\n";
	print OUTGENO "missing_threshold $missingthreshold\n";
	print OUTGENO "estimation_before_clustering no\n";
	print OUTGENO "detect_bad_data yes\n";
	print OUTGENO "objective_function $objectivefunction\n";
	print OUTGENO "number_of_loci $nloci\n";
	print OUTGENO "number_of_individual $nind\n";
	print OUTGENO "\n";
	
	print OUTGENO "locus_name\t";

	for (my $i=0; $i<scalar(@inds)-1; $i++) {
		print OUTGENO "$inds[$i]\t";
	}
	print OUTGENO "$inds[-1]\n";

	while (<BINARY>) {
		chomp;
		my @line = split("\t", $_);
		my ($chr2, $pos, $id, $major, $minor) = ($line[0], $line[1], $line[2], $line[3], $line[4]); ##fill me in
		unless (exists $markerconv{$id}) {
			next;
		}
		my @genos = @line[5..scalar(@line)-1];
		
		my $lg;
		if (exists $markerconv{$id}{orderLG}) {
			$lg = $markerconv{$id}{orderLG};
			$pos = $markerconv{$id}{orderpos};
		}
		else { 
			$lg = $markerconv{$id}{LG};
		}
		my $markerlg = $markerconv{$id}{LG};
		if ($markerlg eq $working) {

			print OUTGENO "$id\t";
			#look up the number for that marker
			$markerGeno{$id} = 0;
			for (my $i=0; $i<scalar(@genos)-1; $i++) {
				if ($genos[$i] eq 1 || $genos[$i] eq 0) {
					print OUTGENO "$markerconv{$id}{$genos[$i]}\t";
					$markerGeno{$id} = $markerGeno{$id} +1;
				}
					else { 
					print OUTGENO "-\t"; 
				}
			}
			if ($genos[-1] eq 1 || $genos[-1] eq 0) {
				print OUTGENO "$markerconv{$id}{$genos[-1]}\n";
				$markerGeno{$id} = $markerGeno{$id} +1;
			}
			else { 
				print OUTGENO "-\n"; 
			}
		}
	}
	close MARKERS;
	close BINARY;
	close OUTGENO;

	#now run MSTMap
	for (my $i=0;$i<=$reorder ;$i++) {
		my $status = `$mstmapexe "$outdir/$working.$i.unordered" "$outdir/$working.$i.ordered" > "$outdir/$working.$i.log"`;
		if ($i == 0) {
			$status = parsemstmap($working, "$outdir/$working.$i.unordered", "$outdir/$working.$i.ordered", "$outdir/$working.$i.log", $nomapsize, $maxPropDXO, $i, $nloci);
		}
		if ($i > 0) {
			my $status = `$mstmapexe "$outdir/$working.$i.unordered" "$outdir/$working.$i.ordered" > "$outdir/$working.$i.log"`;
			$status = parsemstmap($working, "$outdir/$working.$i.unordered", "$outdir/$working.$i.ordered", "$outdir/$working.$i.log", $nomapsize, $maxPropDXO, $i, $nloci);
		}
	}
}

#put them all together in one file for next step

for (my $i=0;$i<=$reorder ;$i++) {
	system "cat $outdir/*.$i.ordered.mstmap > $outdir/all.$i.mstmap";
	system "cat $outdir/*.$i.ordered.markers2drop > $outdir/all.$i.markers2drop";
	system "cat $outdir/*.$i.ordered.genotypes2drop > $outdir/all.$i.genotypes2drop";
}


exit;

sub parsemstmap {

	my ($working, $unordered, $mstmapfile, $mstlogfile, $nomapsize, $maxDXO, $i, $nloci) = @_;
		
	open (MAP, $mstmapfile);
	open (MSTLOG, $mstlogfile);
	open DROPMARK, ">$mstmapfile.markers2drop";
	open KEEP, ">$mstmapfile.mstmap";
	open DROPGENO, ">$mstmapfile.genotypes2drop";
	open FILT, ">$mstmapfile.filters";

	my %dropmarkers;
	my %dropgenos;
	my %markersATpos;

	my $nomapsizefilt = 0;

	while (<MAP>) {
		my $line = $_;
		chomp $line;
		if ($line =~ m/;BEGINOFGROUP/) {
			my $count = 0;
			my @output_drop = ();
			my @output_keep = ();
			until ($line =~ m/;ENDOFGROUP/) {
				$line= <MAP>;
				$count++;
				my @line = split('\s+', $line);
				my $snp = $line[0];
				my $pos = $line[1];
				unless ($line =~ m/;ENDOFGROUP/) {
					push @output_drop, $snp;
				}
				unless ($line =~ m/;ENDOFGROUP/) {
					push @output_keep, "$snp\t$working\t$pos";
				}
			}
			$count = $count -1;
			if ($count>$nomapsize) {
				foreach my $l (@output_keep) {
					print KEEP "$l\n";
					my @line = split('\s+', $l);
					my $snp = $line[0];
					my $pos = $line[2];
					$markersATpos{$pos}{$snp} = 1;
				}
			}
			else {
				$nomapsizefilt = $nomapsizefilt + $count;
				foreach my $l (@output_drop) {
					print DROPMARK "$l\n";
					my @line = split('\s+', $l);
					my $snp = $line[0];
					my $pos = $line[1];
					print FILT "$snp\tnomapsize\t$count\n";
					$dropmarkers{$snp} = 1;
				}
			}
			next;
		}
	}


	my %doublexo;
	my %suspicious;
	my $nind;

	my ($dropped_dxo, $dropped_duplicate, $dropped_genotypes, $killed, $keepmarkers) = (0,0,0,0,0);

	open (MSTLOG, $mstlogfile);
	while (<MSTLOG>) {
		my $line = $_;
		chomp $line;
		if ($line =~ m/^caution! marker:/) {
			$line =~ s/caution! marker://;
			$line =~ s/was killed due to too many missing genotype calls//;
			print DROPMARK "$line\n";
			chomp $line;
			print FILT "$line\tkilled_too_many_missing_genotypes\t1\n"; 
			$killed++;
			$dropmarkers{$line} = 1;
			next;
		}
		if ($line =~ m/number of individuals:/) {
			$nind = $line;
			$nind =~ s/number of individuals://;
			next;
		}
		elsif ($line =~ m/suspicious data detected by our algorithm/) {
			goto SUSPICIOUS;
		}
		elsif ($line =~ m/double cross overs based on the current order/) {
			goto DOUBLE;
		}
		else {
			next;
		}

		SUSPICIOUS: 
		until ($line =~ m/double cross overs based on the current order/) {
			unless ($line =~ m/double/) {
				$line=<MSTLOG>;
				chomp $line;
				unless ($line =~ m/double/ || $line eq "") {
					my ($marker,$ind) = split("\t", $line);
					if (exists $suspicious{$marker}) {
						next;
					}
					else {
						print DROPGENO "$line\n";
						print FILT "$marker\t$ind\tsuspicious_genotype\t1\n";
						$dropgenos{$marker}{$ind} = 1;
						$dropped_genotypes++;
					}
				}
			}
			if ($line =~ m/double/) {
				goto DOUBLE;
			}
		}
		
		DOUBLE: 
		while ($line = <MSTLOG>) {
			chomp $line;
			my ($marker,$ind) = split(",", $line);
			if (exists $doublexo{$marker}) {
				$doublexo{$marker} = $doublexo{$marker} + 1;
			}
			else {
				$doublexo{$marker} = 1;
			}
		}
	}

	foreach my $marker (sort {$a cmp  $b} keys %doublexo) {
		if ($doublexo{$marker}/$nind > $maxDXO) {
			print DROPMARK "$marker\n";
			print FILT "$marker\tdoublexo\t$doublexo{$marker}\n";
			$dropmarkers{$marker} = 1;
			$dropped_dxo++;
		}
	}

	foreach my $pos (sort {$a <=> $b} keys %markersATpos) {
		if (scalar(keys%{$markersATpos{$pos}}) > 1) {
			#figure out which one to drop
			my $maxgeno = 0;
			my $keep;
			foreach my $marker (keys %{$markersATpos{$pos}}) {
				if ($markerGeno{$marker} > $maxgeno) {
					$keep  = $marker;
				}
			}
			#now drop
			foreach my $marker (keys %{$markersATpos{$pos}}) {
				if ($marker eq $keep) {
				}
				else {
					$dropmarkers{$marker} = 1;
					$dropped_duplicate++;
					print DROPMARK "$marker\n";
					print FILT "$marker\tduplicate\t$pos\n";
				}
			}
		}
	}

	if ($i==0) {
		print  LOG "$working\t$i\t$nloci\t";
	}	

	close MAP;
	close MSTLOG;
	close DROPMARK;
	close KEEP;
	close DROPGENO;
	close FILT;

	if ($reorder >0) {
		my $newi = $i+1;
		open FILTERED, ">$outdir/$working.$newi.unordered";
		#create the filtered "unordered" files
		open UNORDERED, $unordered || die "could not open $unordered\n";
		for (my $i=0;$i<=9 ;$i++) {
			my $line = <UNORDERED>;
			print FILTERED $line;
		}
		my $loci = <UNORDERED>;
		chomp $loci;
		my ($tag, $n) = split('\s+', $loci);
		my $newn = $n - scalar(keys %dropmarkers);
		print FILTERED "$tag $newn\n";
		my $l = "no";
		until ($l =~ m/locus_name/) {
			$l = <UNORDERED>;
			print FILTERED $l;
			chomp $l;
		}
		my @indvs = split('\s+', $l);
		# create a lookup table for filtering genotypes
		my %lookup;
		for (my $i=1;$i<scalar(@indvs) ;$i++) {
			$lookup{$indvs[$i]}=$i;
		}
		while (my $line = <UNORDERED>) {
			chomp $line;
			my @line = split('\s+', $line);
			my $locus = $line[0];
			#check if the whole marker gets dropped
			if (exists $dropmarkers{$locus}) {
				next;
			}
			else {
				$keepmarkers = $keepmarkers + 1;
				#check if any genotypes get dropped
				if (exists $dropgenos{$locus}) {
					#figure out which individuals
					foreach my $ind (keys $dropgenos{$locus}) {
						$line[$lookup{$ind}] = "-";
					}
				}
			}
			#print the potentially doctored up line
			my $newline = join("\t", @line);
			print FILTERED "$newline\n";
		}
		close FILTERED;
		close UNORDERED;
		if ($i==0) {
			print LOG "$nomapsizefilt\t$killed\t$dropped_dxo\t$dropped_duplicate\t$dropped_genotypes\n";
			print LOG "$working\t$newi\t$newn\t";
		}
		if ($i > 0) {
			print LOG "$nomapsizefilt\t$killed\t$dropped_dxo\t$dropped_duplicate\t$dropped_genotypes\n";
		}
		#if ($i > 0 && $newi > $reorder) {
		if ($i>0 && $newi <= $reorder) {
			print LOG "$working\t$newi\t$newn\t";
		}

	}
}

