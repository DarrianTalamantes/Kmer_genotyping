
#!/usr/bin/perl
use strict;
use warnings;
use Statistics::R;
use Getopt::Long;
use Scalar::Util qw(reftype);

######################################################################################
#
# synteny_phase_part2.pl separate minor allele contribution by phase 
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


my $usage = "perl synteny_phase_part2.pl [options]\n\n";
$usage = $usage."\t-l lg marker file (required, either the output from reference_based_phase_part1.pl\n\t\tor a custom file with parental designations\n\t\tif output from reference_based_phase_parent2.pl:\n\t\t one of -static, -dynamic, or -hybrid must be specified (defaults to -static)\n\t\tif custom, option -custom must be in effect and\n\t\t the file must be tab delimited and contain the columns \"ID\" \"CHR\" \"POS\" and \"PARENT\"\n"; 
$usage = $usage."\t-b2 binary file of PtMs for progeny (required)\n\t\toutput from reference_based_chr_assignment.pl\n";
$usage = $usage."\t-p parent binary file of PtMs(optional),\n\t\tif available is used for validation; generated from get_pseudotestcross_markers.pl\n";
$usage = $usage."\t-g grandparent binary file of PtMs(optional),\n\t\tif available is used for validation; generated from get_pseudotestcross_markers.pl\n";
$usage = $usage."\t-chr chromosome info file (required),\n\t\ta two column tab delimited file with the columns chromosomeID and length\n\t\tOnly include chromosomes that you want to create linkage groups for\n";
$usage = $usage."\t-od output directory (required)\n";
$usage = $usage."\t-t number of threads to use (defaults to 1)\n";
$usage = $usage."\t-v verbose print figures for: \n\t\thierarchical clustering / cutting of dendrograms\n\t\tthe topological overlap matrix for each chromosome\n\t\tindividual outlier detection plots for each chromosome\n\t\tphased markers vs physical order\n\t\t(all can be used for manual curation)\n";
$usage = $usage."\t-static (must specify either -static, -dynamic, or -hybrid with option -l,\n\t\t indicating which cut method to use for parental assignment. defaults to --static).\n";
$usage = $usage."\t-dynamic (must specify either -static, -dynamic, or -hybrid with option -l,\n\t\t indicating which cut method to use for parental assignment. defaults to --static).\n";
$usage = $usage."\t-hybrid (must specify either -static, -dynamic, or -hybrid with option -l,\n\t\t indicating which cut method to use for parental assignment. defaults to --static).\n";
$usage = $usage."\t-custom (use when the marker file is custom, see option -l)\n";
$usage = $usage."\t-diff (optional, defaults to 2) the minimum value of the ratio of mean correlation\n\t\t of markers on the assigned linkage group to\n\t\t the mean correlation to markeres on an alternate linkage group\n";
$usage = $usage."\t-static_height (optional, defaults to 0.9)\n";
$usage = $usage."\t-static_minsize (optional, deafults to 10)\n";
$usage = $usage."\t-hybrid_height (optional, defaults to 0.998)\n";
$usage = $usage."\t-hybrid_deepsplit (optional, defaults to 2)\n";
$usage = $usage."\t-hybrid_pamrespectsdendro (optional, defaults to FALSE)\n";

my ($markerfile, $binaryfile, $grandparentfile, $parentfile, $chrfile, $outdir, $nthreads, $verbose, $static, $dynamic, $hybrid, $custom, $diff, $height, $minsize, $hybridcut, $hybridsplit, $hybridpam);
my $type;

if (scalar(@ARGV==0)) {
	die "\n\n$usage\n\n";
}

GetOptions(
	'l:s' => \$markerfile,
	'b2:s' => \$binaryfile,
	'p:s' => \$parentfile,
	'g:s' => \$grandparentfile,
	'chr:s' => \$chrfile,
	'od:s' => \$outdir,
	't:i' => \$nthreads,
	'v' => \$verbose,
	'-static' => \$static,
	'-dynamic' => \$dynamic,
	'-hybrid' => \$hybrid,
	'-custom' => \$custom,
	'-diff:f' => \$diff,
	'static_height:f' => \$height,
	'static_minsize:i' => \$minsize,
	'hybrid_height:f' => \$hybridcut,
	'hybrid_deepsplit:i' => \$hybridsplit,
	'hybrid_pamrespectsdendro' => \$hybridpam,
);

if (scalar(@ARGV)>0) {
	die "\n\nunknown options: @ARGV\n\n$usage\n\n";
}

unless (defined $markerfile) {
	die "\n\nmust specify input marker file with option -m\n\n$usage\n\n";
}
unless (defined $binaryfile) {
	die "\n\nmust specify input binary file with option -b\n\n$usage\n\n";
}
unless (defined $chrfile) {
	die "\n\nmust specify chromosome input file with option -i\n\n$usage\n\n";
}

unless (defined $outdir) {
	die "\n\nmust specify output directory with option -o\n\n$usage\n\n";
}
unless (defined $nthreads) {
	$nthreads = 1;
}
unless (defined $height) {
	$height = 0.9
}
unless (defined $minsize) {
	$minsize = 10;
}
unless (defined $hybridcut) {
	$hybridcut = 0.998;
}
unless (defined $hybridsplit) {
	$hybridsplit = 2;
}
unless (defined $hybridpam) {
	$hybridpam = "FALSE";
}
unless (defined $static || defined $dynamic || defined $hybrid || defined $custom) {
	$type = "static";
}
if ((defined $static && defined $dynamic) || (defined $static && defined $hybrid) || (defined $static && defined $custom) || (defined $dynamic && defined $hybrid) || (defined $hybrid && defined $custom)) {
	die "\n\nplease specify only one of -static, -dynamic, -hybrid, or -custom\n";
}
if (defined $static) {
	$type = "static";
}
if (defined $dynamic) {
	$type = "dynamic";
}
if (defined $hybrid) {
	$type = "hybrid";
}
if (defined $custom) {
	$type = "custom";
}
unless (defined $diff) {
	$diff = 2;
}
if ($height > 1) {
	die "\n\nheight must be less than or equal to 1 (not $height)\n\n";
}

unless (-e $markerfile) {
	die "\n\nbinary file $binaryfile does not exist\n\n";
}
unless (-e $binaryfile) {
	die "\n\nbinary file $binaryfile does not exist\n\n";
}

mkdir $outdir;
my ($logfile) = ("$outdir/log");
open LOG, ">$logfile";

print LOG "$0\n\noptions in effect:\n-l $markerfile\n-b2 $binaryfile\n-chr $chrfile\n-od $outdir\n-t $nthreads\n-diff $diff\n-static_height $height\n-static_minsize $minsize\n-hybrid_height $hybridcut\n-hybrid_depsplit $hybridsplit\n-hybrid_pamrespectsdendro $hybridpam\n-$type\n";
if (defined $parentfile) {
	unless (-e $parentfile) {
		die "\n\nparentfile $parentfile does not exist\n\n";
	}
	print LOG "-p $parentfile\n";
}
if (defined $grandparentfile) {
	unless (-e $grandparentfile) {
		die "\n\nparentfile $grandparentfile does not exist\n\n";
	}
	print LOG "-g $grandparentfile\n";
}

if ($verbose) {
	print LOG "-v\n";
}
my $time = localtime();
print LOG "\nscript executed at localtime: $time\n\n";

print LOG "using the following input files:\n";
my $binarymd5 = ` md5sum $binaryfile `;
print LOG $binarymd5;
my $markermd5 = ` md5sum $markerfile `;
print LOG $markermd5;

if (defined $parentfile) {
	my $parentmd5 = ` md5sum $parentfile `;
	print LOG $parentmd5;
}
print LOG "\n\n";


# see WGCNA tutorials, especially http://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf


print LOG "found chromosome info:\n";
my %chrlength;
open CHR, $chrfile || die "\n\ncould not open chromosome info file $chrfile\n\n";
while (<CHR>) {
	my ($chr, $length) = split('\s+', $_);
	$chrlength{$chr} = $length;
	print LOG "$chr\t$length\n";
}
print LOG "\n";

open INMARKERS, $markerfile || die "could not open marker file $markerfile\n\n";

my $header = <INMARKERS>;
if ($type eq "static" || $type eq "dynamic" || $type eq "hybrid") {
	unless ($header =~ m/ID\tCHR\tPOS\tSTATIC_COL\tSTATIC_NUM\tDYNAMIC_COL\tDYNAMIC_NUM\tHYBRID_COL\tHYBRID_NUM/) {
		die "\n\nfile format doesn't look correct for $markerfile using type $type:\n\n$header\n\n";
	}
}
if ($type eq "custom" ) {
	unless ($header =~ m/ID\tCHR\tPOS\tCLUSTER/) {
		die "\n\nfile format doesn't look correct for $markerfile using type $type:\n\n$header\n\n";
	}
}

my %parentassign;
my %id2parent;

while (my $line = <INMARKERS>) {
	chomp $line;
	my @line = split('\t', $line);
	my $parent;
	my ($id, $chr, $pos) = ($line[0], $line[1], $line[2]);
	if ($type eq "static") {
		$parent = $line[4];
	}
	if ($type eq "dynamic") {
		$parent = $line[6];
	}
	if ($type eq "hybrid") {
		$parent = $line[8];
	}
	if ($type eq "custom") {
		$parent = $line[3];
	}
	push (@{$parentassign{$chr}{$parent}{unfilt}}, $id);
	$id2parent{$id} = $parent;
}
close INMARKERS;

my $counter = 0;
my $counter2 = 0;

open OUTMARKERS, ">$outdir/phased.txt";
print OUTMARKERS "ID\tCHR\tPOS\tLG\tPHASE_NUM\n";
close OUTMARKERS;

my $status = Statistics::R->new();
$status->run(qq`library(WGCNA)`);
$status->run(qq`library(cluster)`);
$status->run(qq`options(stringsAsFactors = FALSE)`);
if ($nthreads > 1) {
	$status->run(qq`allowWGCNAThreads($nthreads)`);
}
$status->run(qq`read.table(\"$binaryfile\",header=T)->alldata`);
$status->run(qq`alldatGt = data.frame(t(alldata[,6:dim(alldata)[2]]))`); #transpose it and get rid of the first 5 columns (information columns)
$status->run(qq`colnames(alldatGt)<-alldata[,3] `);

#clean up the dataset
$status->run(qq`gsg=goodSamplesGenes(alldatGt, verbose=3)`);
$status->run(qq`gsg\$allOK`);
$status->run(qq`alldatGtclean = alldatGt`);
$status->run(qq`if (!gsg\$allOK) \n { alldatGtclean = alldatGt[gsg\$goodSamples, gsg\$goodGenes] \n }`);

my $nindall = $status->get( 'dim(alldatGt)[1]');
my $nmarkall = $status->get( 'dim(alldatGt)[2]');
my $nindclean = $status->get( 'dim(alldatGtclean)[1]' );
my $nmarkclean = $status->get( 'dim(alldatGtclean)[2]' );
my $filtind = $nindall - $nindclean;
my $filtmark = $nmarkall - $nmarkclean; 

#report the cleaning results 
my $filt_markers = $status->get( 'colnames(alldatGt)[which(gsg$goodGenes!="TRUE")]' );
my $filt_inds = $status->get( 'rownames(alldatGt)[which(gsg$goodSamples!="TRUE")] ' );
print LOG "\nfiltered the following $filtind individuals due to too much missing data.\n$nindclean individuals out of $nindall remain:\n";
if (reftype($filt_inds)) {
	if (reftype($filt_inds) eq 'ARRAY') {
		foreach my $f (@{$filt_inds}) {
			print LOG "$f\n";
		}
	}
	else {
		print LOG "error parsing individuals";
	}
}
else {
	unless ($filt_inds =~ m/character/) {
		print LOG "$filt_inds\n";
	}
}

print LOG "\nfiltered the following $filtmark markers due to too much missing data.\n$nmarkclean markers out of $nmarkall remain:\n";
if (reftype($filt_markers)) {
	if (reftype($filt_markers) eq 'ARRAY') {
		foreach my $f (@{$filt_markers}) {
			print LOG "$f\n";
		}
	}
	else {
		print LOG "error parsing markers\n";
	}
}
else {
	unless ($filt_markers =~ m/character/) {
		print LOG "$filt_markers\n";
	}
}


foreach my $workingchr (sort {$a <=> $b} keys %chrlength) {
	my @chrmarkers;

	foreach my $workingparent (sort {$a <=> $b} keys %{$parentassign{$workingchr}}) {
		unless ($workingparent == 0) {
			push (@chrmarkers, @{$parentassign{$workingchr}{$workingparent}{unfilt}});
		}
	}

	#clean up mis-assigned markers
	print LOG "\ncleaning up mis-assigned markers on chromosome $workingchr\n";
	if ($verbose) {
		print LOG "\nN\tID\tchr\tassigned_cluster\tmeancorr_on\tother_cluster\tmeancorr_other\tdifference\n";
	}
	$status->set('chrmarkers', \@chrmarkers);
	$status->run(qq`datGtchr <- alldatGtclean[which(colnames(alldatGtclean) %in% chrmarkers)]`);
	$status->run(qq`datGtchr = sapply(datGtchr, as.numeric)`);
	$status->run(qq`cor_chr <- abs(cor(datGtchr, use=\"p\"))`);
	my %filt;
	for (my $z=0; $z<scalar(@chrmarkers);$z++) {
		my $id = $chrmarkers[$z];
		$id =~ s/\"//g;
		my $assigned = $id2parent{$id};
		my $meancorr_on;
		if ($assigned == 0) {
			next;
		}
		else {
			my @workingmarkers = @{$parentassign{$workingchr}{$assigned}{unfilt}};		
			$status->set('workingmarkers', \@workingmarkers);
			$meancorr_on = $status->get("as.numeric(mean(cor_chr[\"$id\", which(colnames(cor_chr) %in% workingmarkers)],use=\"p\"))");
			if ($meancorr_on eq "NaN") {
				$meancorr_on = 0;
			}
		}
		foreach my $otherparent (sort {$a <=> $b} keys %{$parentassign{$workingchr}}) {
			my ($meancorr_off, $calcdiff);
			if ($otherparent == 0 || $otherparent == $assigned) {
				next;
			}
			else {
				my @workingmarkers = ($id, @{$parentassign{$workingchr}{$otherparent}{unfilt}});		
				$status->set('workingmarkers', \@workingmarkers);
				$meancorr_off = $status->get("as.numeric(mean(cor_chr[\"$id\", which(colnames(cor_chr) %in% workingmarkers)],use=\"p\"))");
				if ($meancorr_off eq "NaN") {
					$meancorr_off = 0;
				}
			}
			$calcdiff = $meancorr_on / ($meancorr_off + 0.00001);
			if ($verbose) {
				print LOG "$z\t$id\t$workingchr\t$assigned\t$meancorr_on\t$otherparent\t$meancorr_off\t$calcdiff\n";
			}
			if ($calcdiff < $diff) {
				push @{$filt{$id}}, "filtered $id from chromosome $workingchr cluster $assigned: mean corr to cluster $assigned is $meancorr_on, mean corr to cluster $otherparent is $meancorr_off, difference is $calcdiff";
				$id2parent{$id} = 0;
				push (@{$parentassign{$workingchr}{0}{filt}}, $id);
			}
			else {
				push (@{$parentassign{$workingchr}{$assigned}{filt}}, $id);
			}
		}
	}
	
	my $nfilt = scalar(keys %filt);
	print LOG "\nfiltered $nfilt markers from chromosome $workingchr due to potential mis-assignment (mean correlation to assigned cluster / mean correlation to another cluster was < diff value $diff):\n\n";
	foreach my $f (keys %filt) {
		foreach my $k (@{$filt{$f}}) {
			print LOG "$k\n";
		}
	}


	foreach my $workingparent (sort {$a <=> $b} keys %{$parentassign{$workingchr}}) {
		if ($workingparent == 0) {
			next;
		}
		print LOG "\nworking on chromosome $workingchr cluster $workingparent\n";

		my @workingmarkers = @{$parentassign{$workingchr}{$workingparent}{filt}};
	
		#get the genotype information
		$status->set('workingmarkers', \@workingmarkers);
		$status->run(qq`print(workingmarkers)`);
		$status->run(qq`datGtclean <- alldatGtclean[which(colnames(alldatGtclean) %in% workingmarkers)]`);

		#check for outliers
		if ($verbose) {
			$status->run(qq`indTree <- flashClust(dist(datGtclean),method=\"average\")`);
			$status->run(qq`pdf(file=\"$outdir/$workingchr.$workingparent.indOutlier.pdf\", title=\"$workingchr.$workingparent.outlier\", height=8.5, width=11)`);
			$status->run(qq`plot(indTree, main=\"chr $workingchr parent $workingparent individual clustering to detect outlier\", cex=0.6, mar=c(0,4,2,0))`);
			$status->run(qq`dev.off()`);

		}
		
		$status->run(qq`allmarkers = data.frame(ID=character(0), CHR=numeric(0), POS=numeric(0), CLUSTER_NUM=numeric(0), PHASE_NUM=numeric())`);
			
		#create the dendrogram and cut it up a few different ways;
		$status->run(qq`datGtclean = sapply(datGtclean, as.numeric)`);
		$status->run(qq`ADJ = cor(datGtclean, use = \"p\")`);
		$status->run(qq`ADJ=(ADJ+abs(ADJ))/2`);
			
		$status->run(qq`dissTOM=TOMdist(ADJ)`);
		$status->run(qq`collectGarbage()`);
		$status->run(qq`hierTOM = hclust(as.dist(dissTOM), method= \"average\")`);
		$status->run(qq`hierTOM\$height <- round(hierTOM\$height,6)`); #identical height values were giving a non-sorted error with the cutree function; round the height values to 6 decimal places and it seems to work. Thanks google.
			
		$status->run(qq`colorOrder = c("grey", standardColors(50))`);
		$status->run(qq`cut2 <- cutree(hierTOM, k=2)`);
			
		$status->run(qq`colorStaticTOM=as.character(cutreeStaticColor(hierTOM, cutHeight=$height, minSize=$minsize))`);
		$status->run(qq`branch.number=cutreeDynamic(hierTOM, method = \"tree\")`);
		$status->run(qq`colorDynamicTOM=labels2colors(branch.number)`);
		$status->run(qq`colorDynamicHybridTOM=labels2colors(cutreeDynamic(hierTOM,distM=dissTOM, cutHeight = $hybridcut, deepSplit=$hybridsplit, pamRespectsDendro = $hybridpam))`);
		if ($verbose) {
			$status->run(qq`pdf(file=\"$outdir/$workingchr.$workingparent.hierTOM.pdf\", title=\"$workingchr.$workingparent.hierTOM\", height=8.5, width=11)`);
			$status->run(qq`plotDendroAndColors(dendro=hierTOM, colors=data.frame(cut2, colorStaticTOM, colorDynamicTOM, colorDynamicHybridTOM), dendroLabels=FALSE, marAll = c(1, 8, 3, 1), main = \"chromosome $workingchr cluster $workingparent Marker dendrogram and module colors, TOM dissimilarity\")`);
			$status->run(qq`dev.off()`);
		}
		
		#create the TOMplot for the dynamic cut
		$status->run(qq`plotTOM = dissTOM^3`);
		$status->run(qq`diag(plotTOM) = NA`);
		if ($verbose) {
			$status->run(qq`pdf(file=\"$outdir/$workingchr.$workingparent.TOMplot.pdf\", title=\"$workingchr.$workingparent.TOMplot\", height=11, width=8.5)`);
			$status->run(qq`TOMplot(plotTOM, hierTOM, colorStaticTOM, main = \"chromosome $workingchr cluster $workingparent Network heatmap plot, static cut\")`);
			$status->run(qq`dev.off()`);
		}
		
		$status->run(qq`markers = data.frame(ID=colnames(datGtclean), CHR=rep($workingchr,length(colorDynamicTOM)), LG= rep(paste($workingchr,$workingparent,sep="-"),length(colorDynamicTOM)), POS=as.numeric( alldata\$POS[alldata\$ID %in% colnames(datGtclean)] ), STATIC_COL=as.character(colorStaticTOM), STATIC_NUM = match(colorStaticTOM, colorOrder)-1, DYNAMIC_COL = as.character(colorDynamicTOM), DYNAMIC_NUM = match(colorDynamicTOM, colorOrder)-1, HYBRID_COL = as.character(colorDynamicHybridTOM), HYBRID_NUM = match(colorDynamicHybridTOM, colorOrder)-1)`);
	
		my $nmodules = $status->get( 'max(as.numeric(markers$STATIC_NUM))' );
		if ($nmodules != 2) {
			die "could not resolve 2 phases for chromosome $workingchr cluster $workingparent, found $nmodules\n";
		}

		my $nmark = $status->get("dim(markers)[1]");
		my $unassigned = 0;
		my $test = $status->get("dim(subset(markers,markers\$STATIC_NUM==0))[1]");
		if ($test > 0 ) {
			$unassigned = $test;
		}

		my $assigned = $nmark - $unassigned;

		print LOG "\n$assigned markers (of $nmark total) on chromosome $workingchr cluster $workingparent separated into $nmodules phases\n";
		print LOG "$unassigned markers are unassigned\n";
		for (my $i=1;$i<=$nmodules ;$i++) {
			my $n = $status->get("dim(subset(markers,markers\$STATIC_NUM==$i))[1]");
			print LOG "$n markers are assigned to phase $i\n";
		}
		
		if (defined $parentfile) {
			$status->run(qq`data.parents <- read.table(\"$parentfile\", header=T)`);
			$status->run(qq`datGt.parents <- data.frame(t(data.parents[,6:dim(data.parents)[2]]))`);
			$status->run(qq`colnames(datGt.parents) <- data.parents[,3]`);
			$status->run(qq`markers.overlap <- markers[markers\$ID %in% colnames(datGt.parents),]`);
				$status->run(qq`GENOvPARENT <- data.frame(LG=sort(unique(markers.overlap\$LG)),PHASE=sort(unique(markers.overlap\$STATIC_NUM)))`);
				$status->run(qq`for (i in 6:dim(data.parents)[2]) {
					markers.overlap[colnames(data.parents[i])] <- factor(data.parents[,i][data.parents\$ID %in% markers.overlap\$ID],levels=c('0','1'))
					mtable <- table(paste(markers.overlap\$LG, markers.overlap\$STATIC_NUM), markers.overlap[colnames(data.parents[i])][,1])
					GENOvPARENT[paste(colnames(data.parents[i]),".hom",sep="")] <- c(mtable[1,1],mtable[2,1])
					GENOvPARENT[paste(colnames(data.parents[i]),".het",sep="")] <- c(mtable[1,2],mtable[2,2])
				}`);
				if ($counter ==0) {
					$status->run(qq`write.table(GENOvPARENT,\"$outdir/GENOvPARENT.txt\",sep="\t", quote=FALSE, append=F, col.names=TRUE, row.names=FALSE)`);
				}
				else {
					$status->run(qq`write.table(GENOvPARENT,\"$outdir/GENOvPARENT.txt\",sep="\t", quote=FALSE, append=T, col.names=FALSE, row.names=FALSE)`);
				}

				$counter++;
		}

		if (defined $grandparentfile) {
			$status->run(qq`data.parents <- read.table(\"$grandparentfile\", header=T)`);
			$status->run(qq`datGt.parents <- data.frame(t(data.parents[,6:dim(data.parents)[2]]))`);
			$status->run(qq`colnames(datGt.parents) <- data.parents[,3]`);
			$status->run(qq`markers.overlap <- markers[markers\$ID %in% colnames(datGt.parents),]`);
				$status->run(qq`GENOvPARENT <- data.frame(LG=sort(unique(markers.overlap\$LG)),PHASE=sort(unique(markers.overlap\$STATIC_NUM)))`);
				$status->run(qq`for (i in 6:dim(data.parents)[2]) {
					markers.overlap[colnames(data.parents[i])] <- factor(data.parents[,i][data.parents\$ID %in% markers.overlap\$ID],levels=c('0','1'))
					mtable <- table(paste(markers.overlap\$LG, markers.overlap\$STATIC_NUM), markers.overlap[colnames(data.parents[i])][,1])
					GENOvPARENT[paste(colnames(data.parents[i]),".hom",sep="")] <- c(mtable[1,1],mtable[2,1])
					GENOvPARENT[paste(colnames(data.parents[i]),".het",sep="")] <- c(mtable[1,2],mtable[2,2])
				}`);
				if ($counter2 ==0) {
					$status->run(qq`write.table(GENOvPARENT,\"$outdir/GENOvGRANDPARENT.txt\",sep="\t", quote=FALSE, append=F, col.names=TRUE, row.names=FALSE)`);
				}
				else {
					$status->run(qq`write.table(GENOvPARENT,\"$outdir/GENOvGRANDPARENT.txt\",sep="\t", quote=FALSE, append=T, col.names=FALSE, row.names=FALSE)`);
				}

				$counter2++;
		}


		$status->run(qq'parents <- data.frame(ID = markers\$ID, CLUSTER_NUM= rep($workingparent, length(markers\$POS)))');	
		$status->run(qq'getmarkers <- data.frame(ID = markers\$ID, CHR=markers\$CHR, POS=markers\$POS, LG=rep(paste(markers\$CHR[1],"-",$workingparent, sep=""), length(markers\$POS)), PHASE_NUM=markers\$STATIC_NUM)' );
		$status->run(qq'allmarkers <- rbind(allmarkers, getmarkers)');

		$status->run(qq`colortable = data.frame(reds = c("red", "red4", "magenta", "purple", "darkorchid"), blues = c("blue", "turquoise", "navyblue", "skyblue","slategray3"), greens = c("green", "green4", "greenyellow", "pine", "cactus"))`);
		
		$status->run(qq'subset.phased <- datGtclean');
		for (my $z = 1; $z<=$nmark; $z++) {
			$status->run(qq`mname <- colnames(datGtclean)[$z]`);
			my $parentnum = $status->get("parents[$z,]\$CLUSTER_NUM");
			my $phasenum = $status->get("getmarkers[$z,]\$PHASE_NUM");
			$status->run(qq`subset.phased[,mname] <- replace(subset.phased[,mname], subset.phased[,mname] == 0, "transparent")`);
			$status->run(qq`subset.phased[,mname] <- replace(subset.phased[,mname], subset.phased[,mname] == 1, colortable[$parentnum, $phasenum])`);
		}
		
		my $nprogeny = $status->get( 'dim(datGtclean)[1]' );
		if ($verbose) {
			$status->run(qq`pdf(file="$outdir/$workingchr.$workingparent.phased.pdf", title="$outdir/$workingchr.$workingparent.phased.pdf", height=8.5, width=11)`);
			$status->run(qq`plot(rep(1, $nmark), getmarkers\$POS, col=as.character(subset.phased[1,]), pch=20, main=paste(\"phased progeny\", \"chr $workingchr parent $workingparent\"), xlim = c(0, $nprogeny), ylim=c(0,$chrlength{$workingchr} + 1000000), cex=0.5)`);
			$status->run(qq`abline(h=$chrlength{$workingchr})`);
			for (my $i=2; $i <= $nprogeny; $i++) {
				$status->run(qq`points(rep($i, $nmark), getmarkers\$POS, col=as.character(subset.phased[$i,]),pch=20,cex=0.5)`);
			}
			$status->run(qq`dev.off()`);
		}
		#now write the marker list
		$status->run(qq`allmarkers <- allmarkers[order(allmarkers\$POS),]`);
		my $nmarkersall = $status->get('length(allmarkers[,1])');
		$status->run(qq`write.table(allmarkers, \"$outdir/phased.txt\", sep="\t", quote=FALSE, append=T, col.names=F, row.names=FALSE)`);

	}

}

	
$status->stop();

		
			
