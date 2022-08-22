#!/usr/bin/perl
use strict;
use warnings;
use Statistics::R;
use Getopt::Long;
use Scalar::Util qw(reftype);


######################################################################################
#
# synteny_phase_part1.pl separate minor allele contribution by parent 
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


my $usage = "perl synteny_phase_part1.pl [options]\n\n";
$usage = $usage."\t-b2 binary file of PtMs for progeny (required)\n\t\toutput from synteny_chr_assignment.pl\n";
$usage = $usage."\t-p parent binary file of PtMs(optional),\n\t\tif available is used for validation; generated from get_pseudotestcross_markers.pl\n";
$usage = $usage."\t-chr chromosome info file (required),\n\t\ta two column tab delimited file with the columns chromosomeID and length.\n\t\tOnly include chromosomes that you want to create linkage groups for\n";
$usage = $usage."\t-od output directory (required)\n";
$usage = $usage."\t-t number of threads to use (defaults to 1)\n";
$usage = $usage."\t-v verbose print figures for: \n\t\thierarchical clustering / cutting of dendrograms\n\t\tthe topological overlap matrix for each chromosome\n\t\tindividual outlier detection plots for each chromosome\n\t\t(all can be used for manual curation)\n";
$usage = $usage."\t-static_height (optional, defaults to 0.9)\n";
$usage = $usage."\t-static_minsize (optional, deafults to 30)\n";
$usage = $usage."\t-hybrid_height (optional, defaults to 0.998)\n";
$usage = $usage."\t-hybrid_deepsplit (optional, defaults to 2)\n";
$usage = $usage."\t-hybrid_pamrespectsdendro (optional, defaults to FALSE)\n";

my ($binaryfile, $parentfile, $chrfile, $outdir, $nthreads, $verbose, $height, $minsize, $hybridcut, $hybridsplit, $hybridpam);

if (scalar(@ARGV==0)) {
	die "\n\n$usage\n\n";
}

GetOptions(
	'b2:s' => \$binaryfile,
	'p:s' => \$parentfile,
	'chr:s' => \$chrfile,
	'od:s' => \$outdir,
	't:i' => \$nthreads,
	'v' => \$verbose,
	'static_height:f' => \$height,
	'static_minsize:i' => \$minsize,
	'hybrid_height:f' => \$hybridcut,
	'hybrid_deepsplit:i' => \$hybridsplit,
	'hybrid_pamrespectsdendro' => \$hybridpam,
);

if (scalar(@ARGV)>0) {
	die "\n\nunknown options: @ARGV\n\n$usage\n\n";
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
	$minsize = 30;
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
if ($height > 1) {
	die "\n\nheight must be less than or equal to 1 (not $height)\n\n";
}

unless (-e $binaryfile) {
	die "\n\nbinary file $binaryfile does not exist\n\n";
}

if (-d $outdir) {
	die "\n\noutput directory $outdir already exists\n\n";
}
mkdir $outdir; 

my ($logfile) = ("$outdir/log");
open LOG, ">$logfile";

print LOG "$0\n\noptions in effect:\n-b2 $binaryfile\n-i $chrfile\n-od $outdir\n-t $nthreads\n-static_height $height\n-static_minsize $minsize\n-hybrid_height $hybridcut\n-hybrid_depsplit $hybridsplit\n-hybrid_pamrespectsdendro $hybridpam\n";
if (defined $parentfile) {
	unless (-e $parentfile) {
		die "\n\nparentfile $parentfile does not exist\n\n";
	}
	print LOG "-p $parentfile\n";
}
if ($verbose) {
	print LOG "-v\n";
}
my $time = localtime();
print LOG "\nscript executed at localtime: $time\n\n";

print LOG "using the following input files:\n";
my $binarymd5 = ` md5sum $binaryfile `;
print LOG $binarymd5;
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
close CHR;

my $counter = 0;

open MARKERS, ">$outdir/markers_by_parent.txt";
print MARKERS "ID\tCHR\tPOS\tSTATIC_COL\tSTATIC_NUM\tDYNAMIC_COL\tDYNAMIC_NUM\tHYBRID_COL\tHYBRID_NUM\n";
close MARKERS;

my $status = Statistics::R->new();
$status->run(qq`library(WGCNA)`);
$status->run(qq`library(cluster)`);
$status->run(qq`options(stringsAsFactors = FALSE)`);
if ($nthreads > 1) {
	$status->run(qq`allowWGCNAThreads($nthreads)`);
}
$status->run(qq`read.table(\"$binaryfile\",header=T)->alldata`);

foreach my $workingchr (sort {$a <=> $b} keys %chrlength) {
	
	print LOG "\nworking on chromosome $workingchr\n";
	#separate the parental contribution based on the progeny
	$status->run(qq`data <- subset(alldata, alldata\$CHROM==$workingchr)`);
	$status->run(qq`datGt = data.frame(t(data[,6:dim(data)[2]]))`); #transpose it and get rid of the first 5 columns (information columns)
	$status->run(qq`colnames(datGt)<-data[,3] `);
	
	#clean up the dataset
	$status->run(qq`gsg=goodSamplesGenes(datGt, verbose=3)`);
	$status->run(qq`gsg\$allOK`);
	$status->run(qq`datGtclean = datGt`);
	$status->run(qq`if (!gsg\$allOK) \n { datGtclean = datGt[gsg\$goodSamples, gsg\$goodGenes] \n }`);

	my $nind = $status->get( 'dim(data)[1]');
	my $nmarkers = $status->get( 'dim(data)[2]');
	my $nindclean = $status->get( 'dim(datGtclean)[1]' );
	my $nmarkclean = $status->get( 'dim(datGtclean)[2]' );
	
	#report the cleaning results 
	my $filt_markers = $status->get( 'colnames(datGt)[which(gsg$goodGenes!="TRUE")]' );
	my $filt_inds = $status->get( 'rownames(datGt)[which(gsg$goodSamples!="TRUE")] ' );
	print LOG "\nfiltered the following individuals due to too much missing data:\n";
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

	print LOG "\nfiltered the following markers due to too much missing data:\n";
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


	
	#check for outliers
	if ($verbose) {
		$status->run(qq`indTree <- flashClust(dist(datGtclean),method=\"average\")`);
		$status->run(qq`pdf(file=\"$outdir/$workingchr.indOutlier.pdf\", title=\"$workingchr.outlier\", height=8.5, width=11)`);
		$status->run(qq`plot(indTree, main=\"$workingchr individual clustering to detect outlier\", cex=0.6, mar=c(0,4,2,0))`);
		$status->run(qq`dev.off()`);

	}
	
	#create the dendrogram and cut it up a few different ways;
	$status->run(qq`datGtclean = sapply(datGtclean, as.numeric)`);
	$status->run(qq`ADJ = abs(cor(datGtclean, use=\"p\"))^2`);
	$status->run(qq`dissTOM=TOMdist(ADJ)`);
	$status->run(qq`collectGarbage()`);
	$status->run(qq`hierTOM = hclust(as.dist(dissTOM), method= \"average\")`);
	$status->run(qq`hierTOM\$height <- round(hierTOM\$height,6)`); #identical height values were giving a non-sorted error with the cutree function; round the height values to 6 decimal places and it seems to work. Thanks google.
	$status->run(qq`colorStaticTOM=as.character(cutreeStaticColor(hierTOM, cutHeight=$height, minSize=$minsize))`);
	$status->run(qq`branch.number=cutreeDynamic(hierTOM, method = \"tree\")`);
	$status->run(qq`colorDynamicTOM=labels2colors(branch.number)`);
	$status->run(qq`colorDynamicHybridTOM=labels2colors(cutreeDynamic(hierTOM,distM=dissTOM, cutHeight = $hybridcut, deepSplit=$hybridsplit, pamRespectsDendro = $hybridpam))`);
	if ($verbose) {
		#create the hierTO PDF for all cuts
		$status->run(qq`pdf(file=\"$outdir/$workingchr.hierTOM.pdf\", title=\"$workingchr.hierTOM\", height=8.5, width=11)`);
		$status->run(qq`plotDendroAndColors(dendro=hierTOM, colors=data.frame(colorStaticTOM, colorDynamicTOM, colorDynamicHybridTOM), dendroLabels=FALSE, marAll = c(1, 8, 3, 1), main = \"$workingchr Marker dendrogram and module colors, TOM dissimilarity\")`);
		$status->run(qq`dev.off()`);
		#create the TOMplot for the static cut
		$status->run(qq`plotTOM = dissTOM^3`);
		$status->run(qq`diag(plotTOM) = NA`);
		$status->run(qq`pdf(file=\"$outdir/$workingchr.TOMplot.pdf\", title=\"$workingchr.TOMplot\", height=11, width=8.5)`);
		$status->run(qq`TOMplot(plotTOM, hierTOM, colorStaticTOM, main = \"$workingchr Network heatmap plot, static cut 0.9\")`);
		$status->run(qq`dev.off()`);
	}
	
	$status->run(qq`colorOrder = c("grey", standardColors(50))`);
	$status->run(qq`markers = data.frame(ID=colnames(datGtclean), CHR=rep($workingchr,length(colorStaticTOM)), POS=as.numeric(data\$POS[data\$ID %in% colnames(datGtclean)]), STATIC_COL=as.character(colorStaticTOM), STATIC_NUM = match(colorStaticTOM, colorOrder)-1, DYNAMIC_COL = as.character(colorDynamicTOM), DYNAMIC_NUM = match(colorDynamicTOM, colorOrder)-1, HYBRID_COL = as.character(colorDynamicHybridTOM), HYBRID_NUM = match(colorDynamicHybridTOM, colorOrder)-1)`);
		
	my $nmodules = $status->get( 'max(as.numeric(markers$STATIC_NUM))' );

	my $unassigned = 0;
	my $test = $status->get("dim(subset(markers,markers\$STATIC_NUM==0))[1]");
	if ($test > 0 ) {
		$unassigned = $test;
	}
	my $assigned = $nmarkclean- $unassigned;

	print LOG "\n$assigned markers (of $nmarkclean total) on chromosome $workingchr separated into $nmodules clusters\n";
	print LOG "$unassigned markers are unassigned\n";
	for (my $i=1;$i<=$nmodules ;$i++) {
		my $n = $status->get("dim(subset(markers,markers\$STATIC_NUM==$i))[1]");
		print LOG "$n markers are assigned to cluster $i\n";
	}

	$status->run(qq`write.table(markers,\"$outdir/markers_by_cluster.txt\",sep=\"\t\", quote=FALSE, append=T, col.names=FALSE, row.names=FALSE)`);
	
	#if the parents are sequenced, plot against parents to see how we are doing
	if (-e $parentfile) {
		$status->run(qq`alldata.parents <- read.table(\"$parentfile\", header=T)`);
		$status->run(qq`data.parents <- subset(alldata.parents, alldata.parents\$CHROM==$workingchr)`);
		$status->run(qq`datGt.parents <- data.frame(t(data.parents[,6:dim(data.parents)[2]]))`);
		$status->run(qq`colnames(datGt.parents) <- data.parents[,3]`);
		$status->run(qq`markers.overlap <- markers[markers\$ID %in% colnames(datGt.parents),]`);
		$status->run(qq`GENOvPARENT <- data.frame(CHR = rep($workingchr,times=length(sort(unique(markers.overlap\$STATIC_NUM)))), STATIC_NUM=sort(unique(markers.overlap\$STATIC_NUM)))`);
		$status->run(qq`for (i in 6:dim(data.parents)[2]) {
			markers.overlap[colnames(data.parents[i])] <- factor(data.parents[,i][data.parents\$ID %in% markers.overlap\$ID],levels=c('0','1'))
			GENOvPARENT[paste(colnames(data.parents[i]),".hom",sep="")] <- table(markers.overlap\$STATIC_NUM, markers.overlap[colnames(data.parents[i])][,1])[,1]
			GENOvPARENT[paste(colnames(data.parents[i]),".het",sep="")] <- table(markers.overlap\$STATIC_NUM, markers.overlap[colnames(data.parents[i])][,1])[,2]
		}`);
		if ($counter ==0) {
			$status->run(qq`write.table(GENOvPARENT,\"$outdir/GENOvPARENT.txt\",sep="\t", quote=FALSE, append=F, col.names=TRUE, row.names=FALSE)`);
		}
		else {
			$status->run(qq`write.table(GENOvPARENT,\"$outdir/GENOvPARENT.txt\",sep="\t", quote=FALSE, append=T, col.names=FALSE, row.names=FALSE)`);
		}

		$counter++;
	}

} 	

$status->stop();	
