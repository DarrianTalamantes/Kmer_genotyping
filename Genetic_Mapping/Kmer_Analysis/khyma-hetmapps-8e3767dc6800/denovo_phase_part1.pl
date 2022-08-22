#!/usr/bin/perl
use strict;
use warnings;
use Statistics::R;
use Getopt::Long;
use Scalar::Util qw(reftype);


######################################################################################
#
# denovo_phase_part1.pl separate minor allele contribution by parent 
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


my $usage = "perl denovo_phase_part1.pl [options]\n\n";
$usage = $usage."\t-b binary file of PtMs for progeny (required)\n\t\toutput from get_pseudotestcross_markers.pl\n";
$usage = $usage."\t-p parent binary file of PtMs(optional),\n\t\tif available is used for validation; generated from get_pseudotestcross_markers.pl\n";
$usage = $usage."\t-od output directory (required)\n";
$usage = $usage."\t-t number of threads to use (defaults to 1)\n";
$usage = $usage."\t-v verbose print figures for: \n\t\thierarchical clustering / cutting of dendrograms\n\t\tindividual outlier detection plot\n\t\t(all can be used for manual curation)\n";
$usage = $usage."\t-diff (optional, defaults to 2) the minimum value of the ratio of mean correlation^2 of a marker\n\t\t to other markers on the same assigned linkage group\n\t\t to the mean correlation^2 of markers on an alternate linkage group\n";
$usage = $usage."\t-height_min dendrogram cut height minimum to try (optional, defaults to 0.8)\n";
$usage = $usage."\t-height_max dendrogram cut height maximum to try (optional, defaults to 0.9)\n";
$usage = $usage."\t-height_step dendrogram cut height step size between min and max (optional, defaults to 0.0125)\n";
$usage = $usage."\t-minsize_min dendrogram cutoff LG size minimum to try (optional, deafults to 50)\n";
$usage = $usage."\t-minsize_max dendrogram cutoff LG size maximum to try (optional, defaults to 300)\n";
$usage = $usage."\t-minsize_step dendrogram cutoff LG step size between min and max (optional, defaults to 50)\n";

my ($binaryfile, $parentfile, $outdir, $verbose);
my ($height_min, $height_max, $height_step, $size_min, $size_max, $size_step, $nthreads, $diff) = (0.8, 0.9, 0.0125, 50, 300, 50, 1, 2);

if (scalar(@ARGV==0)) {
	die "\n\n$usage\n\n";
}

GetOptions(
	'b:s' => \$binaryfile,
	'p:s' => \$parentfile,
	'od:s' => \$outdir,
	't:i' => \$nthreads,
	'v' => \$verbose,
	'height_min:f' => \$height_min,
	'height_max:f' => \$height_max,
	'height_step:i' => \$height_max,
	'minsize_min:i' => \$size_min,
	'minsize_max:i' => \$size_max,
	'minsize_step:i' => \$size_step,
	'diff:f' => \$diff,
);

if (scalar(@ARGV)>0) {
	die "\n\nunknown options: @ARGV\n\n$usage\n\n";
}

unless (defined $binaryfile) {
	die "\n\nmust specify input binary file with option -b\n\n$usage\n\n";
}

unless (defined $outdir) {
	die "\n\nmust specify output directory\n\n";
}
if ($height_max > 1 || $height_min < 0) {
	die "\n\nheights must be between 0 and 1\n\n";
}
if ($size_min < 0 || $size_max < 0 ) {
	die "\n\nsizes must be greater than 0\n\n";
}
unless (-e $binaryfile) {
	die "\n\nbinary file $binaryfile does not exist\n\n";
}

mkdir $outdir; 

open LOG, ">$outdir/log.txt";
open SUMMARY, ">$outdir/lg_summary.txt";
print SUMMARY "HEIGHT\tSIZE\tLGs\tTOTAL\tASSIGNED\tUNASSIGNED\tFILTERED\n";

print LOG "$0\n\noptions in effect:\n";
print LOG "-b $binaryfile\n";
print LOG "-od $outdir\n";
print LOG "-t $nthreads\n";
print LOG "-height_min $height_min\n";
print LOG "-height_max $height_max\n";
print LOG "-height_step $height_step\n";
print LOG "-size_min $size_min\n";
print LOG "-size_max $size_max\n";
print LOG "-size_step $size_step\n";
print LOG "-diff $diff\n";


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
print LOG "\nlocaltime: $time\n\n";

print LOG "using the following input files:\n";
my $binarymd5 = ` md5sum $binaryfile `;
print LOG $binarymd5;
if (defined $parentfile) {
	my $parentmd5 = ` md5sum $parentfile `;
	print LOG $parentmd5;
}
print LOG "\n\n";


# see WGCNA tutorials, especially http://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf

#get the steps for heigh and size and make the combos

my @heights;
my @sizes;
print LOG  "using the following cut heights:\n";
for (my $i = $height_max;$i>=$height_min ;$i = $i - $height_step) {
	push @heights, $i;
	print LOG "$i\n";
	if ($height_min < $i && $i - $height_step < $height_min) {
		push @heights, $height_min;
		print LOG "$height_min\n";
	}
}
print LOG "using the following min LG sizes:\n";
for (my $i = $size_max;$i>=$size_min ;$i = $i - $size_step) {
	push @sizes, $i;
	print LOG "$i\n";
	if ($size_min < $i && $i - $size_step < $size_min) {
		push @sizes, $size_min;
		print LOG "$size_min\n";
	}
}
print LOG "\n\n";
print LOG "\nreading in data, clustering and grouping. Hold on, this may take a while\n";
print STDERR "\nreading in data, clustering and grouping. Hold on, this may take a while\n";

my $status = Statistics::R->new();
$status->run(qq`library(WGCNA)`);
$status->run(qq`library(cluster)`);
$status->run(qq`options(stringsAsFactors = FALSE)`);
if ($nthreads > 1) {
	$status->run(qq`allowWGCNAThreads($nthreads)`);
}
$status->run(qq`read.table(\"$binaryfile\",header=T)->data`);

#separate the parental contribution based on the progeny
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
	$status->run(qq`pdf(file=\"$outdir/indOutlier.pdf\", title=\"$binaryfile.indoutlier\", height=8.5, width=11)`);
	$status->run(qq`plot(indTree, main=\"$binaryfile \nindividual clustering to detect outlier\", cex=0.6, mar=c(0,4,2,0))`);
	$status->run(qq`dev.off()`);

}

#create the dendrogram and cut it up a few different ways;
$status->run(qq`datGtclean = sapply(datGtclean, as.numeric)`);
$status->run(qq`ADJ = abs(cor(datGtclean, use=\"p\"))^2`);
$status->run(qq`dissTOM=TOMdist(ADJ)`);
$status->run(qq`collectGarbage()`);
$status->run(qq`hierTOM = hclust(as.dist(dissTOM), method= \"average\")`);
$status->run(qq`hierTOM\$height <- round(hierTOM\$height,6)`); #identical height values were giving a non-sorted error with the cutree function; round the height values to 6 decimal places and it seems to work. Thanks google.

foreach my $cutheight (@heights) {
	foreach my $minsize (@sizes) {
		print LOG "\n\nworking on cutheight $cutheight and minsize $minsize:\n";
		$status->run(qq`colorStaticTOM=as.character(cutreeStaticColor(hierTOM, cutHeight=$cutheight, minSize=$minsize))`);
		if ($verbose) {
			#create the hierTO PDF for all cuts
			$status->run(qq`pdf(file=\"$outdir/h$cutheight.s$minsize.hierTOM.pdf\", title=\"h$cutheight s$minsize.hierTOM\", height=8.5, width=11)`);
			$status->run(qq`plotDendroAndColors(dendro=hierTOM, colors=colorStaticTOM, dendroLabels=FALSE, marAll = c(1, 8, 3, 1), main = \"h$cutheight s$minsize\b\nMarker dendrogram and module colors, TOM dissimilarity\")`);
			$status->run(qq`dev.off()`);
		}

		$status->run(qq`colorOrder = c("grey", standardColors(100))`);
		$status->run(qq`markers = data.frame(ID=colnames(datGtclean), CHR=as.numeric(data\$CHROM[data\$ID %in% colnames(datGtclean)]), POS=as.numeric(data\$POS[data\$ID %in% colnames(datGtclean)]), LG=match(colorStaticTOM, colorOrder)-1, PHASE = rep(\"NA\",times=dim(datGtclean)[2]))`);
			
		my $nmodules = $status->get( 'max(as.numeric(markers$LG),na.rm=T)' );

		my $unassigned = 0;
		my $test = $status->get("dim(subset(markers,markers\$LG==0))[1]");
		if ($test > 0 ) {
			$unassigned = $test;
		}
		my $assigned = $nmarkclean- $unassigned;

		print LOG "\n$assigned markers (of $nmarkclean total) separated into $nmodules linkage groups\n";
		print LOG "$unassigned markers are unassigned\n";
		my $p_assigned = sprintf ("%.2f", ($assigned / $nmarkclean + 0.0000001));
		my $p_unassigned = sprintf ("%.2f", ($unassigned / $nmarkclean + 0.0000001));


		#clean up mis-assigned markers
		print LOG "\ncleaning up mis-assigned markers\n";
		my %filt;
		#loop through each marker
		for (my $z=1; $z<=$assigned;$z++) {
			#if it has already been filtered, just go to the next marker
			$status->run(qq`currentmarker <- markers[$z,]\$ID`);
			$status->run(qq`currentLG <- markers[$z,]\$LG`);
			my $id = $status->get("currentmarker");
			my $assignedLG = $status->get("currentLG");
			my $meancorr_on;
			#get mean on
			if ($assignedLG eq '0') {
				next;
			}
			if ($assignedLG eq 'NA') {
				print LOG "no assignment for marker $assignedLG\n";
			}
			else {
				$status->run(qq`workingmarkers <- subset(markers, markers\$LG == currentLG)\$ID`);
				$meancorr_on = $status->get("as.numeric(mean(ADJ[$z, which(colnames(ADJ) %in% workingmarkers)]))");
				if ($meancorr_on eq "NaN") {
					$meancorr_on = 0;
				}
			}
			#foreach other parent, get mean off

			for (my $i=1;$i<=$nmodules ;$i++) {
				if (exists $filt{$id}) {
					next;
				}
				my ($meancorr_off, $calcdiff);
				if ($i == 0 || $i == $assignedLG) {
					next;
				}
				#if it's already been filtered, don't test again
				else {
					$status->run(qq`workingmarkers <- subset(markers, markers\$LG == $i)\$ID`);
					$meancorr_off = $status->get("as.numeric(mean(ADJ[$z, which(colnames(ADJ) %in% workingmarkers)]))");
					if ($meancorr_off eq "NaN") {
						$meancorr_off = 0;
					}
				}
				$calcdiff = $meancorr_on / ($meancorr_off + 0.00001);
				if ($calcdiff < $diff) {
					push @{$filt{$id}}, "filtered $id from LG $assignedLG: mean corr^2 to LG $assignedLG is $meancorr_on, mean corr^2 to LG $i is $meancorr_off, difference is $calcdiff";
					$status->run(qq`markers[$z,]\$LG = 0`);
				}
			}
		}
		
		my $nfilt = scalar(keys %filt);
		my $pfilt = sprintf ("%.2f", ($nfilt / $nmarkclean + 0.0000001));
		$assigned = $assigned - $nfilt;
		$p_assigned = sprintf ("%.2f", ($assigned / $nmarkclean + 0.0000001));

		print LOG "\nfiltered $nfilt markers due to potential mis-assignment (mean correlation to assigned LG / mean correlation to another LG was < diff value $diff):\n";
		print LOG "$assigned markers remain\n\n";

		foreach my $f (keys %filt) {
			foreach my $k (@{$filt{$f}}) {
				print LOG "$k\n";
			}
		}

		print SUMMARY "$cutheight\t$minsize\t$nmodules\t$nmarkclean\t$assigned ($p_assigned)\t$unassigned ($p_unassigned)\t$nfilt ($pfilt)\n";
		
		unless ($nmodules eq "NA") {
			for (my $i=1;$i<=$nmodules ;$i++) {
				my $n = $status->get("dim(subset(markers,markers\$LG==$i))[1]");
				print LOG "$n markers are assigned to LG $i\n";
			}
		}


		$status->run(qq`write.table(markers,\"$outdir/h$cutheight.s$minsize.LG.txt\",sep=\"\t\", quote=FALSE, append=F, col.names=T, row.names=FALSE)`);

		#Give some information on LGs vs CHRs
		$status->run(qq`table(markers\$CHR, markers\$LG) -> CHRvLG`);
		$status->run(qq`write.table(CHRvLG, "$outdir/h$cutheight.s$minsize.CHRvLG.txt",quote=F,sep="\t", row.names=T,col.names=T)`);

		$status->run(qq`maxn <- function(n) function(x) order(x,decreasing=TRUE)[n]`);
		#get the values and positions of the first two largest
		$status->run(qq`chr1 <- apply(CHRvLG,2,maxn(1))`);
		$status->run(qq`chr2 <- apply(CHRvLG,2,maxn(2))`);
		$status->run(qq`val1 <- apply(CHRvLG,2,function(x)x[maxn(1)(x)])`);
		$status->run(qq`val2 <- apply(CHRvLG,2,function(x)x[maxn(2)(x)])`);
		$status->run(qq`CHRvLGtop <- data.frame(LG = colnames(CHRvLG), chr1 = row.names(CHRvLG)[chr1], val1, pos2 = row.names(CHRvLG)[chr2], val2, diff=val2/val1)`);
		$status->run(qq`write.table(CHRvLGtop,"$outdir/h$cutheight.s$minsize.CHRvLGtop.txt",sep="\t", quote=FALSE, append=F, col.names=TRUE, row.names=FALSE)`);


		#check parents if available
		if (-e $parentfile ) {
			$status->run(qq`data.parents <- read.table(\"$parentfile\", header=T)`);
			$status->run(qq`datGt.parents <- data.frame(t(data.parents[,6:dim(data.parents)[2]]))`);
			$status->run(qq`colnames(datGt.parents) <- data.parents[,3]`);
			$status->run(qq`markers.overlap <- markers[markers\$ID %in% colnames(datGt.parents),]`);
			$status->run(qq`GENOvPARENT <- data.frame(LG=sort(unique(markers.overlap\$LG)))`);
			$status->run(qq`for (i in 6:dim(data.parents)[2]) {
				markers.overlap[colnames(data.parents[i])] <- factor(data.parents[,i][data.parents\$ID %in% markers.overlap\$ID],levels=c('0','1'))
				GENOvPARENT[paste(colnames(data.parents[i]),".hom",sep="")] <- table(markers.overlap\$LG, markers.overlap[colnames(data.parents[i])][,1])[,1]
				GENOvPARENT[paste(colnames(data.parents[i]),".het",sep="")] <- table(markers.overlap\$LG, markers.overlap[colnames(data.parents[i])][,1])[,2]
			}`);
			$status->run(qq`write.table(GENOvPARENT,"$outdir/h$cutheight.s$minsize.GENOvPARENT.txt",sep="\t", quote=FALSE, append=F, col.names=TRUE, row.names=FALSE)`);
		}
	}
}

$status->stop();	

close LOG;
close SUMMARY;
exit;
