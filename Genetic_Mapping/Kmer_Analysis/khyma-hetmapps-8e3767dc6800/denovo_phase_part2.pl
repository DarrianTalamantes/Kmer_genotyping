
#!/usr/bin/perl
use strict;
use warnings;
use Statistics::R;
use Getopt::Long;
use Scalar::Util qw(reftype);


######################################################################################
#
# denovo_phase_part2.pl separate minor allele contribution by parent 
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


my $usage = "perl denovo_phase_part2.pl [options]\n\n";
$usage = $usage."\t-b binary file of PtMs for progeny (required)\n\t\toutput from get_pseudotestcross_markers.pl\n";
$usage = $usage."\t-l LG marker file (required, output from denovo_phase_part1.pl)\n";
$usage = $usage."\t-p parent binary file of PtMs(optional),\n\t\tif available is used for validation; generated from get_pseudotestcross_markers.pl\n";
$usage = $usage."\t-g grandparent binary file ot PtMs (optional),\n\t\tif available is used for validation; generated from get_pseudotestcross_markers.pl\n";
$usage = $usage."\t-od output directory (required)\n";
$usage = $usage."\t-t number of threads to use (defaults to 1)\n";
$usage = $usage."\t-v verbose print figures for: \n\t\thierarchical clustering / cutting of dendrograms\n\t\tthe topological overlap matrix for each linkage group\n\t\tindividual outlier detection plots for each linkage group\n\t\t(all can be used for manual curation)\n";
$usage = $usage."\t-static_height (optional, defaults to 0.9)\n";
$usage = $usage."\t-static_minsize (optional, deafults to 10)\n";
$usage = $usage."\t-hybrid_height (optional, defaults to 0.998)\n";
$usage = $usage."\t-hybrid_deepsplit (optional, defaults to 2)\n";
$usage = $usage."\t-hybrid_pamrespectsdendro (optional, defaults to FALSE)\n";

my ($binaryfile, $lgfile, $parentfile, $grandparentfile, $outdir, $verbose);
my ($nthreads, $height, $minsize, $hybridcut, $hybridsplit, $hybridpam) = (1, 0.9, 10, 0.998, 2, "FALSE");

if (scalar(@ARGV==0)) {
	die "\n\n$usage\n\n";
}

GetOptions(
	'b:s' => \$binaryfile,
	'l:s' => \$lgfile,
	'p:s' => \$parentfile,
	'g:s' => \$grandparentfile,
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
unless (defined $lgfile) {
	die "\n\nmust specify input LG file with option -l\n\n$usage\n\n";
}
unless (defined $outdir) {
	die "\n\nmust specify output directory\n\n";
}
if ($height > 1 || $height < 0) {
	die "\n\nheight must be between 0 and 1\n\n";
}
if ($minsize < 0 || $minsize < 0 ) {
	die "\n\nsizes must be greater than 0\n\n";
}
if ($hybridcut <0 || $hybridcut > 1) {
	die "\n\nhybrid_height must be between 0 and 1\n\n";
}
unless (-e $binaryfile) {
	die "\n\nbinary file $binaryfile does not exist\n\n";
}

mkdir $outdir; 

open LOG, ">$outdir/log.txt";

print LOG "$0\n\noptions in effect:\n";
print LOG "-b $binaryfile\n";
print LOG "-l $lgfile\n";
print LOG "-od $outdir\n";
print LOG "-t $nthreads\n";
print LOG "-static_height $height\n";
print LOG "-static_minsize $minsize\n";
print LOG "-hybrid_height $hybridcut\n";
print LOG "-hybrid_deepsplit $hybridsplit\n";
print LOG "-hybrid_pamrespectsdendro $hybridpam\n";

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
print LOG "\nlocaltime: $time\n\n";

print LOG "using the following input files:\n";
my $binarymd5 = ` md5sum $binaryfile `;
print LOG $binarymd5;
if (defined $parentfile) {
	my $parentmd5 = ` md5sum $parentfile `;
	print LOG $parentmd5;
}
print LOG "\n\n";

open INLG, $lgfile || die "could not open LG file $lgfile\n\n";

my $header = <INLG>;
unless ($header =~ m/^ID\tCHR\tPOS\tLG/) {
	die "doesn't look look like your LG file is in the right format.\nexpecting tab delimited file with header starting ID CHR POS LG\n";
}

my %lgassign;
while (my $line = <INLG>) {
	chomp $line;
	my @line = split('\s+', $line);
	my ($id, $chr, $pos, $lg) = ($line[0], $line[1], $line[2], $line[3]);
	unless ($lg eq "NA") {
		push @{$lgassign{$lg}}, $id;
	}
}
close INLG;
if (scalar (keys %lgassign)==0 ) {
	die "\n\nno valid linkage groups found in lgfile $lgfile\n\n";
}
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

foreach my $lg (sort {$a <=> $b} keys %lgassign) {
	if ($lg == 0) {
		next;
	}
	print LOG "\nworking on linkage group $lg\n";
	
	my @workingmarkers = @{$lgassign{$lg}};

	$status->set('workingmarkers', \@workingmarkers);
	$status->run(qq`print(workingmarkers)`);
	$status->run(qq`data <- alldata[which(alldata\$ID %in% workingmarkers),]`);
	$status->run(qq`datGt = data.frame(t(data[,6:dim(data)[2]]))`); #transpose it and get rid of the first 5 columns (information columns)
	$status->run(qq`colnames(datGt)<-data[,3] `);
	
	#clean up the dataset
	$status->run(qq`gsg=goodSamplesGenes(datGt, verbose=3)`);
	$status->run(qq`gsg\$allOK`);
	$status->run(qq`datGtclean = datGt`);
	$status->run(qq`if (!gsg\$allOK) \n { datGtclean = datGt[gsg\$goodSamples, gsg\$goodGenes] \n }`);

	#my $nind = $status->get( 'dim(data)[1]');
	#my $nmarkers = $status->get( 'dim(data)[2]');
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
		$status->run(qq`pdf(file=\"$outdir/$lg.indOutlier.pdf\", title=\"$lg.outlier\", height=8.5, width=11)`);
		$status->run(qq`plot(indTree, main=\"lg $lg\b\nindividual clustering to detect outlier\", cex=0.6, mar=c(0,4,2,0))`);
		$status->run(qq`dev.off()`);

	}

	$status->run(qq`allmarkers = data.frame(ID=character(0), CHR=numeric(0), POS=numeric(0), PARENT_NUM=numeric(0), PHASE_NUM=numeric())`);
			
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
		$status->run(qq`pdf(file=\"$outdir/$lg.hierTOM.pdf\", title=\"$lg.hierTOM\", height=8.5, width=11)`);
		$status->run(qq`plotDendroAndColors(dendro=hierTOM, colors=data.frame(cut2, colorStaticTOM, colorDynamicTOM, colorDynamicHybridTOM), dendroLabels=FALSE, marAll = c(1, 8, 3, 1), main = \"linkage group $lg parent\b\nMarker dendrogram and module colors, TOM dissimilarity\")`);
		$status->run(qq`dev.off()`);
	}
	
	#create the TOMplot for the dynamic cut
	$status->run(qq`plotTOM = dissTOM^3`);
	$status->run(qq`diag(plotTOM) = NA`);
	if ($verbose) {
		$status->run(qq`pdf(file=\"$outdir/$lg.TOMplot.pdf\", title=\"$lg.TOMplot\", height=11, width=8.5)`);
		$status->run(qq`TOMplot(plotTOM, hierTOM, colorStaticTOM, main = \"linkage group $lg\b\nNetwork heatmap plot, static cut\")`);
		$status->run(qq`dev.off()`);
	}
	
	#$status->run(qq`markers = data.frame(ID=colnames(datGtclean), CHR=as.numeric(data\$CHROM[data\$ID %in% colnames(datGtclean)]), LG=rep($lg,times=length(data\$CHROM[data\$ID %in% colnames(datGtclean)]), POS=as.numeric(data\$POS[data\$ID %in% colnames(datGtclean)]), STATIC_COL=as.character(colorStaticTOM), STATIC_NUM = match(colorStaticTOM, colorOrder)-1, DYNAMIC_COL = as.character(colorDynamicTOM), DYNAMIC_NUM = match(colorDynamicTOM, colorOrder)-1, HYBRID_COL = as.character(colorDynamicHybridTOM), HYBRID_NUM = match(colorDynamicHybridTOM, colorOrder)-1)`);
	$status->run(qq`markers = data.frame(ID=colnames(datGtclean), CHR=as.numeric(data\$CHROM[data\$ID %in% colnames(datGtclean)]), LG=rep($lg, times=length(colnames(datGtclean))), POS=as.numeric(data\$POS[data\$ID %in% colnames(datGtclean)]), STATIC_COL=as.character(colorStaticTOM), STATIC_NUM = match(colorStaticTOM, colorOrder)-1, DYNAMIC_COL = as.character(colorDynamicTOM), DYNAMIC_NUM = match(colorDynamicTOM, colorOrder)-1, HYBRID_COL = as.character(colorDynamicHybridTOM), HYBRID_NUM = match(colorDynamicHybridTOM, colorOrder)-1)`);
	
	my $nmodules = $status->get( 'max(as.numeric(markers$STATIC_NUM),na.rm=T)' );
	if ($nmodules != 2) {
		die "could not resolve 2 phases for linakge group $lg, found $nmodules\n";
	}

	my $unassigned = 0;
	my $test = $status->get("dim(subset(markers,markers\$STATIC_NUM==0))[1]");
	if ($test > 0 ) {
		$unassigned = $test;
	}
	my $assigned = $nmarkclean- $unassigned;

	print LOG "\n$assigned markers (of $nmarkclean total) on linkage group $lg separated into $nmodules phases\n";
	print LOG "$unassigned markers are unassigned\n";
	unless ($nmodules eq "NA") {
		for (my $i=1;$i<=$nmodules ;$i++) {
			my $n = $status->get("dim(subset(markers,markers\$STATIC_NUM==$i))[1]");
			print LOG "$n markers are assigned to phase $i\n";
		}
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
				GENOvPARENT[paste(colnames(data.parents[i]),".0",sep="")] <- c(mtable[1,1],mtable[2,1])
				GENOvPARENT[paste(colnames(data.parents[i]),".1",sep="")] <- c(mtable[1,2],mtable[2,2])
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


	$status->run(qq'parents <- data.frame(ID = markers\$ID, LG=rep($lg, times=length( markers\$ID)))');	
	$status->run(qq'getmarkers <- data.frame(ID = markers\$ID, CHR=markers\$CHR, POS=markers\$POS, LG=rep($lg, times=length(markers\$ID)), PHASE_NUM=markers\$STATIC_NUM)' );
	$status->run(qq'allmarkers <- rbind(allmarkers, getmarkers)');

	$status->run(qq`colortable = data.frame(reds = c("red", "red4", "magenta"), blues = c("blue", "turquoise", "navyblue"), greens = c("green", "green4", "greenyellow"))`);
	
	$status->run(qq'subset.phased <- datGtclean');
	for (my $z = 1; $z<=$nmarkclean; $z++) {
		$status->run(qq`mname <- colnames(datGtclean)[$z]`);
		my $lg = $status->get("parents[$z,]\$LG");
		my $phasenum = $status->get("getmarkers[$z,]\$PHASE_NUM");
		$status->run(qq`subset.phased[,mname] <- replace(subset.phased[,mname], subset.phased[,mname] == 0, "transparent")`);
		$status->run(qq`subset.phased[,mname] <- replace(subset.phased[,mname], subset.phased[,mname] == 1, colortable[$lg, $phasenum])`);
	}
	
	my $nprogeny = $status->get( 'dim(datGtclean)[1]' );
	#now write the marker list
	$status->run(qq`allmarkers <- allmarkers[order(allmarkers\$POS),]`);
	my $nmarkersall = $status->get('length(allmarkers[,1])');
	$status->run(qq`write.table(allmarkers, \"$outdir/phased.txt\", sep="\t", quote=FALSE, append=T, col.names=F, row.names=FALSE)`);

}

$status->stop();
exit;

close LOG;

