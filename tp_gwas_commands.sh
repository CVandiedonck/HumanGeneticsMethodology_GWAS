# Merge cases and controls files
plink --bfile input/ic_cases --make-bed --out output/ic_cases
plink --bfile input/ic_controls --make-bed --out output/ic_controls
plink --bfile output/ic_cases --exclude input/toexclude.txt --make-bed --out output/ic_cases_autosomes
plink --bfile output/ic_controls --bmerge output/ic_cases_autosomes.bed output/ic_cases_autosomes.bim output/ic_cases_autosomes.fam --make-bed --out output/t1dcc
plink --bfile output/ic_controls --exclude input/rs3754055.txt --make-bed --out input/ic_controls2
plink --bfile output/ic_cases_autosomes --exclude input/rs3754055.txt --make-bed --out input/ic_cases_autosomes2
plink --bfile output/ic_controls --flip output/t1dcc-merge.missnp --make-bed --out output/ic_controls_flipped
plink --bfile output/ic_controls_flipped --bmerge input/ic_cases_autosomes2.bed input/ic_cases_autosomes2.bim input/ic_cases_autosomes2.fam --make-bed --out output/t1dcc

# clean SNPs
plink --bfile output/t1dcc --geno 0.05 --make-bed --out output/filter1
plink --bfile output/filter1 --test-missing --pfilter 0.00001 --out output/missingcc
plink --bfile output/filter1 --exclude output/missingcc.missing --make-bed --out output/filter2
plink --bfile output/filter2 --hwe 0.0001 --make-bed --out output/filter3
plink --bfile output/filter3 --maf 0.01 --make-bed --out output/filter4

# clean subjects
plink --bfile output/filter4 --mind 0.05 --make-bed --out output/filter5
# look for relatedness with KING
king -b output/filter5.bed --kinship --related --prefix output/filter5
king -b output/filter5.bed --kinship --unrelated --prefix output/filter5

plink --bfile output/filter5 --keep output/filter5unrelated.txt --make-bed --out output/filter6
# ou NEW : enlevez les individus apparentes avec:
plink --bfile output/filter5 --remove output/filter5unrelated_toberemoved.txt --make-bed --out output/filter6

# First association test:
plink --bfile output/filter6 --assoc --adjust --out output/allelictest

# Start R
if (!requireNamespace("qqman", quietly = TRUE)){
    install.packages("qqman")}
library(qqman)
sessionInfo()    
gwasoutputults <- read.table("output/allelictest.assoc",header=T)
head(gwasResults)
png("output/manhattan.png")
manhattan(gwasResults)
dev.off()
png("output/qqplot.png")
qq(gwasResults$P, main="Q-Qplot of P-values")
dev.off()
quit('no')

# perform MDS analysis
plink --bfile output/filter6 --exclude input/mhc_range.txt --range --make-bed --out output/withoutmhc
plink --bfile output/withoutmhc --indep 50 5 2 --out output/indep
plink --bfile output/withoutmhc --extract output/indep.prune.in --make-bed --out output/withoutld
king -b output/withoutld.bed --mds --ibs --prefix output/stratifwithoutld

# Start R
mds <- read.table("output/stratifwithoutldpc.txt", header=F)
png("output/mds.png")
par(mfrow=c(2,2))
plot(mds$V7, mds$V8, type="n", main="MDS vectors 1 and 2 in Controls")
points(mds$V7[1:357], mds$V8[1:357], col="black")
plot(mds$V7, mds$V8, type="n",main="MDS vectors 1 and 2 in Cases" )
points(mds$V7[358:577], mds$V8[358:577], col="red")
plot(mds$V8, mds$V9, type="n", main="MDS vectors 2 and 3 in Controls")
points(mds$V8[1:357], mds$V9[1:357], col="black")
plot(mds$V8, mds$V9, type="n", main="MDS vectors 2 and 3 in Cases")
points(mds$V8[358:577], mds$V9[358:577], col="red")
dev.off()
quit("no")

# 2nd association test, after MDS:
plink --bfile output/filter6 --logistic hide-covar --adjust --covar output/stratifwithoutldpc.txt --covar-number 7-11 --ci 0.95 --out output/logistic_mds_PC5

# Start R
gwasResultsAfterMDS <- read.table("output/logistic_mds_PC5.assoc.logistic", header=T)
head(gwasResultsAfterMDS)
library(qqman)
png("output/manhattan_mds.png")
manhattan(gwasResultsAfterMDS[! is.na(gwasResultsAfterMDS$P),])
dev.off()
png("output/qqplot_mds.png")
qq(gwasResultsAfterMDS$P, main="Q-Qplot of P-values")
dev.off()
quit('no')

# trios analysis
plink --bfile input/ic_trio --tdt --adjust --out output/tdttrio

# Start R
library(qqman)
tdt <- read.table("output/tdttrio.tdt", header=T)
head(tdt)
png("output/manhattan_trios.png")
manhattan(tdt[! is.na(tdt$P),], p="P")
dev.off()
png("output/qqplot_trios.png")
qq(tdt$P)
dev.off()
quit('no')
