library(docopt)
## adapted from https://github.com/oicr-gsi/wf_cfmedip/blob/master/workflow/runMedips/runMedips.r

doc <- "Get MEDIPS QC metrics.
Usage:
    QC_MEDIPS.R --bamFile <FILE> --outputDir <DIR> --windowSize <SIZE> [ --genome <GENOME> ]

Options:
    --bamFile FILE       Aligned, sorted, filtered reads (bam)
    --outputDir DIR      Path to output folder
    --windowSize SIZE    Size of genomic windows (bp, e.g. 300)
    --genome GENOME      Path to a folder containing a custom BSgenome as a package, 
                             which will be loaded using devtools::load_all();
                             or the name of BSgenome (usually BSgenome.Hsapiens.UCSC.hg38
                             or BSgenome.Athaliana.TAIR.TAIR9)
    --help               show this help text
"
opt <- docopt(doc)

library(tidyverse)
library(gtools)
library(MEDIPS)
library(BSgenome.Hsapiens.UCSC.hg38)

#if (file.exists(paste(opt[['genome']], 'DESCRIPTION', sep='/'))) {
#    devtools::load_all(opt[['genome']])
#    bsgenome <- getBSgenome(basename(opt[['genome']]))
#} else {
#    bsgenome <- getBSgenome(opt[['genome']])
#}

if (!file.exists(opt$bamFile)){
  stop(paste0("ERROR: bam file not found ", opt$bamFile), call.=FALSE)
}

## get user parameters
bam_file = opt$bamFile
ws = as.numeric(opt$windowSize)
out_dir = paste0(opt$outputDir, "/")
chr.select=paste0("chr", c(1:22,"X","Y"))
#chr.select=c("chr1","chr22")

BSgenome="BSgenome.Hsapiens.UCSC.hg38"
uniq = 0 ## WARNING: default settings normalize the data, must be set to 0 to disable this transformation
extend = 0 ## relevant for single-end: https://support.bioconductor.org/p/81098/
shift = 0
paired = TRUE

# disables the scientific notation to avoid powers in genomic coordinates (i.e. 1e+10)
options(scipen = 999)

## create MEDIPS set ###################################
message("Processing MEDIPS QC metrics: window size: ", ws)
MeDIPset =
  MEDIPS.createSet(file = bam_file,
                   BSgenome = BSgenome,
                   uniq = uniq,
                   extend = extend,
                   shift = shift,
                   paired = paired,
                   window_size = ws,
                   chr.select = chr.select)

fname <- unlist(strsplit(basename(bam_file),split="\\."))[1]
# fname

## coupling set: maps CG densities across the genome
CS = MEDIPS.couplingVector(pattern="CG", refObj=MeDIPset)

## saturation analysis #################################
## whether a given set of mapped reads is sufficient to generate a saturated and reproducible coverage profile
## calculates Pearson correlation of coverage profile between exclusive sets A and B from a sample

sr =
  MEDIPS.saturation(file = bam_file,
                    BSgenome = BSgenome,
                    uniq = uniq,
                    extend = extend,
                    shift = shift,
                    paired = paired,
                    window_size = ws,
                    chr.select = chr.select,
                    nit = 10, nrit = 1, empty_bins = TRUE, rank = FALSE)
print(paste0("Estimated correlation is: ", round(sr$maxEstCor[2], 5)))
print(paste0("True correlation is: ", round(sr$maxTruCor[2],5)))

pdf(paste0(out_dir, fname, ".MEDIPS.SaturationPlot.pdf"), width = 5, height = 4)
MEDIPS.plotSaturation(sr)
dev.off()

## sequence coverage analysis ##########################
## outputs #of CpGs covered/not covered in a sample

cr =
  MEDIPS.seqCoverage(file = bam_file,
                     pattern = "CG",
                     BSgenome = BSgenome,
                     uniq = uniq,
                     extend = extend,
                     shift = shift,
                     paired = paired,
                     chr.select = chr.select)
print(paste0("Total number of reads: ", cr$numberReads))
print(paste0("Number of reads NOT covering a CpG: ", cr$numberReadsWO))
print(paste0("Fraction of reads NOT covering a CpG: ", round(cr$numberReadsWO / cr$numberReads, 5)))

print(paste0("Number of CpGs in reference: ", length(cr$cov.res)))
print(paste0("Number of CpG not covered by a read: ", length(cr$cov.res[cr$cov.res < 1])))
print(paste0("Number of CpG covered by 1 read: ", length(cr$cov.res[cr$cov.res == 1])))
print(paste0("Number of CpG covered by 2 reads: ", length(cr$cov.res[cr$cov.res == 2])))
print(paste0("Number of CpG covered by 3 reads: ", length(cr$cov.res[cr$cov.res == 3])))
print(paste0("Number of CpG covered by 4 reads: ", length(cr$cov.res[cr$cov.res == 4])))
print(paste0("Number of CpG covered by 5 reads: ", length(cr$cov.res[cr$cov.res == 5])))
print(paste0("Number of CpG covered by >5 reads: ", length(cr$cov.res[cr$cov.res > 5])))

print(paste0("Fraction of CpG not covered by a read: ", round(length(cr$cov.res[cr$cov.res < 1]) / length(cr$cov.res),5)))
print(paste0("Fraction of CpG covered by 1 read: ", round(length(cr$cov.res[cr$cov.res == 1]) / length(cr$cov.res),5)))
print(paste0("Fraction of CpG covered by 2 reads: ", round(length(cr$cov.res[cr$cov.res == 2]) / length(cr$cov.res),5)))
print(paste0("Fraction of CpG covered by 3 reads: ", round(length(cr$cov.res[cr$cov.res == 3]) / length(cr$cov.res),5)))
print(paste0("Fraction of CpG covered by 4 reads: ", round(length(cr$cov.res[cr$cov.res == 4]) / length(cr$cov.res),5)))
print(paste0("Fraction of CpG covered by 5 reads: ", round(length(cr$cov.res[cr$cov.res == 5]) / length(cr$cov.res),5)))
print(paste0("Fraction of CpG covered by >5 reads: ", round(length(cr$cov.res[cr$cov.res > 5]) / length(cr$cov.res),5)))


pdf(paste0(out_dir, fname, ".MEDIPS.seqCovPie.pdf"), width = 5, height = 4)
MEDIPS.plotSeqCoverage(seqCoverageObj=cr,
                       type="pie",
                       cov.level = c(0,1,2,3,4,5),
                       main="Sequence pattern coverage, pie chart")
dev.off()

pdf(paste0(out_dir, fname, ".MEDIPS.seqCovHist.pdf"), width = 5, height = 4)
MEDIPS.plotSeqCoverage(seqCoverageObj=cr,
                       type="hist",
                       t = 15,
                       main="Sequence pattern coverage, histogram")
dev.off()

## CpG enrichment #####################################
## test CpG enrichment of given set of short reads covering a set of genomic regions vs reference genome
## regions.relH - relative freq of CpGs within a sample's immunoprecipitated regions
## genome.relH - relative freq of CpGs within reference genome
## enrichment.score.relH - regions.relH/genome.relH
## regions.GoGe - obs/exp ratio of CpGs within a sample's immunoprecipitated regions
## genome.GoGe - obs/exp ratio of CpGs within genomic regions
## enrichment.score.GoGe - regions.GoGe/genome.GoGe
## (relH and GoGe = 2 different ways of calculating enrichment)


## original MEDIPS.CpGenrich has IRanges issue
## this is adapted from script by Nick Cheng
MEDIPS.CpGenrichNew <-
  function(file=NULL, BSgenome=NULL, extend=0, shift=0, uniq=1e-3, chr.select=NULL, paired=F){

    ## Proof correctness....
    if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}

    ## Read region file
    fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
    path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1], collapse="/")
    if(path==""){path=getwd()}
    if(!fileName%in%dir(path)){stop(paste("File", fileName, " not found in", path, sep =" "))}

    dataset = get(ls(paste("package:", BSgenome, sep = ""))[1])

    if(!paired){GRange.Reads = getGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}
    else{GRange.Reads = getPairedGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}

    ## Sort chromosomes
    if(length(unique(seqlevels(GRange.Reads)))>1){chromosomes=mixedsort(unique(seqlevels(GRange.Reads)))}
    if(length(unique(seqlevels(GRange.Reads)))==1){chromosomes=unique(seqlevels(GRange.Reads))}

    ## Get chromosome lengths for all chromosomes within data set.
    cat(paste("Loading chromosome lengths for ",BSgenome, "...\n", sep=""))

    chr_lengths=as.numeric(seqlengths(dataset)[chromosomes])

    ranges(GRange.Reads) <- restrict(ranges(GRange.Reads),+1)

    ##Calculate CpG density for regions
    total=length(chromosomes)
    cat("Calculating CpG density for given regions...\n")

    ## new code ##################################
    readsChars <- unlist(getSeq(dataset, GRange.Reads, as.character=TRUE))

    regions.CG = sum(vcountPattern("CG",readsChars))
    regions.C  = sum(vcountPattern("C",readsChars))
    regions.G  = sum(vcountPattern("G",readsChars))
    all.genomic= sum(width(readsChars))

    nReads <- length(readsChars)
    ###############################################

    regions.relH=as.numeric(regions.CG)/as.numeric(all.genomic)*100
    regions.GoGe=(as.numeric(regions.CG)*as.numeric(all.genomic))/(as.numeric(regions.C)*as.numeric(regions.G))

    cat(paste("Calculating CpG density for the reference genome",
              BSgenome, "...\n", sep = " "))

    CG <- DNAStringSet("CG")
    pdict0 <- PDict(CG)
    params <- new("BSParams", X = dataset, FUN = countPDict, simplify = TRUE, exclude = c("rand", "chrUn"))
    genome.CG=sum(bsapply(params, pdict = pdict0))
    params <- new("BSParams", X = dataset, FUN = alphabetFrequency, exclude = c("rand", "chrUn"), simplify=TRUE)
    alphabet=bsapply(params)
    genome.l=sum(as.numeric(alphabet))
    genome.C=as.numeric(sum(alphabet[2,]))
    genome.G=as.numeric(sum(alphabet[3,]))
    genome.relH=genome.CG/genome.l*100
    genome.GoGe=(genome.CG*genome.l)/(genome.C*genome.G);

    ##Calculate CpG density for reference genome

    enrichment.score.relH=regions.relH/genome.relH
    enrichment.score.GoGe=regions.GoGe/genome.GoGe

    gc()
    return(list(genome=BSgenome,
                regions.CG=regions.CG,
                regions.C=regions.C,
                regions.G=regions.G,
                regions.relH=regions.relH,
                regions.GoGe=regions.GoGe,
                genome.C=genome.C,
                genome.G=genome.G,
                genome.CG=genome.CG,
                genome.relH=genome.relH,
                genome.GoGe=genome.GoGe,
                enrichment.score.relH=enrichment.score.relH,
                enrichment.score.GoGe=enrichment.score.GoGe))
  }

er =
  MEDIPS.CpGenrichNew(file = bam_file,
                   BSgenome = BSgenome,
                   uniq = uniq,
                   extend = extend,
                   shift = shift,
                   paired = paired,
                   chr.select = chr.select)

## medips.satr.est_cor and medips.satr.tru_cor involve randomness, will not give identical results on repeat runs
## rest of metrics should be identical on repeat runs

message("Writing out MEDIPS QC metrics: saturation, CpG coverage and CpG enrichment.")
QC_MEDIPS.df =
  data.frame(QC_type = rep("medips_QC", 33),
             metrics = c("ref_genome",
                         "satr.est_cor",
                         "satr.tru_cor",
                         "CpG_cov.totalNumReads",
                         "CpG_cov.numReadsWoCpG",
                         "CpG_cov.fracReadsWoCpG",
                         "CpG_cov.numCpGinRef",
                         "CpG_cov.numCpGwoReads",
                         "CpG_cov.numCpGw1read",
                         "CpG_cov.numCpGw2Reads",
                         "CpG_cov.numCpGw3Reads",
                         "CpG_cov.numCpGw4Reads",
                         "CpG_cov.numCpGw5Reads",
                         "CpG_cov.numCpGgt5Reads",
                         "CpG_cov.fracCpGwoReads",
                         "CpG_cov.fracCpGw1read",
                         "CpG_cov.fracCpGw2Reads",
                         "CpG_cov.fracCpGw3Reads",
                         "CpG_cov.fracCpGw4Reads",
                         "CpG_cov.fracCpGw5Reads",
                         "CpG_cov.fracCpGgt5Reads",
                         "enrich.regions.C",
                         "enrich.regions.G",
                         "enrich.regions.CG",
                         "enrich.genome.C",
                         "enrich.genome.G",
                         "enrich.genome.CG",
                         "enrich.regions.relH",
                         "enrich.genome.relH",
                         "enrich.regions.GoGe",
                         "enrich.genome.GoGe",
                         "enrich.enrichment.score.relH",
                         "enrich.enrichment.score.GoGe"),
             values = c(er$genome,
                        round(sr$maxEstCor[2], 5),
                        round(sr$maxTruCor[2], 5),
                        cr$numberReads,
                        cr$numberReadsWO,
                        round(cr$numberReadsWO / cr$numberReads, 5),
                        length(cr$cov.res),
                        length(cr$cov.res[cr$cov.res < 1]),
                        length(cr$cov.res[cr$cov.res == 1]),
                        length(cr$cov.res[cr$cov.res == 2]),
                        length(cr$cov.res[cr$cov.res == 3]),
                        length(cr$cov.res[cr$cov.res == 4]),
                        length(cr$cov.res[cr$cov.res == 5]),
                        length(cr$cov.res[cr$cov.res > 5]),
                        round(length(cr$cov.res[cr$cov.res < 1]) / length(cr$cov.res), 5),
                        round(length(cr$cov.res[cr$cov.res == 1]) / length(cr$cov.res), 5),
                        round(length(cr$cov.res[cr$cov.res == 2]) / length(cr$cov.res), 5),
                        round(length(cr$cov.res[cr$cov.res == 3]) / length(cr$cov.res), 5),
                        round(length(cr$cov.res[cr$cov.res == 4]) / length(cr$cov.res), 5),
                        round(length(cr$cov.res[cr$cov.res == 5]) / length(cr$cov.res), 5),
                        round(length(cr$cov.res[cr$cov.res > 5]) / length(cr$cov.res), 5),
                        er$regions.C,
                        er$regions.G,
                        er$regions.CG,
                        er$genome.C,
                        er$genome.G,
                        er$genome.CG,
                        round(er$regions.relH, 5),
                        round(er$genome.relH, 5),
                        round(er$regions.GoGe, 5),
                        round(er$genome.GoGe, 5),
                        round(er$enrichment.score.relH, 5),
                        round(er$enrichment.score.GoGe, 5)))
names(QC_MEDIPS.df) = c("QC_type", "metrics", fname)
# QC_MEDIPS.df

write.table(QC_MEDIPS.df, paste0(out_dir, fname, "_QC_MEDIPS.csv"), row.names=F, quote=F, sep = '\t')

## EOF
