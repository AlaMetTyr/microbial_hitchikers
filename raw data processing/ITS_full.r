if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

#Load Libraries
library(dada2)

#Setting Path for input files
setwd('/home/a.vaughan/nobackup_ga03488/Amy/tents/Demultiplexed_ITS')
path='/home/a.vaughan/nobackup_ga03488/Amy/tents/Demultiplexed_ITS'
#list.files(path)

#Specify forward & reverse read fastqs
fnFs <- sort(list.files(path, pattern="__ITS1f_KYO1.forward.fastq.gz$", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="__ITS1f_KYO1.reverse.fastq.gz$", full.names=TRUE))

sampleF <- sapply(strsplit(basename(fnFs), "__"), `[`, 2)
sampleR <- sapply(strsplit(basename(fnRs), "__"), `[`, 2)
fnRs <- fnRs[match(sampleF, sampleR)]

data.frame(
  sampleF = sampleF,
  sampleR = sapply(strsplit(basename(fnRs), "__"), `[`, 2)
)

sample.names <- sampleF


#Filtering reads based on sequence quality scores and primers
#Set file path
filtFs <- file.path(path, "filtered_ITS", paste0(sample.names, "_F_filt.fastq.gz")) #creating filtered folder
filtRs <- file.path(path, "filtered_ITS", paste0(sample.names, "_R_filt.fastq.gz"))
#Pulling sample names from filtered fasta files
names(filtFs) <- sample.names 
names(filtRs) <- sample.names
#Trimming primers - set your primer sequence
FWD <- "NNNNNNCTHGGTCATTTAGAGGAASTAA"  #ITS1F
REV <- "NNNNNNTTYRCTRCGTTCTTCATC"   #ITS2
trimLeft = c(FWD,REV)

#Use known primer sequences to trim from your amplicon sequences
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     maxN=0, maxEE=c(2,2),
                     truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,
                     trimLeft=c(25,26))  # remove primers

# rm.phix is default and removes reads that match the phiX genome
# truncQ=2 is deafult and truncates reads at first instance of a quality score less than or equal to 2

#exists <- file.exists(filtFs) & file.exists(filtRs)
#filtFs <- filtFs[exists]
#filtRs <- filtRs[exists]

# Learning Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plotting out the errors
png(filename="Error_ITS_F.png")
plotErrors(errF, nominalQ=TRUE)
dev.off()
png(filename="Error_ITS_R.png")
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)  # if throws an error, this will be because some reads didn't pass filter.
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names:
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference from forward reads
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)  # instead of filtFs
# Sample inference from reverse reads
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Making your ASV table. This is synonymous to OTU table
seqtab <- makeSequenceTable(mergers)
#dim(seqtab)

# Length of ASV's and their frequencies
table(nchar(getSequences(seqtab)))
# The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by)

# Identifying and removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",multithread=TRUE, verbose=TRUE)
# Identified 2437 bimeras out of 8312 input sequences.
#dim(seqtab.nochim)
# 138 5875
#sum(seqtab.nochim)/sum(seqtab)
# 0.7872123

# Track reads through the pipeline
# Good checkpoint to ensure you did not lose too many reads

getN <- function(x) sum(getUniques(x))
track<- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,getN), rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names
#head(track)

write.csv(track, "SequencingStatistics_ITS.csv")

taxa <- assignTaxonomy(seqtab.nochim, "/home/a.vaughan/nobackup_ga03488/Amy/tents/sh_general_release_dynamic_s_19.02.2025.fasta", 
                       tryRC=TRUE, multithread=TRUE)

# Check your taxonomy assignment

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
#dim(taxa.print)
#head(taxa.print)

# Giving our seq headers more manageable names (ASV_1, ASV_2...)
   
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
    
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs_ITS_bf.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts_ITS_bf.tsv", sep="\t", quote=F, col.names=NA)

##  Giving taxonomy table corresponding names as above (ASV_1, ASV_2...)

row.names(taxa.print) <- sub(">", "", asv_headers)

write.table(taxa.print, "ASVs_assigned_ITS_bf2.tsv")
write.table(taxa.print, "ASVs_assigned_ITS_bf.tsv", sep="\t", quote=F, col.names=NA)