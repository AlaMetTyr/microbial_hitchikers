#Load Libraries
library(dada2)

#Setting Path for input files
setwd('/home/a.vaughan/nobackup_ga03488/Amy/tents/Demultiplexed_16S')
path='/home/a.vaughan/nobackup_ga03488/Amy/tents/Demultiplexed_16S'

#list.files(path, full.names = TRUE)

#Specify forward & reverse read fastqs
fnFs <- sort(list.files(path, pattern="_515f.forward.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_515f.reverse.fastq.gz", full.names = TRUE))

sample.names <- sub(".*(16S_\\d+)_.*", "\\1", basename(fnFs))

sample.names

#Filtering reads based on sequence quality scores and primers
#Set file path
filtFs <- file.path(path, "filtered_16S", paste0(sample.names, "_F_filt.fastq.gz")) #creating filtered folder
filtRs <- file.path(path, "filtered_16S", paste0(sample.names, "_R_filt.fastq.gz"))

#Pulling sample names from filtered fasta files
names(filtFs) <- sample.names 
names(filtRs) <- sample.names

length(fnFs)
length(fnRs)
length(filtFs)
length(filtRs)

#Trimming primers - set your primer sequence
FWD <- "NNNNNNGTGYCAGCMGCCGCGGTAA"  #515f
REV <- "NNNNNNGGACTACNVGGGTWTCTAAT"   #806r
trimLeft = c(FWD,REV)

#Use known primer sequences to trim from your amplicon sequences
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,200),   
            maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,             # filter out all reads with > maxN=0 ambiguous nucleotides and >2 expected errors
              compress=TRUE, multithread=TRUE,trimLeft = c(25,26))         # remove first 26 nucleotides of F/R reads (length of primers?)


# rm.phix is default and removes reads that match the phiX genome
# truncQ=2 is deafult and truncates reads at first instance of a quality score less than or equal to 2

#exists <- file.exists(filtFs) & file.exists(filtRs)
#filtFs <- filtFs[exists]
#filtRs <- filtRs[exists]

# Learning Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plotting out the errors
png(filename="Error_16S_F.png")
plotErrors(errF, nominalQ=TRUE)
dev.off()
png(filename="Error_16S_R.png")
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
# 304 bimeras out of 7501 input sequences
dim(seqtab.nochim)
# 65 7197
sum(seqtab.nochim)/sum(seqtab)
# 0.96132

# Track reads through the pipeline
# Good checkpoint to ensure you did not lose too many reads

getN <- function(x) sum(getUniques(x))
track<- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,getN), rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "SequencingStatistics_16S.csv")

taxa <- assignTaxonomy(seqtab.nochim, "/home/a.vaughan/nobackup_ga03488/Amy/tents/silva_nr99_v138.2_toSpecies_trainset.fa.gz", 
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
write(asv_fasta, "ASVs_16S_bf.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts_16S_bf.tsv", sep="\t", quote=F, col.names=NA)

##  Giving taxonomy table corresponding names as above (ASV_1, ASV_2...)

row.names(taxa.print) <- sub(">", "", asv_headers)

write.table(taxa.print, "ASVs_assigned_16S_bf2.tsv")
write.table(taxa.print, "ASVs_assigned_16S_bf.tsv", sep="\t", quote=F, col.names=NA)