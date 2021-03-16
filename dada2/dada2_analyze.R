library(dada2)
library(here)
library(stringr)
library(readr)
library(tibble)
library(tictoc)

##################################
## RETRIEVE SEQUENCE FILE PATHS ##
# Note: If files are stored in per-sample folders, use the 'recursive = T' parameter to list.files()
fwdFPs <- sort(list.files(here("data"), pattern = "*.R1_001.fastq.gz", full.names = TRUE))
revFPs <- sort(list.files(here("data"), pattern = "*.R2_001.fastq.gz", full.names = TRUE))

# Extract Sample IDs from the file names by splitting on "_"
sample_ids <- str_split_n(basename(fwdFPs), "_", n = 1)


#####################################
## SEQUENCE TRIMMING AND FILTERING ##

## Check average per-base quality for the first 4 samples ##
plotQualityProfile(fwdFPs[1:4])
plotQualityProfile(revFPs[1:4])

## Set up output file paths for filtering sample data ##
filtFs <- here("filtered", paste0(sample_ids, "_F_filt.fastq.gz"))
filtRs <- here("filtered", paste0(sample_ids, "_R_filt.fastq.gz"))
names(filtFs) <- sample_ids
names(filtRs) <- sample_ids

tic()
out <- filterAndTrim(fwdFPs, filtFs, revFPs, filtRs, truncLen = c(240,160),
                     compress = TRUE, multithread = TRUE)
toc()
filt_trim_res <- tibble(SampleIDs = sample_ids, Input = out[,1], Filtered = out[,2])
# display the # of reads pre- and post-filtering
filt_trim_res


#######################
## LEARN ERROR RATES ##

# forward reads
tic()
errF <- learnErrors(filtFs, multithread=TRUE)
toc()
# reverse reads
errR <- learnErrors(filtRs, multithread=TRUE)

# Optional
#plotErrors(errF, nominalQ=TRUE)


#########################################
## Sequence Correction (ASV Inference) ##

# forward reads
tic()
dadaFs <- dada(filtFs, err=errF, pool = "pseudo", multithread=TRUE)
toc()
# reverse reads
dadaRs <- dada(filtRs, err=errR, pool = "pseudo", multithread=TRUE)

# inspect 'dada' objects
dadaFs[[1]]
help("dada-class")


#######################
## MERGE PAIRED ENDS ##
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample
View(head(mergers[[1]]))


###########################
## CREATE SEQUENCE TABLE ##
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove ASVs significantly shorter/longer than expected
seqtab <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]

#####################
## REMOVE CHIMERAS ##
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

############################################
## OVERVIEW OF READS THROUGH THE PIPELINE ##
getN <- function(x) sum(getUniques(x))
track <- tibble(filt_trim_res, DenoisedFwd = sapply(dadaFs, getN), DenoisedRev = sapply(dadaRs, getN),
                Merged = sapply(mergers, getN), NonChimeric = rowSums(seqtab.nochim))
head(track)
write_tsv(track, here("dada2_pipeline_reads.tsv"))

###################################################
# EXPORT DADA2 RESULTS FOR DOWNSTREAM PROCESSING ##

# Write ASVs to FASTA file 
asv_ids <- paste0("ASV13_", seq(length(getSequences(seqtab.nochim))))
uniquesToFasta(seqtab.nochim, "dada2_asvs.fasta", ids=asv_ids)

# export counts table to biom
library(biomformat)
seqtab.asvids <- t(seqtab.nochim)
rownames(seqtab.asvids) <- asv_ids
st.biom <- make_biom(seqtab.asvids)
write_biom(st.biom, "asv13.biom")

#####################
## ASSIGN TAXONOMY ##
tic
taxa <- assignTaxonomy(seqtab.nochim, here("train", "silva_nr99_v138.1_wSpecies_train_set.fa.gz"), multithread=TRUE)
taxa <- addSpecies(taxa, here("train", "silva_species_assignment_v138.1.fa.gz"))

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


###########################################
## PREPARE DATA FOR PHYLOSEQ FOR ANALYIS ##
library(phyloseq)
subject <- str_split_n(sample_ids, "D", n = 1)
gender <- str_sub(subject, 1, 1)
day <- as.integer(str_split_n(sample_ids, "D", n = 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- sample_ids
write_delim(samdf, "mg-metadata.tsv", delim = "\t")

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")
