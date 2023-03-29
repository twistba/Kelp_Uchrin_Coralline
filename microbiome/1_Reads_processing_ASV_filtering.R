#######################################
# Prepare dataset with dada2 pipeline #
#######################################

# 0. Prep environment 
# 1. Remove Primers
# 2. Filter and Trim 
# 3. Error model and Run dada2
# 4. Merge reads and remove chimera
# 5. Assign TAXONOMY
# 6. Filter ASVs 
# 7. Rarefy 

#########################
# 0. Prep environment 
#########################

library(dada2)
library(tidyverse)
library(DECIPHER)
library(vegan)
library(phyloseq)

cutadapt <- "/Users/fmazel/Softwares/miniconda3/envs/cutadaptenv/bin/cutadapt" # path to cutadapt on the computer 


path <- "Data/raw_sequencing_output/Fastq_files" # location of data on computer (after downloadign from repo)
list.files(path)

# Define other paths
path_dada2 <- "Data/dada2_processed"
path_to_trimmed_reads = "Data/Primers_trimmed_Sequences/"

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#########################
# 1. Remove Primers
#########################

# summary of pathways 
reads_path = tibble(ForwardReadsPath = fnFs,
                    ReverseReadsPath = fnRs,
                    ForwardReadsName = sort(list.files(path, pattern="_R1_001.fastq", full.names = F)),
                    ReverseReadsName = sort(list.files(path, pattern="_R2_001.fastq", full.names = F)))


reads_path = reads_path %>% 
  mutate(Trimmed_ForwardReadsPath = paste0(path_to_trimmed_reads,ForwardReadsName),
         Trimmed_ReverseReadsPath = paste0(path_to_trimmed_reads,ReverseReadsName))


system2(cutadapt, args = "--help")
F_primers_Freads = "GTGYCAGCMGCCGCGGTAA" #515F
F_primers_Rreads = "CCGYCAATTYMTTTRAGTTT" #926R

error=0.1
#error 0 is the baseline 

for (i in 1:52)
{
  params_cutadapt <- c('-e', error, # 
                       '-m', 100, # discard all reads < 100 pb
                       '-g', F_primers_Freads, # Forward barcode 
                       '-G', F_primers_Rreads,  # Reverse barcode 
                       '-o',reads_path$Trimmed_ForwardReadsPath[i], # R1 output
                       '-p',reads_path$Trimmed_ReverseReadsPath[i], # R2 output 
                       reads_path$ForwardReadsPath[i], # Input R1
                       reads_path$ReverseReadsPath[i] # Input R2 
  )
  
  
  system2(cutadapt, args = params_cutadapt)
  
}


#########################
# 2. Filter and Trim 
#########################

fnFs <- sort(list.files(path_to_trimmed_reads, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path_to_trimmed_reads, pattern="_R2_001.fastq", full.names = TRUE))

# Quality profiles 
plotQualityProfile(fnFs[1:2]) #trim at 270pb 
plotQualityProfile(fnRs[1:5]) #trim at 210pb


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path_dada2 , "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_dada2 , "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#########################
# 2. Filter and Trim 
#########################

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(270,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)


#########################
# 3. Error model and Run
#########################

errF <- learnErrors(filtFs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Sample inferences
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)



#########################
# 4. Merge reads and remove chimera
#########################

# Merge Forwards and Reverses reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Create ASV*sample table 
seqtab <- makeSequenceTable(mergers)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Remove chimera
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


write.table(seqtab.nochim, "Data/dada2_processed/unfiltered_ASV_table.txt")


dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Check out reads 
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign Taxonomy 
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/Data/SILVA/SILVA_SSU_r138_2019.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors

ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

write.table(taxid, "Data/dada2_processed/unfiltered_ASV_taxonomy.txt")


#################################################
# 6. Remove ASVs based on taxonomy and length 
#################################################


taxid = read.table("Data/dada2_processed/unfiltered_ASV_taxonomy.txt")
seqtab.nochim = read.table("Data/dada2_processed/unfiltered_ASV_table.txt")


# Normalize
depth <- apply(seqtab.nochim,1,sum)
seqtab.nochim_norm = seqtab.nochim / depth

rel_read_counts <- tibble(rel_read_counts=apply(seqtab.nochim_norm,2,sum, na.rm=T),
                          seq=colnames(seqtab.nochim_norm))

taxonomy = as_tibble(taxid, rownames = "seq") %>% 
  left_join(rel_read_counts)

# Check
taxonomy %>% group_by(domain) %>% 
  summarise(n=n(),
            sum_rrc=sum(rel_read_counts))

#domain       n sum_rrc
#Archaea     19  0.0577
#Bacteria  4484 50.0   
#NA        1377  1.95  


# Remove ASVs assigned to chloroplast, mitochondria, or any family 
filtered_taxonomy = taxonomy %>% 
  subset(!order == "Chloroplast") %>% 
  subset(!family == "Mitochondria") %>% 
  subset(!is.na(domain))


# Remove Long ASVs (non cutted primer)
expected_length = 926-515-nchar(F_primers_Freads)-nchar(F_primers_Rreads)
expected_length # 372

filtered_taxonomy = filtered_taxonomy %>% 
  subset(length_ASV<380) 

write.table(filtered_taxonomy,"Data/dada2_processed/filtered_ASV_taxonomy.txt") # 1988 ASVs

filtered_taxonomy = read.table("Data/dada2_processed/filtered_ASV_taxonomy.txt")
seqtab.nochim.filtASVs= seqtab.nochim[,filtered_taxonomy$seq]
write.table(seqtab.nochim.filtASVs, "Data/dada2_processed/filtered_ASV_table.txt")


#############
# 7. Rarefy 
#############

depth <- tibble(depth=apply(seqtab.nochim.filtASVs,1,sum),
                swab.ID=rownames(seqtab.nochim.filtASVs))

sample_info <-  read_csv("Data/metadata/Metadata_updated_BT.csv") %>% 
  left_join(depth)

sample_info %>% 
  group_by(Substrate) %>% 
  summarise(n=n(),depth=median(depth))

hist(sample_info$depth)

set.seed(21)

# Rarefying 1000
seqtab.nochim.filtASVs_1000reads_rar = rrarefy(seqtab.nochim.filtASVs, 1000)
ASV_present = apply(seqtab.nochim.filtASVs_1000reads_rar,2,sum)>0
seqtab.nochim.filtASVs_1000reads_rar = seqtab.nochim.filtASVs_1000reads_rar[,ASV_present]
dim(seqtab.nochim.filtASVs_1000reads_rar)
write.table(seqtab.nochim.filtASVs_1000reads_rar, "Data/dada2_processed/filtered_ASV_table_rarefied.to.1000.txt")


#############
# 8.Phyloseq objects
#############


# Format Phyloseq object
MT <- as.data.frame(sample_info); rownames(MT)=MT$swab.ID
PT <- as.matrix(filtered_taxonomy); rownames(PT)=PT[,'seq']

phyloseq_obj=phyloseq(tax_table(PT),otu_table(seqtab.nochim.filtASVs_1000reads_rar,taxa_are_rows=F),sample_data(MT))
phyloseq_obj
saveRDS(phyloseq_obj,"Data/dada2_processed/Phyloseq_object_rarefied.1000.RDS")















