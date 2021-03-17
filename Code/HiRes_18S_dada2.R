
#########################################
### 18S + STDs processing, MiSeq plate #2
### Lou LaMartina, Mar 16 2021
#########################################


library(dada2)
library(decontam)
library(ggplot2)
library(phyloseq)
library(vegan)
library(ggplot2)


###################################
### prepare data for processing ###
###################################

# set working directory
setwd("~/Desktop/Lab/Projects/Kazuaki/MiSeq2/18S")


# set file paths
path <- "./cutadapt"
pathF <- "./cutadapt/fastqF"
pathR <- "./cutadapt/fastqR"


# set paths for filtered reads
filtered_pathF <- "./cutadapt/fastqF/Filtered"
filtered_pathR <- "./cutadapt/fastqR/Filtered"


# sort forward and reverse reads into file paths
fastqFs <- sort(list.files(pathF, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fastqRs <- sort(list.files(pathR, pattern = "_R2_001.fastq.gz", full.names = TRUE))


# extract file names
sample_names <- gsub("-18S", "", sapply(strsplit(basename(fastqFs), "_"), '[', 1))




################################
### Inspect sequence quality ###
################################

# visualize quality of reads
qualityF.plot <- plotQualityProfile(fastqFs[1:4]); qualityF.plot
qualityR.plot <- plotQualityProfile(fastqRs[1:4]); qualityR.plot


# check: if there are not the same number of F and R files, stop.
if(length(fastqFs) != length(fastqRs)) 
  stop("Forward and reverse files do not match.")


# save quality profiles
#ggsave("./Plots/qualityF.pdf", plot = qualityF.plot, device = "pdf", width = 12, height = 8, units = "in")
#ggsave("./Plots/qualityR.pdf", plot = qualityR.plot, device = "pdf", width = 12, height = 8, units = "in")




#######################
### Quality control ###
#######################

# give filtered files new names and paths
filteredFs <- file.path(filtered_pathF, paste0(sample_names, "_F_filt.fastq.gz"))
filteredRs <- file.path(filtered_pathR, paste0(sample_names, "_R_filt.fastq.gz"))


# filter based on quality and read length
filtered_out <- filterAndTrim(fastqFs, filteredFs, fastqRs, filteredRs, 
                              maxEE = 2, maxN = 0, truncQ = 10, # not truncating
                              rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = TRUE)


# inspect how many reads were filtered out of each sample
filtered_out
(1 - (filtered_out[,2] / filtered_out[,1])) * 100
cat(mean((1 - (filtered_out[,2] / filtered_out[,1])) * 100), "% removed\n")
# 12.22103 % removed


# set sample names to the ID only
names(filteredFs) <- sample_names
names(filteredRs) <- sample_names


# plot quality profiles of filtered reads
filtF.plot <- plotQualityProfile(filteredFs[1:4]); filtF.plot
filtR.plot <- plotQualityProfile(filteredRs[1:4]); filtR.plot


# save quality profiles
#ggsave("./Plots/filt_qualityF.pdf", plot = filtF.plot, device = "pdf", width = 12, height = 8, units = "in")
#ggsave("./Plots/filt_qualityR.pdf", plot = filtR.plot, device = "pdf", width = 12, height = 8, units = "in")


# save env
save.image("./RData/MiSeq2_18S_dada2_env.RData")




############################
### Learning error rates ###
############################

# learn and visualize error rates of F reads
errorF <- learnErrors(filteredFs, multithread = TRUE)
errorF.plot <- plotErrors(errorF, nominalQ = TRUE, err_in = TRUE, err_out = TRUE); errorF.plot


# learn and visualize error rates of R reads
errorR <- learnErrors(filteredRs, multithread = TRUE)
errorR.plot <- plotErrors(errorR, nominalQ = TRUE, err_in = TRUE, err_out = TRUE); errorR.plot


# save error plots
#ggsave("./Plots/errorF.pdf", plot = errorF.plot, device = "pdf", width = 12, height = 8, units = "in")
#ggsave("./Plots/errorR.pdf", plot = errorR.plot, device = "pdf", width = 12, height = 8, units = "in")




################################
### Merging paired-end reads ###
################################

# create list of merged reads
mergers <- vector("list", length(sample_names))
names(mergers) <- sample_names


# sample inference and merging paired-end reads
for(i in sample_names) {
  cat("\nProcessing", i, "(", match(i, sample_names),"/", 
      length(sample_names), ") -", format(Sys.time(), "%I:%M %p"), "\n")
  derepF <- derepFastq(filteredFs[[i]])
  dadaF <- dada(derepF, err = errorF, multithread = TRUE)
  derepR <- derepFastq(filteredRs[[i]])
  dadaR <- dada(derepR, err = errorR, multithread = TRUE)
  merger <- mergePairs(dadaF, derepF, dadaR, derepR)
  mergers[[i]] <- merger
}


# remove dereps to save memory
rm(derepF, derepR)


# construct a sequence table: number of unique sequences (ASVs) in each sample
counts_table <- makeSequenceTable(mergers)
cat(ncol(counts_table), "ASVs in", nrow(counts_table), "samples\n")
# 1140 ASVs in 33 samples




##################################
### Quality control: processed ###
##################################


########
### trim

# ~~~ not doing - we are not familiar with the 
#     distribution of this 18S gene region lengths ~~~



###################
### Remove chimeras

# removing chimeras with denovo screening
counts_noChim <- removeBimeraDenovo(counts_table, method = "consensus",
                                    multithread = TRUE, verbose = TRUE)


# how many unique sequences were moved?
cat(ncol(counts_noChim), "out of", ncol(counts_table), "ASVs passed (", 
    (1 - ncol(counts_noChim) / ncol(counts_table)) * 100, "% removed )\n")
# without trimming: 712 out of 1177 ASVs passed ( 39.50722 % removed )
# with trimming: 311 out of 730 ASVs passed ( 57.39726 % removed )
# second time without trimming: 684 out of 1140 ASVs passed ( 40 % removed )


# what percentage of reads were identified as chimeras?
cat((1 - sum(counts_noChim) / sum(counts_table)) * 100, "% reads removed\n")
# without trimming: 32.64434 % reads removed
# with trimming: 32.58731 % reads removed
# second time without trimming: 32.64817 % reads removed
# decided to not trim



######################
### identify standards

# save as fasta
uniquesToFasta(counts_noChim,
               ids = paste0("NOCHIM", sprintf("%06d", 1:ncol(counts_noChim))), 
               fout = "./RData/Kazuaki_18S_nochim.fasta")


# in terminal:
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# $ makeblastdb -dbtype nucl -in Kazuaki_18S_nochim.fasta -input_type fasta
#
# $ blastn -db Kazuaki_18S_nochim.fasta -query 18S_STDs_trim.fasta -task blastn -perc_identity 100 -outfmt 6 -out STD_align_nochim.txt
#
# $ echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n$(cat STD_align_nochim.txt)" > STD_align_nochim.txt
#
# $ sed 's/\t/,/g' STD_align_nochim.txt > STD_align_nochim.csv
#
# $ rm Kazuaki_18S_nochim.fasta.n*
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# load alignment results
blastn_nochim <- read.csv("./RData/STD_align_nochim.csv")


# create vector of fastas associated with SEQ IDs
fastas_nochim <- data.frame(FASTA = colnames(counts_noChim))
fastas_nochim$sseqid <- paste0("NOCHIM", sprintf("%06d", 1:ncol(counts_noChim)))


# locate those that match standards
blastn_nochim[grep("pro1", blastn_nochim$qseqid), "STD"] <- "STD1"
blastn_nochim[grep("pro2", blastn_nochim$qseqid), "STD"] <- "STD2"
blastn_nochim[grep("pro3", blastn_nochim$qseqid), "STD"] <- "STD3"
blastn_nochim[grep("pro4", blastn_nochim$qseqid), "STD"] <- "STD4"


# add FASTA seqs
blastn_nochim <- merge(blastn_nochim, fastas_nochim, by = "sseqid")


# remove STDs
STDs_counts <- counts_noChim[,colnames(counts_noChim) %in% blastn_nochim$FASTA]
ncol(STDs_counts) == length(unique(blastn_nochim$FASTA))

counts_noChimSTDs <- counts_noChim[, ! colnames(counts_noChim) %in% blastn_nochim$FASTA]
ncol(counts_noChimSTDs) == ncol(counts_noChim) - ncol(STDs_counts)

cat(ncol(counts_noChimSTDs), "ASVs in", nrow(counts_noChimSTDs), "samples\n")
# 627 ASVs in 33 samples


# combine STDs; multiple matches for each but derived from single sequence
STD_counts.glom <- data.frame(STD1 = rowSums(STDs_counts[, colnames(STDs_counts) %in% subset(blastn_nochim, STD == "STD1")$FASTA]),
                         STD2 = rowSums(STDs_counts[, colnames(STDs_counts) %in% subset(blastn_nochim, STD == "STD2")$FASTA]),
                         STD3 = rowSums(STDs_counts[, colnames(STDs_counts) %in% subset(blastn_nochim, STD == "STD3")$FASTA]),
                         STD4 = rowSums(STDs_counts[, colnames(STDs_counts) %in% subset(blastn_nochim, STD == "STD4")$FASTA]))




#######################
### Assign taxonomy ###
#######################

# yes, you can use the silva database for fungi!
taxa_table <- assignTaxonomy(counts_noChimSTDs, 
                             "~/Desktop/Lab/Projects/Misc/silva_nr_v132_train_set.fa.gz", 
                             multithread = TRUE)


# save
save.image("./RData/MiSeq2_18S_dada2_env.RData")




##################################
### remove mock community ASVs ###
##################################

# create fasta
uniquesToFasta(counts_noChimSTDs, ids = paste0("nochimstd", sprintf("%06d", 1:ncol(counts_noChimSTDs))),
               fout = "./RData/SeqtabnochimSTDs.fasta")


# create taxonomy data frame
taxa_table.df <- data.frame(taxa_table)
taxa_table.df$sseqid <- paste0("nochimstd", sprintf("%06d", 1:ncol(counts_noChimSTDs)))
taxa_table.df$FASTA <- rownames(taxa_table.df)


# in terminal:
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# $ makeblastdb -dbtype nucl -in SeqtabnochimSTDs.fasta -input_type fasta
#
# $ blastn -db SeqtabnochimSTDs.fasta -query ~/Desktop/Lab/Projects/Misc/zymo_mock_all.fasta -task blastn -perc_identity 100 -outfmt 6 -out mock_align.txt
#
# $ echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n$(cat mock_align.txt)" > mock_align.txt
#
# $ sed 's/\t/,/g' mock_align.txt > mock_align.csv
#
# $ rm SeqtabnochimSTDs.fasta.n*
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# load
mock_align <- read.csv("RData/mock_align.csv")


# keep only best matches
mock_align <- subset(mock_align, bitscore > 400)


# compare
mock_ASVs <- unique(merge(taxa_table.df, mock_align[1:2], by = "sseqid"))


# remove
counts_noMock <- counts_noChimSTDs[,! colnames(counts_noChimSTDs) %in% mock_ASVs$FASTA]
ncol(counts_noChimSTDs) - ncol(counts_noMock) == length(unique(mock_ASVs$FASTA))


# save
mock_ASVs <- subset(taxa_table.df, FASTA %in% mock_ASVs$FASTA)
mock_ASVs$Source <- "Mock"




###############################
### remove nonspecific ASVs ###
###############################

# subset eukaryotes, mitochondria, chloroplasts
bac_ASVs <- subset(taxa_table.df, Kingdom == "Bacteria")
chloro_ASVs <- subset(taxa_table.df, Order == "Chloroplast")
mito_ASVs <- subset(taxa_table.df, Family == "Mitochondria")


# combine them
bac_ASVs$Source <- "Bacteria"
chloro_ASVs$Source <- "Chloroplasts"
mito_ASVs$Source <- "Mitochondria"
contaminants <- rbind(mock_ASVs, bac_ASVs, mito_ASVs, chloro_ASVs)


# remove from data
counts_noContam <- counts_noMock[,! colnames(counts_noMock) %in% contaminants$FASTA]
ncol(counts_noChimSTDs) - ncol(counts_noContam) == length(unique(contaminants$FASTA))

taxa_noContam.df <- subset(taxa_table.df, ! FASTA %in% contaminants$FASTA)
taxa_noContam.df <- data.frame(FASTA = taxa_noContam.df$FASTA,
                               ASV = paste0("ASV", sprintf("%06d", 1:ncol(counts_noContam))), 
                               taxa_noContam.df[1:6])
identical(rownames(taxa_noContam.df), colnames(counts_noContam))




###################
### rarefaction ###
###################
# https://stat.ethz.ch/pipermail/r-sig-ecology/2018-December/005867.html


# export for ggplot
raredata <- rarecurve(counts_noContam, step = 20, sample = min(rowSums(counts_noContam)), col = "blue", cex = 0.6)
names(raredata) <- rownames(counts_noContam)


# coerce data into "long" form.
protox <- mapply(FUN = function(x, y) {
  mydf <- as.data.frame(x)
  colnames(mydf) <- "value"
  mydf$species <- y
  mydf$subsample <- attr(x, "Subsample")
  mydf
}, x = raredata, y = as.list(names(raredata)), SIMPLIFY = FALSE)

raredata <- do.call(rbind, protox)
rownames(raredata) <- NULL


# plot
rare.plot <-
  ggplot(raredata, aes(x = subsample, y = value, color = species)) +
  geom_vline(xintercept = min(rowSums(counts_noContam))) +
  theme_bw() +
  scale_color_discrete() +
  geom_line() +
  theme(legend.text = element_text(size = 6)) +
  guides(color = guide_legend(keyheight = 0.5, keywidth = 0.2, units = "in", ncol = 1)) +
  labs(x = "Sample size", y = "No. of ASVs", color = "Sample")
rare.plot

#ggsave("./Plots/rarefaction.pdf", plot = rare.plot, device = "pdf", width = 8, height = 5)




#########################
### Organize and save ### 
#########################

# save FASTA
uniquesToFasta(counts_noContam,
               ids = paste0(taxa_noContam.df$ASV, "__",
                            taxa_noContam.df$Kingdom, "__",
                            taxa_noContam.df$Phylum, "__",
                            taxa_noContam.df$Class, "__",
                            taxa_noContam.df$Order, "__",
                            taxa_noContam.df$Family, "__",
                            taxa_noContam.df$Genus),
               fout = "./RData/MiSeq2_18S_dada2.fasta")


# add STDs
counts_noContam.df <- data.frame(Sample_name = rownames(counts_noContam), counts_noContam)
STD_counts.glom <- data.frame(Sample_name = rownames(STD_counts.glom), STD_counts.glom)
identical(rownames(counts_noContam.df), rownames(STD_counts.glom))
counts.df <- cbind(STD_counts.glom[-1], counts_noContam.df[-1])


# transpose
counts.df <- data.frame(t(counts.df))
counts.df <- data.frame(FASTA = rownames(counts.df), counts.df)


# combine
all <- merge(taxa_noContam.df, counts.df, by = "FASTA", all.y = TRUE)
all <- all[order(all$ASV),]


# bring STDs to top
all$ASV[grep("STD", all$FASTA)] <- "STD"
all <- rbind(subset(all, ASV == "STD"), subset(all, ASV != "STD"))


# save
write.csv(all, "./RData/MiSeq2_18S_dada2.csv", row.names = FALSE, na = "")



# # # # # # # # # # #
# save R environment
save.image("./RData/MiSeq2_18S_dada2_env.RData")



