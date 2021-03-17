
################################################
### Weekly sewage time series standards analysis
### Lou LaMartina, started Mar 16, 2021
################################################


# https://github.com/ong8181/eDNA-qmiseq


# there are 110 samples from 2 16S sequencing runs:

# PLATE 1 had 94 samples,
# but only 71 were kept.               71

# PLATE 2 had 39 samples -
#   - 23 were redo's from PLATE 1    + 23
#   - 16 were remaining in the set   + 16  
#                                   = 110


setwd("~/Desktop/Lab/Projects/Kazuaki/Final/16S")
load("./RData/conversions_env.RData")

library(reshape2)
library(ggplot2)



########################
### calculate slopes ###
########################


stds_copies <- c(300000, 40000, 7000, 1000)
names(stds_copies) <- c("STD1", "STD2", "STD3", "STD4")


###########
### PLATE 1

r2.vec <- vector() # r squareds
m.vec <- vector()  # slopes


# calculate for each sample
for(smp in rownames(plate1_counts)){
  r2.vec[[smp]] <- summary(lm(stds_copies ~ 
                                as.numeric(plate1_counts[smp, names(stds_copies)]) + 0))$r.squared
  m.vec[[smp]] <- summary(lm(stds_copies ~ 
                               as.numeric(plate1_counts[smp, names(stds_copies)]) + 0))$coefficient[1]
}


# combine, add sample info, add std counts
plate1_cor.df <- data.frame(Sample_name = names(r2.vec), r2 = r2.vec, m = m.vec, plate1_stds)



###########
### PLATE 2

r2.vec <- vector() # r squareds
m.vec <- vector()  # slopes


# calculate for each sample
for(smp in rownames(plate2_counts)){
  r2.vec[[smp]] <- summary(lm(stds_copies ~ 
                                as.numeric(plate2_counts[smp, names(stds_copies)]) + 0))$r.squared
  m.vec[[smp]] <- summary(lm(stds_copies ~ 
                               as.numeric(plate2_counts[smp, names(stds_copies)]) + 0))$coefficient[1]
}


# combine, add sample info, add std counts
plate2_cor.df <- data.frame(Sample_name = names(r2.vec), r2 = r2.vec, m = m.vec, plate2_stds)




#####################################
### convert to absolute abundance ###
#####################################


# number of eDNA copies (copies/μL) = 
#    MiSeq sequence reads / 
#    a sample-specific regression slope

# in order words...

# ASV absolute abundance (copies/μL) = 
#    ASV count / 
#    slope(STD count ~ STD copies spiked in)


###########
### PLATE 1

# function
ConvertReads <- function(sample.row = 1){
  x <- plate1_counts[sample.row, 1:4]      # std counts
  y <- plate1_counts[sample.row,]          # all (std & asv) counts
  slope <- lm(as.numeric(x) ~              # std counts ~ std copies
                as.numeric(stds_copies) + 0)$coefficients 
  calc.copy <- y / slope
  return(calc.copy)
}


# convert
plate1_convert <- data.frame()
for(i in 1:nrow(plate1_counts)){
  plate1_convert <- rbind(plate1_convert, ConvertReads(sample.row = i))
}



###########
### PLATE 2

# function
ConvertReads <- function(sample.row = 1){
  x <- plate2_counts[sample.row, 1:4]      # std counts
  y <- plate2_counts[sample.row,]          # all (std & asv) counts
  slope <- lm(as.numeric(x) ~              # std counts ~ std copies
                as.numeric(stds_copies) + 0)$coefficients 
  calc.copy <- y / slope
  return(calc.copy)
}


# convert (takes a few mins)
plate2_convert <- data.frame()
for(i in 1:nrow(plate2_counts)){
  plate2_convert <- rbind(plate2_convert, ConvertReads(sample.row = i))
}




########################
### combine datasets ###
########################


#############
### originals

# transpose
plate1_counts.t <- data.frame(t(plate1_counts))
plate2_counts.t <- data.frame(t(plate2_counts))


# add fasta variable
plate1_counts.t$FASTA <- rownames(plate1_counts.t)
plate2_counts.t$FASTA <- rownames(plate2_counts.t)


# how many are there?
length(unique(c(plate2_counts.t$FASTA, plate1_counts.t$FASTA))) # [1] 14958


# merge
total_counts.t <- merge(plate1_counts.t, plate2_counts.t, by = "FASTA", all = TRUE)
length(unique(total_counts.t$FASTA)) # [1] 14958


# change NA to zero
total_counts.t[is.na(total_counts.t)] <- 0


# remove FASTA column
rownames(total_counts.t) <- total_counts.t$FASTA
total_counts.t <- total_counts.t[-1]


# order by most abundant ASVs
total_counts.t <- total_counts.t[order(rowSums(total_counts.t), decreasing = TRUE),]


# bring stds to top
total_stds.t <- total_counts.t[grep("STD", rownames(total_counts.t)),]
total_counts.t <- total_counts.t[! rownames(total_counts.t) %in% rownames(total_stds.t),]
total_counts.t <- rbind(total_stds.t, total_counts.t)


# transpose back
total_counts <- data.frame(t(total_counts.t))
total_mini <- total_counts[,1:100]


# remove empty ASVs
total_counts <- total_counts[rowSums(total_counts) > 0, colSums(total_counts) > 0]
cat(nrow(total_counts.t) - ncol(total_counts), "empty ASVs out of", nrow(total_counts.t), "\n")
# 2218 empty ASVs out of 14958 



#############
### converted

# transpose
plate1_convert.t <- data.frame(t(plate1_convert))
plate2_convert.t <- data.frame(t(plate2_convert))


# add fasta variable
plate1_convert.t$FASTA <- rownames(plate1_convert.t)
plate2_convert.t$FASTA <- rownames(plate2_convert.t)


# how many are there?
length(unique(c(plate2_convert.t$FASTA, plate1_convert.t$FASTA))) # [1] 14958


# merge
total_convert.t <- merge(plate1_convert.t, plate2_convert.t, by = "FASTA", all = TRUE)
length(unique(total_convert.t$FASTA)) # [1] 14958


# change NA to zero
total_convert.t[is.na(total_convert.t)] <- 0


# remove FASTA column
rownames(total_convert.t) <- total_convert.t$FASTA
total_convert.t <- total_convert.t[-1]


# order by most abundant ASVs
total_convert.t <- total_convert.t[order(rowSums(total_convert.t), decreasing = TRUE),]


# bring stds to top
total_stds.t <- total_convert.t[grep("STD", rownames(total_convert.t)),]
total_convert.t <- total_convert.t[! rownames(total_convert.t) %in% rownames(total_stds.t),]
total_convert.t <- rbind(total_stds.t, total_convert.t)


# transpose back
total_convert <- data.frame(t(total_convert.t))
total_mini2 <- total_convert[,1:100]


# remove empty ASVs
total_convert <- total_convert[rowSums(total_convert) > 0, colSums(total_convert) > 0]
cat(nrow(total_convert.t) - ncol(total_convert), "empty ASVs out of", nrow(total_convert.t), "\n")
# 2218 empty ASVs out of 14958 




########################
### assess standards ###
########################

# subset standards
stds_convert <- total_convert[1:4]
stds_counts <- total_counts[1:4]


# melt for plotting
stds_convert.m <- melt(stds_convert, variable.name = "Standard", value.name = "Abundance")
stds_counts.m <- melt(stds_counts, variable.name = "Standard", value.name = "Abundance")


# change zeros to ones, so can see them on plot x axis
stds_convert.m$Abundance[stds_convert.m$Abundance == 0] <- 1
stds_counts.m$Abundance[stds_counts.m$Abundance == 0] <- 1


# add data points from expected std copies
expect_data <- data.frame(STD = c("STD1", "STD2", "STD3", "STD4"),
                          Abundance = c(300000, 40000, 7000, 1000),
                          Label = c("300,000", "40,000", "7,000", "1,000"))


# plot converted data
convert.plot <-
  ggplot(stds_convert.m, aes(x = Standard, y = log10(Abundance))) +
  geom_jitter(size = 1, width = 0.1, alpha = 0.5, shape = 1) +
  geom_point(data = expect_data, aes(x = STD, y = log10(Abundance)), color = "red") +
  geom_text(data = expect_data, aes(x = STD, y = log10(Abundance), label = Label), 
            color = "red", hjust = -0.7, size = 2.5) +
  scale_y_continuous(breaks = c(1,2,3,4,5),
                     labels = c("10", "100", "1,000", "10,000", "100,000")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.title.x = element_text(size = 8, color = "black", face = "bold"),
        legend.title = element_text(size = 8, color = "black", face = "bold"),
        legend.text = element_text(size = 7, color = "black", face = "italic"),
        plot.title = element_text(size = 9, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 7, color = "red"),
        panel.border = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.25)) +
  labs(y = "Absolute abundance (copies/µL)", title = "Standard reads after conversion",
       subtitle = ". Standard copies/µL spiked into PCR")
convert.plot

#ggsave("./Plots/convertjitter.pdf", convert.plot, width = 4, height = 3, device = "pdf")


# plot count data
counts.plot <-
  ggplot(stds_counts.m, aes(x = Standard, y = log10(Abundance))) +
  geom_jitter(size = 1, width = 0.1, alpha = 0.5, shape = 1) +
  geom_point(data = expect_data, aes(x = STD, y = log10(Abundance)), color = "red") +
  geom_text(data = expect_data, aes(x = STD, y = log10(Abundance), label = Label), 
            color = "red", hjust = -0.7, size = 2.5) +
  scale_y_continuous(breaks = c(1,2,3,4,5),
                     labels = c("10", "100", "1,000", "10,000", "100,000")) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.y = element_text(size = 8, color = "black", face = "bold"),
        axis.title.x = element_text(size = 8, color = "black", face = "bold"),
        legend.title = element_text(size = 8, color = "black", face = "bold"),
        legend.text = element_text(size = 7, color = "black", face = "italic"),
        plot.title = element_text(size = 9, color = "black", face = "bold"),
        plot.subtitle = element_text(size = 7, color = "red"),
        panel.border = element_blank(), 
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.25)) +
  labs(y = "Illumina read count", title = "Standard reads before conversion",
       subtitle = ". Standard copies/µL spiked into PCR")
counts.plot

#ggsave("./Plots/countjitter.pdf", counts.plot, width = 4, height = 3, device = "pdf")


# # save
# save.image("./RData/conversions_env.RData")


