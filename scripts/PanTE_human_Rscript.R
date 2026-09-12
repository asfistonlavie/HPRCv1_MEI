#### Characterization of Transposable Elements Variations in the Human Pangenome - v2026
#### Shadi Shahatit, Master's Thesis, ISEM, UM 2023-2024
#### Last Update: 12 September 2026
# Libraries ---------------------------------------------------------------

library("org.Hs.eg.db")
library("clusterProfiler")
library("ggforce")
library("fuzzyjoin")
library("qqman")
library("dplyr")
library("readxl")
library("ggplot2")
library("ggpattern")
library("diveRsity")
library("rehh")
library("vcfR")
library("rehh.data")
library("writexl")
library("data.table")
library("RColorBrewer")
library("tidyverse")
library("forcats")
library("wesanderson")
library("stringr")
library("tidyr")
library("readr")
library("splitstackshape")
library("ggridges")
library("karyoploteR")
library("BiocManager")
library("GenomicRanges")
library("cowplot")
library("scales")
library("palmerpenguins")
library("ggbeeswarm")
library("systemfonts")
library("ggforce")
library("ggpubr")
library("rtracklayer")
library("ggcorrplot")
library("scMethrix")
library("BRGenomics")
library("minpack.lm")
library("sjPlot")
library("DHARMa")
library("RRPP")
library("devtools")
library("lme4")
library("betareg")
library("nnet")
library("MASS")
library("ordinal")
library("brant")
library("ggeffects")
library("lvplot")
library("ggthemes")
library("caret")
library("patchwork")
library("ggrepel")
library("MuMIn")
library("lmtest")
library("fitdistrplus")
library("VennDiagram")
library("viridis")
library("grid")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library("GenometriCorr")
library("rGREAT")
library("regioneR")
library("plyranges")
library("flexplot")

## note: replace your system's directory in sys_dir
sys_dir <- "/home/shadi/Desktop/S3_Project"

# Figure 1 ----------------------------------------------------------------



# Define AFR, outofAFR, global allele freq, and TE annotation ----------------------------

## outofAFR freq

outofAFR_freq <- read.table(file = file.path(sys_dir,"PanTE_human/frequency_distribution/outofAFR/freqoutofAFR.frq"), header = T,
                            sep = " ",
                            dec = ".",
                            fill = T) 
colnames(outofAFR_freq)= c("V1")
outofAFR_freq$V1 <- gsub(":", "\t", outofAFR_freq$V1)
outofAFR_freq_2 <- str_split(outofAFR_freq$V1, "\t", simplify=TRUE) %>% data.frame(.)
outofAFR_freq_3 <- outofAFR_freq_2 %>% dplyr::select(c(1,2,5,6,7,8))
colnames(outofAFR_freq_3) = c("CHR","POS","REF.o","FREQ_REF.o","ALT.o","TE_FREQ.o")
outofAFR_freq_3[,6] <- as.numeric(unlist(outofAFR_freq_3[,6]))
outofAFR_freq_3[,4] <- as.numeric(unlist(outofAFR_freq_3[,4]))
outofAFR_freq_3 <- data.frame(append(outofAFR_freq_3, c(superpop.o="outofAFR"), after=6))

## AFR freq

AFR_freq <- read.table(file = file.path(sys_dir,"PanTE_human/frequency_distribution/AFR/freqAFR.frq"), header = T,
                       sep = " ",
                       dec = ".",
                       fill = T) 
colnames(AFR_freq)= c("V1")
AFR_freq$V1 <- gsub(":", "\t", AFR_freq$V1)
AFR_freq_2 <- str_split(AFR_freq$V1, "\t", simplify=TRUE) %>% data.frame(.)
AFR_freq_3 <- AFR_freq_2 %>% dplyr::select(c(1,2,5,6,7,8))
colnames(AFR_freq_3) = c("CHR","POS","REF","FREQ_REF","ALT","FREQ_ALT")
AFR_freq_3[,6] <- as.numeric(unlist(AFR_freq_3[,6]))
AFR_freq_3[,4] <- as.numeric(unlist(AFR_freq_3[,4]))
AFR_freq_3 <- data.frame(append(AFR_freq_3, c(superpop="AFR"), after=6))

## global freq

super_freq <- read.table(file = file.path(sys_dir,"PanTE_human/frequency_distribution/superpop/freqSuper.frq"), header = T,
                         sep = " ",
                         dec = ".",
                         fill = T) 

colnames(super_freq)= c("V1")

super_freq$V1 <- gsub(":", "\t", super_freq$V1)

super_freq_2 <- str_split(super_freq$V1, "\t", simplify=TRUE) %>% data.frame(.)

super_freq_3 <- super_freq_2 %>% dplyr::select(c(1,2,5,6,7,8))
colnames(super_freq_3) = c("CHR","POS","REF","FREQ_REF","ALT","FREQ_ALT")

super_freq_3[,6] <- as.numeric(unlist(super_freq_3[,6]))
super_freq_3[,4] <- as.numeric(unlist(super_freq_3[,4]))

super_freq_3[,2] <- as.numeric(unlist(super_freq_3[,2]))
super_freq_7 <- super_freq_3 %>% mutate(end=POS+nchar(ALT)-1)
super_freq_7 <- super_freq_7 %>% dplyr::select(c(1,2,6,7))

## TE_anno

TE_anno <- read.table(file = file.path(sys_dir,"PanTE_human/frequency_distribution/RTE_norm_mm90_1alt.vcf"), header = F,
                      sep = " ",
                      dec = ".")
colnames(TE_anno) = c("CHR","POS","svid","ALT","ANNO")
TE_anno$ANNO <- gsub(";", "", TE_anno$ANNO)
TE_anno$ANNO <- gsub("=", "", TE_anno$ANNO)
TE_anno$ANNO <- gsub("rC", "\t", TE_anno$ANNO)
TE_anno[,6:7] <- stringr::str_split_fixed(TE_anno$ANNO, "\t", 2)
TE_anno <- TE_anno %>% dplyr::select(c(1,2,4,7))
TE_anno[,2] <- as.numeric(unlist(TE_anno[,2]))
TE_anno <- TE_anno[order(TE_anno[,1],TE_anno[,2]),]
colnames(TE_anno) = c("CHR","POS","ALT","ANNO")
TE_anno <- TE_anno %>% mutate(END=POS+nchar(ALT)-1)
super_freq_anno <- merge(super_freq_3,TE_anno)

TE_anno_Gr <- makeGRangesFromDataFrame(df=TE_anno,
                                       keep.extra.columns=T,
                                       ignore.strand=T,
                                       seqnames.field="CHR",
                                       start.field="POS",
                                       end.field="end",
                                       starts.in.df.are.0based=FALSE)

# Freq & size assignment - v2026 ------------------------------------------

super_freq_3
super_freq_anno

nrow(super_freq_3)
nrow(super_freq_anno)
colnames(super_freq_3)
colnames(super_freq_anno)

super_freq_anno$ALT_size <- nchar(super_freq_anno$ALT)
super_freq_anno$REF_size <- nchar(super_freq_anno$REF)

nrow(super_freq_anno[super_freq_anno$ALT_size < 100,])
nrow(super_freq_anno[super_freq_anno$REF_size > 100,])

super_freq_anno_v2026 <- super_freq_anno %>%
  mutate(
    TE_allele = case_when(
      ALT_size >= 100 & REF_size <  100 ~ "ALT_TE",        
      ALT_size <  100 & REF_size >= 100 ~ "REF_TE",        
      ALT_size >= 100 & REF_size >= 100 ~ "Questioning_both_large",
      ALT_size <  100 & REF_size <  100 ~ "Questioning_both_small", # should be 0
      TRUE ~ NA_character_),
    
    TE_size = case_when(
      TE_allele == "ALT_TE" | TE_allele == "Questioning_both_large" ~ ALT_size,
      TE_allele == "REF_TE" ~ REF_size,
      TRUE ~ NA_real_),
    
    TE_FREQ = case_when(
      TE_allele == "ALT_TE" | TE_allele == "Questioning_both_large" ~ FREQ_ALT,
      TE_allele == "REF_TE" ~ FREQ_REF,
      TRUE ~ NA_real_)) %>%
  mutate(
    new_class = case_when(
      TE_FREQ >= 0.95                   ~ "Fixed",
      TE_FREQ >= 0.05 & TE_FREQ < 0.95  ~ "Polymorphic",
      TE_FREQ < 0.05                    ~ "Rare",
      TRUE ~ NA_character_))

nrow(super_freq_anno_v2026)
table(super_freq_anno_v2026$TE_allele)

table(super_freq_anno_v2026[super_freq_anno_v2026$TE_allele == "ALT_TE",]$REF_size)
table(super_freq_anno_v2026[super_freq_anno_v2026$TE_allele == "REF_TE",]$ALT_size)
table(super_freq_anno_v2026[super_freq_anno_v2026$TE_allele == "Questioning_both_large",]$REF_size)
table(super_freq_anno_v2026[super_freq_anno_v2026$TE_allele == "Questioning_both_large",]$ALT_size)

table(super_freq_anno_v2026$new_class)

nrow(super_freq_anno_v2026[super_freq_anno_v2026$TE_allele == "ALT_TE" &
                             # super_freq_anno_v2026$REF_size < 100 &
                             super_freq_anno_v2026$new_class == "Fixed",])

## for AFR_freq_3

AFR_freq_3$ALT_size <- nchar(AFR_freq_3$ALT)
AFR_freq_3$REF_size <- nchar(AFR_freq_3$REF)

nrow(AFR_freq_3[AFR_freq_3$ALT_size < 100,])
nrow(AFR_freq_3[AFR_freq_3$REF_size > 100,])

AFR_freq_3_v2026 <- AFR_freq_3 %>%
  mutate(
    TE_allele = case_when(
      ALT_size >= 100 & REF_size <  100 ~ "ALT_TE",        
      ALT_size <  100 & REF_size >= 100 ~ "REF_TE",        
      ALT_size >= 100 & REF_size >= 100 ~ "Questioning_both_large",
      ALT_size <  100 & REF_size <  100 ~ "Questioning_both_small", # should be 0
      TRUE ~ NA_character_),
    
    TE_size = case_when(
      TE_allele == "ALT_TE" | TE_allele == "Questioning_both_large" ~ ALT_size,
      TE_allele == "REF_TE" ~ REF_size,
      TRUE ~ NA_real_),
    
    TE_FREQ = case_when(
      TE_allele == "ALT_TE" | TE_allele == "Questioning_both_large" ~ FREQ_ALT,
      TE_allele == "REF_TE" ~ FREQ_REF,
      TRUE ~ NA_real_)) %>%
  mutate(
    new_class = case_when(
      TE_FREQ >= 0.95                   ~ "Fixed",
      TE_FREQ >= 0.05 & TE_FREQ < 0.95  ~ "Polymorphic",
      TE_FREQ < 0.05                    ~ "Rare",
      TRUE ~ NA_character_))

nrow(AFR_freq_3_v2026)
table(AFR_freq_3_v2026$TE_allele) 
table(AFR_freq_3_v2026$new_class) - table(super_freq_anno_v2026$new_class)
plot(density(AFR_freq_3_v2026$TE_size))
plot(density(super_freq_anno_v2026$TE_size))
plot(density(AFR_freq_3_v2026$TE_FREQ))
plot(density(super_freq_anno_v2026$TE_FREQ))

identical(AFR_freq_3_v2026$TE_allele, super_freq_anno_v2026$TE_allele)

## for outofAFR_freq_3

outofAFR_freq_3$ALT_size <- nchar(outofAFR_freq_3$ALT)
outofAFR_freq_3$REF_size <- nchar(outofAFR_freq_3$REF)

nrow(outofAFR_freq_3[outofAFR_freq_3$ALT_size < 100,])
nrow(outofAFR_freq_3[outofAFR_freq_3$REF_size > 100,])

outofAFR_freq_3_v2026 <- outofAFR_freq_3 %>%
  mutate(
    TE_allele = case_when(
      ALT_size >= 100 & REF_size <  100 ~ "ALT_TE",        
      ALT_size <  100 & REF_size >= 100 ~ "REF_TE",        
      ALT_size >= 100 & REF_size >= 100 ~ "Questioning_both_large",
      ALT_size <  100 & REF_size <  100 ~ "Questioning_both_small", # should be 0
      TRUE ~ NA_character_),
    
    TE_size = case_when(
      TE_allele == "ALT_TE" | TE_allele == "Questioning_both_large" ~ ALT_size,
      TE_allele == "REF_TE" ~ REF_size,
      TRUE ~ NA_real_),
    
    TE_FREQ = case_when(
      TE_allele == "ALT_TE" | TE_allele == "Questioning_both_large" ~ FREQ_ALT.o,
      TE_allele == "REF_TE" ~ FREQ_REF.o,
      TRUE ~ NA_real_)) %>%
  mutate(
    new_class = case_when(
      TE_FREQ >= 0.95                   ~ "Fixed",
      TE_FREQ >= 0.05 & TE_FREQ < 0.95  ~ "Polymorphic",
      TE_FREQ < 0.05                    ~ "Rare",
      TRUE ~ NA_character_))

nrow(outofAFR_freq_3_v2026)
table(outofAFR_freq_3_v2026$TE_allele) 
table(outofAFR_freq_3_v2026$new_class) - table(super_freq_anno_v2026$new_class)
plot(density(outofAFR_freq_3_v2026$TE_size))
plot(density(super_freq_anno_v2026$TE_size))
plot(density(outofAFR_freq_3_v2026$TE_FREQ))
plot(density(super_freq_anno_v2026$TE_FREQ))

identical(outofAFR_freq_3_v2026$TE_allele, super_freq_anno_v2026$TE_allele)

# TE size and Karyoplote ---------------------------------------------------

## TE size

chr_length <- read.table(file = file.path(sys_dir,"PanTE_human/karyoplote_and_TE_size/chr_lenght_hg28"), header = FALSE,
                         sep = "",
                         dec = ".")
colnames(chr_length) <- c("CHR","chr_length")
chr_length$CHR <- gsub("chr01", "chr1", chr_length$CHR)
chr_length$CHR <- gsub("chr02", "chr2", chr_length$CHR)
chr_length$CHR <- gsub("chr03", "chr3", chr_length$CHR)
chr_length$CHR <- gsub("chr04", "chr4", chr_length$CHR)
chr_length$CHR <- gsub("chr05", "chr5", chr_length$CHR)
chr_length$CHR <- gsub("chr06", "chr6", chr_length$CHR)
chr_length$CHR <- gsub("chr07", "chr7", chr_length$CHR)
chr_length$CHR <- gsub("chr08", "chr8", chr_length$CHR)
chr_length$CHR <- gsub("chr09", "chr9", chr_length$CHR)
# super_freq_3$TE_size <- nchar(super_freq_3$ALT)
# super_freq_3_chr_length <- merge(super_freq_3,chr_length) 
super_freq_anno_v2026_chr_length <- merge(super_freq_anno_v2026,chr_length) 

size_length <- super_freq_anno_v2026_chr_length %>% dplyr::select(c("CHR","POS","chr_length","TE_size"))

size_length$TE_percen <- (size_length$TE_size/size_length$chr_length)*100
(sum(size_length$TE_percen))

size_length <- size_length %>% filter(TE_size >= 100)
nrow(size_length)

TE_size_sum <- sum(size_length$TE_size)
chr_length_sum <- sum(unique(size_length$chr_length))
TE_genome_percentage <- ((TE_size_sum)/(chr_length_sum))*100
median(size_length$TE_size)
mean(size_length$TE_size)

size_length_counted <- size_length %>% dplyr::count(CHR,chr_length)
cor <- cor(size_length_counted$chr_length,size_length_counted$n, method = "pearson")
# cor <- cor(size_length_counted$chr_length,size_length_counted$n, method = "spearman")
qqnorm(size_length_counted$chr_length)
qqline(size_length_counted$chr_length, col = "red")
qqnorm(size_length_counted$n)
qqline(size_length_counted$n, col = "red")
shapiro.test(size_length_counted$chr_length)
shapiro.test(size_length_counted$n)

figsupp_TEcountchrlength <- ggplot(data=size_length_counted, mapping=aes(x=chr_length, y=n)) +
  geom_point(col="steelblue") +
  geom_smooth(method="lm", col="#EE353E") +
  labs(x="Chromosome length", y="TE count")+
  theme_classic() +
  annotate("text", x=min(size_length_counted$chr_length), y=max(size_length_counted$n), 
           label=sprintf("R2: %.2f", cor), hjust=0, vjust=1, size=5, color="black")

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_TEcountchrlength.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_TEcountchrlength.png"),
       plot = figsupp_TEcountchrlength, width = 6, height = 4, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

TE_size_chr <- size_length %>%
  group_by(CHR) %>%
  summarise(total_TE_size = sum(TE_size),
            mean_chr_length = mean(chr_length)) %>% as.data.frame()
cor <- cor(TE_size_chr$mean_chr_length,TE_size_chr$total_TE_size, method = "pearson")
# cor <- cor(TE_size_chr$mean_chr_length,TE_size_chr$total_TE_size, method = "spearman")
qqnorm(TE_size_chr$mean_chr_length)
qqline(TE_size_chr$mean_chr_length, col = "red")
qqnorm(TE_size_chr$total_TE_size)
qqline(TE_size_chr$total_TE_size, col = "red")
shapiro.test(TE_size_chr$mean_chr_length)
shapiro.test(TE_size_chr$total_TE_size)

figsupp_TEsizechrlength <- ggplot(data=TE_size_chr, mapping=aes(x=mean_chr_length, y=total_TE_size)) +
  geom_point(col="steelblue") +
  geom_smooth(method="lm", col="#EE353E") +
  labs(x="Chromosome length", y="TE size")+
  theme_classic() +
  annotate("text", x=min(TE_size_chr$mean_chr_length), y=max(TE_size_chr$total_TE_size), 
           label=sprintf("R2: %.2f", cor), hjust=0, vjust=1, size=5, color="black")

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_TEsizechrlength.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_TEsizechrlength.png"),
       plot = figsupp_TEsizechrlength, width = 6, height = 4, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

TE_percen_chr <- size_length %>%
  group_by(CHR) %>%
  summarise(total_TE_percen = sum(TE_percen),
            mean_chr_length = mean(chr_length)) %>% as.data.frame()
cor <- cor(TE_percen_chr$mean_chr_length,TE_percen_chr$total_TE_percen)

cor <- cor(TE_percen_chr$mean_chr_length,TE_percen_chr$total_TE_percen, method = "pearson")
cor <- cor(TE_percen_chr$mean_chr_length,TE_percen_chr$total_TE_percen, method = "spearman")
qqnorm(TE_percen_chr$mean_chr_length)
qqline(TE_percen_chr$mean_chr_length, col = "red")
qqnorm(TE_percen_chr$total_TE_percen)
qqline(TE_percen_chr$total_TE_percen, col = "red")
shapiro.test(TE_percen_chr$mean_chr_length)
shapiro.test(TE_percen_chr$total_TE_percen)

ggplot(data=TE_percen_chr, mapping=aes(x=mean_chr_length, y=total_TE_percen)) +
  geom_point(col="steelblue") +
  geom_smooth(method="lm", col="#EE353E") +
  labs(x="Chromosome length", y="TE percentage")+
  theme_classic() +
  annotate("text", x=min(TE_percen_chr$mean_chr_length), y=max(TE_percen_chr$total_TE_percen), 
           label=sprintf("R2: %.2f", cor), hjust=0, vjust=1, size=5, color="black")

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_TEsizepercenchrlength.png"), plot=ggplot2::last_plot())

kruskal_test <- kruskal.test(total_TE_percen ~ CHR, data = TE_percen_chr)
print(kruskal_test)

# size_length$CHR <- gsub("chr1", "chr01", size_length$CHR)
# size_length$CHR <- gsub("chr2", "chr02", size_length$CHR)
# size_length$CHR <- gsub("chr3", "chr03", size_length$CHR)
# size_length$CHR <- gsub("chr4", "chr04", size_length$CHR)
# size_length$CHR <- gsub("chr5", "chr05", size_length$CHR)
# size_length$CHR <- gsub("chr6", "chr06", size_length$CHR)
# size_length$CHR <- gsub("chr7", "chr07", size_length$CHR)
# size_length$CHR <- gsub("chr8", "chr08", size_length$CHR)
# size_length$CHR <- gsub("chr9", "chr09", size_length$CHR)

# hist(size_length$TE_size,
#      col = "skyblue",
#      border = "black",
#      main = "Histogram of Size",
#      xlab = "Size",
#      ylab = "Frequency",
#      xlim = c(min(size_length$TE_size), max(size_length$TE_size)),
#      breaks = 15)

TE_size_anno <- merge(size_length,TE_anno)

figsupp_TEsize <- ggplot(size_length)+
  geom_histogram(binwidth = 100, color = "black", alpha = 0.7 ,aes(x = TE_size))+
  geom_vline(aes(xintercept=median(TE_size)),color="#EE353E", linetype="dashed", size=1)+
  labs(
    title = "TE Size Distribution",
    x = "Size (bp)",
    y = "TE Count")+
  scale_color_grey() +
  theme_classic()

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_TEsize.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_TEsize.png"),
       plot = figsupp_TEsize, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

L1_size <- ggplot(data=subset(TE_size_anno, ANNO == "LINE/L1"))+
  geom_histogram(binwidth = 100, color = "black", alpha = 0.7 ,aes(x = TE_size))+
  geom_vline(aes(xintercept=median(TE_size)),color="#EE353E", linetype="dashed", size=1)+
  labs(
    title = "L1 Size Distribution",
    x = "Size (bp)",
    y = "TE Count")+
  scale_color_grey()+
  theme_classic()
Alu_size <- ggplot(data=subset(TE_size_anno, ANNO == "SINE/Alu"))+
  geom_histogram(binwidth = 100, color = "black", alpha = 0.7 ,aes(x = TE_size))+
  geom_vline(aes(xintercept=median(TE_size)),color="#EE353E", linetype="dashed", size=1)+
  labs(
    title = "Alu Size Distribution",
    x = "Size (bp)",
    y = "TE Count")+
  scale_color_grey()+
  theme_classic()
SVA_size <- ggplot(data=subset(TE_size_anno, ANNO == "Retroposon/SVA"))+
  geom_histogram(binwidth = 100, color = "black", alpha = 0.7 ,aes(x = TE_size))+
  geom_vline(aes(xintercept=median(TE_size)),color="#EE353E", linetype="dashed", size=1)+
  labs(
    title = "SVA Size Distribution",
    x = "Size (bp)",
    y = "TE Count")+
  scale_color_grey()+
  theme_classic()
ERV_size <- ggplot(data=subset(TE_size_anno, ANNO == "LTR/ERV"))+
  geom_histogram(binwidth = 100, color = "black", alpha = 0.7 ,aes(x = TE_size))+
  geom_vline(aes(xintercept=median(TE_size)),color="#EE353E", linetype="dashed", size=1)+
  labs(
    title = "ERV Size Distribution",
    x = "Size (bp)",
    y = "TE Count")+
  scale_color_grey()+
  theme_classic()
figsupp_TEsize_perTE <- L1_size+Alu_size+SVA_size+ERV_size

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_TEsize_perTE.png"),
       plot = figsupp_TEsize_perTE, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

## Karyoplote

super_freq_7 <- super_freq_7[order(super_freq_7[,1],super_freq_7[,2]),]

super_freq_Gr <- makeGRangesFromDataFrame(df=super_freq_7,
                                          keep.extra.columns=T,
                                          ignore.strand=T,
                                          seqnames.field="CHR",
                                          start.field="POS",
                                          end.field="end",
                                          starts.in.df.are.0based=FALSE)

kpden <- plotKaryotype(genome="hg38", plot.type = 2)
kpPlotDensity(kpden, data=super_freq_Gr,col="goldenrod")
kpAddMainTitle(kpden,"TE density distribution across the chromosomes")

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_karyotype.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_karyotype.png"),
       plot = figsupp_karyotype, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

kpli <- plotKaryotype(genome="hg38", plot.type = 2)
kpLines(kpli, data=super_freq_Gr, y=super_freq_Gr$TE_FREQ,data.panel=1, col="darkgoldenrod",r0=0.25, r1=1)
kpAxis(kpli, ymax=1, r0=0.25, r1=1,numticks = 2, col="#666666", cex=0.5, text.col="black")
kpAddMainTitle(kpli,"Global Allele freqeuncy")

## check the size of REF and ALT alleles - HG review

shadyTE <- (super_freq_3_chr_length[(super_freq_3_chr_length$TE_size < 100) & (super_freq_3_chr_length$REF_length > 99),])$POS

TE_freq_anno_selscan_master <- TE_freq_anno_selscan_master %>%
  mutate(shadyTE = case_when(
    POS %in% shadyTE ~ T,
    TRUE ~ F))

table(TE_freq_anno_selscan_master$shadyTE)

TE_freq_anno_selscan_master_shadyTE <- TE_freq_anno_selscan_master %>%
  filter(shadyTE == T)

table(TE_freq_anno_selscan_master_shadyTE$ANNO)
table(TE_freq_anno_selscan_master_shadyTE$LOC)
table(TE_freq_anno_selscan_master_shadyTE$selscans)
table(TE_freq_anno_selscan_master_shadyTE$new_class)

# TE counts per chr -------------------------------------------------------

# super_freq_anno_v2026_chr_length instead of super_freq_3_chr_length

biTE_count <- super_freq_anno_v2026_chr_length %>% dplyr::select(c("CHR","POS","TE_size","chr_length"))
biTE_count <- merge(biTE_count,TE_anno)
biTE_count <- biTE_count[order(biTE_count[,1],biTE_count[,2]),]

# no duplicates are present
# dup_indices <- duplicated(biTE_count$POS) &
#   duplicated(biTE_count$POS, fromLast = TRUE) &
#   duplicated(biTE_count$CHR[biTE_count$CHR == "chrX"])
# biTE_count <- biTE_count[!dup_indices, ]

biTE_count_2 <- biTE_count %>%
  dplyr::group_by(CHR, ANNO) %>%
  dplyr::summarize(count = dplyr::n(),
            mean_TE_size = mean(TE_size),
            mean_chr_length = mean(chr_length)) %>% as.data.frame()

biTE_count_2$norm_count <- biTE_count_2$count/biTE_count_2$mean_chr_length

biTE_count_2[biTE_count_2 == "chr1"] <- "chr01"
biTE_count_2[biTE_count_2 == "chr2"] <- "chr02"
biTE_count_2[biTE_count_2 == "chr3"] <- "chr03"
biTE_count_2[biTE_count_2 == "chr4"] <- "chr04"
biTE_count_2[biTE_count_2 == "chr5"] <- "chr05"
biTE_count_2[biTE_count_2 == "chr6"] <- "chr06"
biTE_count_2[biTE_count_2 == "chr7"] <- "chr07"
biTE_count_2[biTE_count_2 == "chr8"] <- "chr08"
biTE_count_2[biTE_count_2 == "chr9"] <- "chr09"

anno_order <- fct_relevel(biTE_count_2$ANNO,"DNA/misc","Retroposon/SVA","LTR/misc","LTR/ERV",
                          "LINE/misc","LINE/L1","SINE/MIR","SINE/Alu")
new_colors <- c("SINE/Alu"="#440154FF", "SINE/MIR"="#453781FF",
                "LINE/L1"="#39558CFF", "LINE/misc"="#238A8DFF",
                "LTR/ERV"="#29AF7FFF","LTR/misc"="#74D055FF",
                "Retroposon/SVA"= "#B8DE29FF","DNA/misc"="#FDE725FF")

figsupp_TEcountperchr <- ggplot(biTE_count_2, aes(x = (CHR), y=(norm_count), fill = anno_order)) + 
  geom_bar(position="fill", stat="identity")+
  theme_minimal()+
  xlab("Chromosomes")+
  ylab("Normalized TE counts")+
  labs(fill ="TE type")+
  scale_fill_manual(values=c(new_colors))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_TEcountperchr.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_TEcountperchr.png"),
       plot = figsupp_TEcountperchr, width = 14, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

contingency_table <- table(biTE_count$CHR)
c_test_result <- chisq.test(contingency_table)
print(c_test_result)
chisq.test(contingency_table)$expected

biTE_count_3 <- aggregate(cbind(count, norm_count) ~ CHR, data = biTE_count_2, sum)

kruskal.test(count ~ CHR, data = biTE_count_3)
kruskal.test(norm_count ~ CHR, data = biTE_count_3)

# Allele freq distribution and pie charts ----------------------------------

## define freq ranges - v1

# super_freq_anno$new_class <- lapply(super_freq_anno$FREQ_ALT, function(x) if(x>0.901){
#   super_freq_anno$new_class="Fixed"
# }else if (0.634 < x && x <= 0.901){
#   super_freq_anno$new_class="Major"
# }else if (0.367 < x && x <= 0.634){
#   super_freq_anno$new_class="Common"
# }else if ((0.10) < x && x <= 0.367){
#   super_freq_anno$new_class="Rare"
# }else 
#   super_freq_anno$new_class="Very Rare")
# super_freq_anno$new_class <- as.character(unlist(super_freq_anno$new_class))
# new_class_table <- table(super_freq_anno$new_class) %>% as.data.frame()
# colnames(new_class_table) <- c("new_class","count")
# new_class_table$count_sum <- sum(new_class_table$count)
# new_class_table$freq <- (new_class_table$count/new_class_table$count_sum)
# 
# super_freq_anno$ANNO <- as.character(unlist(super_freq_anno$ANNO))
# super_freq_anno_counted <- super_freq_anno %>% dplyr::count(new_class,ANNO)
# colnames(super_freq_anno_counted) = c("CLASS","ANNO","NUMBER")
# super_freq_anno_counted$CLASS <- as.factor(unlist(super_freq_anno_counted$CLASS))
# 
# ## freq distribution
# 
# anno_order <- fct_relevel(super_freq_anno_counted$ANNO,"SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")
# color_class <- c("#fc4e2a","#fd8d3c","#feb24c", "#fed976", "#fff3d1")
# new_color_class <- c("#fff3d1","#fed976","#feb24c", "#fd8d3c", "#fc4e2a")
# new_class_order <- fct_relevel(super_freq_anno$new_class,"Fixed","Major","Common","Rare","Very Rare")
# super_freq_anno_counted$ANNO <- reorder(super_freq_anno_counted$ANNO, -super_freq_anno_counted$NUMBER)
# 
# ggplot(super_freq_anno,mapping=aes(fill = new_class_order))+
#   geom_histogram(binwidth = 0.07,aes(x = FREQ_ALT),color="black")+
#   labs(
#     # title = "Allele frequency distribution",
#     x = "TE frequency",
#     y = "TE count",
#     fill="Frequency ranges")+
#   scale_fill_manual(values = color_class) +
#   theme_classic()+
#   theme(
#     plot.title = element_text(size = 16, face = "bold"),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     axis.title.y = element_text(size = 14, face = "bold"))
# 
# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_freqspec.png"), plot=ggplot2::last_plot())
# 
# p <- ggplot(super_freq_anno_counted, mapping = aes(x=ANNO,
#                                                    y=NUMBER,
#                                                    fill=anno_order))+
#   geom_bar(position="stack", stat="identity")+
#   scale_fill_viridis(discrete = TRUE) +
#   scale_fill_manual(values=c("#440154FF","#453781FF","#39558CFF","#2D718EFF","#29AF7FFF","#56C667FF","#B8DE29FF","#FDE725FF"))+
#   # scale_color_manual(values = COLS,breaks=legend_order)+
#   labs(
#     x = "TE type",
#     y = "TE count",
#     fill=("TE type")
#   ) +
#   theme_classic()+
#   theme(plot.title = element_text( size = 16, face = "bold", hjust = 0.5))+
#   scale_y_continuous(breaks = seq(0,max(super_freq_anno_counted$NUMBER), by = 500))+
#   # theme(axis.text.y=element_blank())+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#   facet_grid(~factor(CLASS,levels=c("Very Rare","Rare","Common","Major","Fixed")))+
#   theme(
#     plot.title = element_text(size = 16, face = "bold"),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     axis.title.y = element_text(size = 14, face = "bold"))
# g <- ggplot_gtable(ggplot_build(p))
# striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) )
# fills <- new_color_class
# k <- 1
# for (i in striprt) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid.draw(g)
# 
# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_freqdist.png"), plot=ggplot2::last_plot())

## define freq ranges - v2 - already done in "Freq & size assignment - v2026" 

## use super_freq_anno_v2026

# super_freq_anno$new_class <- lapply(super_freq_anno$FREQ_ALT, function(x) if(x>=0.95){
#   super_freq_anno$new_class="Fixed"
# }else if (0.05 <= x && x < 0.95){
#   super_freq_anno$new_class="Polymorphic"
# }else 
#   super_freq_anno$new_class="Rare")
# super_freq_anno$new_class <- as.character(unlist(super_freq_anno$new_class))
new_class_table <- table(super_freq_anno_v2026$new_class) %>% as.data.frame()
colnames(new_class_table) <- c("new_class","count")
new_class_table$count_sum <- sum(new_class_table$count)
new_class_table$freq <- (new_class_table$count/new_class_table$count_sum)
# new_class count count_sum       freq
# 1       Fixed   234      6407 0.03652255
# 2 Polymorphic   857      6407 0.13375995
# 3        Rare  5316      6407 0.82971750

table(super_freq_anno_v2026$new_class,super_freq_anno_v2026$ANNO)

super_freq_anno_v2026 %>%
  count(new_class, ANNO) %>%
  arrange(new_class, n)

## Fixed but not in REF? lets highlight

super_freq_anno_v2026 %>%
  filter(new_class == "Fixed") %>%
  ggplot(aes(x = POS, y = TE_FREQ, col = TE_allele)) +
  geom_point() +
  theme_classic() +
  labs(x = "Position", y = "TE frequency")

nrow(super_freq_anno_v2026[super_freq_anno_v2026$TE_allele == "ALT_TE" &
                             super_freq_anno_v2026$REF_size < 100 &
                             super_freq_anno_v2026$new_class == "Fixed",])
table(super_freq_anno_v2026[super_freq_anno_v2026$TE_allele == "ALT_TE" &
                             super_freq_anno_v2026$REF_size < 100 &
                             super_freq_anno_v2026$new_class == "Fixed",]$ANNO)

super_freq_anno_v2026$ANNO <- as.character(unlist(super_freq_anno_v2026$ANNO))
super_freq_anno_v2026_counted <- super_freq_anno_v2026 %>% dplyr::count(new_class,ANNO)
colnames(super_freq_anno_v2026_counted) = c("CLASS","ANNO","NUMBER")
super_freq_anno_v2026_counted$CLASS <- as.factor(unlist(super_freq_anno_v2026_counted$CLASS))

## freq distribution

anno_order <- fct_relevel(super_freq_anno_v2026_counted$ANNO,"SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")
color_class <- c("#fc4e2a","#feb24c", "#fff3d1")
# new_color_class <- c("#fff3d1","#feb24c", "#fc4e2a")
new_color_class <- c("#fc4e2a","#feb24c", "#fff3d1")
new_class_order <- fct_relevel(super_freq_anno_v2026$new_class,"Fixed","Polymorphic","Rare")
super_freq_anno_v2026_counted$ANNO <- reorder(super_freq_anno_v2026_counted$ANNO, -super_freq_anno_v2026_counted$NUMBER)

figsupp_freqspec <- ggplot(super_freq_anno_v2026,mapping=aes(fill = new_class_order))+
  geom_histogram(binwidth = 0.07,aes(x = TE_FREQ),color="black")+
  labs(
    # title = "Allele frequency distribution",
    x = "TE frequency",
    y = "TE count",
    fill="Frequency ranges")+
  scale_fill_manual(values = new_color_class) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_freqspec.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_freqspec.png"),
       plot = figsupp_freqspec, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

ggplot(super_freq_anno_v2026_counted, mapping = aes(x=ANNO,
                                              y=NUMBER,
                                              fill=anno_order))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#440154FF","#453781FF","#39558CFF","#2D718EFF","#29AF7FFF","#56C667FF","#B8DE29FF","#FDE725FF"))+
  # scale_color_manual(values = COLS,breaks=legend_order)+
  labs(
    x = "TE type",
    y = "TE count",
    fill=("TE type")
  ) +
  theme_classic()+
  theme(plot.title = element_text( size = 16, face = "bold", hjust = 0.5))+
  scale_y_continuous(breaks = seq(0,max(super_freq_anno_v2026_counted$NUMBER), by = 500))+
  # theme(axis.text.y=element_blank())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  facet_grid(~factor(CLASS,levels=c("Rare","Polymorphic","Fixed")))+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# Global folded SFS from super_freq_3, conflict-proof

bin_width <- 0.02
breaks <- seq(0, 0.50, by = bin_width)
mids <- head(breaks, -1) + bin_width/2

global_maf <- super_freq_3 %>%
  dplyr::mutate(FREQ_ALT = as.numeric(FREQ_ALT)) %>%
  dplyr::filter(!is.na(FREQ_ALT), FREQ_ALT >= 0, FREQ_ALT <= 1) %>%
  dplyr::transmute(maf = pmin(FREQ_ALT, 1 - FREQ_ALT))

counts_tbl <- global_maf %>%
  mutate(
    bin = findInterval(maf, breaks, rightmost.closed = TRUE, all.inside = TRUE),
    bin = ifelse(bin == length(breaks), length(breaks) - 1, bin)
  ) %>%
  group_by(bin) %>%
  summarise(N = dplyr::n(), .groups = "drop") %>%
  right_join(tibble::tibble(bin = seq_len(length(breaks) - 1)), by = "bin") %>%
  mutate(
    N         = dplyr::coalesce(N, 0L),
    bin_start = breaks[bin],
    bin_end   = breaks[bin + 1],
    bin_mid   = mids[bin],
    density   = N / sum(N) / bin_width
  )

ggplot(counts_tbl, aes(x = bin_mid, y = N)) +
  geom_col(color = "black", fill = "grey70", width = bin_width) +
  geom_line(aes(y = density * max(N) * bin_width), linewidth = 0.7) +
  scale_x_continuous(breaks = seq(0, 0.5, by = 0.05), limits = c(0, 0.5)) +
  labs(title = "Global Folded Site-Frequency Spectrum",
                x = "Minor-allele frequency (folded)",
                y = "Variant count") +
  theme_classic()

# Unfolded SFS from super_freq_3 with clash-safe n()

bin_width <- 0.02
breaks <- seq(0, 1, by = bin_width)
mids <- head(breaks, -1) + bin_width/2

global_p <- super_freq_3 %>%
  mutate(FREQ_ALT = as.numeric(FREQ_ALT)) %>%
  filter(!is.na(FREQ_ALT), FREQ_ALT >= 0, FREQ_ALT <= 1) %>%
  transmute(p = FREQ_ALT)

counts_tbl <- global_p %>%
  mutate(
    bin = findInterval(p, breaks, rightmost.closed = TRUE, all.inside = TRUE),
    bin = ifelse(bin == length(breaks), length(breaks) - 1, bin)
  ) %>%
  group_by(bin) %>%
  summarise(N = dplyr::n(), .groups = "drop") %>%
  right_join(tibble(bin = seq_len(length(breaks) - 1)), by = "bin") %>%
  mutate(
    N         = coalesce(N, 0L),
    bin_start = breaks[bin],
    bin_end   = breaks[bin + 1],
    bin_mid   = mids[bin],
    density   = N / sum(N) / bin_width
  )

ggplot(counts_tbl, aes(x = bin_mid, y = N)) +
  geom_col(color = "black", fill = "grey70", width = bin_width) +
  geom_line(aes(y = density * max(N) * bin_width), linewidth = 0.7) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1)) +
  labs(title = "Global Unfolded Site-Frequency Spectrum",
       x = "Derived allele frequency (p)",
       y = "Variant count") +
  theme_classic()

## define freq ranges - DISCARD

# # super_freq_anno$new_class <- lapply(super_freq_anno$FREQ_ALT, function(x) if(x>=0.95){
# #   super_freq_anno$new_class="fixed"
# # }else if (0.05 < x && x < 0.95){
# #   super_freq_anno$new_class="common"
# # }else if (0.04 < x && x < 0.05){
# #   super_freq_anno$new_class="rare5"
# # }else if (0.03 < x && x < 0.04){
# #   super_freq_anno$new_class="rare4"
# # }else if (0.02 < x && x < 0.03){
# #   super_freq_anno$new_class="rare3"
# # }else if (0.01 < x && x < 0.02){
# #   super_freq_anno$new_class="rare2"
# # }else 
# #   super_freq_anno$new_class="rare1")
# super_freq_anno$new_class <- as.character(unlist(super_freq_anno$new_class))
# new_class_table <- table(super_freq_anno$new_class) %>% as.data.frame()
# colnames(new_class_table) <- c("new_class","count")
# new_class_table$count_sum <- sum(new_class_table$count)
# new_class_table$freq <- (new_class_table$count/new_class_table$count_sum)
# 
# super_freq_anno$ANNO <- as.character(unlist(super_freq_anno$ANNO))
# super_freq_anno_counted <- super_freq_anno %>% dplyr::count(new_class,ANNO)
# colnames(super_freq_anno_counted) = c("CLASS","ANNO","NUMBER")
# super_freq_anno_counted$CLASS <- as.factor(unlist(super_freq_anno_counted$CLASS))
# 
# ## freq distribution
# 
# anno_order <- fct_relevel(super_freq_anno_counted$ANNO,"SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")
# color_class <- c("#fc4e2a","#fd8d3c","#fff3d1")
# color_class <- c("#fff3d1","#fed976","#feb24c", "#fd8d3c", "#fc4e2a","#fc4e2c")
# new_class_order <- fct_relevel(super_freq_anno$new_class,"fixed","common","rare5","rare4","rare3","rare2")
# super_freq_anno_counted$ANNO <- reorder(super_freq_anno_counted$ANNO, -super_freq_anno_counted$NUMBER)
# 
# ggplot(super_freq_anno,mapping=aes(fill = new_class_order))+
#   geom_histogram(binwidth = 0.07,aes(x = FREQ_ALT),color="black")+
#   labs(
#     # title = "Allele frequency distribution",
#     x = "TE frequency",
#     y = "TE count",
#     fill="Frequency ranges")+
#   scale_fill_manual(values = color_class) +
#   theme_classic()+
#   theme(
#     plot.title = element_text(size = 16, face = "bold"),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     axis.title.y = element_text(size = 14, face = "bold"))
# 
# ggplot(super_freq_anno_counted, mapping = aes(x=ANNO,
#                                               y=NUMBER,
#                                               fill=anno_order))+
#   geom_bar(position="stack", stat="identity")+
#   scale_fill_viridis(discrete = TRUE) +
#   scale_fill_manual(values=c("#440154FF","#453781FF","#39558CFF","#2D718EFF","#29AF7FFF","#56C667FF","#B8DE29FF","#FDE725FF"))+
#   # scale_color_manual(values = COLS,breaks=legend_order)+
#   labs(
#     x = "TE type",
#     y = "TE count",
#     fill=("TE type")
#   ) +
#   theme_classic()+
#   theme(plot.title = element_text( size = 16, face = "bold", hjust = 0.5))+
#   scale_y_continuous(breaks = seq(0,max(super_freq_anno_counted$NUMBER), by = 500))+
#   # theme(axis.text.y=element_blank())+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
#   facet_grid(~factor(CLASS,levels=c("rare2","rare3","rare4","rare5","common","fixed")))+
#   theme(
#     plot.title = element_text(size = 16, face = "bold"),
#     axis.title.x = element_text(size = 14, face = "bold"),
#     axis.title.y = element_text(size = 14, face = "bold"))

## pie charts - v1

# super_freq_anno_counted_Q4 <- super_freq_anno_counted %>% filter(CLASS == "Fixed")
# super_freq_anno_counted_Q3 <- super_freq_anno_counted %>% filter(CLASS == "Major")
# super_freq_anno_counted_Q2 <- super_freq_anno_counted %>% filter(CLASS == "Common")
# super_freq_anno_counted_Q1 <- super_freq_anno_counted %>% filter(CLASS == "Rare")
# super_freq_anno_counted_Q1_rare <- super_freq_anno_counted %>% filter(CLASS == "Very Rare")
# 
# A <- ggplot(super_freq_anno_counted_Q4, aes(x="", y=NUMBER,
#                                             fill=fct_relevel(ANNO,"LINE/L1","SINE/Alu","SINE/MIR","LTR/ERV","Retroposon/SVA","DNA/misc")))+
#   geom_bar(stat="identity", width=1,color = "black") +
#   coord_polar("y", start=0,)+
#   theme_void()+
#   scale_fill_viridis(discrete = TRUE) +
#   scale_fill_manual(values=c("#39558CFF","#440154FF","#453781FF","#29AF7FFF","#B8DE29FF","#FDE725FF"))+
#   labs(fill="TE type")+
#   ggtitle("Fixed")
# 
# B <- ggplot(super_freq_anno_counted_Q3, aes(x="", y=NUMBER,
#                                             fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","Retroposon/SVA","DNA/misc")))+
#   geom_bar(stat="identity", width=1,color = "black") +
#   coord_polar("y", start=0,)+
#   theme_void()+
#   scale_fill_viridis(discrete = TRUE) +
#   scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#B8DE29FF","#FDE725FF"))+
#   labs(fill="TE type")+
#   ggtitle("Major")
# 
# C <- ggplot(super_freq_anno_counted_Q2, aes(x="", y=NUMBER,
#                                             fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","Retroposon/SVA")))+
#   geom_bar(stat="identity", width=1,color = "black") +
#   coord_polar("y", start=0,)+
#   theme_void()+
#   scale_fill_viridis(discrete = TRUE) +
#   scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#B8DE29FF"))+
#   labs(fill="TE type")+
#   ggtitle("Common")
# 
# D <- ggplot(super_freq_anno_counted_Q1, aes(x="", y=NUMBER,
#                                             fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","Retroposon/SVA","DNA/misc")))+
#   geom_bar(stat="identity", width=1,color = "black") +
#   coord_polar("y", start=0,)+
#   theme_void()+
#   scale_fill_viridis(discrete = TRUE) +
#   scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#B8DE29FF","#FDE725FF"))+
#   labs(fill="TE type")+
#   ggtitle("Rare")
# 
# E <- ggplot(super_freq_anno_counted_Q1_rare, aes(x="", y=NUMBER,
#                                                  fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")))+
#   geom_bar(stat="identity", width=1,color = "black") +
#   coord_polar("y", start=0,)+
#   theme_void()+
#   scale_fill_viridis(discrete = TRUE) +
#   scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF"))+
#   labs(fill="TE type")+
#   ggtitle("Very Rare")
# 
# combined <- E+D+C+B+A & theme(legend.position = "none")
# combined + plot_layout(ncol=5, nrow=1, guides = "collect", axis_titles = "collect", tag_level="keep")

## pie charts v2

super_freq_anno_v2026_counted_Q3 <- super_freq_anno_v2026_counted %>% filter(CLASS == "Fixed")
super_freq_anno_v2026_counted_Q2 <- super_freq_anno_v2026_counted %>% filter(CLASS == "Polymorphic")
super_freq_anno_v2026_counted_Q1 <- super_freq_anno_v2026_counted %>% filter(CLASS == "Rare")

unique(super_freq_anno_v2026_counted_Q1$ANNO)
length(unique(super_freq_anno_v2026_counted_Q1$ANNO))
unique(super_freq_anno_v2026_counted_Q2$ANNO)
length(unique(super_freq_anno_v2026_counted_Q2$ANNO))
unique(super_freq_anno_v2026_counted_Q3$ANNO)
length(unique(super_freq_anno_v2026_counted_Q3$ANNO))

A <- ggplot(super_freq_anno_v2026_counted_Q3, aes(x="", y=NUMBER,
                                            fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("Fixed")

B <- ggplot(super_freq_anno_v2026_counted_Q2, aes(x="", y=NUMBER,
                                            fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("Polymorphic")

C <- ggplot(super_freq_anno_v2026_counted_Q1, aes(x="", y=NUMBER,
                                            fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("Rare")

pie_charts_combined <- C+B+A & theme(legend.position = "none")
fig1_pie_charts_combined <- pie_charts_combined + 
  plot_layout(ncol=3, nrow=1, guides = "collect", axis_titles = "collect", tag_level="keep")

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "fig1_pie_charts_combined.png"),
       plot = fig1_pie_charts_combined, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

# Annovar -----------------------------------------------------------------

## annovar - genomic location count

## out put of cut -f 4 anno_g.hg38_multianno.txt | sort | uniq -c 

annovar_out <- data.frame(
  class = c("downstream","Exons","Intergenes","Introns","ncRNA exons","ncRNA introns","upstream","upstream;downstream","3'UTR","5'UTR"),
  counts = c(39 ,3,3452 ,2390 ,10 ,430 ,34 ,2,41 ,6))
new_class <- annovar_out %>%
  filter(class == 'upstream'| class == "upstream;downstream" | class == "downstream") 
sum(new_class$counts)
upstream_downstream <- data.frame(class = "Upstream/Downstream", counts = sum(new_class$counts))
updated_annovar_out <- rbind(annovar_out, upstream_downstream)
updated_annovar_out <- updated_annovar_out %>%
  filter(class != 'upstream', class != "upstream;downstream" , class != "downstream") 
updated_annovar_out$counts <- as.numeric(updated_annovar_out$counts)
updated_annovar_out$TE <- as.character("TE")
region_order <- fct_relevel(updated_annovar_out$class,"Intergenes","Introns","Exons","5'UTR","3'UTR","Upstream/Downstream","ncRNA introns","ncRNA exons")

ggplot(updated_annovar_out, aes(x=region_order,y=counts,fill=region_order)) +
  geom_bar(stat = "identity", position = "dodge",color = "black") +
  geom_text(aes(label = counts), vjust = -0.5, color = "black") +
  labs(title = "",
       x = "Genomic region",
       y = "TE counts",
       fill = "Genomic region" ) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2")+
  # + theme(axis.text.x = element_text(angle = 360, hjust = 1))
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

## annovar - location per family

annovar_res <- read.table(file = file.path(sys_dir,"PanTE_human/annovar/R_1alt_annovar_out.txt"),
                          header=T)
colnames(annovar_res) <- c("CHR","POS","LOC")
annovar_res$LOC <- gsub("upstream;downstream", "Upstream/Downstream", annovar_res$LOC)
annovar_res$LOC <- gsub("upstream", "Upstream/Downstream", annovar_res$LOC)
annovar_res$LOC <- gsub("downstream", "Upstream/Downstream", annovar_res$LOC)
annovar_res$LOC <- gsub("intergenic", "Intergenes", annovar_res$LOC)
annovar_res$LOC <- gsub("intronic", "Introns", annovar_res$LOC)
annovar_res$LOC <- gsub("exonic", "Exons", annovar_res$LOC)
annovar_res$LOC <- gsub("UTR5", "5'UTR", annovar_res$LOC)
annovar_res$LOC <- gsub("UTR3", "3'UTR", annovar_res$LOC)
annovar_res$LOC <- gsub("ncRNA_intronic", "ncRNA_Introns", annovar_res$LOC)
annovar_res$LOC <- gsub("ncRNA_exonic", "ncRNA_Exons", annovar_res$LOC)

# super_freq_anno_v2026 instead of super_freq_anno
TE_loc <- merge(super_freq_anno_v2026,annovar_res)
nrow(TE_loc)
TE_loc_c <- table(TE_loc$LOC, TE_loc$new_class, TE_loc$ANNO) %>% data.frame()
colnames(TE_loc_c) <- c("LOC","Freq","Anno","counts")

region_order <- fct_relevel(TE_loc_c$LOC,"Intergenes","Introns","Exons","5'UTR","3'UTR","Upstream/Downstream","ncRNA_Introns","ncRNA_Exons")
anno_order <- fct_relevel(TE_loc_c$Anno,"SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")
anno_colors <- c("#440154FF","#453781FF","#39558CFF","#238A8DFF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF")
color_class <- c("#fff3d1","#feb24c", "#fc4e2a")
freq_order <- c("Rare","Polymorphic","Fixed")
# color_class <- c("#fff3d1","#fed976","#feb24c", "#fd8d3c", "#fc4e2a")
# freq_order <- c("Very Rare","Rare","Common","Major","Fixed")
TE_loc_c$Freq <- factor(TE_loc_c$Freq, levels = freq_order)
TE_loc_c$LOC <- reorder(TE_loc_c$LOC, -TE_loc_c$counts)
TE_loc_c$Anno <- reorder(TE_loc_c$Anno, -TE_loc_c$counts)

P <- ggplot(TE_loc_c, aes(x=Anno,y=(log10(counts)),fill=region_order)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "",
       x = "Genomic region",
       y = "TE counts (log10)",
       fill = "Genomic region" ) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2")+
  facet_wrap(~Freq,nrow=1,ncol=5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g <- ggplot_gtable(ggplot_build(P))
striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) )
fills <- color_class
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

p <- ggplot(TE_loc_c, aes(x=Anno,y=(counts),fill=region_order)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "",
       x = "Genomic region",
       y = "TE counts",
       fill = "TE types" ) +
  theme_classic() +
  scale_fill_brewer(palette = "Dark2")+
  # scale_fill_manual(values = anno_colors)+
  facet_wrap(~Freq,nrow=1,ncol=5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
g <- ggplot_gtable(ggplot_build(p))
striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) )
fills <- color_class
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

## combine freq distribution and annovar

TE_newloc <- TE_loc_c
TE_newloc$LOC <- gsub("5'UTR", "Regulatory", TE_newloc$LOC)
TE_newloc$LOC <- gsub("3'UTR", "Regulatory", TE_newloc$LOC)
TE_newloc$LOC <- gsub("Upstream/Downstream", "Regulatory", TE_newloc$LOC)
TE_newloc$LOC <- gsub("ncRNA_Exons", "Exons", TE_newloc$LOC)
TE_newloc$LOC <- gsub("ncRNA_Introns","Introns", TE_newloc$LOC)

region_order <- fct_relevel(TE_newloc$LOC,"Intergenes","Regulatory","Introns","Exons")
anno_order <- fct_relevel(TE_newloc$Anno,"SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")
anno_order_by_type <- fct_relevel(TE_newloc$Anno, "SINE/Alu", "SINE/MIR", "LINE/L1", "LINE/misc", "LTR/ERV", "LTR/misc", "Retroposon/SVA", "DNA/misc")
anno_colors <- c("#440154FF","#453781FF","#39558CFF","#238A8DFF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF")
color_class <- c("#fff3d1","#feb24c", "#fc4e2a")
freq_order <- c("Rare","Polymorphic","Fixed")
# color_class <- c("#fff3d1","#fed976","#feb24c", "#fd8d3c", "#fc4e2a")
# freq_order <- c("Very Rare","Rare","Common","Major","Fixed")
TE_newloc$Freq <- factor(TE_newloc$Freq, levels = freq_order)
TE_newloc$LOC <- reorder(TE_newloc$LOC, -TE_newloc$counts)
TE_newloc$Anno <- reorder(TE_newloc$Anno, -TE_newloc$counts)

sum(TE_loc_c[TE_loc_c$LOC == "3'UTR",]$counts)+
  sum(TE_loc_c[TE_loc_c$LOC == "Upstream/Downstream",]$counts)+
  sum(TE_loc_c[TE_loc_c$LOC == "5'UTR",]$counts)

sum(TE_newloc[TE_newloc$LOC == "Regulatory",]$counts)

anno_levels <- c("SINE/Alu", "SINE/MIR", "LINE/L1", "LINE/misc", "LTR/ERV", "LTR/misc", "Retroposon/SVA", "DNA/misc")
color_gradients <- lapply(anno_colors, function(color) {
  colorRampPalette(c(color, "#F4F3EE"))(4)
})
names(color_gradients) <- anno_levels
color_mapping <- unlist(lapply(seq_along(color_gradients), function(i) {
  setNames(color_gradients[[i]], paste0(names(color_gradients)[i], "_", 1:4))
}))
TE_newloc$Color_Group <- paste0(TE_newloc$Anno, "_", as.numeric(region_order))

p <- ggplot(TE_newloc, aes(x=Anno, y= log10(counts+1), fill = Color_Group)) +
  geom_bar(stat = "identity", position = "stack",) +
  labs(title = "",
       x = "TE types",
       y = "log10 (counts+1)",
       fill = "") +
  theme_classic() +
  scale_fill_manual(values = color_mapping) +
  facet_wrap(~Freq, nrow = 1, ncol = 5) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"),
        strip.background = element_rect(linewidth = 1.2))
g <- ggplot_gtable(ggplot_build(p))
striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) )
fills <- c("white","white", "white")
k <- 1
for (i in striprt) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
grid.draw(g)

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_freqanno_mod.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "fig1_freqanno_mod.png"),
       plot = fig1_freqanno_mod, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

# ggplot(TE_newloc, aes(x=Anno, y= log10(counts+1), fill = Color_Group)) +
#   geom_bar(stat = "identity", position = "stack",) +
#   labs(title = "",
#        x = "TE types",
#        y = "log10 (counts+1)",
#        fill = "") +
#   theme_classic() +
#   scale_fill_manual(values = color_mapping) +
#   facet_wrap(~Freq, nrow = 1, ncol = 5) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.position = "none") +
#   theme(strip.background = element_blank(),
#         strip.text.x = element_text(colour = "black", margin = margin(b = 4))) +
#   geom_segment(data = data.frame(x = -Inf, xend = Inf, y = Inf, yend = Inf),
#                aes(x = x, xend = xend, y = y, yend = yend),
#                inherit.aes = FALSE)

fig1_freqannolegand <- ggplot(TE_newloc, aes(x= Anno, y= log10(counts+1), fill = region_order)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(fill = "Genomic region") +
  theme_classic() +
  scale_fill_manual(values = colorRampPalette(c("#463f3a","#F4F3EE"))(4)) +
  theme_void()

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_freqannolegand.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "fig1_freqannolegand.png"),
       plot = fig1_freqannolegand, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

table(TE_newloc$Freq)
table(super_freq_anno_v2026$new_class)
# Fixed Polymorphic        Rare 
# 234         857        5316 

Rare   <- 5316
Polymorphic <- 857
Fixed  <- 234
total <- Rare+Fixed+Polymorphic
100-((total-Rare)/total)*100 
100-((total-Polymorphic)/total)*100 
100-((total-Fixed)/total)*100

## old
# Rare   <- 5455
# Polymorphic <- 857
# Fixed  <- 95
# total <- Rare+Fixed+Polymorphic
# 100-((total-Rare)/total)*100 
# 100-((total-Polymorphic)/total)*100 
# 100-((total-Fixed)/total)*100

## test whether bar-stacking and log transformation works together

TE_newloc_veryrare_alu <- TE_newloc %>% filter(Freq == "Very Rare") %>% filter(Anno == "SINE/Alu")
(sum(log10(TE_newloc_veryrare_alu$counts+1)))
region_order_test <- fct_relevel(TE_newloc_veryrare_alu$LOC,"Exons","Regulatory","Introns","Intergenes")
TE_newloc_veryrare_alu$Color_Group_test <- paste0(TE_newloc_veryrare_alu$Anno, "_", as.numeric(region_order_test))

ggplot(TE_newloc_veryrare_alu, aes(x=Anno, y= log10(counts+1),fill = region_order_test)) +
  geom_bar(stat = "identity", position = "stack")+
  scale_fill_manual(values = c("Exons"="red","Regulatory"="blue","Introns"="green","Intergenes"="black"))

# p <- ggplot(TE_newloc, aes(x= Anno, y= log10(counts + 1), fill= anno_order))+
#   geom_bar_pattern(stat = "identity",
#                    pattern_color = "white",
#                    pattern_fill = "white",
#                    pattern_density = 0.35,
#                    aes(pattern = LOC, pattern_angle = LOC, pattern_spacing = LOC))+
#   labs(x = "TE types",
#        y = "log10 (counts+1)",
#        fill = "TE types",
#        pattern = "Genomic location")+
#   theme_classic()+
#   scale_pattern_manual(
#     values = c("none","stripe","circle","wave"),
#     guide = guide_legend(override.aes = list(pattern_fill = 'grey')))+
#   scale_fill_manual(values = anno_colors) +
#   facet_wrap(~Freq, nrow = 1, ncol = 5) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# g <- ggplot_gtable(ggplot_build(p))
# striprt <- which( grepl('strip-r', g$layout$name) | grepl('strip-t', g$layout$name) )
# fills <- color_class
# k <- 1
# for (i in striprt) {
#   j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
#   g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
#   k <- k+1
# }
# grid.draw(g)

## attempt to diagnosis the missing NA in AnnoVar tool - 2026 analysis [in progress]

TE_loc <- merge(
  super_freq_anno,
  annovar_res,
  by = c("CHR", "POS"),
  all.x = TRUE)

NA_df <- TE_loc %>%
  filter(is.na(LOC))

NA_df_anno <- NA_df %>%
  left_join(
    annovar_res_allinfo %>% select(CHR, POS, LOC),
    by = c("CHR", "POS"),
    suffix = c("", "_AnnoVar"))

NA_df_anno %>%
  filter(is.na(LOC_AnnoVar))

loc_counts <- NA_df_anno %>%
  count(LOC_AnnoVar, sort = TRUE)

ggplot(loc_counts, aes(x = reorder(LOC_AnnoVar, n), y = n)) +
  geom_col() +
  coord_flip() +
  labs(x = "LOC", y = "Count") +
  theme_classic()

# rGREAT per freq - v2 ----------------------------------------------------

TE_anno_freq_Gr <- makeGRangesFromDataFrame(df=super_freq_anno_v2026,
                                            keep.extra.columns=T,
                                            ignore.strand=T,
                                            seqnames.field="CHR",
                                            start.field="POS",
                                            end.field="end",
                                            starts.in.df.are.0based=FALSE)

super_freq_anno_v2026 <- super_freq_anno_v2026 %>%
  dplyr::select(c("CHR","POS","ALT","ANNO","END","new_class"))
super_freq_anno_v2026_rare <- super_freq_anno_v2026 %>%
  filter(new_class == "Rare")
super_freq_anno_v2026_polymorphic <- super_freq_anno_v2026 %>%
  filter(new_class == "Polymorphic")
super_freq_anno_v2026_fixed <- super_freq_anno_v2026 %>%
  filter(new_class == "Fixed")

nrow(super_freq_anno_v2026_polymorphic)+nrow(super_freq_anno_v2026_rare)+nrow(super_freq_anno_v2026_fixed)

TE_anno_freq_rare_Gr <- makeGRangesFromDataFrame(df=super_freq_anno_v2026_rare,
                                            keep.extra.columns=T,
                                            ignore.strand=T,
                                            seqnames.field="CHR",
                                            start.field="POS",
                                            end.field="end",
                                            starts.in.df.are.0based=FALSE)
TE_anno_freq_polymorphic_Gr <- makeGRangesFromDataFrame(df=super_freq_anno_v2026_polymorphic,
                                            keep.extra.columns=T,
                                            ignore.strand=T,
                                            seqnames.field="CHR",
                                            start.field="POS",
                                            end.field="end",
                                            starts.in.df.are.0based=FALSE)
TE_anno_freq_fixed_Gr <- makeGRangesFromDataFrame(df=super_freq_anno_v2026_fixed,
                                            keep.extra.columns=T,
                                            ignore.strand=T,
                                            seqnames.field="CHR",
                                            start.field="POS",
                                            end.field="end",
                                            starts.in.df.are.0based=FALSE)

# set.seed(123)
# TE_anno_freq_rare_job = submitGreatJob(TE_anno_freq_rare_Gr,species="hg38")
# plotRegionGeneAssociations(TE_anno_freq_rare_job)
# TE_anno_freq_rare_res = great(TE_anno_freq_rare_job, "GO:BP", "txdb:hg38") ## GO:BP or MSigDB:H
# plotVolcano(TE_anno_freq_rare_res)
# set.seed(1234)
# TE_anno_freq_common_job = submitGreatJob(TE_anno_freq_common_Gr,species="hg38")
# plotRegionGeneAssociations(TE_anno_freq_common_job)
# TE_anno_freq_common_res = great(TE_anno_freq_common_job, "GO:BP", "txdb:hg38")
# plotVolcano(TE_anno_freq_common_res)
# set.seed(12345)
# TE_anno_freq_fixed_job = submitGreatJob(TE_anno_freq_fixed_Gr,species="hg38")
# plotRegionGeneAssociations(TE_anno_freq_fixed_job)
# TE_anno_freq_fixed_res = great(TE_anno_freq_fixed_job, "GO:BP", "txdb:hg38")
# plotVolcano(TE_anno_freq_fixed_res)

## OR run locally

TE_anno_freq_rare_res = great(TE_anno_freq_rare_Gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene")
plotVolcano(TE_anno_freq_rare_res)
TE_anno_freq_polymorphic_res = great(TE_anno_freq_polymorphic_Gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene")
plotVolcano(TE_anno_freq_polymorphic_res)
TE_anno_freq_fixed_res = great(TE_anno_freq_fixed_Gr, "GO:BP", "TxDb.Hsapiens.UCSC.hg38.knownGene")
plotVolcano(TE_anno_freq_fixed_res)

tb_rare         <- getEnrichmentTable(TE_anno_freq_rare_res)   %>% mutate(class="rare")
tb_polymorphic  <- getEnrichmentTable(TE_anno_freq_polymorphic_res) %>% mutate(class="polymorphic")
tb_fixed        <- getEnrichmentTable(TE_anno_freq_fixed_res)  %>% mutate(class="fixed")

tb_allfreq <- bind_rows(tb_rare, tb_polymorphic, tb_fixed) %>%
  # filter(p_adjust <= 0.05) %>%
  filter(p_adjust_hyper <= 0.05 | p_adjust <= 0.05)
  # filter(mean_tss_dist < 300000)
  # filter(observed_region_hits >= 500)
  # group_by(class) %>%
  # slice_min(p_adjust_hyper, n = 10)
  # filter(class == "rare")
  # select(class, id, description, mean_tss_dist,fold_enrichment ,p_adjust, fold_enrichment_hyper,p_adjust_hyper)

# write_xlsx(tb_allfreq, file.path(sys_dir,'PanTE_human_2025_v2/scripts/rGREAT_tb_res_allfreq_v2.xlsx'))
# write_xlsx(tb_allfreq, file.path(sys_dir,'PanTE_human_2025_v2/scripts/rGREAT_tb_res_allfreq_v3.xlsx'))

mean(tb_allfreq$mean_tss_dist)

## Human Genomics Responses

# Top enriched GO:BP terms for each TE frequency class

tb_allfreq_master_top10 <- bind_rows(
  getEnrichmentTable(TE_anno_freq_rare_res) %>% mutate(class = "Rare"),
  getEnrichmentTable(TE_anno_freq_polymorphic_res) %>% mutate(class = "Polymorphic"),
  getEnrichmentTable(TE_anno_freq_fixed_res) %>% mutate(class = "Fixed")) %>%
  filter(p_adjust_hyper <= 0.05 | p_adjust <= 0.05) %>%
  filter(mean_tss_dist < 400000) %>%
  mutate(
    stats_sig = case_when(
      p_adjust_hyper <= 0.05 & p_adjust <= 0.05 ~ "Both",
      p_adjust_hyper <= 0.05 ~ "Hypergeometric",
      p_adjust <= 0.05 ~ "Binomial",
      TRUE ~ "Not_significant")) %>%
  group_by(class) %>%
  filter(
    (stats_sig %in% c("Hypergeometric", "Both") &
       rank(p_adjust_hyper, ties.method = "first", na.last = "keep") <= 10) |
      (stats_sig %in% c("Binomial", "Both") &
         rank(p_adjust, ties.method = "first", na.last = "keep") <= 10)
  ) %>%
  ungroup() %>%
  dplyr::select(
    class, id, description,
    mean_tss_dist, stats_sig,
    fold_enrichment_hyper, p_value_hyper, p_adjust_hyper,
    fold_enrichment, p_value, p_adjust,
    genome_fraction, observed_region_hits, observed_gene_hits, gene_set_size)

table(tb_allfreq_master_top10$stats_sig)

hyper_rGREAT_fc_plot <- ggplot(subset(tb_allfreq_master_top10, stats_sig %in% c("Hypergeometric", "Both")),
                               aes(x = fold_enrichment_hyper, y = reorder(description, -log10(p_adjust_hyper)), fill = p_adjust_hyper)) +
  geom_col() +
  facet_wrap(~class, scales = "free_y") +
  scale_fill_gradient(low = "firebrick", high = "darkgreen") +
  labs(x = "Hypergeometric fold enrichment", y = "GO Biological Processes", fill = "Adjusted p-value") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

binomial_rGREAT_fc_plot <- ggplot(subset(tb_allfreq_master_top10, stats_sig %in% c("Binomial", "Both")),
                                  aes(x = fold_enrichment, y = reorder(description, -log10(p_adjust)), fill = p_adjust)) +
  geom_col() +
  facet_wrap(~class, scales = "free_y") +
  scale_fill_gradient(low = "firebrick", high = "darkgreen") +
  labs(x = "Binomial fold enrichment", y = "GO Biological Processes", fill = "Adjusted p-value") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

hyper_rGREAT_fc_plot / binomial_rGREAT_fc_plot

write_xlsx(tb_allfreq_master_top10, file.path(sys_dir,'PanTE_human_2025_v2/scripts/rGREAT_tb_res_allfreq_master_top10_v4_final.xlsx'))

## get genes of sig enrich

top10_terms_rGREAT <- tb_allfreq_master_top10 %>%
  dplyr::select(class, id, description)

RegionGeneAss_Rare <- as.data.frame(getRegionGeneAssociations(TE_anno_freq_rare_res))
RegionGeneAss_Rare_map <- unlist(RegionGeneAss_Rare$annotated_genes)
genes_top10_terms_rGREAT_rare <- top10_terms_rGREAT %>%
  filter(class == "Rare") %>%
  rowwise() %>%
  mutate(
    Entrez_ID = list(intersect(
      names(RegionGeneAss_Rare_map),
      TE_anno_freq_rare_res@gene_sets[[id]])),
    annotated_genes = paste(unname(RegionGeneAss_Rare_map[Entrez_ID]), collapse = ", ")) %>%
  ungroup() %>%
  filter(!description %in% c("cellular process", "biological_process"))

RegionGeneAss_Polymorphic <- as.data.frame(getRegionGeneAssociations(TE_anno_freq_polymorphic_res))
RegionGeneAss_Polymorphic_map <- unlist(RegionGeneAss_Polymorphic$annotated_genes)
genes_top10_terms_rGREAT_polymorphic <- top10_terms_rGREAT %>%
  filter(class == "Polymorphic") %>%
  rowwise() %>%
  mutate(
    Entrez_ID = list(intersect(
      names(RegionGeneAss_Polymorphic_map),
      TE_anno_freq_polymorphic_res@gene_sets[[id]])),
    annotated_genes = paste(unname(RegionGeneAss_Polymorphic_map[Entrez_ID]), collapse = ", ")) %>%
  ungroup()

RegionGeneAss_Fixed <- as.data.frame(getRegionGeneAssociations(TE_anno_freq_fixed_res))
RegionGeneAss_Fixed_map <- unlist(RegionGeneAss_Fixed$annotated_genes)
genes_top10_terms_rGREAT_fixed <- top10_terms_rGREAT %>%
  filter(class == "Fixed") %>%
  rowwise() %>%
  mutate(
    Entrez_ID = list(intersect(
      names(RegionGeneAss_Fixed_map),
      TE_anno_freq_fixed_res@gene_sets[[id]])),
    annotated_genes = paste(unname(RegionGeneAss_Fixed_map[Entrez_ID]), collapse = ", ")) %>%
  ungroup()

sum(is.na(RegionGeneAss_Rare_map))
sum(is.na(RegionGeneAss_Polymorphic_map))
sum(is.na(RegionGeneAss_Fixed_map))

table(lengths(RegionGeneAss_Rare$annotated_genes))
table(lengths(genes_top10_terms_rGREAT_rare$annotated_genes))
table(lengths(RegionGeneAss_Polymorphic$annotated_genes))
table(lengths(genes_top10_terms_rGREAT_polymorphic$annotated_genes))
sum(length(unique(unlist(strsplit(genes_top10_terms_rGREAT_rare$annotated_genes, ", ")))))
sum(length(unique(unlist(strsplit(genes_top10_terms_rGREAT_polymorphic$annotated_genes, ", ")))))
nrow(genes_top10_terms_rGREAT_rare)
nrow(genes_top10_terms_rGREAT_polymorphic)

genes_rGREAT_rare_vec <- unique(unlist(strsplit(genes_top10_terms_rGREAT_rare$annotated_genes, ", ")))
genes_rGREAT_polymorphic_vec <- unique(unlist(strsplit(genes_top10_terms_rGREAT_polymorphic$annotated_genes, ", ")))

genes_rGREAT_vec <- list(Rare = genes_rGREAT_rare_vec, Polymorphic = genes_rGREAT_polymorphic_vec)

sum(length(unique(unlist(genes_rGREAT_vec))))

genes_top10_terms_rGREAT_rare %>%
  filter(grepl("PLCB1", annotated_genes)) %>%
  pull(description)
genes_top10_terms_rGREAT_rare %>%
  filter(grepl("SULF1", annotated_genes)) %>%
  pull(description)

genes_top10_terms_rGREAT_polymorphic %>%
  filter(grepl("PLCB1", annotated_genes)) %>%
  pull(description)
genes_top10_terms_rGREAT_polymorphic %>%
  filter(grepl("SULF1", annotated_genes)) %>%
  pull(description)

# Recombination rate, TE size, and gene distance ------------------------------------------------------

## define recombination rate date

rec_rate <- read.table(file = file.path(sys_dir,"PanTE_human/recombination_rate/recombAvg.bed"), header = F,
                       sep = "",
                       dec = ".") %>% dplyr::select(c(-4))
colnames(rec_rate) = c("chr","start","end","rate")

## create 10 000 bins for recombination rate

chromo <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7", "chr8", "chr9", "chr10",
            "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
            "chr21", "chr22","chrX")
bin          <- 10000

df2          <- rec_rate %>% 
  mutate(window = start %/% bin) %>% 
  group_by(window,chr) %>%
  summarise(mean = mean(rate)) %>%
  mutate(start = (window*bin)+1, stop=(window+1)*bin) 

df_out       <- c()
for(i in unique(df2$chr)) {
  df2_tmp      <- df2 %>% filter(chr == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(chr    = i,
           mean = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(chr,window)
  df_out       <- bind_rows(df_out,df_tmp)
}

rec_rate_binned <- df_out %>% arrange(chr, window) 

fc.ymax <- ceiling(max(abs(range(rec_rate_binned$mean))))
fc.ymin <- -fc.ymax
col.over <- "tomato"
col.under <- "gray"
sign.col <- rep(col.over, length(rec_rate_binned))
sign.col[rec_rate_binned$mean>150] <- col.under

ggplot(rec_rate_binned, mapping=aes(x=window,y=mean,col=sign.col))+
  geom_point()+
  facet_wrap(~factor(chr, levels=chromo))

rec_rate_binned_Gr <- makeGRangesFromDataFrame(df=rec_rate_binned,
                                               keep.extra.columns=T,
                                               ignore.strand=T,
                                               seqnames.field="chr",
                                               start.field="start",
                                               end.field="stop",
                                               starts.in.df.are.0based=FALSE)

rec_rate_binned_noNA <- rec_rate_binned %>% drop_na()

ggplot()+
  geom_boxplot(rec_rate_binned_noNA, mapping=aes(x=chr,y=log10(mean)))

## define gene distance and size for super_freq_anno

## size

colnames(super_freq_anno)
colnames(super_freq_anno_v2026)

# super_freq_anno$SIZE <- nchar(super_freq_anno$ALT)
super_freq_anno_v2026$SIZE <- super_freq_anno_v2026$TE_size
nrow(super_freq_anno_v2026)

## gene distance rGREAT

set.seed(123)
TE_anno_job = submitGreatJob(TE_anno_Gr,species="hg38")
## Error in (function (type, msg, asError = TRUE) : Could not resolve host: great.stanford.edu
## Run locally 
TE_anno_job = great(TE_anno_Gr, "GO:BP", "hg38")
## or use old data
table = getEnrichmentTables(TE_anno_job)
str(table)
plotRegionGeneAssociations(TE_anno_job)
plotRegionGeneAssociations(TE_anno_job,ontology = "GO Biological Process",term_id="GO:0003251")
availableOntologies(TE_anno_job)
res = great(TE_anno_Gr, "MSigDB:H", "txdb:hg38")
plotRegionGeneAssociations(res)
plotVolcano(res)

TE_distTSS <- data.frame(TE_anno_job@association_tables[["all"]][["chr"]],
                           (TE_anno_job@association_tables[["all"]][["start"]]+1),
                           abs(TE_anno_job@association_tables[["all"]][["distTSS"]]),
                           TE_anno_job@association_tables[["all"]][["name"]])
colnames(TE_distTSS) <- c("CHR","POS","distTSS","names")

min_indices <- tapply(seq_along(TE_distTSS$distTSS), TE_distTSS$names, function(x) x[which.min(TE_distTSS$distTSS[x])])
TE_distTSS_filtered <- TE_distTSS[unlist(min_indices), ]

# TE_distTSS <- TE_distTSS[duplicated(TE_distTSS$names),] ######## alawys keep the first one (need to get new gene dis)

# rm(list = setdiff(ls(), c("TE_anno_job", "TE_distTSS", "TE_distTSS_filtered")))

nrow(TE_distTSS_filtered)
nrow(TE_distTSS)

# FREQ_ALT should be replaced with TE_FREQ

TE_info_distTSS <- merge(super_freq_anno_v2026,TE_distTSS_filtered)
TE_info_distTSS$ANNO <- as.factor(unlist(TE_info_distTSS$ANNO))
TE_info_distTSS <- TE_info_distTSS[order(TE_info_distTSS[,1],TE_info_distTSS[,2]),]
TE_size_distTSS <- TE_info_distTSS %>%
  dplyr::select(c("CHR","POS","END","TE_FREQ","distTSS","SIZE","ANNO"))

## create 10 000 bins for TE_size_distTSS

bin          <- 10000

df2          <- TE_size_distTSS %>% 
  mutate(window = POS %/% bin)
df2$window_CHR <- paste(df2$window, df2$CHR, sep = "_") 
df2          <- df2 %>%
  group_by(window_CHR) %>%
  summarise(sum_anno = list(ANNO),across(c(TE_FREQ, distTSS,SIZE), mean))
df2[,6:7] <- stringr::str_split_fixed(df2$window_CHR, "_", 2)
df2 <- df2[-c(1)] 
colnames(df2) <- c("sum_anno","mean_TEfreq","mean_disttss","mean_size","window","CHR")
df2$window <- as.numeric(df2$window)
df2          <- df2 %>%
  mutate(start = (window*bin)+1, stop=(window+1)*bin) %>% 
  relocate(c("window","CHR","start","stop","sum_anno","mean_TEfreq","mean_disttss","mean_size"))

df_out       <- c()
for(i in unique(df2$CHR)) {
  df2_tmp      <- df2 %>% filter(CHR == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(CHR    = i,
           mean_TEfreq = NA,
           mean_disttss = NA,
           mean_size = NA,
           sum_anno = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(CHR,window)
  df_out       <- bind_rows(df_out,df_tmp)
}

TE_size_distTSS_binned <- df_out %>% arrange(CHR, window) %>% drop_na()

fc.ymax <- ceiling(max(abs(range(TE_size_distTSS_binned$mean_TEfreq))))
fc.ymin <- -fc.ymax
col.over <- "tomato"
col.under <- "gray"
sign.col <- rep(col.over, length(TE_size_distTSS_binned))
sign.col[TE_size_distTSS_binned$mean_TEfreq>0.2] <- col.under

ggplot(TE_size_distTSS_binned, mapping=aes(x=window,y=mean_TEfreq,col=sign.col))+
  geom_point()+
  facet_wrap(~factor(CHR, levels=chromo))

colnames(TE_size_distTSS_binned)
colnames(rec_rate_binned_noNA) = c("window","CHR","mean_recrate","pos","end")
TE_size_distTSS_rec_rate <- merge(rec_rate_binned_noNA,TE_size_distTSS_binned)
TE_size_distTSS_rec_rate$new_class <- lapply(TE_size_distTSS_rec_rate$mean_TEfreq, function(x) if(x>=0.95){
  TE_size_distTSS_rec_rate$new_class="Fixed"
}else if (0.05 <= x && x < 0.95){
  TE_size_distTSS_rec_rate$new_class="Polymorphic"
}else 
  TE_size_distTSS_rec_rate$new_class="Rare")
# TE_size_distTSS_rec_rate$new_class <- lapply(TE_size_distTSS_rec_rate$mean_TEfreq, function(x) if(x>0.901){
#   TE_size_distTSS_rec_rate$new_class="Fixed"
# }else if (0.634 < x && x <= 0.901){
#   TE_size_distTSS_rec_rate$new_class="Major"
# }else if (0.367 < x && x <= 0.634){
#   TE_size_distTSS_rec_rate$new_class="Common"
# }else if ((0.10) < x && x <= 0.367){
#   TE_size_distTSS_rec_rate$new_class="Rare"
# }else 
#   TE_size_distTSS_rec_rate$new_class="Very Rare")
TE_size_distTSS_rec_rate$new_class <- as.character(unlist(TE_size_distTSS_rec_rate$new_class))
table(TE_size_distTSS_rec_rate$new_class)
TE_genomic_windows <- TE_size_distTSS_rec_rate %>% arrange(CHR,window) %>% 
  dplyr::select(c("CHR","start","stop","window","sum_anno","mean_recrate","mean_TEfreq","mean_disttss","mean_size","new_class"))
TE_genomic_windows <- TE_genomic_windows[order(TE_genomic_windows[,1],TE_genomic_windows[,2]),]

nrow(TE_genomic_windows)

## TE_genomic_windows is a data frame that has recombination rate, gene distance, size, Global alt allele freq, annotation, and freq ranges information for all TE vars 

## rec rate boxplots, scatter plots, and stats

## filter outliers (SD > 3)

threshold <- 3
z_scores <- scale(TE_genomic_windows$mean_recrate)
outliers <- which(abs(z_scores) > threshold)
TE_dataframe_clean <- TE_genomic_windows[-outliers, ]

z_scores <- scale(TE_dataframe_clean$mean_disttss)
outliers <- which(abs(z_scores) > threshold)
TE_dataframe_clean <- TE_dataframe_clean[-outliers, ]

z_scores <- scale(TE_dataframe_clean$mean_size)
outliers <- which(abs(z_scores) > threshold)
TE_dataframe_clean <- TE_dataframe_clean[-outliers, ]

TE_dataframe <- TE_dataframe_clean %>%
  # filter(mean_recrate < 500, mean_recrate != 0) %>%  # %>% filter(mean_TEfreq > 0.1, mean_TEfreq < 0.90) # %>% filter(mean_size < 0.98) %>% filter(mean_disttss < 0.98) %>% 
  # mutate(freqeuncy = fct_relevel(new_class, "Fixed", "Major","Common","Rare","Very Rare")) %>% 
  mutate(freqeuncy = fct_relevel(new_class, "Fixed","Polymorphic","Rare")) %>% 
  filter(mean_recrate > 0) %>%
  filter(mean_TEfreq > 0) %>% 
  filter(mean_size >= 100) %>% 
  filter(CHR != "chrX")

## log transformation

TE_dataframe_log <- TE_dataframe
TE_dataframe_log$mean_size <- log10(TE_dataframe_log$mean_size)
TE_dataframe_log$mean_recrate <- log10(TE_dataframe_log$mean_recrate)
TE_dataframe_log$mean_disttss <- log10(TE_dataframe_log$mean_disttss)

color_class <- c("#fc4e2a","#feb24c", "#fff3d1")

# Boxplots for recombination rate, TE size, and gene distance -----------------------------------------------------------

## all TE boxplot

allTE_size <-  ggplot()+
  geom_boxplot(TE_dataframe_log,mapping=aes(x= freqeuncy,y= mean_size,fill=freqeuncy))+
  labs(y = "TE size (log10)",
       x = "Allele frequency",
       fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

allTE_disttss <- ggplot()+
  geom_boxplot(TE_dataframe_log,mapping=aes(x= freqeuncy,y= mean_disttss,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Gene distance (log10)",
       x = "Allele frequency",
       fill = "Frequency ranges"
  )+
scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

color_class <- c("#fc4e2a","#feb24c", "#fff3d1")
freq_order <- fct_relevel(TE_dataframe_log$new_class ,"Rare","Polymorphic","Fixed")
# color_class <- c("#fc4e2a","#fd8d3c","#feb24c", "#fed976", "#fff3d1")
# freq_order <- fct_relevel(TE_dataframe_log$new_class ,"Very Rare","Rare","Common","Major","Fixed")

allTE_recrate <- ggplot() +
  geom_boxplot(TE_dataframe_log,mapping=aes(x= freqeuncy,y= mean_recrate,fill=freqeuncy),outlier.shape = "|",outlier.size = 3)+
  labs(y="log10(Recombination rate)",
       x= "TE frequency",
       fill= "Frequency ranges")+
  # geom_jitter(color="black", size=0.2, alpha=0.2)
  scale_fill_manual(values = color_class)+
  theme_classic()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

TE_dataframe_mod <- TE_dataframe %>%
  filter(mean_recrate < 10)
freq_order <- fct_relevel(TE_dataframe_mod$new_class ,"Rare","Polymorphic","Fixed")

fig1_recrate_boxplot_2 <- ggplot() +
  geom_boxplot(TE_dataframe_mod,mapping=aes(x=freq_order,y= mean_recrate,fill=freqeuncy),outlier.size = 1)+
  labs(y="Recombination rate",
       x= "TE frequency",
       fill= "Frequency ranges")+
  # geom_jitter(color="black", size=0.2, alpha=0.2)
  scale_fill_manual(values = color_class)+
  theme_classic()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_recrate_boxplot_2.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "fig1_recrate_boxplot_2.png"),
       plot = fig1_recrate_boxplot_2, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

# allTE_size <- ggplot(TE_dataframe_log,aes(x = freqeuncy, y=mean_size,fill=freqeuncy)) +
#   geom_lv(color='black', size=0.75) +
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "TE length (log10 scaled)",
#        x = "Allele frequency"
#        )+
#   # geom_jitter(color="black", size=0.2, alpha=0.2)+
#   ggtitle("All TE families")+
#   scale_fill_viridis(discrete = TRUE,option="H")+
#   theme_minimal() +
#   theme(
#     panel.grid.major.y = element_line(color = "white", linetype = "dashed"),
#     axis.line = element_line(color = "white"),
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#     axis.title = element_text(size = 12, face = "bold"),
#     axis.text = element_text(size = 10))
# allTE_disttss <- ggplot(TE_dataframe_log,aes(x = freqeuncy, y=mean_disttss,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "Distance to closest gene (log10 scaled)",
#      x = "Allele frequency"
#   )+ # +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
#   scale_fill_viridis(discrete = TRUE,option="H")+
#   theme_minimal() +
#   theme(
#     panel.grid.major.y = element_line(color = "white", linetype = "dashed"),
#     axis.line = element_line(color = "white"),
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
#     axis.title = element_text(size = 12, face = "bold"),
#     axis.text = element_text(size = 10))

# color_class <- c("#fc4e2a","#fd8d3c","#feb24c", "#fed976", "#fff3d1")
# freq_order <- fct_relevel(TE_dataframe_log$new_class ,"Very Rare","Rare","Common","Major","Fixed")
# 
# allTE_recrate <-ggplot(TE_dataframe_log,aes(x = freq_order, y=mean_recrate,fill=freqeuncy)) +
#   geom_lv(color='black', size=0.75,outlier.shape = "|",outlier.size = 3) +
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000",outlier.shape = "|") +
#   labs(y="log10(Recombination rate)",
#        x= "TE frequency",
#        fill= "Frequency ranges")+
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
#   scale_fill_manual(values = color_class)+
#   theme_classic()+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_recrate_boxplot.png"), plot=ggplot2::last_plot())

## L1, Alu, SVA boxplot

TE_dataframe_log$new_class <- (unlist(TE_dataframe_log$new_class))
TE_dataframe_log_nolist <- TE_dataframe_log %>%  unnest(sum_anno)

TE_dataframe_log_Alu <- TE_dataframe_log_nolist %>% filter(sum_anno == "SINE/Alu")
TE_dataframe_log_L1 <- TE_dataframe_log_nolist %>% filter(sum_anno == "LINE/L1")
TE_dataframe_log_SVA <- TE_dataframe_log_nolist %>% filter(sum_anno == "Retroposon/SVA")

## SVA

# SVA_size <- ggplot(TE_dataframe_log_SVA, aes(x = freqeuncy, y=mean_size,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "TE length",
#        x = "Allele frequency"
#   )+
#   # geom_jitter(color="black", size=0.2, alpha=0.2)+
#   ggtitle("SVA")
SVA_size <- ggplot()+
  geom_boxplot(TE_dataframe_log_SVA,mapping=aes(x= freqeuncy,y= mean_size,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "TE size (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# SVA_disttss <- ggplot(TE_dataframe_log_SVA,aes(x = freqeuncy, y=mean_disttss,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "Distance to closest gene",
#      x = "Allele frequency"
#   )# +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
SVA_disttss <- ggplot()+
  geom_boxplot(TE_dataframe_log_SVA,mapping=aes(x= freqeuncy,y= mean_disttss,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Gene distance (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# SVA_recrate <-ggplot(TE_dataframe_log_SVA,aes(x = freqeuncy, y=mean_recrate,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette= "Dark1")+
#   labs(y = "Recombination rate",
#      x = "Allele frequency"
#   ) # +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
SVA_recrate <- ggplot()+
  geom_boxplot(TE_dataframe_log_SVA,mapping=aes(x= freqeuncy,y= mean_recrate,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Recombination rate (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+   ggtitle("SVA")+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

## Alu

# Alu_size <- ggplot(TE_dataframe_log_Alu, aes(x = freqeuncy, y=mean_size,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "TE length",
#        x = "Allele frequency"
#   )+
#   # geom_jitter(color="black", size=0.2, alpha=0.2)+
#   ggtitle("Alu")
Alu_size <- ggplot()+
  geom_boxplot(TE_dataframe_log_Alu,mapping=aes(x= freqeuncy,y= mean_size,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "TE size (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# Alu_disttss <- ggplot(TE_dataframe_log_Alu,aes(x = freqeuncy, y=mean_disttss,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "Distance to closest gene",
#      x = "Allele frequency"
#   ) # +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
Alu_disttss <- ggplot()+
  geom_boxplot(TE_dataframe_log_Alu,mapping=aes(x= freqeuncy,y= mean_disttss,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Gene distance (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# Alu_recrate <-ggplot(TE_dataframe_log_Alu,aes(x = freqeuncy, y=mean_recrate,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette= "Dark1")+
#   labs(y = "Recombination rate",
#      x = "Allele frequency"
#   )# +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
Alu_recrate <- ggplot()+
  geom_boxplot(TE_dataframe_log_Alu,mapping=aes(x= freqeuncy,y= mean_recrate,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Recombination rate (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  ggtitle("Alu")+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

## L1

# L1_size <- ggplot(TE_dataframe_log_L1, aes(x = freqeuncy, y=mean_size,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "TE length",
#        x = "Allele frequency"
#   )+
#   # geom_jitter(color="black", size=0.2, alpha=0.2)+
#   ggtitle("L1")
L1_size <- ggplot()+
  geom_boxplot(TE_dataframe_log_L1,mapping=aes(x= freqeuncy,y= mean_size,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "TE size (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# L1_disttss <- ggplot(TE_dataframe_log_L1,aes(x = freqeuncy, y=mean_disttss,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette = "Dark1")+
#   labs(y = "Distance to closest gene",
#      x = "Allele frequency"
#   )# +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
L1_disttss <- ggplot()+
  geom_boxplot(TE_dataframe_log_L1,mapping=aes(x= freqeuncy,y= mean_disttss,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Gene distance (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

# L1_recrate <-ggplot(TE_dataframe_log_L1,aes(x = freqeuncy, y=mean_recrate,fill=freqeuncy)) + 
#   geom_lv(color='black', size=0.75) + 
#   geom_boxplot(outlier.alpha = 0, coef=0, fill="#00000000") +
#   scale_fill_brewer(palette= "Dark1")+
#   labs(y = "Recombination rate",
#      x = "Allele frequency"
#   )# +
#   # geom_jitter(color="black", size=0.2, alpha=0.2)
L1_recrate <- ggplot()+
  geom_boxplot(TE_dataframe_log_L1,mapping=aes(x= freqeuncy,y= mean_recrate,fill=freqeuncy))+
  scale_fill_brewer(palette = "Dark1")+
  labs(y = "Recombination rate (log10)",
       x = "Allele frequency",fill = "Frequency ranges"
  )+   ggtitle("L1")+
  scale_fill_manual(values=color_class)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

combined <- allTE_recrate + Alu_recrate + L1_recrate + SVA_recrate +
  allTE_disttss + Alu_disttss + L1_disttss + SVA_disttss +
  allTE_size + Alu_size + L1_size + SVA_size &
  theme(legend.position = "bottom")
combined + plot_layout(ncol=4, nrow=3, guides = "collect",
                       axis_titles = "collect")

# Scatter plots for recombination rate, TE size, and gene distance -------------------------------------------------

# color_class <- c("#fc4e2a","#fd8d3c","#feb24c", "#fed976", "#fff3d1")

a <- ggplot()+
  geom_point(TE_dataframe,mapping=aes(y=mean_TEfreq,x=mean_size,fill=freqeuncy),shape=21,color = "black")+
  scale_fill_manual(values = color_class)+
  labs(y="TE size",
       x= "Allele frequency",
       fill= "Frequency ranges")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))
b <- ggplot()+
  geom_point(TE_dataframe,mapping=aes(y=mean_TEfreq,x=mean_disttss,fill=freqeuncy),shape=21,color = "black")+
  scale_fill_manual(values = color_class)+
  labs(y="Gene distance",
       x= "Allele frequency",
       fill= "Frequency ranges")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))
c <- ggplot()+
  geom_point(TE_dataframe,mapping=aes(y=mean_TEfreq,x=mean_recrate,fill=freqeuncy),shape=21,color = "black")+
  scale_fill_manual(values = color_class)+
  labs(x="Recombination rate",
       y= "TE frequency",
       fill= "Frequency ranges")+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  # theme(
  #   plot.title = element_text(size = 16, face = "bold"),
  #   axis.title.x = element_text(size = 14, face = "bold"),
  #   axis.title.y = element_text(size = 14, face = "bold"))

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig1_recrate.png"), plot=ggplot2::last_plot())

combined <- a + b + c & theme(legend.position = "bottom")
combined + plot_layout(ncol=3, nrow=3, guides = "collect",
                       axis_titles = "collect")

combined <- c + b + a +
  allTE_recrate + allTE_disttss + allTE_size +
  Alu_recrate + Alu_disttss + Alu_size + 
  L1_recrate + L1_disttss + L1_size + 
  SVA_recrate + SVA_disttss + SVA_size & theme(legend.position = "bottom")
combined + plot_layout(ncol=3, nrow=5, guides = "collect",
                       axis_titles = "collect")

# Stats and background intervals for recombination rate, TE size, and gene distance ----------------

## correlations

plot(TE_dataframe$mean_size,TE_dataframe$mean_TEfreq)
plot(TE_dataframe$mean_disttss,TE_dataframe$mean_TEfreq)
plot(TE_dataframe$mean_recrate,TE_dataframe$mean_TEfreq)

cor(TE_dataframe$mean_size,TE_dataframe$mean_TEfreq, method = "spearman")
cor(TE_dataframe$mean_disttss,TE_dataframe$mean_TEfreq, method = "spearman")
cor(TE_dataframe$mean_recrate,TE_dataframe$mean_TEfreq, method = "spearman")

cor.test(TE_dataframe$mean_size,TE_dataframe$mean_TEfreq, method = "spearman")
cor.test(TE_dataframe$mean_disttss,TE_dataframe$mean_TEfreq, method = "spearman")
cor.test(TE_dataframe$mean_recrate,TE_dataframe$mean_TEfreq, method = "spearman")

qqnorm(TE_dataframe$mean_TEfreq)
qqline(TE_dataframe$mean_TEfreq, col = "red")
shapiro.test(TE_dataframe$mean_TEfreq)

pairs(~mean_TEfreq+mean_disttss+mean_recrate+mean_size,data=TE_dataframe,
      pch = 21,
      diag.panel=NULL,
      upper.panel = panel.smooth,
      bg = TE_dataframe$sum_anno,
      main = "All TE types"
)

## stats for low/high rec rates regions

TE_dataframe 

add_recrate_category <- function(data, threshold = 5) {
  data <- mutate(data,
                 recrate_category = ifelse(mean_recrate > threshold, "high_recrate", "low_recrate"))
  return(data)
}
TE_dataframe_recratecat <- add_recrate_category(TE_dataframe, threshold = 5)

color_class <- c("#fc4e2a","#feb24c", "#fff3d1")
color_class <- c("#fff3d1","#feb24c", "#fc4e2a")
freq_order <- fct_relevel(TE_dataframe_recratecat$new_class ,"Rare","Polymorphic","Fixed")
# color_class <- c("#fff3d1","#fed976","#feb24c","#fd8d3c", "#fc4e2a")
# freq_order <- fct_relevel(TE_dataframe_recratecat$new_class ,"Very Rare","Rare","Common","Major","Fixed")
recrate_order <- fct_relevel(TE_dataframe_recratecat$recrate_category ,"low_recrate","high_recrate")

ggplot()+
  geom_boxplot(TE_dataframe_recratecat,mapping=aes(x= freq_order,y=log10(mean_recrate),fill=recrate_order),outlier.shape = "|",outlier.size = 3)+
  labs(y="log10(Recombination rate)",
       x= "TE frequency",
       fill= "Recombination rate")+
  scale_fill_manual(values = c("#4B878BFF","#D01C1FFF"))+
  theme_classic()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_recrate_cat.png"), plot=ggplot2::last_plot())

TE_recrate <- table(TE_dataframe_recratecat$recrate_category,TE_dataframe_recratecat$new_class) %>% as.data.frame()
colnames(TE_recrate) <- c("recrate_category","new_class","count")
freq_order <- fct_relevel(TE_recrate$new_class ,"Rare","Polymorphic","Fixed")
# freq_order <- fct_relevel(TE_recrate$new_class ,"Very Rare","Rare","Common","Major","Fixed")
recrate_order <- fct_relevel(TE_recrate$recrate_category ,"low_recrate","high_recrate")
figsupp_recrate_cat_2 <- ggplot(TE_recrate,mapping=aes(x=freq_order,y=count,fill=recrate_category))+
  geom_bar(position="fill", stat="identity")+
  labs(y="Count",
       x= "TE frequency",
       fill= "Recombination rate")+
  scale_fill_manual(values = c("#D01C1FFF","#4B878BFF"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_recrate_cat_2.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_recrate_cat_2.png"),
       plot = figsupp_recrate_cat_2, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

low_df <- TE_recrate %>% filter(recrate_category == "low_recrate")
high_df <- TE_recrate %>% filter(recrate_category == "high_recrate")
categories <- unique(TE_recrate$new_class)
results <- lapply(categories, function(category) {
  low_count <- low_df$count[low_df$new_class == category]
  high_count <- high_df$count[high_df$new_class == category]
  chi_squared <- chisq.test(c(low_count, high_count))
  list(
    category = category,
    chi_squared = chi_squared$statistic,
    df = chi_squared$parameter,
    p_value = chi_squared$p.value
  )
})
for (result in results) {
  cat("Category:", result$category, "\n")
  cat("Chi-squared value:", result$chi_squared, "\n")
  cat("Degrees of freedom:", result$df, "\n")
  cat("P-value:", result$p_value, "\n\n")}

# Category: 1 
# Chi-squared value: 57.81752 
# Degrees of freedom: 1 
# P-value: 2.87597e-14 
# Category: 2 
# Chi-squared value: 250.6738 
# Degrees of freedom: 1 
# P-value: 1.851562e-56 
# Category: 3 
# Chi-squared value: 1690.425 
# Degrees of freedom: 1 
# P-value: 0 

## percentage of counts in low rec region:

nrow(TE_dataframe_recratecat)
sum(TE_recrate$count)

((sum(low_df$count))/(sum(TE_recrate$count)))*100
((sum(high_df$count))/(sum(TE_recrate$count)))*100

# (3166/(3166+547))*100

## background intervals

background_hg38 <- read.table(file = file.path(sys_dir,"PanTE_human/recombination_rate/background_3713_hg38.bed"),header=F)
colnames(background_hg38) <- c("CHR","POS","END","id","size","strand")
background_hg38 <- background_hg38[order(background_hg38[,1],background_hg38[,2]),]
background_hg38_Gr <- makeGRangesFromDataFrame(df=background_hg38,
                                               keep.extra.columns=T,
                                               ignore.strand=F,
                                               seqnames.field="CHR",
                                               start.field="POS",
                                               end.field="END",
                                               starts.in.df.are.0based=FALSE)
set.seed(123)
background_hg38_job = submitGreatJob(background_hg38_Gr,species="hg38")
## Error in (function (type, msg, asError = TRUE) : Could not resolve host: great.stanford.edu
## Use old data
# rm(list = setdiff(ls(), c("background_hg38_job","background_hg38_distTSS","background_hg38_distTSS_filtered","background_hg38_distTSS_binned")))
plotRegionGeneAssociations(background_hg38_job)
background_hg38_distTSS <- data.frame(background_hg38_job@association_tables[["all"]][["chr"]],
                                      (background_hg38_job@association_tables[["all"]][["start"]]+1),
                                      abs(background_hg38_job@association_tables[["all"]][["distTSS"]]),
                                      background_hg38_job@association_tables[["all"]][["name"]])
colnames(background_hg38_distTSS) <- c("CHR","POS","distTSS","names")
min_indices <- tapply(seq_along(background_hg38_distTSS$distTSS), background_hg38_distTSS$names, function(x) x[which.min(background_hg38_distTSS$distTSS[x])])
background_hg38_distTSS_filtered <- background_hg38_distTSS[unlist(min_indices), ]

bin          <- 10000

df2          <- background_hg38_distTSS_filtered %>% 
  mutate(window = POS %/% bin) %>% 
  group_by(window,CHR) %>%
  summarise(across(c(distTSS), mean)) %>% 
  mutate(start = (window*bin)+1, stop=(window+1)*bin)
colnames(df2) <- c("window","CHR","mean_distTSS","start","stop")
df_out       <- c()
for(i in unique(df2$CHR)) {
  df2_tmp      <- df2 %>% filter(CHR == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(CHR    = i,
           mean_distTSS = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(CHR,window)
  df_out       <- bind_rows(df_out,df_tmp)
}
background_hg38_distTSS_binned <- df_out %>% arrange(CHR, window) %>% drop_na()

threshold <- 3
z_scores <- scale(rec_rate_binned_noNA$mean_recrate)
outliers <- which(abs(z_scores) > threshold)
rec_rate_binned_noNA_filterd <- rec_rate_binned_noNA[-outliers, ]
background_rec <- merge(rec_rate_binned_noNA,background_hg38_distTSS_binned)
background_rec <- merge(rec_rate_binned_noNA_filterd,background_hg38_distTSS_binned) %>% drop_na() %>% 
  dplyr::select(c("CHR","start","stop","window","mean_recrate","mean_distTSS"))
colnames(background_rec) <- (c("CHR","start","stop","window","mean_recrate_background","mean_distTSS"))

table(TE_dataframe$new_class)
# background_rec_less100 <- background_rec %>% filter(mean_recrate_background < 100)
# TE_dataframe_less100 <- TE_dataframe %>% filter(mean_recrate < 100)

plot(data = background_rec, mean_recrate_background ~ window, col="lightgray")
plot(data = TE_dataframe, mean_recrate ~ window, col="#EE353E")

ggplot()+
  geom_point(data=subset(background_rec,mean_recrate_background < 1000),mapping=aes(x=window,y=mean_recrate_background),shape=25,col="lightgray")+
  geom_point(data=subset(TE_dataframe,mean_recrate < 1000),mapping=aes(x=window,y=mean_recrate),shape=25,col="#EE353E",fill="#EE353E")+
  theme_classic()+
  labs(
    x="Genomic window",
    y="Recombination rate")

figsupp_recratebackground <- ggplot()+
  geom_point(data=subset(background_rec,mean_recrate_background < 50),mapping=aes(x=window,y=mean_recrate_background),shape=25,col="lightgray")+
  geom_point(data=subset(TE_dataframe,mean_recrate < 50),mapping=aes(x=window,y=mean_recrate),shape=25,col="#EE353E",fill="#EE353E")+
  theme_classic()+
  labs(
    x="Genomic window",
    y="Recombination rate")

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_recratebackground.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_recratebackground.png"),
       plot = figsupp_recratebackground, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

background_rec <- background_rec %>% mutate(Source = "Background")
TE_dataframe <- TE_dataframe %>% mutate(Source = "TE")

ggplot()+
  # geom_jitter(data = background_rec, aes(x = CHR, y = mean_recrate_background,group=Source), color = "#A2A2A1FF", width = 0.1, alpha = 0.1)+
  geom_boxplot(background_rec,mapping=aes(x=CHR,y=(mean_recrate_background),group=Source),fill=alpha('black',0.5),col="black",outlier.shape = "|")+
  labs(x="Genomic window",y="Recombination rate")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14, face = "bold",color="white"))
  
ggplot()+
  # geom_jitter(data = TE_dataframe, aes(x = CHR, y = mean_recrate,group=Source), color = "#A2A2A1FF", width = 0.1, alpha = 0.1)+
  geom_boxplot(TE_dataframe,mapping=aes(x=CHR,y=(mean_recrate),group=Source),fill=alpha('#EE353E',0.5),col="#EE353E",outlier.shape = "|")+
  labs(x="Genomic window",y="Recombination rate")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14, face = "bold",color="white"))

ggplot()+
  # geom_jitter(data = background_rec, aes(x = CHR, y = mean_recrate_background,group=Source), color = "#A2A2A1FF", width = 0.1, alpha = 0.1)+
  geom_boxplot(background_rec,mapping=aes(x=CHR,y=(mean_recrate_background),group=Source),fill=alpha('black',0.5),col="black",outlier.shape = "|")+
  geom_boxplot(TE_dataframe,mapping=aes(x=CHR,y=(mean_recrate),group=Source),fill=alpha('#EE353E',0.5),col="#EE353E",outlier.shape = "|")+
  labs(x="Genomic window",y="Recombination rate")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14, face = "bold",color="white"))

ggplot()+
  # geom_jitter(data = TE_dataframe, aes(x = CHR, y = mean_recrate,group=Source), color = "#A2A2A1FF", width = 0.1, alpha = 0.1)+
  geom_boxplot(TE_dataframe,mapping=aes(x=CHR,y=(mean_recrate),group=Source),fill=alpha('#EE353E',0.5),col="#EE353E",outlier.shape = "|")+
  labs(x="Genomic window",y="Recombination rate")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14, face = "bold",color="white"))

wilcox_test <- wilcox.test(TE_dataframe$mean_recrate,background_rec$mean_recrate_background)
print(wilcox_test)
length(TE_dataframe$mean_recrate)
length(background_rec$mean_recrate_background)

nrow(TE_dataframe)
nrow(background_rec)
mean(TE_dataframe$mean_recrate)
mean(background_rec$mean_recrate_background)
median(TE_dataframe$mean_recrate)
median(background_rec$mean_recrate_background)
plot(density(TE_dataframe$mean_recrate))
plot(density(background_rec$mean_recrate_background))

plot_data <- data.frame(
  recombination_rate = c(TE_dataframe$mean_recrate,
                         background_rec$mean_recrate_background),
  group = c(rep("TE windows", length(TE_dataframe$mean_recrate)),
            rep("Background windows", length(background_rec$mean_recrate_background))))

ggplot(plot_data, aes(x = group, y = log(recombination_rate))) +
  geom_boxplot() +
  geom_hline(yintercept = log(5), linetype = "dashed") +
  annotate("text", x = 1.5, y = max(plot_data$recombination_rate),
           label = paste0("Wilcoxon p = ", format.pval(wilcox_test$p.value))) +
  labs(x = "", y = "Mean recombination rate") +
  theme_classic()

# plot_data$group <- factor(plot_data$group)
# lm_res <- lm(recombination_rate ~ group, data = plot_data)
# summary(lm_res)

aggregate(recombination_rate ~ group, plot_data, median)
aggregate(recombination_rate ~ group, plot_data, mean)

# TE divergence in Global, AFR, and outofAFR -----------------------------------------------------------

## define kimura distance from RepeatMasker

panic_KD <- read.csv(file.path(sys_dir,"PanTE_human/evo_age/Kdistance_filtered.txt"),
                     header = F,
                     sep="",
                     dec = ".")
extract_last_part <- function(x) {
  parts <- strsplit(x, "_")[[1]]
  if (length(parts) > 1 && grepl("Kimura", x)) {
    return(as.numeric(parts[length(parts)]))
  } else {
    return(x)
  }
}
panic_KD$chr_div <- sapply(panic_KD$V1, extract_last_part)
chromosome_indices <- grep("chromosome", panic_KD$chr_div)
for (i in 1:(length(chromosome_indices) - 1)) {
  start_index <- chromosome_indices[i] + 1
  end_index <- chromosome_indices[i + 1] - 1
  if (end_index > start_index) {
    values <- as.numeric(panic_KD$chr_div[start_index:end_index])
    mean_value <- mean(values, na.rm = TRUE)
    panic_KD$chr_div[start_index:end_index] <- mean_value
  }
}
panic_KD <- panic_KD %>% dplyr::select(c(-1))
for (i in 1:(nrow(panic_KD) - 1)) {
  if (panic_KD$chr_div[i] == panic_KD$chr_div[i + 1]) {
    panic_KD$chr_div[i] <- NA
  }
}
panic_KD <- panic_KD %>% drop_na()
panic_KD <- panic_KD[-10361, , drop = FALSE]
panic_KD$div <- NA
panic_KD$chr <- NA
panic_KD$div[seq(2, nrow(panic_KD), by = 2)] <- panic_KD$chr_div[seq(2, nrow(panic_KD), by = 2)]
panic_KD$chr[seq(1, nrow(panic_KD), by = 2)] <- panic_KD$chr_div[seq(1, nrow(panic_KD), by = 2)]
panic_KD_div <- data.frame(panic_KD$div) %>% drop_na()
panic_KD_chr <- data.frame(panic_KD$chr) %>% drop_na()
length(panic_KD_div$panic_KD.div)
length(panic_KD_chr$panic_KD.chr)
panic_KD_2 <- data.frame(cbind(panic_KD_div$panic_KD.div,panic_KD_chr$panic_KD.chr))
colnames(panic_KD_2) = c("div","chr")
panic_KD_2$chr <- gsub("_", "\t", panic_KD_2$chr)
panic_KD_2[,3:5] <- stringr::str_split_fixed(panic_KD_2$chr, "\t", 3) 
panic_KD_2 <- panic_KD_2 %>% dplyr::select(c(-"chr"))
colnames(panic_KD_2) = c("div","CHR","POS","svid")
panic_KD_2$svid <- gsub("\t", ">", panic_KD_2$svid)
panic_KD_2$div <- as.numeric(unlist(panic_KD_2$div))
panic_KD_2$POS <- as.numeric(unlist(panic_KD_2$POS))
panic_KD_2$CHR <- gsub("chromosome","chr",panic_KD_2$CHR)
panic_KD_2 <- panic_KD_2 %>% dplyr::select(c("CHR","POS",-"svid","div"))
panic_KD_2 <- panic_KD_2[order(panic_KD_2[,1],panic_KD_2[,2]),]

TE_size_distTSS_div <- merge(TE_size_distTSS,panic_KD_2)
TE_size_distTSS_div$ANNO <- as.character(unlist(TE_size_distTSS_div$ANNO))
TE_size_distTSS_div$div <- as.numeric(TE_size_distTSS_div$div)

nrow(TE_size_distTSS_div)
nrow(TE_size_distTSS)
nrow(panic_KD_2)

## include AFR and outofAFR info
## use AFR_freq_3_v2026 and outofAFR_freq_3_v2026 instead of AFR_freq_3 and outofAFR_freq_3

colnames(AFR_freq_3)
colnames(AFR_freq_3_v2026)
colnames(outofAFR_freq_3)
colnames(outofAFR_freq_3_v2026)

AFR_freq_edited <- AFR_freq_3_v2026 %>% dplyr::select(c("CHR","POS","TE_FREQ"))
colnames(AFR_freq_edited) = c("CHR","POS","TE_FREQ.a")
outofAFR_freq_edited <- outofAFR_freq_3_v2026 %>% dplyr::select(c("CHR","POS","TE_FREQ"))
colnames(outofAFR_freq_edited) = c("CHR","POS","TE_FREQ.o")

TE_africa_outofafrica <- merge(AFR_freq_edited,outofAFR_freq_edited)
TE_size_distTSS_div_AO <- merge(TE_africa_outofafrica,TE_size_distTSS_div)

TE_size_distTSS_div_AO$ANNO <- as.character(unlist(TE_size_distTSS_div_AO$ANNO))
TE_size_distTSS_div_AO$div <- as.numeric(TE_size_distTSS_div_AO$div)

nrow(TE_africa_outofafrica)
nrow(TE_size_distTSS_div_AO)

## neutral evo age calculation for TE_size_distTSS_div_AO

TE_div_Ao_neu <- TE_size_distTSS_div_AO %>%
  filter(TE_FREQ < 1)

meo <- (0.016 + 0.025 + 0.016)/3
TE_div_Ao_neu$div_age <- ((TE_div_Ao_neu$div)/(2*(meo)))
Ne <- 12500
x <- TE_div_Ao_neu$TE_FREQ
TE_div_Ao_neu$exp_age_ancestral <- (-4*Ne*(x/(1-x))*log(x))
Ne <- 24000
x <- TE_div_Ao_neu$TE_FREQ.a
TE_div_Ao_neu$exp_age_africa <- (-4*Ne*(x/(1-x))*log(x))
Ne <- 7700
x <- TE_div_Ao_neu$TE_FREQ.o
TE_div_Ao_neu$exp_age_outofafrica <- (-4*Ne*(x/(1-x))*log(x))

get_custom_color <- function(value) {
  if (value > +3) {
    return("Old")
  } else if (value < -3) {
    return("Young")
  } else {
    return("Neutral")  
  }
}

TE_neu_L1AluSVA <- TE_div_Ao_neu %>% 
  filter(ANNO == "SINE/Alu" | ANNO == "LINE/L1" | ANNO == "Retroposon/SVA") # %>% filter(div > 0)
TE_neu_L1AluSVA$diff_age <- scale(scale(TE_neu_L1AluSVA$div_age))-(scale(TE_neu_L1AluSVA$exp_age_ancestral))
percentile_99 <- quantile(TE_neu_L1AluSVA$diff_age, 0.99)
percentile_01 <- quantile(TE_neu_L1AluSVA$diff_age, 0.01)
TE_neu_L1AluSVA$custom_color <- sapply(TE_neu_L1AluSVA$diff_age, get_custom_color)

ggplot()+
  geom_point(data= TE_neu_L1AluSVA,mapping=aes(y=diff_age,x=POS,col=custom_color),size = 1)+
  labs(
    y="Scaled age estimate difference",
    col = "TE age") +
  scale_color_manual(values = c("Old" = "orange", "Young" = "steelblue", "Neutral" = "gray")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_evoagediff.png"), plot=ggplot2::last_plot())

## TE divergence/freq plots - similarity based AFR and outofAFR

anno_order <- c("SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")
new_colors <- c("#440154FF","#453781FF","#39558CFF","#238A8DFF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF")

TE_div_Ao_neu$anno_order <- ordered(TE_div_Ao_neu$ANNO,levels = anno_order)
TE_div_Ao_neu$diff_age <- scale(scale(TE_div_Ao_neu$div_age))-(scale(TE_div_Ao_neu$exp_age_ancestral))
TE_div_Ao_neu$custom_color <- sapply(TE_div_Ao_neu$diff_age, get_custom_color)

plot_all <- ggplot()+
  geom_point(TE_div_Ao_neu, mapping=aes(x=scale(div_age),y=TE_FREQ,col=custom_color))+
  geom_line(TE_div_Ao_neu, mapping=aes(x=scale(exp_age_ancestral),y=TE_FREQ,col="black"))+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))
plot_all <- plot_all + scale_colour_manual(values=c("black","grey","orange","steelblue"))
plot_all <- plot_all + labs(x="Divergence age (scaled)",y="Global TE frequency",col="TE type")
plot_all <- plot_all + theme_linedraw()
plot_all +  scale_size(trans ="reverse")
  # theme(
  #   plot.title = element_text(size = 16, face = "bold"),
  #   axis.title.x = element_text(size = 14, face = "bold"),
  #   axis.title.y = element_text(size = 14, face = "bold"))

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_freqdiv.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_freqdiv.png"),
       plot = figsupp_freqdiv, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

## include TE freq ranges to TE_div_Ao_neu

TE_div_Ao_neu_freq <- TE_div_Ao_neu

# TE_div_Ao_neu_freq$new_class <- lapply(TE_div_Ao_neu_freq$TE_FREQ, function(x) if(x>0.901){
#   TE_div_Ao_neu_freq$new_class="Fixed"
# }else if (0.634 < x && x <= 0.901){
#   TE_div_Ao_neu_freq$new_class="Major"
# }else if (0.367 < x && x <= 0.634){
#   TE_div_Ao_neu_freq$new_class="Common"
# }else if ((0.10) < x && x <= 0.367){
#   TE_div_Ao_neu_freq$new_class="Rare"
# }else 
#   TE_div_Ao_neu_freq$new_class="Very Rare")

TE_div_Ao_neu_freq <- TE_div_Ao_neu_freq %>%
  mutate(
    new_class = case_when(
      TE_FREQ >= 0.95                   ~ "Fixed",
      TE_FREQ >= 0.05 & TE_FREQ < 0.95  ~ "Polymorphic",
      TE_FREQ < 0.05                    ~ "Rare",
      TRUE ~ NA_character_))

TE_div_Ao_neu_freq$new_class <- as.character(unlist(TE_div_Ao_neu_freq$new_class))
table(TE_div_Ao_neu_freq$new_class)
# color_class <- c("#fff3d1","#fed976","#feb24c", "#fd8d3c", "#fc4e2a")
# class_order <- fct_relevel(TE_div_Ao_neu_freq$new_class,"Very Rare","Rare","Common","Major","Fixed")
class_order <- fct_relevel(TE_div_Ao_neu_freq$new_class,"Rare","Polymorphic","Fixed")

ggplot()+
  geom_point(TE_div_Ao_neu_freq, mapping=aes(x=scale(div_age),y=TE_FREQ,col=class_order))+
  geom_line(TE_div_Ao_neu_freq, mapping=aes(x=scale(exp_age_ancestral),y=TE_FREQ,col="black"))+
  scale_colour_manual(values=c(color_class,"black"))+
  labs(x="Divergence age (scaled)",y="Global TE frequency",col="TE type")+
  theme_linedraw()+
  scale_size(trans ="reverse")

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_freqdiv_freqrange.png"), plot=ggplot2::last_plot())

## include candi_all_2 from the following selection scans sections to label candidate TE loci under positive selection
## include TE_Global_shared from the following venn diagram sections to extract shared TE loci

## candi_all_2 v2026 should be used with 154 sites

TE_candi_all <- TE_div_Ao_neu

colnames(TE_candi_all)
nrow(TE_candi_all)

TE_candi_all$candilabel[TE_candi_all$POS %in% candi_all_2$POS] <- "candidate"
TE_candi_all$candilabel[is.na(TE_candi_all$candilabel)] <- "not_candidate"

TE_candi_all$candilabel_2[TE_candi_all$POS %in% merged_ohana_fst_xpehh_filtered$pos] <- "candidate_v2"
TE_candi_all$candilabel_2[is.na(TE_candi_all$candilabel_2)] <- "not_candidate_v2"

TE_candi_all$globalsharedlabel[TE_candi_all$POS %in% TE_Global_shared$POS] <- "globalshared"
TE_candi_all$globalsharedlabel[is.na(TE_candi_all$globalsharedlabel)] <- "not_globalshared"

plot_africa <- ggplot()+
  geom_point(data=subset(TE_candi_all,globalsharedlabel=="globalshared"), mapping=aes(x=scale(div_age),y=TE_FREQ.a,col=candilabel),size=1.5) +
  geom_line(data=subset(TE_candi_all,globalsharedlabel=="globalshared"), mapping=aes(x=scale(exp_age_ancestral),y=TE_FREQ,col="black"))+
  # scale_colour_manual(values=c(new_colors,"black"))+
  scale_color_manual(values=c("black","#EE353E","#F5F5F5"))+
  labs(x="Divergence age (scaled)",y="African TE frequency",col="TE type")+
  theme_linedraw()+
  theme_few()+
  scale_size(trans ="reverse")
  # theme(panel.border = element_rect(color = "#00539CFF", size = 1))

plot_outofafrica <- ggplot()+
  geom_point(data=subset(TE_candi_all,globalsharedlabel=="globalshared"), mapping=aes(x=scale(div_age),y=TE_FREQ.o,col=candilabel),size=1.5) +
  geom_line(data=subset(TE_candi_all,globalsharedlabel=="globalshared"), mapping=aes(x=scale(exp_age_ancestral),y=TE_FREQ,col="black"))+
  # scale_colour_manual(values=c(new_colors,"black"))+
  scale_color_manual(values=c("black","#EE353E","#F5F5F5"))+
  labs(x="Divergence age (scaled)",y="Non-African TE frequency",col="TE type")+
  theme_linedraw()+
  theme_few()+
  scale_size(trans ="reverse")
  # theme(panel.border = element_rect(color = "#EEA47FFF", size = 1))

combined <- plot_africa + plot_outofafrica &
  theme(legend.position = "bottom")
combined + plot_layout(ncol=2, nrow=1, guides = "collect",
                       axis_titles = "collect")

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_freqdiv_perpop.png"), plot=ggplot2::last_plot())

## HumGenom Review

TE_candi_all_v2026 <- TE_div_Ao_neu

TE_candi_all_v2026 <- merge(
  TE_candi_all_v2026,
  TE_freq_anno_selscan_rec_rGREAT_master_154[, c("CHR", "POS", "selscans", "selection_res_comb","LOC", "Gene.refGene")],
  by = c("CHR", "POS"),
  all.x = TRUE)

nrow(TE_div_Ao_neu)
nrow(TE_candi_all_v2026)

# candidate     --> selscans == T
# candidate_v2  --> selection_res_comb == T

TE_candi_all_v2026$selscans[is.na(TE_candi_all_v2026$selscans)] <- "none"
TE_candi_all_v2026$selection_res_comb[is.na(TE_candi_all_v2026$selection_res_comb)] <- "none"
TE_candi_all_v2026$LOC[is.na(TE_candi_all_v2026$LOC)] <- "none"
TE_candi_all_v2026$Gene.refGene[is.na(TE_candi_all_v2026$Gene.refGene)] <- "none"

TE_candi_all_v2026$globalsharedlabel <- ifelse(
  TE_candi_all_v2026$POS %in% TE_Global_shared$POS,
  "globalshared", "not_globalshared")

TE_candi_all_v2026$candilabel <- ifelse(
  TE_candi_all_v2026$selscans != "none",
  "candidate", "not_candidate")

colnames(TE_candi_all_v2026)
table(TE_candi_all_v2026$globalsharedlabel)
table(TE_candi_all_v2026$selscans)
table(TE_candi_all_v2026$candilabel)
table(TE_candi_all_v2026$selection_res_comb)
table(TE_candi_all_v2026$LOC)

plot_africa_v2026 <- ggplot()+
  geom_point(data=subset(TE_candi_all_v2026,globalsharedlabel == "globalshared" & candilabel == "not_candidate"),
             mapping=aes(x=scale(div_age),y=TE_FREQ.a,col=candilabel),size=1.5) +
  geom_point(data=subset(TE_candi_all_v2026, globalsharedlabel == "globalshared" & candilabel == "candidate"),
             mapping=aes(x=scale(div_age),y=TE_FREQ.a,col=candilabel),size=1.5) +
  geom_line(data=subset(TE_candi_all_v2026,globalsharedlabel=="globalshared"),
            mapping=aes(x=scale(exp_age_ancestral),y=TE_FREQ), color = "black")+
  geom_text_repel(data=subset(TE_candi_all_v2026, globalsharedlabel == "globalshared" &
                                candilabel == "candidate" & LOC == "Introns"
  ), aes(x=scale(div_age), y=TE_FREQ.a,label=Gene.refGene),
  color = "black", size = 3, box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf) +
  scale_color_manual(values=c("#EE353E","grey70"))+
  labs(x="Divergence age (scaled)",y="African TE frequency",col="TE type")+
  theme_linedraw()+
  theme_few()+
  scale_size(trans ="reverse") + 
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5))

plot_outofafrica_v2026 <- ggplot()+
  geom_point(data=subset(TE_candi_all_v2026,globalsharedlabel == "globalshared" & candilabel == "not_candidate"),
             mapping=aes(x=scale(div_age),y=TE_FREQ.o,col=candilabel),size=1.5) +
  geom_point(data=subset(TE_candi_all_v2026, globalsharedlabel == "globalshared" & candilabel == "candidate"),
             mapping=aes(x=scale(div_age),y=TE_FREQ.o,col=candilabel),size=1.5) +
  geom_line(data=subset(TE_candi_all_v2026,globalsharedlabel=="globalshared"),
            mapping=aes(x=scale(exp_age_ancestral),y=TE_FREQ), color = "black")+
  geom_text_repel(data=subset(TE_candi_all_v2026, globalsharedlabel == "globalshared" &
                                candilabel == "candidate" & LOC == "Introns"
  ), aes(x=scale(div_age), y=TE_FREQ.o,label=Gene.refGene),
  color = "black", size = 3, box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf) +
  scale_color_manual(values=c("#EE353E","grey70"))+
  labs(x="Divergence age (scaled)",y="Non-African TE frequency",col="TE type")+
  theme_linedraw()+
  theme_few()+
  scale_size(trans ="reverse") + 
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5))

combined_v2026 <- plot_africa_v2026 + plot_outofafrica_v2026 &
  theme(legend.position = "bottom")
fig2e_freqdiv_perpop <- combined_v2026 +
  plot_layout(ncol=2, nrow=1, guides = "collect", axis_titles = "collect")

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_freqdiv_perpop.png"), plot=ggplot2::last_plot())

# ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "fig2e_freqdiv_perpop_nolabel.png"),
#        plot = fig2e_freqdiv_perpop, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "fig2e_freqdiv_perpop.png"),
       plot = fig2e_freqdiv_perpop, width = 20, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

## TE divergence/fst plots

Fst <- read.csv(file.path(sys_dir,"PanTE_human/selection_scans/fst/out.weir.fst"), header = T,sep = "")   
colnames(Fst) <- c("CHR","POS","Fst")
# Fst_q <- Fst %>% filter(Fst > 0)
# quantile(Fst_q$Fst, 0.99)
TE_div_Fst <- merge(TE_div_Ao_neu,Fst)

TE_div_Fst_filtered <- TE_div_Fst %>% filter(Fst > 0)
Fst_percentile_99 <- quantile(TE_div_Fst_filtered$Fst, 0.99)

get_custom_fst <- function(value_1,value_2) {
  if (value_1 > Fst_percentile_99 && value_2 == "Young")
  {return("Sig")}
  else {
    return("not_sig")
  }}
TE_div_Fst_filtered$custom_fst <- apply(TE_div_Fst_filtered, 1, function(row) get_custom_fst(row["Fst"],row["custom_color"]))

plot_fst <- ggplot()+
  geom_point(data=(TE_div_Fst_filtered), mapping=aes(x=scale(div_age),y=Fst,col=custom_fst))
  # geom_line(data=subset(TE_div_Fst_filtered, TE_FREQ < 0.5), mapping=aes(x=scale(exp_age_ancestral),y=TE_FREQ,col="black"))
plot_fst <- plot_fst + scale_colour_manual(values=c("grey","orange","steelblue"))
plot_fst <- plot_fst + scale_colour_manual(values=c("grey","#EE353E"))
plot_fst <- plot_fst + labs(x="Divergence age (scaled)",y="Fst",col="TE type")
plot_fst <- plot_fst + theme_linedraw()
plot_fst +  scale_size(trans ="reverse") 
  # theme(
  #   plot.title = element_text(size = 16, face = "bold"),
  #   axis.title.x = element_text(size = 14, face = "bold"),
  #   axis.title.y = element_text(size = 14, face = "bold"))

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_fstdivplot.png"), plot=ggplot2::last_plot())

## TE_div_Fst_filtered is a data frame that has gene distance, size, global, AFR, and outofAFR alt allele freq, annotation, divergence, neutral evo calculation, and Fst information for all TE vars 

candidate_TE_div_Fst <- TE_div_Fst_filtered %>%
  filter(custom_color == "Young") %>%
  filter(Fst > Fst_percentile_99)

# TE_FULL_INFO ------------------------------------------------------------

## create 10 000 bins for TE_div_Fst_filtered

TE_div_Fst_filtered$POS <- as.numeric(TE_div_Fst_filtered$POS)
TE_div_Fst$POS <- as.numeric(TE_div_Fst$POS)

bin          <- 10000

df2          <- TE_div_Fst %>% 
  mutate(window = POS %/% bin)
df2$window_CHR <- paste(df2$window, df2$CHR, sep = "_") 
df2          <- df2 %>%
  group_by(window_CHR) %>%
  summarise(sum_anno = list(ANNO),across(c(TE_FREQ,TE_FREQ.a,TE_FREQ.o,distTSS,SIZE,div,div_age,exp_age_ancestral,diff_age,Fst), mean))
df2[,13:14] <- stringr::str_split_fixed(df2$window_CHR, "_", 2)
df2 <- df2[-c(1)] 
colnames(df2) <- c("sum_anno","mean_TEfreq","TE_FREQ.a","TE_FREQ.o","mean_disttss","mean_size","div","div_age","exp_age_ancestral","diff_age","Fst","window","CHR")
df2$window <- as.numeric(df2$window)
df2          <- df2 %>%
  mutate(start = (window*bin)+1, stop=(window+1)*bin) %>% 
  relocate(c("window","CHR","start","stop","sum_anno","mean_TEfreq","TE_FREQ.a","TE_FREQ.o","mean_disttss","mean_size","div","div_age","exp_age_ancestral","diff_age","Fst"))
df_out       <- c()
for(i in unique(df2$CHR)) {
  df2_tmp      <- df2 %>% filter(CHR == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(CHR    = i,
           mean_TEfreq = NA,
           TE_FREQ.a = NA,
           TE_FREQ.o = NA,
           mean_disttss = NA,
           mean_size = NA,
           div = NA,
           div_age = NA,
           exp_age_ancestral = NA,
           diff_age = NA,
           Fst = NA,
           sum_anno = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(CHR,window)
  df_out       <- bind_rows(df_out,df_tmp)
}

TE_div_Fst_binned <- df_out %>% arrange(CHR, window) %>% drop_na()
colnames(TE_div_Fst_binned)
colnames(rec_rate_binned_noNA) = c("window","CHR","mean_recrate","pos","end")
TE_div_Fst_rec_rate <- merge(rec_rate_binned_noNA,TE_div_Fst_binned)
# TE_div_Fst_rec_rate$new_class <- lapply(TE_div_Fst_rec_rate$mean_TEfreq, function(x) if(x>0.901){
#   TE_div_Fst_rec_rate$new_class="Fixed"
# }else if (0.634 < x && x <= 0.901){
#   TE_div_Fst_rec_rate$new_class="Major"
# }else if (0.367 < x && x <= 0.634){
#   TE_div_Fst_rec_rate$new_class="Common"
# }else if ((0.10) < x && x <= 0.367){
#   TE_div_Fst_rec_rate$new_class="Rare"
# }else
#   TE_div_Fst_rec_rate$new_class="Very Rare")
# TE_div_Fst_rec_rate$new_class <- as.character(unlist(TE_div_Fst_rec_rate$new_class))
# table(TE_div_Fst_rec_rate$new_class)
TE_FULL_INFO <- TE_div_Fst_rec_rate %>% arrange(CHR,window) %>% 
  dplyr::select(c("CHR","start","stop","window","mean_recrate","mean_TEfreq","TE_FREQ.a","TE_FREQ.o","mean_disttss","mean_size",
                  "div","div_age","exp_age_ancestral","diff_age","Fst"))
TE_FULL_INFO <- TE_FULL_INFO[order(TE_FULL_INFO[,1],TE_FULL_INFO[,2]),]

## TE_FULL_INFO is a data frame that has ALL information for all TE vars ## check duplicated entries

# Figure 2 ----------------------------------------------------------------



# TE counts per individual ---------------------------------------------------------------

## define population info, major, singletons, and shared TEs

## major TE

population_info <- read.table(file = file.path(sys_dir,"PanTE_human/population_comp/hprc_year1_sample_metadata.txt"), header = F,
                              sep = "",
                              dec = ".")

major_TE_filtered_counted <- read.table(file = file.path(sys_dir,"PanTE_human/population_comp/major_TE_stat"), header = F,
                                        sep = "",
                                        dec = ".")
major_filtered_counted_1 <- subset(major_TE_filtered_counted, select=(c(3,9)))
colnames(major_filtered_counted_1) <- c("V1","V2")
maj_pop <- right_join(major_filtered_counted_1, population_info, by='V1') %>%
  arrange(desc(V3),(V2.x)) %>%
  mutate(V1 = factor(V1, levels = V1))
ggplot(maj_pop,aes(x = V1, y = V2.x, fill = V3)) +
  geom_col()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  xlab("ind")+
  ylab("variant counts")+
  ggtitle("major_TE_filtered_counted")+
  labs(fill = "superpop")

## singleton TE

singletons_TE_filtered_counted <- read.table(file = file.path(sys_dir,"PanTE_human/population_comp/R_singletons_TE_filtered.vcf.singletons"), header = F,
                                             sep = "",
                                             dec = ".")
col_order <- c("V2","V1")
singletons_TE_filtered_counted <- singletons_TE_filtered_counted[,col_order]
colnames(singletons_TE_filtered_counted) <- c("V1","V2")
sin_pop <- right_join(singletons_TE_filtered_counted, population_info, by='V1') %>%
  arrange(desc(V3),(V2.x)) %>%
  mutate(V1 = factor(V1, levels = V1))
ggplot(sin_pop,aes(x = V1, y = V2.x, fill = V3)) +
  geom_col()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  xlab("ind")+
  ylab("variant counts")+
  ggtitle("singletons_TE_filtered_counted")+
  labs(fill = "superpop")
ggplot(sin_pop,aes(x = reorder(V2.y,V2.x), y = V2.x, fill = V3)) +
  geom_boxplot()+
  theme_minimal()+
  geom_jitter(color="black", size=0.5, alpha=0.9)+
  xlab("subpop")+
  ylab("variant counts")+
  ggtitle("singletons_TE_filtered_counted")+
  labs(fill = "superpop")+
  scale_fill_manual(values=wes_palette(n=4, name="GrandBudapest1"))

sin_pop_AFR <- sin_pop %>%
  filter(V3 == "AFR") %>%
  drop_na()
nrow(sin_pop_AFR)
summary(sin_pop_AFR$V2.x)
sd(sin_pop_AFR$V2.x)

sin_pop_OUTofAFR <- sin_pop %>%
  filter(V3 != "AFR")
nrow(sin_pop_OUTofAFR)
summary(sin_pop_OUTofAFR$V2.x)
sd(sin_pop_OUTofAFR$V2.x)

## major + singletons + shared

colnames(sin_pop)=c("ind","singletons","subpop","superpop")
colnames(maj_pop)=c("ind","major","subpop","superpop")
maj_sin <- merge(maj_pop,sin_pop)
maj_sin_shar <- data.frame(append(maj_sin, c(V5=105), after=5))
## counts of shared var = 105. The number might change depending on the nonref TE var calculated in source

colnames(maj_sin_shar)=c("ind","subpop","superpop","Major allele frequency > 0.5","Singletons","Shared in all but GRCh38")
maj_sin_shar <- cbind(maj_sin_shar, Total = rowSums(maj_sin_shar[,c(4,5,6)]))
maj_sin_shar_tidy <- pivot_longer(maj_sin_shar, cols = (c(4,5,6)), values_drop_na = TRUE)
colnames(maj_sin_shar_tidy)=c("ind","subpop","superpop","Total","filter_type","filter_type_count")
maj_sin_shar_tidy <- maj_sin_shar_tidy %>% drop_na() %>% as.data.frame() 
maj_sin_shar_tidy <- maj_sin_shar_tidy %>% filter(maj_sin_shar_tidy$ind != "NA19240")
maj_sin_shar_tidy$filter_type <- as.factor(maj_sin_shar_tidy$filter_type)
maj_sin_shar_tidy$filter_type <- fct_relevel(maj_sin_shar_tidy$filter_type, "Major allele frequency > 0.5","Singletons","Shared in all but GRCh38")

ggplot(maj_sin_shar_tidy,aes(x = reorder(ind,filter_type_count), y = filter_type_count, fill = filter_type)) +
  geom_bar(position="stack", stat="identity",width=0.7)+
  # theme_set(theme_minimal(base_family = "Georgia"))+
  scale_y_continuous(breaks = seq(0, 465, by = 50)) +
  labs(fill = "TEs",x="Individuals",y="TE counts")+
  scale_fill_manual(values =  c("#5F9E9E","#272772","#B22224"))+
  theme_classic()+
  theme()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, vjust = 0.5))

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig2_countperind.png"), plot=ggplot2::last_plot())

# TE counts per population and Venn diagram ------------------------------------------------

## TE counts per population

AFR <- c("HG01891","HG02257","HG02486","HG02559","HG02572","HG03516","HG02622","HG02630","HG02717","HG02886","HG03453","HG03540","HG03579","HG02109","HG02145","HG02723","HG02818","HG03486","NA18906","NA20129","NA21309","HG02055","HG03098")
outofAFR <- c("HG01123","HG01258","HG01358","HG01361","HG00735","HG00741","HG01071","HG01106","HG01175","HG01928","HG01952","HG01978","HG02148","HG00733","HG01109","HG01243","HG00438","HG00621","HG00673","HG02080","HG03492")

TE_indv <- read.table(file = file.path(sys_dir,"PanTE_human/population_comp/TE_counts_1"),header = F,sep = "") %>%
  dplyr::select(c(3,9)) %>%
  filter(V3 != "CHM13")
colnames(TE_indv) <- c("indv","counts")

pop <- rep(NA, length(TE_indv$indv))
for (value in AFR) {pop[grep(value, TE_indv$indv)] <- "African"}
for (value in outofAFR) {pop[grep(value, TE_indv$indv)] <- "Non-African"}
TE_indv <- as.tibble(data.frame(TE_indv, pop)) %>% as.data.frame()

ggplot(TE_indv,aes(x = reorder(indv,counts), y = counts,fill=pop)) +
  geom_bar(position="stack", stat="identity",width=0.7)+
  scale_y_continuous(breaks = seq(0,700, by = 50)) +
  labs(fill="Population",x="Individuals",y="TE counts")+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
  )+
  scale_fill_manual(values =  c("#00539CFF","#EEA47FFF"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig2_countperindperpop.png"), plot=ggplot2::last_plot())

## TE counts + singletons

sin_pop_v2 <- sin_pop %>%
  mutate(super_pop_label = case_when(
    superpop %in% c("AMR","EAS", "SAS") ~ "Non-African",
    superpop %in% c("AFR")              ~ "African")) %>%
  dplyr::select(c("ind","singletons","super_pop_label"))
colnames(sin_pop_v2) <- c("indv","singletons","pop")

TE_indv_sin <- merge(sin_pop_v2,TE_indv)
TE_indv_sin$NOTsingletons <- TE_indv_sin$counts-TE_indv_sin$singletons

TE_indv_sin_tidy <- TE_indv_sin %>%
  mutate(NOTsingletons = counts - singletons) %>%
  pivot_longer(cols = c(NOTsingletons, singletons),
               names_to = "var_types", values_to = "n")

ggplot(TE_indv_sin_tidy, aes(x = reorder(indv, counts), y = n, fill = var_types)) +
  geom_col(width = 0.7) +
  scale_y_continuous(breaks = seq(0, 700, by = 50)) +
  scale_fill_manual(values = c("NOTsingletons" = "#5F9E9E", "singletons"="#B22224")) +
  labs(fill = "TEs", x = "Individuals", y = "TE counts") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(angle = 90, vjust = 0.5))

TE_indv_sin_tidy <- TE_indv_sin %>%
  mutate(NOTsingletons = counts - singletons) %>%
  pivot_longer(c(NOTsingletons, singletons),
               names_to = "var_types", values_to = "n")

fig2_TEcountsperindv_sing <- ggplot(TE_indv_sin_tidy, aes(x = reorder(indv, counts), y = n, fill = interaction(pop, var_types))) +
  geom_col(width = 0.7) +
  scale_y_continuous(breaks = seq(0, 700, by = 50)) +
  scale_fill_manual(
    breaks = c("African.NOTsingletons","African.singletons",
               "Non-African.NOTsingletons","Non-African.singletons"),
    labels = c("African – non-singletons","African – singletons",
               "Non-African – non-singletons","Non-African – singletons"),
    values = c(
      "African.NOTsingletons"      = "#00539CFF",
      "African.singletons"         = "#00539C80",  
      "Non-African.NOTsingletons"  = "#EEA47FFF",  
      "Non-African.singletons"     = "#EEA47F80")) +
  labs(fill = "TE components", x = "Individuals", y = "TE counts") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(angle = 90, vjust = 0.5),
    legend.position = "none")

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "fig2_TEcountsperindv_sing.png"),
       plot = fig2_TEcountsperindv_sing, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

ggplot(TE_indv_sin, aes(x = reorder(indv, counts), y = counts, fill = pop)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = singletons, y = counts - 10, colour  = pop), size = 2.5, fontface = "bold") +
  scale_y_continuous(breaks = seq(0,700,50)) +
  scale_fill_manual(values = c("African" = "#00539CFF", "Non-African" = "#EEA47FFF")) +
  scale_colour_manual(values = c("African" = "white", "Non-African" = "black")) +
  labs(fill = "Super Population", x = "Individuals", y = "TE counts") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x  = element_text(angle = 90, vjust = 0.5))

TE_indv_sin_v2 <- TE_indv_sin %>%
  mutate(xpos = as.numeric(reorder(indv, counts)),
         ysplit = singletons)

ggplot(TE_indv_sin_v2, aes(x = reorder(indv, counts), y = counts, fill = pop)) +
  geom_col(width = 0.7) +
  geom_segment(
    data = TE_indv_sin_v2,
    aes(x = xpos - 0.35, xend = xpos + 0.35, y = ysplit, yend = ysplit),
    inherit.aes = FALSE,
    colour = "white",
    linewidth = 1.2) +
  scale_fill_manual(values = c("African" = "#00539CFF", "Non-African" = "#EEA47FFF")) +
  scale_y_continuous(breaks = seq(0, 700, 50)) +
  labs(x = "Individuals", y = "TE counts", fill = "Population") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"))

## stats

TE_indv_AFR <- TE_indv %>% filter(pop=="African")
TE_indv_outofAFR <- TE_indv %>% filter(pop=="Non-African")
summary(TE_indv_AFR$counts)
sd(TE_indv_AFR$counts)
summary(TE_indv_outofAFR$counts)
sd(TE_indv_outofAFR$counts)

shapiro.test(TE_indv$counts[TE_indv$pop == "African"])
shapiro.test(TE_indv$counts[TE_indv$pop == "Non-African"])

wilcox_test_result <- wilcox.test(counts ~ pop, data = TE_indv)
print(wilcox_test_result)

p <- ggplot(TE_indv,aes(x = pop, y = counts, fill=pop)) +
  geom_boxplot()+
  geom_jitter(shape=16, position=position_jitter(0.2),color="black",alpha=0.3)+
  # geom_text(aes(x = 1.5, y = max(TE_indv$counts)+10, label = significance_level), size = 6)+
  scale_y_continuous(breaks = seq(0,700, by = 50)) +
  labs(fill="Population",x="Population",y="TE counts")+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )+
  # annotate("text", x = 1.5, y = 650, label = paste("Wilcoxon, p-value=", round(wilcox_test_result$p.value,8)), size = 4)+
  scale_fill_manual(values =  c("#00539CFF","#EEA47FFF"))+
  theme_classic()
p + stat_compare_means(method = "wilcox", label.y = 680)

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig2_countperpopboxplot.png"), plot=ggplot2::last_plot())

## Venn diagram 

# TE_africa_outofafrica_no_dup <- TE_africa_outofafrica %>% filter(!duplicated(.))
# nrow(TE_africa_outofafrica)
# nrow(TE_africa_outofafrica_no_dup)

TE_africa_outofafrica <- merge(AFR_freq_edited,outofAFR_freq_edited)
TE_africa_outofafrica$id <- seq(1,6407)

TE_OO <- TE_africa_outofafrica %>% filter(TE_FREQ.o == 0)
TE_AA <- TE_africa_outofafrica %>% filter(TE_FREQ.a == 0)
TE_O  <- TE_africa_outofafrica %>% filter(TE_FREQ.o > 0) 
TE_A  <- TE_africa_outofafrica %>% filter(TE_FREQ.a > 0) 
TE_Global <- TE_africa_outofafrica
TE_Global_shared <- TE_africa_outofafrica %>%
  filter(TE_FREQ.a != 0) %>%
  filter(TE_FREQ.o != 0) 
TE_Global_zero <- TE_africa_outofafrica %>%
  filter(TE_FREQ.a == 0 & TE_FREQ.o == 0) 

a <- length(TE_OO$id)
b <- length(TE_AA$id)
c <- length(TE_O$id)
d <- length(TE_A$id)
e <- length(TE_Global$id)
f <- length(TE_Global_shared$id)
g <- length(TE_Global_zero$id)

# diff_TE <- dplyr::anti_join(TE_Global_shared, TE_Global_shared_v2, by = "id")

length(TE_O$id) - length(TE_Global_shared$id)
length(TE_A$id) - length(TE_Global_shared$id)

length(TE_Global$id) - length(TE_Global_zero$id) ==
length(TE_O$id) - length(TE_Global_shared$id) +
  length(TE_A$id) - length(TE_Global_shared$id) +
  length(TE_Global_shared$id)

venn.plot <- venn.diagram(x= list(TE_A$id, TE_O$id),
                          filename = "fig2_vendia.png",
                          category.names = c("African","Non-African"),
                          output=TRUE,
                          col=c("#00539CFF", '#EEA47FFF'),
                          fill = c(alpha("#00539CFF",0.3), alpha('#EEA47FFF',0.3)),
                          fontface = "bold",
                          cat.fontface = "bold",
                          cat.cex = 0.4,
                          cat.pos = c(-27, 27),
                          cat.default.pos = "outer",
                          cex = 0.6,
                          # cat.pos = c(0, 0),
                          cat.dist = c(0.05, 0.05),
                          imagetype="png" ,
                          height = 600 , 
                          width = 600 , 
                          resolution = 300,
                          compression = "lzw",
                          print.mode = c("raw", "percent"),
                          sigdigs = 2)

# ggvenn(
#   list(African = TE_A$id, `Non-African` = TE_O$id),
#   fill_color = c("#00539CFF", "#EEA47FFF"),
#   stroke_color = "black",
#   stroke_size = 1.5,
#   set_name_size = 6,
#   text_size = 5,
#   show_percentage = TRUE)

# Allele frequency per population -----------------------------------------

TE_africa_outofafrica <- merge(AFR_freq_edited,outofAFR_freq_edited)
# TE_africa_outofafrica_no_dup <- do.call(rbind, lapply(split(TE_africa_outofafrica, TE_africa_outofafrica$CHR), function(group) group[!duplicated(group), ]))
# unique(length(TE_africa_outofafrica_no_dup$POS))
TE_africa_outofafrica$id <- seq(1,6407)
TE_africa_outofafrica$pop <- c("population")
TE_africa <- TE_africa_outofafrica %>% dplyr::select(c(1,2,3,5)) %>% filter(TE_FREQ.a > 0)
TE_africa$pop <- c("African")
TE_outofafrica <- TE_africa_outofafrica %>% dplyr::select(c(1,2,4,5)) %>% filter(TE_FREQ.o > 0)
TE_outofafrica$pop <- c("Non-African")

fig2_freqperpopbox <- ggplot()+
  geom_jitter(data = TE_outofafrica, aes(x = pop, y = TE_FREQ.o), color = "#A2A2A1FF", width = 0.1, alpha = 0.1)+
  geom_jitter(data = TE_africa, aes(x = pop, y = TE_FREQ.a), color = "#A2A2A1FF", width = 0.1, alpha = 0.1)+
  geom_boxplot(TE_outofafrica,mapping=aes(x=pop,y=(TE_FREQ.o)),fill=alpha('#EEA47FFF',0.5),col="#EEA47FFF",outlier.shape = "|")+
  geom_boxplot(TE_africa,mapping=aes(x=pop,y=(TE_FREQ.a)),fill=alpha('#00539CFF',0.5),col="#00539CFF",outlier.shape = "|")+
  labs(x="",y="TE frequency")+
  theme_classic()+
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold",color="black"))

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig2_freqperpopbox.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "fig2_freqperpopbox.png"),
       plot = fig2_freqperpopbox, width = 12, height = 10, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

ggplot() +
  geom_violin(data = TE_africa,
              aes(pop, TE_FREQ.a, fill = "African", colour = "African"),
              alpha = 0.25, width = 0.95, size = 0.5, trim = TRUE) +
  geom_violin(data = TE_outofafrica,
              aes(pop, TE_FREQ.o, fill = "Non-African", colour = "Non-African"),
              alpha = 0.25, width = 0.95, size = 0.5, trim = TRUE) +
  # geom_boxplot(data = TE_africa,
  #              aes(pop, TE_FREQ.a, colour = "African"),
  #              width = 0.12, outlier.shape = NA, fill = NA, size = 0.6) +
  # geom_boxplot(data = TE_outofafrica,
  #              aes(pop, TE_FREQ.o, colour = "Non-African"),
  #              width = 0.12, outlier.shape = NA, fill = NA, size = 0.6) +
  geom_sina(data = TE_africa,
            aes(pop, TE_FREQ.a, colour = "African"),
            maxwidth = 0.7, alpha = 0.5, size = 0.8) +
  geom_sina(data = TE_outofafrica,
            aes(pop, TE_FREQ.o, colour = "Non-African"),
            maxwidth = 0.7, alpha = 0.5, size = 0.8) +
  scale_fill_manual(values = c("African" = "#00539CFF", "Non-African" = "#EEA47FFF")) +
  scale_colour_manual(values = c("African" = "#00539CFF", "Non-African" = "#EEA47FFF")) +
  labs(x = "", y = "TE frequency") +
  # scale_y_log10() +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 14, face = "bold", colour = "black"))

# TE family distribution per population - pie charts ----------------------

TE_africa_outofafrica <- merge(AFR_freq_edited,outofAFR_freq_edited)
TE_africa_outofafrica$id <- seq(1,6407)
TE_africa_outofafrica$pop <- c("population")
TE_africa_outofafrica_anno <- merge(TE_africa_outofafrica,TE_anno) %>%
  dplyr::select(c("CHR","POS","TE_FREQ.a","TE_FREQ.o","id","ANNO"))
TE_africa_anno <- TE_africa_outofafrica_anno %>% dplyr::select(c(1,2,3,5,6)) %>% filter(TE_FREQ.a > 0)
TE_africa_anno$pop <- c("Africa")
TE_outofafrica_anno <- TE_africa_outofafrica_anno %>% dplyr::select(c(1,2,4,5,6)) %>% filter(TE_FREQ.o > 0)
TE_outofafrica_anno$pop <- c("Out of Africa")

TE_africa_anno_counted <- TE_africa_anno %>% dplyr::count(pop,ANNO)
colnames(TE_africa_anno_counted) = c("pop","ANNO","NUMBER")
TE_outofafrica_anno_counted <- TE_outofafrica_anno %>% dplyr::count(pop,ANNO)
colnames(TE_outofafrica_anno_counted) = c("pop","ANNO","NUMBER")

AFRpie <- ggplot(TE_africa_anno_counted, aes(x="", y=NUMBER,
                                             fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("AFR")

outofAFRpie <- ggplot(TE_outofafrica_anno_counted, aes(x="", y=NUMBER,
                                                       fill=fct_relevel(ANNO,"LINE/L1","LINE/misc","SINE/Alu","SINE/MIR","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")))+
  geom_bar(stat="identity", width=1,color = "black") +
  coord_polar("y", start=0,)+
  theme_void()+
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#39558CFF","#238A8DFF","#440154FF","#453781FF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  ggtitle("outofAFR")

combined <- AFRpie+outofAFRpie & theme(legend.position = "right")
combined + plot_layout(ncol=2, nrow=1, guides = "collect", axis_titles = "collect", tag_level="keep")

TE_africaoutofafrica_anno_counted <- rbind(TE_outofafrica_anno_counted,TE_africa_anno_counted)

ggplot(TE_africaoutofafrica_anno_counted, mapping=aes(x=pop, y=NUMBER,
                                        fill=fct_relevel(ANNO,"SINE/Alu","SINE/MIR","LINE/L1","LINE/misc","LTR/ERV","LTR/misc","Retroposon/SVA","DNA/misc")))+
  geom_bar(position="fill", stat="identity") +
  scale_fill_viridis(discrete = TRUE) +
  scale_fill_manual(values=c("#440154FF","#453781FF","#39558CFF","#238A8DFF","#29AF7FFF","#74D055FF","#B8DE29FF","#FDE725FF"))+
  labs(fill="TE type")+
  theme_classic()

# PCA ---------------------------------------------------------------------

## define variables and plink output files

pca <- read_table(file.path(sys_dir,"PanTE_human/pca/TE_pca.eigenvec"), col_names = FALSE)
eigenval <- scan(file.path(sys_dir,"PanTE_human/pca/TE_pca.eigenval"))

AFR <- c("HG01891","HG02257","HG02486","HG02559","HG02572","HG03516","HG02622","HG02630","HG02717","HG02886","HG03453","HG03540","HG03579","HG02109","HG02145","HG02723","HG02818","HG03486","NA18906","NA20129","NA21309","HG02055","HG03098")
outofAFR <- c("HG01123","HG01258","HG01358","HG01361","HG00735","HG00741","HG01071","HG01106","HG01175","HG01928","HG01952","HG01978","HG02148","HG00733","HG01109","HG01243","HG00438","HG00621","HG00673","HG02080","HG03492")
pop_info <- read.table(file = file.path(sys_dir,"PanTE_human/frequency_distribution/hprc_year1_sample_metadata.txt"), header = F)
colnames(pop_info) <- c("ind","subpop","pop")
subpop_vectors <- list()
unique_subpops <- unique(pop_info$subpop)
for (subpop in unique_subpops) {
  subpop_inds <- pop_info[pop_info$subpop == subpop, "ind"]
  subpop_vectors[[subpop]] <- subpop_inds}

## housekeeping for pca variables

pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))
# spp <- rep(NA, length(pca$ind))
# spp[grep("", pca$ind)] <- "human"
pop <- rep(NA, length(pca$ind))
pop[grep("CHM13", pca$ind)] <- "CHM13"
for (value in AFR) {pop[grep(value, pca$ind)] <- "African"}            ## replace with African
for (value in outofAFR) {pop[grep(value, pca$ind)] <- "Non-African"}   ## replace with Non-African
# spp_pop <- paste0(spp, "_", pop)
subpop <- rep(NA, length(pca$ind))
subpop[grep("CHM13", pca$ind)] <- "CHM13"
for (value in  subpop_vectors[["CLM"]] ) {subpop[grep(value, pca$ind)] <- "CLM"}
for (value in  subpop_vectors[["ACB"]] ) {subpop[grep(value, pca$ind)] <- "ACB"}
for (value in  subpop_vectors[["GWD"]] ) {subpop[grep(value, pca$ind)] <- "GWD"}
for (value in  subpop_vectors[["ESN"]] ) {subpop[grep(value, pca$ind)] <- "ESN"}
for (value in  subpop_vectors[["CHS"]] ) {subpop[grep(value, pca$ind)] <- "CHS"}
for (value in  subpop_vectors[["PUR"]] ) {subpop[grep(value, pca$ind)] <- "PUR"}
for (value in  subpop_vectors[["PEL"]] ) {subpop[grep(value, pca$ind)] <- "PEL"}
for (value in  subpop_vectors[["MSL"]] ) {subpop[grep(value, pca$ind)] <- "MSL"}
for (value in  subpop_vectors[["KHV"]] ) {subpop[grep(value, pca$ind)] <- "KHV"}
for (value in  subpop_vectors[["PJL"]] ) {subpop[grep(value, pca$ind)] <- "PJL"}
for (value in  subpop_vectors[["YRI"]] ) {subpop[grep(value, pca$ind)] <- "YRI"}
for (value in  subpop_vectors[["ASW"]] ) {subpop[grep(value, pca$ind)] <- "ASW"}
for (value in  subpop_vectors[["KENYA"]] ) {subpop[grep(value, pca$ind)] <- "KENYA"}

pca <- as.tibble(data.frame(pca, pop, subpop)) # removed spp and spp_pop because all ind are humans

## plot the % of explained variance and the PCA

pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)

pop_order <- fct_relevel(pca$pop,"African","Non-African","CHM13")

ggplot(pca, aes(PC1, PC2, col = pop_order, shape = pop_order))+
  geom_point(size = 3)+
  scale_colour_manual(values = c("#00539CFF","#EEA47FFF","#A2A2A1FF"))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  labs(color="Population",shape="Population")+
  # coord_equal()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pca_nochm13 <- pca %>% filter(pca$ind != "CHM13")
pop_order <- fct_relevel(pca_nochm13$pop,"African","Non-African")
ggplot(pca_nochm13, aes(PC1, PC2, col = pop_order, shape = pop_order))+
  geom_point(size = 3)+
  scale_colour_manual(values = c("#00539CFF","#EEA47FFF"))+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  labs(color="Population",shape="Population")+
  # coord_equal()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/fig2_pca.png"), plot=ggplot2::last_plot())

subpop_order <- fct_relevel(pca$subpop,"ACB","GWD","MSL","ESN","YRI","ASW","KENYA",
                            "CHS","PUR","CLM","PEL","KHV","PJL",
                            "CHM13")
pop_order <- fct_relevel(pca$pop,"African","Non-African","CHM13")

fig2_pca <- ggplot(pca, aes(PC1, PC2, col = subpop_order, shape = pop_order))+
  geom_point(size = 3)+
  scale_color_viridis(discrete = TRUE, option = "C")+
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))+
  # labs(color="Subpopulation",shape="Population")+
  labs(color="Population",shape="Superpopulation")+
  # coord_equal()+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(angle = 90, vjust = 0.5))

# ggsave(file.path(sys_dir,"PanTE_human/paper_fig/figsupp_pca.png"), plot=ggplot2::last_plot())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "fig2_pca.png"),
       plot = fig2_pca, width = 9, height = 8, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

# Frequency based selection scans -----------------------------------------

Fst <- read.csv(file.path(sys_dir,"PanTE_human/selection_scans/fst/out.weir.fst"), header = T,sep = "")
Fst <- Fst %>% filter(Fst$WEIR_AND_COCKERHAM_FST >= 0)
Fst$CHROM <- gsub("chr","",Fst$CHROM)
Fst$CHROM <- as.numeric(Fst$CHROM)
Fstsubset <- Fst[complete.cases(Fst),]
TE <- c(1:(nrow(Fstsubset)))
mydf <- data.frame(TE,Fstsubset)
gray_chr_colors <- c("#4D4D4D","#AEAEAE","#E6E6E6")
mydf_99 <- quantile(mydf$WEIR_AND_COCKERHAM_FST, 0.99)
mydf_sig <- mydf[mydf$WEIR_AND_COCKERHAM_FST > mydf_99,] 
manhattan(mydf,chr="CHROM",bp="POS",p="WEIR_AND_COCKERHAM_FST",snp="TE",logp=FALSE,ylab="Fst", cex=2, xaxt="n",
          col=c(gray_chr_colors))
abline(h = mydf_99, col = "#EE353E", lty = 2)

# Haplotype based selection scans -----------------------------------------

## rehh - define data files

rehh_wrapper <- function(vcf_file, chr_name) {
  result <- data2haplohh(hap_file = vcf_file,
                         chr.name = chr_name,
                         polarize_vcf = F,
                         verbose = FALSE,
                         recode.allele = FALSE,
                         vcf_reader = "vcfR")
  scan <- scan_hh(result)
  return(scan)
  # assign(paste0("AFR_", chr_name), result, envir = .GlobalEnv)
  cat("Converting and scanning data:", vcf_file, "\n")
}

## outofAFR

directory_path <- file.path(sys_dir,"PanTE_human/selection_scans/ihs_xpehh/outofAFR")
file_pattern <- "outofAFR_1alt_phased_"
file_list <- list.files(directory_path, pattern = file_pattern, full.names = TRUE)
chr_id = as.factor(c("chr1", "chr10" ,"chr11" ,"chr12" ,"chr13" ,"chr14" ,"chr15" ,"chr16", "chr17", "chr18", "chr19" ,"chr2", "chr20","chr21" ,"chr22" ,"chr3" ,"chr4", "chr5" ,"chr6" ,"chr7" ,"chr8" ,"chr9"))
rehh_outofAFR <- data.frame() 

for (vcf_file in file_list) {
  for (chr_level in levels(chr_id)) {
    chr_name <- chr_level
    chr_name <- as.character(chr_level)
    if (grepl(paste0("phased_", chr_name, "\\.vcf"), vcf_file)) {
      # results_list_outofAFR[[paste0("outofAFR_", chr_name)]] <- get(paste0("outofAFR_", chr_name), envir = .GlobalEnv)
      scan <- rehh_wrapper(vcf_file, chr_name)
      rehh_outofAFR <- rbind(rehh_outofAFR, scan)  }}}

## AFR

directory_path <- file.path(sys_dir,"PanTE_human/selection_scans/ihs_xpehh/AFR")
file_pattern <- "AFR_1alt_phased_"
file_list <- list.files(directory_path, pattern = file_pattern, full.names = TRUE)
chr_id = as.factor(c("chr1", "chr10" ,"chr11" ,"chr12" ,"chr13" ,"chr14" ,"chr15" ,"chr16", "chr17", "chr18", "chr19" ,"chr2", "chr20","chr21" ,"chr22" ,"chr3" ,"chr4", "chr5" ,"chr6" ,"chr7" ,"chr8" ,"chr9"))
rehh_AFR <- data.frame() 

for (vcf_file in file_list) {
  for (chr_level in levels(chr_id)) {
    chr_name <- chr_level
    chr_name <- as.character(chr_level)
    if (grepl(paste0("phased_", chr_name, "\\.vcf"), vcf_file)) {
      scan <- rehh_wrapper(vcf_file, chr_name) 
      # results_list_AFR[[paste0("AFR_", chr_name)]] <- get(paste0("AFR_", chr_name), envir = .GlobalEnv)
      rehh_AFR <- rbind(rehh_AFR, scan)  }}}

rehh_AFR
rehh_outofAFR

rehh_AFR$CHR <- gsub("chr1([^0-9]|$)", "chr01\\1", rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr2([^0-9]|$)", "chr02\\1", rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr3","chr03",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr4","chr04",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr5","chr05",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr6","chr06",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr7","chr07",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr8","chr08",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr9","chr09",rehh_AFR$CHR)
rehh_AFR$CHR <- gsub("chr","",rehh_AFR$CHR)
rehh_AFR$CHR <- as.numeric(rehh_AFR$CHR)

rehh_outofAFR$CHR <- gsub("chr1([^0-9]|$)", "chr01\\1", rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr2([^0-9]|$)", "chr02\\1", rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr3","chr03",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr4","chr04",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr5","chr05",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr6","chr06",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr7","chr07",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr8","chr08",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr9","chr09",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- gsub("chr","",rehh_outofAFR$CHR)
rehh_outofAFR$CHR <- as.numeric(rehh_outofAFR$CHR)

chr_colors <- brewer.pal(n = 8, name = "Dark2")

## ihs - calculate ihs and manhattan plots

ihs_AFR <- ihh2ihs(rehh_AFR, min_maf = -1)
palette(chr_colors)
ihs_AFR_2 <- ihs_AFR$ihs %>% drop_na
ihs_AFR_2$IHS_abs <- abs(ihs_AFR_2$IHS)
per_99 <- quantile(ihs_AFR_2$IHS, 0.99)
per_01 <- quantile(ihs_AFR_2$IHS, 0.01)
manhattanplot(ihs_AFR,
              cex = 1.5,
              threshold = c(per_01,per_99),
              chr.name = c("1", "2","3" , "4","5", "6" , '7',"8","9","10","11" ,"12" ,"13", "14","15","16", "17" ,"18" ,"19","20","21","22"))
freqbinplot(ihs_AFR)
# manhattanplot(ihs_AFR,
#               pval = TRUE,
#               threshold = 2)

ihs_outofAFR <- ihh2ihs(rehh_outofAFR, min_maf = -1) 
palette(chr_colors)
ihs_outofAFR_2 <- ihs_outofAFR$ihs %>% drop_na
ihs_outofAFR_2$IHS_abs <- abs(ihs_outofAFR_2$IHS)
per_99 <- quantile(ihs_outofAFR_2$IHS, 0.99)
per_01 <- quantile(ihs_outofAFR_2$IHS, 0.01)
manhattanplot(ihs_outofAFR,
              cex = 1.5,
              threshold = c(per_01,per_99),
              chr.name = c("1", "2","3" , "4","5", "6" , '7',"8","9","10","11" ,"12" ,"13", "14","15","16", "17" ,"18" ,"19","20","21","22"))
freqbinplot(ihs_outofAFR)
# manhattanplot(ihs_outofAFR,
#               pval = TRUE,
#               threshold = 2)

chr_colors <- palette("default")
chr_colors <- brewer.pal(n = 8, name = "Dark2")
gray_chr_colors <- c("#4D4D4D","#AEAEAE","#E6E6E6")

ihs_AFR_outofAFR <- rbind(ihs_outofAFR_2,ihs_AFR_2)
ihs_AFR_outofAFR$TE <- seq(1,nrow(ihs_AFR_outofAFR))
# ihs_AFR_outofAFR$IHS_abs <- abs(ihs_AFR_outofAFR$IHS) 

manhattan(ihs_AFR_outofAFR,chr="CHR",bp="POSITION",p="IHS_abs",snp="TE",logp=FALSE,ylab="iHS", cex=2, xaxt="n",
          col=c(gray_chr_colors),suggestiveline=F,genomewideline=F)
per_99 <- quantile(ihs_AFR_2$IHS_abs, 0.99)
abline(h = per_99, col = "#EE353E", lty = 2)
per_99 <- quantile(ihs_outofAFR_2$IHS_abs, 0.99)
abline(h = per_99, col = "steelblue", lty = 2)

## ihs - extract candidate TE loci with highest 99% absolute ihs

ihs_outofAFR <- ihs_outofAFR$ihs %>% drop_na
ihs_outofAFR$IHS_abs <- abs(ihs_outofAFR$IHS)
per_99 <- quantile(ihs_outofAFR$IHS_abs, 0.99)
ihs_outofAFR_per_99 <- ihs_outofAFR[ihs_outofAFR$IHS_abs > per_99,]
ihs_outofAFR_candidate_TE <- ihs_outofAFR_per_99
# per_99 <- quantile(ihs_outofAFR$IHS, 0.99)
# per_01 <- quantile(ihs_outofAFR$IHS, 0.01)
# ihs_outofAFR_per_01 <- ihs_outofAFR[ihs_outofAFR$IHS < per_01,]
# ihs_outofAFR_per_99 <- ihs_outofAFR[ihs_outofAFR$IHS > per_99,]
# ihs_outofAFR_candidate_TE <- ihs_outofAFR %>% 
#   filter(LOGPVALUE > 2 | LOGPVALUE < -2) %>% 
#   dplyr::select(c("CHR","POSITION"))

ihs_AFR <- ihs_AFR$ihs %>% drop_na
ihs_AFR$IHS_abs <- abs(ihs_AFR$IHS)
per_99 <- quantile(ihs_AFR$IHS_abs, 0.99)
ihs_AFR_per_99 <- ihs_AFR[ihs_AFR$IHS_abs > per_99,]
ihs_AFR_candidate_TE <- ihs_AFR_per_99
# per_99 <- quantile(ihs_AFR$IHS, 0.99)
# per_01 <- quantile(ihs_AFR$IHS, 0.01)
# ihs_AFR_per_01 <- ihs_AFR[ihs_AFR$IHS < per_01,]
# ihs_AFR_per_99 <- ihs_AFR[ihs_AFR$IHS > per_99,]
# ihs_AFR_candidate_TE <- ihs_AFR %>% 
#   filter(LOGPVALUE > 2 | LOGPVALUE < -2) %>% 
#   dplyr::select(c("CHR","POSITION"))

ihs_AFR_outofAFR_candidate_TE <- rbind(ihs_AFR_candidate_TE,ihs_outofAFR_candidate_TE)
ihs_AFR_outofAFR_candidate_TE <- ihs_AFR_outofAFR_candidate_TE[order(ihs_AFR_outofAFR_candidate_TE[,1],ihs_AFR_outofAFR_candidate_TE[,2]),]
ihs_AFR_outofAFR_candidate_TE$POSITION <- as.numeric(ihs_AFR_outofAFR_candidate_TE$POSITION)
# write.table(ihs_AFR_outofAFR_candidate_TE, file = file.path(sys_dir,'PanTE_human/ihs_AFR_outofAFR_candidate_TE.bed'), sep = '\t', row.names = F,quote=F)

## xpehh - calculate xpehh and manhattan plots

xpehh <- ies2xpehh(scan_pop1 =  rehh_AFR,
                   scan_pop2 =  rehh_outofAFR,
                   popname1 = "AFR",
                   popname2 = "outofAFR") %>% drop_na()

manhattanplot(xpehh,
              cex = 1.5,
              # threshold = c(per_01,per_99),
              chr.name = c("1", "2","3" , "4","5", "6" , '7',"8","9","10","11" ,"12" ,"13", "14","15","16", "17" ,"18" ,"19","20"))

chr_colors <- palette("default")
chr_colors <- brewer.pal(n = 8, name = "Dark2")
gray_chr_colors <- c("#4D4D4D","#AEAEAE","#E6E6E6")

xpehh$TE <- seq(1,nrow(xpehh))

manhattan(xpehh,chr="CHR",bp="POSITION",p="XPEHH_AFR_outofAFR",snp="TE",logp=FALSE,ylab="xp-EHH", cex=2, xaxt="n",
          col=c(gray_chr_colors))
abline(h = per_99, col = "#EE353E", lty = 2)

## xpehh - extract candidate TE loci with highest 99% absolute ihs

per_99 <- quantile(xpehh$XPEHH_AFR_outofAFR, 0.99)
per_01 <- quantile(xpehh$XPEHH_AFR_outofAFR, 0.01)
xpehh_per_01 <- xpehh[xpehh$XPEHH_AFR_outofAFR < per_01,]
xpehh_per_99 <- xpehh[xpehh$XPEHH_AFR_outofAFR > per_99,]
xpehh_per_0199 <- rbind(xpehh_per_99,xpehh_per_01)
xpehh_per_0199_new_candidate_TE <- xpehh_per_0199 %>% dplyr::select(c("CHR","POSITION"))
xpehh_per_0199_new_candidate_TE <- xpehh_per_0199_new_candidate_TE[order(xpehh_per_0199_new_candidate_TE[,1],xpehh_per_0199_new_candidate_TE[,2]),]
xpehh_per_0199_new_candidate_TE$POSITION <- as.numeric(xpehh_per_0199_new_candidate_TE$POSITION)
# write.table(xpehh_per_0199_new_candidate_TE, file = file.path(sys_dir,'PanTE_human/xpehh_per_0199_new_candidate_TE.bed'), sep = '\t', row.names = F,quote=F)

# Ohana demography based selection scans ----------------------------------

ohana <- read.table(file = file.path(sys_dir,"PanTE_human/selection_scans/ohana/scansel.txt"), header=T)
pos <-  read.table(file = file.path(sys_dir,"PanTE_human/selection_scans/ohana/TE.map"), header=F) %>% dplyr::select(c(4))
all_pos <-  read.table(file = file.path(sys_dir,"PanTE_human/selection_scans/ohana/input_pos_ohana.txt"), header=F)
colnames(all_pos) <- c("chr","pos")
ohana_pos <- cbind(ohana,pos) %>% dplyr::select(c(-1)) %>% data.frame()
colnames(ohana_pos) <- c("global.lle","local.lle","lle.ratio","pos")
ohana_postions <- merge(ohana_pos,all_pos)

percentile_99 <- quantile(ohana_postions$lle.ratio, 0.99)
ohana_pos_sig <- ohana_postions %>% filter(ohana_postions$lle.ratio > percentile_99)
ohana_pos_sig_pos <- ohana_pos_sig %>% dplyr::select("chr","pos")
# write.table(ohana_pos_sig_pos, file = file.path(sys_dir,'PanTE_human/ohana_candidate_TE.bed'), sep = '\t', row.names = F,quote=F)

sig <-  ohana_postions %>% filter(ohana_postions$lle.ratio > percentile_99)
nrow(sig)

ohana_postions$chr <- gsub("chr1","1",ohana_postions$chr)
ohana_postions$chr <- gsub("chr2","2",ohana_postions$chr)
ohana_postions$chr <- gsub("chr(1[1-9])([^0-9]|$)", "\\1\\2", ohana_postions$chr)
ohana_postions$chr <- gsub("chr3","3",ohana_postions$chr)
ohana_postions$chr <- gsub("chr4","4",ohana_postions$chr)
ohana_postions$chr <- gsub("chr5","5",ohana_postions$chr)
ohana_postions$chr <- gsub("chr6","6",ohana_postions$chr)
ohana_postions$chr <- gsub("chr7","7",ohana_postions$chr)
ohana_postions$chr <- gsub("chr8","8",ohana_postions$chr)
ohana_postions$chr <- gsub("chr9","9",ohana_postions$chr)
ohana_postions$chr <- gsub("chr20","20",ohana_postions$chr)
ohana_postions$chr <- gsub("chr21","21",ohana_postions$chr)
ohana_postions$chr <- gsub("chr22","22",ohana_postions$chr)
ohana_postions$chr <- gsub("chrX","22",ohana_postions$chr)
ohana_postions$chr <- gsub("chrY","23",ohana_postions$chr)
unique(ohana_postions$chr)
ohana_postions$chr <- as.numeric(ohana_postions$chr)
ohana_postions$TE <- seq(1,nrow(ohana_postions))
ohana_postions$lle.ratio_log <- log10(ohana_postions$lle.ratio)
percentile_99 <- quantile(ohana_postions$lle.ratio_log, 0.99)

chr_colors <- palette("default")
chr_colors <- brewer.pal(n = 8, name = "Dark2")
gray_chr_colors <- c("#4D4D4D","#AEAEAE","#E6E6E6")

manhattan(ohana_postions,chr="chr",bp="pos",p="lle.ratio_log",snp="TE",logp=F,ylab="ohana (log scale)", cex=2, xaxt="n",
          col=c(gray_chr_colors),suggestiveline=F,genomewideline=F)
abline(h = (percentile_99), col = "#EE353E", lty = 2)

# Candidate loci - All selection scans ------------------------------------

## extract and merge candidate TE loci from all selection scans 

candidate_fst <- candidate_TE_div_Fst %>% dplyr::select("CHR","POS")                    # top 1% and young
candidate_ihs <- ihs_AFR_outofAFR_candidate_TE %>% dplyr::select("CHR","POSITION")      # top 1% of ihs in AFR and non AFR 
candidate_xpehh <- xpehh_per_0199_new_candidate_TE %>% dplyr::select("CHR","POSITION")  # top 1% and 99% of xpehh 
candidate_ohana <- ohana_pos_sig_pos %>% dplyr::select("chr","pos")                     # top 1% of lrt
colnames(candidate_fst) <- c("chr","pos")
colnames(candidate_ihs) <- c("chr","pos")
colnames(candidate_xpehh) <- c("chr","pos")
colnames(candidate_ohana) <- c("chr","pos")
candidate_fst$selscan <- "fst"
candidate_ihs$selscan <- "ihs"
candidate_xpehh$selscan <- "xpehh"
candidate_ohana$selscan <- "ohana"

candi_all <- rbind(candidate_fst,candidate_ihs,candidate_xpehh,candidate_ohana) %>% data.frame() 
candi_all$chr <- gsub("chr1","1",candi_all$chr)
candi_all$chr <- gsub("chr2","2",candi_all$chr)
candi_all$chr <- gsub("chr(1[1-9])([^0-9]|$)", "\\1\\2", candi_all$chr)
candi_all$chr <- gsub("chr3","3",candi_all$chr)
candi_all$chr <- gsub("chr4","4",candi_all$chr)
candi_all$chr <- gsub("chr5","5",candi_all$chr)
candi_all$chr <- gsub("chr6","6",candi_all$chr)
candi_all$chr <- gsub("chr7","7",candi_all$chr)
candi_all$chr <- gsub("chr8","8",candi_all$chr)
candi_all$chr <- gsub("chr9","9",candi_all$chr)
candi_all$chr <- gsub("chr20","20",candi_all$chr)
candi_all$chr <- gsub("chr21","21",candi_all$chr)
candi_all$chr <- gsub("chr22","22",candi_all$chr)
candi_all <- candi_all %>% mutate(chr = paste0("chr", chr))
candi_all <- candi_all[order(candi_all[,1],candi_all[,2]),]
rownames(candi_all) <- 1:nrow(candi_all)

# duplicated_entries <- duplicated(candi_all[c("chr", "pos")]) | duplicated(candi_all[c("chr", "pos")], fromLast = TRUE)
# candi_all_2 <- candi_all[!duplicated_entries, ]

nrow(candi_all)
nrow(candi_all_2) ## candi_all_2 is wrong as it removes the duplicated entries completely

candi_all_2 <- candi_all %>%
  group_by(chr, pos) %>%
  summarise(selscan = paste(sort(unique(selscan)), collapse = "_"), .groups = "drop")

nrow(candi_all)
nrow(candi_all_2) ## candi_all_2 is now correct

# write_xlsx(candi_all_2, file.path(sys_dir,'PanTE_human/candi_all.xlsx'))
# write.table(candi_all_2, file.path(sys_dir,'PanTE_human/candi_all.bed'),sep = '\t', row.names = F,quote=F)

## NOTE: consider HumGenom Review master data set and summary table instead of the following section

## merge with TE_div_Fst_filtered and TE_FULL_INFO to get full info on candidate TE loci

# > nrow(TE_div_Fst_filtered)
# [1] 2145
# > nrow(TE_FULL_INFO)
# [1] 4232

colnames(TE_div_Fst_filtered)
colnames(candi_all_2) <- c("CHR","POS","selscan")
colnames(candi_all) <- c("CHR","POS","selscan")

candi_all_2$POS <- as.numeric(candi_all_2$POS)
candi_all$POS <- as.numeric(candi_all$POS)

bin          <- 10000

df2          <- candi_all_2 %>% 
  mutate(window = POS %/% bin)
df2$window_CHR <- paste(df2$window, df2$CHR, sep = "_") 
df2          <- df2 %>%
  group_by(window_CHR) %>%
  summarise(sum_selscan = list(selscan))
df2[,3:4] <- stringr::str_split_fixed(df2$window_CHR, "_", 2)
df2 <- df2[-c(1)] 
colnames(df2) <- c("sum_selscan","window","CHR")
df2$window <- as.numeric(df2$window)
df2          <- df2 %>%
  mutate(start = (window*bin)+1, stop=(window+1)*bin) %>% 
  relocate(c("window","CHR","start","stop","sum_selscan"))
df_out       <- c()
for(i in unique(df2$CHR)) {
  df2_tmp      <- df2 %>% filter(CHR == i)
  missing_win  <- setdiff(0:max(df2_tmp$window),unique(df2_tmp$window))
  df_tmp       <- data.frame(window = missing_win) %>% 
    mutate(CHR    = i,
           sum_selscan = NA,
           start  = (window*bin)+1, 
           stop   = (window+1)*bin) %>% 
    bind_rows(df2_tmp) %>% 
    arrange(CHR,window)
  df_out       <- bind_rows(df_out,df_tmp)
}
candi_all_binned <- df_out %>% arrange(CHR, window) %>% drop_na()

colnames(TE_FULL_INFO)
colnames(TE_div_Fst_filtered)
colnames(candi_all_binned)
colnames(TE_size_distTSS_div_AO)
colnames(candi_all)

nrow(TE_FULL_INFO)
nrow(candi_all_binned)
nrow(TE_size_distTSS_div_AO)
nrow(candi_all)

candi_div_Fst <- merge(TE_div_Fst_filtered,candi_all)
candi_size_distTSS_div_AO <- merge(TE_size_distTSS_div_AO,candi_all)
candi_FULL_INFO <- merge(TE_FULL_INFO,candi_all_binned)

nrow(candi_div_Fst)
nrow(candi_size_distTSS_div_AO)
nrow(candi_FULL_INFO)

summary(candi_div_Fst$TE_FREQ)
summary(candi_div_Fst$SIZE)
table(candi_div_Fst$anno_order)
table(candi_div_Fst$selscan)

summary(candi_size_distTSS_div_AO$TE_FREQ)
summary(candi_size_distTSS_div_AO$SIZE)
table(candi_size_distTSS_div_AO$ANNO)
table(candi_size_distTSS_div_AO$selscan)

candi_FULL_INFO_HR <- candi_FULL_INFO %>% filter(mean_recrate > 5)

summary(candi_FULL_INFO_HR$mean_recrate)
summary(candi_FULL_INFO_HR$mean_disttss)

## run AnnotSV & VEP on candi_all.bed and extract summary stats - OLD

# ## remove the # from VEP_file: sed 's/#//g' TE_candi_all_VEP.txt > R_TE_candi_all_VEP.txt
# 
# VEP_res <- read.table(file = file.path(sys_dir,"PanTE_human/candidate_TE/R_TE_candi_all_VEP.txt"), header = T)
# 
# VEP_res$Consequence <- gsub("upstream_gene_variant", "Regulatory", VEP_res$Consequence)
# VEP_res$Consequence <- gsub("downstream_gene_variant", "Regulatory", VEP_res$Consequence)
# VEP_res$Consequence <- gsub("regulatory_region_ablation,regulatory_region_variant", "Regulatory", VEP_res$Consequence)
# VEP_res$Consequence <- gsub("TF_binding_site_variant", "Regulatory", VEP_res$Consequence)
# VEP_res$Consequence <- gsub("regulatory_region_variant", "Regulatory", VEP_res$Consequence)
# VEP_res$Consequence <- gsub("intron_variant,non_coding_transcript_variant", "Introns", VEP_res$Consequence)
# VEP_res$Consequence <- gsub("intron_variant,NMD_transcript_variant", "Introns", VEP_res$Consequence)
# VEP_res$Consequence <- gsub("intron_variant", "Introns", VEP_res$Consequence)
# VEP_res$Consequence <- gsub("intergenic_variant", "Intergenes", VEP_res$Consequence)
# VEP_res$Consequence <- gsub("non_coding_transcript_exon_variant", "non_coding_transcript_exon_variant", VEP_res$Consequence)
# 
# unique(VEP_res$Consequence)
# 
# VEP_res_filtered <- VEP_res %>%
#   group_by(Uploaded_variation) %>%
#   summarize(Consequence_summary = paste(unique(Consequence), collapse = ", "))
# 
# table(VEP_res_filtered$Consequence_summary)
# 
# # only Intergenes = 47
# # only Introns = 90
# # only Regulatory = 13 

# Selection scans stats and new combine -----------------------------------

## NOTE: consider HumGenom Review master data set and summary table instead of the following section

# Compute p-values, FDR, and ACAT combination for Fst, iHS, XP-EHH, and Ohana

ACAT_function <- function(p_values) {
  p_values <- p_values[is.finite(p_values) & !is.na(p_values) & p_values > 0 & p_values < 1]
  if (!length(p_values)) return(NA_real_)
  weights <- rep(1 / length(p_values), length(p_values))
  t_value <- sum(weights * tan((0.5 - p_values) * pi))
  pmax(pmin(0.5 - atan(t_value) / pi, 1), 0)
}

## Fst

fst_dt <- as.data.table(Fst)[is.finite(Fst$WEIR_AND_COCKERHAM_FST)]
fst_dt <- fst_dt[order(-Fst$WEIR_AND_COCKERHAM_FST)]
fst_dt[, p_fst := (seq_len(.N)) / (.N + 1)]
fst_dt[, q_fst := p.adjust(p_fst, "BH")]
fst_tab <- fst_dt[, .(chr = as.numeric(Fst$CHROM), pos = as.numeric(Fst$POS), p_fst, q_fst)]

plot(density(fst_dt$WEIR_AND_COCKERHAM_FST))
plot(density(fst_dt$p_fst))

nrow(fst_dt[])
nrow(fst_dt[fst_dt$p_fst < 0.05,])
nrow(fst_dt[fst_dt$q_fst < 0.05,])
nrow(fst_dt[fst_dt$WEIR_AND_COCKERHAM_FST > 0.25,])
nrow(fst_dt[fst_dt$WEIR_AND_COCKERHAM_FST > 0.4,])
nrow(fst_dt[fst_dt$WEIR_AND_COCKERHAM_FST > quantile(Fst$WEIR_AND_COCKERHAM_FST, 0.99),])

ggplot(fst_dt, aes(x = WEIR_AND_COCKERHAM_FST, y = -log10(p_fst))) +
  geom_point(alpha = 0.6) +
  # geom_point(data = fst_dt[WEIR_AND_COCKERHAM_FST > quantile(Fst$WEIR_AND_COCKERHAM_FST, 0.99), ], color = "red") +
  geom_point(data = fst_dt[fst_dt$p_fst < 0.05, ], color = "red") +
  labs(x = "Fst", y = "-log10(p-value)") +
  theme_minimal()

## iHS

## v1

ihs_A <- as.data.table(ihs_AFR$ihs)[, .(chr = as.numeric(ihs_AFR$CHR),
                                        pos = as.numeric(ihs_AFR$POSITION),
                                        z = ihs_AFR$IHS)]
ihs_O <- as.data.table(ihs_outofAFR$ihs)[, .(chr = as.numeric(ihs_outofAFR$CHR),
                                             pos = as.numeric(ihs_outofAFR$POSITION),
                                             z = ihs_outofAFR$IHS)]
ihs_all <- rbind(ihs_A, ihs_O, use.names = TRUE, fill = TRUE)[is.finite(z)]
ihs_all <- ihs_all[, .(z = z[which.max(abs(z))]), by = .(chr, pos)]
ihs_all[, p_ihs := 2 * pnorm(-abs(z))]
ihs_all[, q_ihs := p.adjust(p_ihs, "BH")]
ihs_tab <- ihs_all[, .(chr, pos, p_ihs, q_ihs)]

## v2 

ihs_A <- as.data.table(ihs_AFR$ihs)[
  , .(pop = "AFR",
      chr = as.numeric(ihs_AFR$CHR),
      pos = as.numeric(ihs_AFR$POSITION),
      z   = ihs_AFR$IHS)]

ihs_O <- as.data.table(ihs_outofAFR$ihs)[
  , .(pop = "OOA",
      chr = as.numeric(ihs_outofAFR$CHR),
      pos = as.numeric(ihs_outofAFR$POSITION),
      z   = ihs_outofAFR$IHS)]

ihs_all <- rbind(ihs_A, ihs_O, use.names = TRUE, fill = TRUE)[is.finite(z)]
ihs_all[, p_ihs := 2 * pnorm(-abs(z))]
ihs_all[, q_ihs := p.adjust(p_ihs, "BH")]
ihs_tab <- ihs_all[, .(pop, chr, pos, z, p_ihs, q_ihs)]

ihs_all <- rbind(ihs_A, ihs_O, use.names = TRUE, fill = TRUE)[is.finite(z)]
ihs_all[, p_ihs := 2 * pnorm(-abs(z)), by = pop]
ihs_all[, q_ihs := p.adjust(p_ihs, method = "BH"), by = pop]
ihs_tab <- ihs_all[, .(pop, chr, pos, z, p_ihs, q_ihs)]

plot(density(ihs_all$p_ihs))
plot(density(ihs_all$q_ihs))

sig <- ihs_all[q_ihs < 0.05]
sig_tab <- sig[, .(pop, chr, pos, z, p_ihs, q_ihs)]
nrow(sig_tab)
nrow(candidate_ihs)

## v3 

nrow(ihs_AFR[ihs_AFR$LOGPVALUE > -log10(0.05),])
nrow(ihs_outofAFR[ihs_outofAFR$LOGPVALUE > -log10(0.05),])
nrow(ihs_AFR_outofAFR_candidate_TE) ## top 1%

## XP-EHH

xdt <- as.data.table(xpehh)[, .(chr = as.numeric(xpehh$CHR),
                                pos = as.numeric(xpehh$POSITION),
                                z = xpehh$XPEHH_AFR_outofAFR)]
xdt <- xdt[is.finite(z)]
xdt[, p_xpehh := 2 * pnorm(-abs(z))]
xdt[, q_xpehh := p.adjust(p_xpehh, "BH")]
xpehh_tab <- xdt[, .(chr, pos, p_xpehh, q_xpehh)]

plot(density(xpehh$XPEHH_AFR_outofAFR))
plot(density(xpehh$LOGPVALUE))
plot(density(xdt$p_xpehh))

nrow(xdt[])
nrow(xdt[xdt$p_xpehh < 0.05,])
nrow(xdt[xdt$q_xpehh < 0.05,])
nrow(xdt[xpehh$LOGPVALUE > -log10(0.05),])
nrow(xdt[xpehh$XPEHH_AFR_outofAFR > quantile(xpehh$XPEHH_AFR_outofAFR, 0.99),])

ggplot(xpehh, aes(x = XPEHH_AFR_outofAFR, y = 10^LOGPVALUE)) +
  geom_point(alpha = 0.6) +
  geom_point(data = subset(xpehh, LOGPVALUE > 1.30103, na.rm = TRUE), color = "red", alpha = 0.8) +
  # geom_point(data = subset(xpehh, XPEHH_AFR_outofAFR > quantile(XPEHH_AFR_outofAFR, 0.99, na.rm = TRUE)), color = "red", alpha = 0.8) +
  # geom_point(data = subset(xpehh, XPEHH_AFR_outofAFR < quantile(XPEHH_AFR_outofAFR, 0.01, na.rm = TRUE)), color = "red", alpha = 0.8) +
  labs(x = "XP-EHH (AFR vs out-of-AFR)", y = "2^LOGPVALUE") +
  theme_minimal()

## Ohana

odt <- as.data.table(ohana_postions)[, .(chr = as.numeric(ohana_postions$chr),
                                         pos = as.numeric(ohana_postions$pos),
                                         lrt = as.numeric(ohana_postions$lle.ratio))]
odt <- odt[is.finite(lrt) & lrt >= 0]
odt[, p_ohana := pchisq(lrt, df = 1, lower.tail = FALSE)]
odt[, q_ohana := p.adjust(p_ohana, "BH")]
ohana_tab <- odt[, .(chr, pos, p_ohana, q_ohana)]

plot(density(odt$lrt))
plot(density(odt$p_ohana))

nrow(odt[])
nrow(odt[odt$p_ohana < 0.05,])
nrow(odt[odt$q_ohana < 0.05,])
nrow(odt[odt$lrt > quantile(odt$lrt, 0.99),])

ggplot(odt, aes(x = lrt, y = -log10(q_ohana))) +
  geom_point(alpha = 0.6) +
  geom_point(data = odt[lrt > quantile(odt$lrt, 0.99), ], color = "red") +
  # geom_point(data = odt[q_ohana < 0.05, ], color = "red") +
  labs(x = "lrt", y = "-log10(p-value)") +
  theme_minimal()

ohana_postions_stats <- merge(ohana_postions,odt)

## combine

selection_combined <- Reduce(function(a, b) full_join(a, b, by = c("chr", "pos")),
                             list(fst_tab, ihs_tab, xpehh_tab, ohana_tab)) %>% as.data.table()
selection_combined <- selection_combined[order(chr, pos)]
selection_combined[, p_acat := apply(cbind(p_fst, p_ihs, p_xpehh, p_ohana), 1, ACAT_function)]
selection_combined[, q_acat := p.adjust(p_acat, "BH")]

print(list(
  fst_significant = sum(selection_combined$q_fst < 0.05, na.rm = TRUE),
  ihs_significant = sum(selection_combined$q_ihs < 0.05, na.rm = TRUE),
  xpehh_significant = sum(selection_combined$q_xpehh < 0.05, na.rm = TRUE),
  ohana_significant = sum(selection_combined$q_ohana < 0.05, na.rm = TRUE),
  acat_significant = sum(selection_combined$q_acat < 0.05, na.rm = TRUE)))

# intersect Ohana and Fst loci

colnames(fst_dt) <- c("chr","pos","WEIR_AND_COCKERHAM_FST","p_fst","q_fst")    
colnames(xpehh) <- c("chr","pos","XPEHH_AFR_outofAFR","LOGPVALUE","q_fst")
colnames(odt)
summary(odt$lrt)

merged_ohana_fst <- merge(fst_dt, odt, by = c("chr", "pos"), all = FALSE)
merged_ohana_fst_xpehh <- merge(merged_ohana_fst, xpehh, by = c("chr", "pos"), all = FALSE)

ggplot(merged_ohana_fst_xpehh, aes(y = WEIR_AND_COCKERHAM_FST, x = lrt)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_point(data = merged_ohana_fst_xpehh[q_ohana < 0.05, ], color = "red") +
  labs(y = "WEIR_AND_COCKERHAM_FST", x = "lrt") +
  theme_minimal()

ggplot(merged_ohana_fst_xpehh, aes(y = XPEHH_AFR_outofAFR, x = lrt)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_point(data = merged_ohana_fst_xpehh[q_ohana < 0.05, ], color = "red") +
  labs(y = "XPEHH_AFR_outofAFR", x = "lrt") +
  theme_minimal()

ggplot(merged_ohana_fst_xpehh, aes(y = WEIR_AND_COCKERHAM_FST, x = XPEHH_AFR_outofAFR)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_point(data = merged_ohana_fst_xpehh[q_ohana < 0.05, ], color = "red") +
  labs(y = "WEIR_AND_COCKERHAM_FST", x = "XPEHH_AFR_outofAFR") +
  theme_minimal()

ggplot(merged_ohana_fst_xpehh, aes(y = WEIR_AND_COCKERHAM_FST, x = XPEHH_AFR_outofAFR)) +
  geom_point(alpha = 0.6, color = "steelblue") +
  geom_point(data = merged_ohana_fst_xpehh[
    q_ohana < 0.05 &
      WEIR_AND_COCKERHAM_FST > 0.25 &
      # abs(XPEHH_AFR_outofAFR) > 2.3199,
      LOGPVALUE > 1.30103,
  ], color = "red") +
  labs(y = "Fst", x = "XP-EHH") +
  theme_minimal()

ggplot(merged_ohana_fst_xpehh, aes(x = XPEHH_AFR_outofAFR, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(alpha = 0.5, color = "darkgray", size = 2.5, stroke = 0) +
  geom_point(data = merged_ohana_fst_xpehh[
    merged_ohana_fst_xpehh$p_ohana < 0.05 &
      merged_ohana_fst_xpehh$WEIR_AND_COCKERHAM_FST > Fst_percentile_99 &
      # abs(XPEHH_AFR_outofAFR) > per_99
      merged_ohana_fst_xpehh$LOGPVALUE > 1.30103
  ], color = "firebrick", size = 2.5, stroke = 0) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1.30103, 1.30103), linetype = "dashed", color = "black") +
  labs(x = "XP-EHH", y = "Fst") +
  scale_x_continuous(limits = c(-4.5, 4.5),breaks = seq(-4, 4, 1),expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank())

## final plot for NAR paper v2 on integrated selection scans
## Fst > 0.25; xp-EHH and Ohana likelihood ratio p-value < 0.05

merged_ohana_fst_xpehh_filtered <- as.data.frame(merged_ohana_fst_xpehh %>%
  dplyr::filter(p_ohana < 0.05, 
                WEIR_AND_COCKERHAM_FST > 0.25, 
                LOGPVALUE > 1.30103) %>%
  dplyr::select(c("chr","pos","WEIR_AND_COCKERHAM_FST","XPEHH_AFR_outofAFR","LOGPVALUE","lrt","p_ohana")))

ggplot(merged_ohana_fst_xpehh, aes(x = XPEHH_AFR_outofAFR, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(alpha = 0.5, color = "darkgray", size = 2.5, stroke = 0) +
  geom_point(merged_ohana_fst_xpehh_filtered, mapping=aes(x = XPEHH_AFR_outofAFR, y = WEIR_AND_COCKERHAM_FST), color = "firebrick", size = 2.5, stroke = 0) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1.30103, 1.30103), linetype = "dashed", color = "black") +
  labs(x = "XP-EHH", y = "Fst") +
  scale_x_continuous(limits = c(-4.5, 4.5),breaks = seq(-4, 4, 1),expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank())

colnames(candidate_fst) <- c("chr","pos","sel")
candidate_fst$chr <- as.character(candidate_fst$chr)
candidate_fst$pos <- as.character(candidate_fst$pos)
merged_ohana_fst_xpehh$chr <- as.character(merged_ohana_fst_xpehh$chr)
merged_ohana_fst_xpehh$pos <- as.character(merged_ohana_fst_xpehh$pos)

merged_ohana_fst_xpehh_mod <- merged_ohana_fst_xpehh %>%
  semi_join(candidate_fst, by = c("chr", "pos"))

ggplot(merged_ohana_fst_xpehh, aes(x = XPEHH_AFR_outofAFR, y = WEIR_AND_COCKERHAM_FST)) +
  geom_point(alpha = 0.5, color = "darkgray", size = 2.5, stroke = 0) +
  geom_point(data = merged_ohana_fst_xpehh[
    q_ohana < 0.05 &
      # WEIR_AND_COCKERHAM_FST > 0.25 &
      LOGPVALUE > 1.30103 & 
      pos %in% c(candidate_fst$pos)
  ], color = "firebrick", size = 2.5, stroke = 0) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1.30103, 1.30103), linetype = "dashed", color = "black") +
  labs(x = "XP-EHH", y = "Fst") +
  scale_x_continuous(limits = c(-4.5, 4.5),breaks = seq(-4, 4, 1),expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank())

## final plot for HumGenom Review on integrated selection scans
## Fst > 0.25; xp-EHH and Ohana likelihood ratio p-value < 0.05

nrow(TE_freq_anno_selscan_rec_master)
nrow(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod)
table(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod$selection_res)
table(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod$selection_res_comb)

figsupp_Fst_xpEHH_Ohana_labled <- ggplot()+
  geom_point(TE_freq_anno_selscan_rec_master, mapping=aes(x = xpehh, y = fst), color = "darkgray", size = 2.5, stroke = 0) +
  geom_point(subset(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod, selection_res_comb == "positive_selection_3meth"),
             mapping=aes(x = xpehh, y = fst), color = "firebrick", size = 2.5, stroke = 0) +
  geom_text_repel(subset(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod, selection_res_comb == "positive_selection_3meth"),
                  mapping=aes(x = xpehh, y = fst, label = Gene.refGene),
                  color = "black", size = 3, box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf) +
  geom_hline(yintercept = 0.25, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1.30103, 1.30103), linetype = "dashed", color = "black") +
  labs(x = "XP-EHH", y = "Fst") +
  scale_x_continuous(limits = c(-4.5, 4.5),breaks = seq(-4, 4, 1),expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank())

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_Fst_xpEHH_Ohana_labled.png"),
       plot = figsupp_Fst_xpEHH_Ohana_labled, width = 9, height = 8, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

# merged_ohana_fst_xpehh_filtered <- merged_ohana_fst_xpehh %>%
#   dplyr::filter(p_ohana < 0.05, 
#                 WEIR_AND_COCKERHAM_FST > 0.25, 
#                 # abs(XPEHH_AFR_outofAFR) > per_99
#                 LOGPVALUE > 1.30103
#                 ) %>%
#   # filter(pos %in% c(candidate_fst$pos)) %>%
#   dplyr::select(c("chr","pos"))
# 
# anti_join(
#   merged_ohana_fst_xpehh_filtered,
#   candidate_fst[, c("CHR", "POS")],
#   by = c("CHR", "POS"))

## NOTE: consider HumGenom Review master data set and summary table instead of the following section

## to merge with other df 

annovar_res$CHR <- gsub("chr","",annovar_res$CHR)
annovar_res$CHR <- as.character(annovar_res$CHR)
annovar_res$POS <- as.character(annovar_res$POS)
TE_div_Fst_filtered$CHR <- gsub("chr","",TE_div_Fst_filtered$CHR)
TE_div_Fst_filtered$CHR <- as.character(TE_div_Fst_filtered$CHR)
TE_div_Fst_filtered$POS <- as.character(TE_div_Fst_filtered$POS)
super_freq_anno$CHR <- gsub("chr","",super_freq_anno$CHR)
super_freq_anno$CHR <- as.character(super_freq_anno$CHR)
super_freq_anno$POS <- as.character(super_freq_anno$POS)

colnames(merged_ohana_fst_xpehh_filtered) <- c("CHR", "POS")
merged_ohana_fst_xpehh_filtered$CHR <- as.character(merged_ohana_fst_xpehh_filtered$CHR)
merged_ohana_fst_xpehh_filtered$POS <- as.character(merged_ohana_fst_xpehh_filtered$POS)

tmp <- copy(merged_ohana_fst_xpehh_filtered)
tmp_shift <- copy(tmp)
tmp_shift$POS <- as.integer(tmp_shift$POS) + 1
tmp_all <- rbind(tmp, tmp_shift, use.names = TRUE)

merged_ohana_fst_xpehh_filtered_anno <- merge(tmp_all, annovar_res,
                                              by = c("CHR", "POS"), all = FALSE)

merged_ohana_fst_xpehh_filtered_TEdiv <- merge(tmp_all, TE_div_Fst_filtered,
                                               by = c("CHR", "POS"), all = FALSE)

anti_join(
  merged_ohana_fst_xpehh_filtered,
  TE_div_Fst_filtered[, c("CHR", "POS")],
  by = c("CHR", "POS"))

merged_ohana_fst_xpehh_filtered_freq <- merge(tmp_all, super_freq_anno,
                                               by = c("CHR", "POS"), all = FALSE)

candidate_fst$chr <- gsub("chr","",candidate_fst$chr)
candidate_fst$chr <- as.character(candidate_fst$chr)
candidate_fst$pos <- as.character(candidate_fst$pos)

colnames(candidate_fst) <- c("CHR", "POS", "scan")

anti_join(
  merged_ohana_fst_xpehh_filtered,
  candidate_fst[, c("CHR", "POS")],
  by = c("CHR", "POS"))

TE_div_Fst_filtered$CHR <- gsub("chr","",TE_div_Fst_filtered$CHR)
TE_div_Fst_filtered$CHR <- as.character(TE_div_Fst_filtered$CHR)
TE_div_Fst_filtered$POS <- as.character(TE_div_Fst_filtered$POS)

merged_ohana_fst_xpehh_filtered_anno <- merge(tmp_all, TE_div_Fst_filtered,
                                              by = c("CHR", "POS"), all = FALSE)
merged_ohana_fst_xpehh_filtered_anno <- merge(merged_ohana_fst_xpehh_filtered_anno, annovar_res,
                                              by = c("CHR", "POS"), all = FALSE)
## EHH decay

outofAFR_haplohh <- data2haplohh(hap_file = "/home/shadi/Desktop/S3_Project/PanTE_human/selection_scans/ihs_xpehh/outofAFR/outofAFR_1alt_phased_chr11.vcf.gz",
                                 chr.name = "chr11",
                                 polarize_vcf = F,
                                 verbose = FALSE,
                                 recode.allele = FALSE,
                                 vcf_reader = "vcfR")
AFR_haplohh <- data2haplohh(hap_file = "/home/shadi/Desktop/S3_Project/PanTE_human/selection_scans/ihs_xpehh/AFR/AFR_1alt_phased_chr11.vcf.gz",
                            chr.name = "chr11",
                            polarize_vcf = F,
                            verbose = FALSE,
                            recode.allele = FALSE,
                            vcf_reader = "vcfR")

ehh_1 <- calc_ehh(AFR_haplohh, mrk = ">52081166>52081164")
plot(ehh_1)

ehh_2 <- calc_ehh(outofAFR_haplohh, mrk = ">52081166>52081164")
plot(ehh_2)

ehh_1 <- calc_ehh(AFR_haplohh, mrk = ">18963811>18963814")
plot(ehh_1)

ehh_2 <- calc_ehh(outofAFR_haplohh, mrk = ">18963811>18963814")
plot(ehh_2)

ehh_1 <- calc_ehh(AFR_haplohh, mrk = ">21356704>21356702")
plot(ehh_1)

ehh_2 <- calc_ehh(outofAFR_haplohh, mrk = ">21356704>21356702")
plot(ehh_2)

pos_chr8 <- ">40258312>40258309"
pos_chr3 <- ">13234392>13234390"

ehh_1 <- calc_ehh(AFR_haplohh, mrk = pos_chr3)
plot(ehh_1)
ehh_2 <- calc_ehh(outofAFR_haplohh, mrk = pos_chr3)
plot(ehh_2)

(outofAFR_haplohh@positions[4])-(outofAFR_haplohh@positions[1])

# Master TE dataset - HumGenom Review --------------------------------------

key_dfs <- c(
  "super_freq_anno","annovar_res","TE_div_Fst_filtered",
  "merged_ohana_fst_xpehh","merged_ohana_fst_xpehh_filtered",
  "selection_combined","candi_all","candi_all_binned",
  "TE_FULL_INFO","TE_size_distTSS_div_AO",
  "candi_div_Fst","candi_size_distTSS_div_AO","candi_FULL_INFO"
  )

for (nm in key_dfs) {
  cat("\n==========", nm, "==========\n")
  df <- get(nm)
  print(colnames(df))
  str(df)
  print(head(df))
}

## old - had some dup issues
# TE_freq_anno_selscan_master <- super_freq_anno %>%
#   mutate(POS = as.numeric(POS)) %>%
#   dplyr::select(CHR, POS, ANNO, new_class) %>%
#   # genomic location
#   left_join(annovar_res %>%
#               mutate(POS = as.numeric(POS)) %>% dplyr::select(CHR, POS, LOC), by = c("CHR", "POS")) %>%
#   # selection scan hits (candi_all uses lowercase chr/pos, chr numeric-free e.g. "chr1")
#   left_join(candi_all %>%
#               mutate(POS = as.numeric(pos)) %>%
#               mutate(CHR = chr) %>%
#               group_by(CHR, POS) %>%
#               summarise(selscans = paste(sort(unique(selscan)), collapse = ";"),
#                         n_selscans = n_distinct(selscan), .groups = "drop"),
#             by = c("CHR", "POS")) %>%
#   # per-scan p-/q-values (selection_combined uses numeric chr with NO "chr" prefix)
#   left_join(selection_combined %>%
#               mutate(CHR = paste0("chr", selection_combined$chr), POS = as.numeric(pos)) %>%
#               dplyr::select(CHR, POS, p_fst, q_fst, p_ihs, q_ihs, p_xpehh, q_xpehh, p_ohana, q_ohana),
#             by = c("CHR", "POS")) %>%
#   mutate(n_selscans = tidyr::replace_na(n_selscans, 0),
#          high_confidence = !is.na(selscans) & CHR %in% paste0("chr", merged_ohana_fst_xpehh_filtered$chr) &
#            POS %in% merged_ohana_fst_xpehh_filtered$pos)
# dups <- TE_freq_anno_selscan_rec_master[duplicated(TE_freq_anno_selscan_rec_master[, c("CHR", "POS")]) | duplicated(TE_freq_anno_selscan_rec_master[, c("CHR", "POS")], fromLast = TRUE), ]
# dups

## the old selection_combined do not retain the acutal values of Fst, ihs, and xphh
## also, no need to calculate p-value for ihs and xphh; the text already pops out a p-value
## final analysis should include a better calling of the final values for sel scans
# rm(odt,ohana_tab,xpehh_tab,xdt,fst_dt,fst_tab,ihs_A,ihs_tab,ihs_O)

## Fst
fst_dt <- as.data.table(Fst)[is.finite(Fst$WEIR_AND_COCKERHAM_FST)]
fst_dt <- fst_dt[order(-Fst$WEIR_AND_COCKERHAM_FST)]
fst_dt[, p_fst := seq_len(.N) / (.N + 1)]
fst_dt[, q_fst := p.adjust(p_fst, "BH")]
fst_tab <- fst_dt[, .(chr = as.numeric(CHROM), pos = as.numeric(POS), fst = WEIR_AND_COCKERHAM_FST, p_fst, q_fst)]

## iHS
ihs_A <- as.data.table(ihs_AFR$ihs)[, .(pop = "AFR", chr = as.numeric(ihs_AFR$CHR), pos = as.numeric(ihs_AFR$POSITION), ihs = ihs_AFR$IHS, LOGPVALUE = ihs_AFR$LOGPVALUE)]
ihs_O <- as.data.table(ihs_outofAFR$ihs)[, .(pop = "OOA", chr = as.numeric(ihs_outofAFR$CHR), pos = as.numeric(ihs_outofAFR$POSITION), ihs = ihs_outofAFR$IHS, LOGPVALUE = ihs_outofAFR$LOGPVALUE)]
ihs_tab <- rbind(ihs_A, ihs_O, use.names = TRUE, fill = TRUE)[is.finite(ihs) & is.finite(LOGPVALUE)]
# pivot iHS wide by population BEFORE rebuilding selection_combined
ihs_tab <- ihs_tab %>%
  mutate(chr = as.integer(chr), pos = as.integer(round(pos))) %>%
  tidyr::pivot_wider(
    id_cols = c(chr, pos),
    names_from = pop,
    values_from = c(ihs, LOGPVALUE),
    names_sep = "_")

## XP-EHH
xdt <- as.data.table(xpehh)[, .(chr = as.numeric(chr), pos = as.numeric(pos), xpehh = XPEHH_AFR_outofAFR, LOGPVALUE = LOGPVALUE)]
xdt <- xdt[is.finite(xpehh) & is.finite(LOGPVALUE)]
xpehh_tab <- xdt[, .(chr, pos, xpehh, LOGPVALUE)]

## Ohana
odt <- as.data.table(ohana_postions)
odt[, `:=`(chr = as.numeric(chr), pos = as.numeric(pos), lrt = as.numeric(lle.ratio))]
odt[, lle.ratio := NULL]
odt <- odt[is.finite(lrt) & lrt >= 0]
odt[, p_ohana := pchisq(lrt, df = 1, lower.tail = FALSE)]
odt[, q_ohana := p.adjust(p_ohana, "BH")]
ohana_tab <- odt

## combine
selection_combined <- Reduce(function(a, b) full_join(a, b, by = c("chr", "pos")), list(fst_tab, ihs_tab, xpehh_tab, ohana_tab)) %>% as.data.table()
selection_combined <- selection_combined[order(chr, pos)]

## prep AFR_freq_3 and outofAFR_freq_3 to merge with the master TE df 
## use AFR_freq_3_v2026 and outofAFR_freq_3_v2026 instead of AFR_freq_3 and outofAFR_freq_3

AFR_freq_3_v2026 <- AFR_freq_3_v2026 %>% dplyr::select(c(-"new_class"))
colnames(AFR_freq_3_v2026) <- c("CHR","POS","REF","FREQ_REF.AFR","ALT","FREQ_ALT.AFR","superpop.AFR", "ALT_size.AFR", "REF_size.AFR",
                                "TE_allele.AFR", "TE_size.AFR","TE_FREQ.AFR")

outofAFR_freq_3_v2026 <- outofAFR_freq_3_v2026 %>% dplyr::select(c(-"new_class"))
colnames(outofAFR_freq_3_v2026) <- c("CHR","POS","REF.o","FREQ_REF.OOA","ALT.o","FREQ_ALT.OOA" ,"superpop.OOA" ,"ALT_size.OOA","REF_size.OOA",
                                     "TE_allele.OOA","TE_size.OOA", "TE_FREQ.OOA")

AFR_freq_3_v2026$POS <- as.numeric(AFR_freq_3_v2026$POS)
outofAFR_freq_3_v2026$POS <- as.numeric(outofAFR_freq_3_v2026$POS)

## create the master table
# super_freq_anno_v2026 instead of super_freq_anno

TE_freq_anno_selscan_master <- super_freq_anno_v2026 %>%
  mutate(POS = as.numeric(POS)) %>%
  dplyr::select(CHR, POS, ANNO, new_class, TE_allele, TE_size, TE_FREQ) %>%
  left_join(AFR_freq_3_v2026 %>%
              dplyr::select(CHR, POS, FREQ_REF.AFR, FREQ_ALT.AFR, TE_FREQ.AFR),
            by = c("CHR", "POS")) %>%
  left_join(outofAFR_freq_3_v2026 %>%
              dplyr::select(CHR, POS, FREQ_REF.OOA, FREQ_ALT.OOA, TE_FREQ.OOA),
            by = c("CHR", "POS")) %>%
  left_join(TE_size_distTSS %>%
              mutate(POS = as.numeric(POS)) %>% dplyr::select(CHR, POS, distTSS), by = c("CHR", "POS")) %>%
  left_join(annovar_res %>%
              mutate(POS = as.numeric(POS)) %>% dplyr::select(CHR, POS, LOC), by = c("CHR", "POS")) %>%
  left_join(candi_all %>%
              # mutate(POS = as.numeric(pos), CHR = chr) %>%
              group_by(CHR, POS) %>%
              summarise(selscans = paste(sort(unique(selscan)), collapse = ";"),
                        n_selscans = n_distinct(selscan), .groups = "drop"),
            by = c("CHR", "POS")) %>%
  left_join(selection_combined %>%
              mutate(CHR = paste0("chr", chr), POS = as.numeric(pos)) %>%
              dplyr::select(CHR, POS,
                            fst, p_fst, q_fst,
                            ihs_AFR, LOGPVALUE_AFR,
                            ihs_OOA, LOGPVALUE_OOA,
                            xpehh, LOGPVALUE,
                            lrt, p_ohana, q_ohana),
            by = c("CHR", "POS")) %>%
  mutate(n_selscans = tidyr::replace_na(n_selscans, 0)) %>%
  as.data.table()

nrow(super_freq_anno_v2026)
nrow(TE_freq_anno_selscan_master)

# te  <- as.data.table(TE_freq_anno_selscan_master)[, `:=`(POS_start = POS, POS_end = POS)]
# bin <- as.data.table(TE_FULL_INFO)[, .(CHR, start, stop, mean_recrate, mean_disttss)]
# setkey(te,  CHR, POS_start, POS_end)
# setkey(bin, CHR, start, stop)
# TE_freq_anno_selscan_master_tmp <- foverlaps(te, bin, by.x = c("CHR", "POS_start", "POS_end"),
#                        by.y = c("CHR", "start", "stop"), type = "within") %>%
#   dplyr::select(-POS_start, -POS_end, -start, -stop) %>%
#   mutate(recomb_category = ifelse(mean_recrate > 5, "High_Rec_Rate", "Low_Rec_Rate"))
# TE_freq_anno_selscan_rec_master <- TE_freq_anno_selscan_master_tmp

TE_recrate <- TE_FULL_INFO %>%
  dplyr::select(c("CHR","start","stop","mean_recrate")) %>%
  as.data.table()

TE_freq_anno_selscan_rec_master <- TE_freq_anno_selscan_master[
  TE_recrate,
  mean_recrate := i.mean_recrate,
  on = .(CHR, POS >= start, POS <= stop)] %>%
  mutate(recomb_category = ifelse(mean_recrate > 5, "High_Rec_Rate", "Low_Rec_Rate"))

sum(is.na(TE_recrate$mean_recrate))
sum(is.na(TE_freq_anno_selscan_rec_master$mean_recrate))

nrow(TE_recrate) - nrow(TE_freq_anno_selscan_master)

nrow(TE_freq_anno_selscan_rec_master) # should be 6407
length(unique(TE_freq_anno_selscan_rec_master$POS)) # should be 6407

sum(is.na(TE_freq_anno_selscan_rec_master$LOC)) # 743
sum(is.na(annovar_res$LOC)) # 0
sum(TE_newloc$counts) # 5664
nrow(annovar_res) # 6407

sum(TE_freq_anno_selscan_rec_master$n_selscans > 0) # 154 

# TE_freq_anno_selscan_rec_master_filtered <- TE_freq_anno_selscan_rec_master %>%
#   filter(n_selscans > 0) %>%
#   mutate(selection_res = "positive_selection") %>%
#   filter(p_ohana < 0.05 & 
#          fst > 0.25 &
#          LOGPVALUE > 1.30103)

TE_freq_anno_selscan_rec_master_filtered <- TE_freq_anno_selscan_rec_master %>%
  # filter(n_selscans > 0) %>%
  filter(n_selscans > 0 | fst > 0.25) %>%
  mutate(
    selection_res = case_when(
      n_selscans > 0 ~ "positive_selection",
      TRUE ~ NA),
    selection_res_comb = case_when(
      p_ohana < 0.05 & fst > 0.25 & LOGPVALUE > 1.30103 ~ "positive_selection_3meth", # or if POS %in% c(merged_ohana_fst_xpehh_filtered$pos)
      TRUE ~ NA)
    # LOC = case_when(
    #   CHR == "chr19" & POS == 36555290   ~ "Introns",
    #   CHR == "chr2"  & POS == 179200614  ~ "Introns",
    #   CHR == "chr8"  & POS == 69626129   ~ "Introns",
    #   CHR == "chr11" & POS == 62143187   ~ "Introns",
    #   CHR == "chr20" & POS == 8835847   ~ "Introns",
    #   POS == 66711937   ~ "Intergenes",
    #   POS == 113714286   ~ "Intergenes",
    #   POS == 116387572   ~ "Intergenes",
    #   POS == 123236804   ~ "Intergenes",
    #   TRUE ~ LOC)
    )

## add proper anno info for genomic location

annovar_res_allinfo <- read.table(file = file.path(sys_dir,"PanTE_human_2025_v2/annovar/R_1alt_annovar_out_allinfo_v2.txt"),
                                  header=T,
                                  sep = "\t",
                                  stringsAsFactors = FALSE,
                                  quote = "")
colnames(annovar_res_allinfo)

colnames(annovar_res_allinfo) <- c("CHR","Start","End","LOC","Gene.refGene","GeneDetail.refGene","ExonicFunc.refGene","AAChange.refGene")
annovar_res_allinfo$LOC <- gsub("upstream;downstream", "Upstream/Downstream", annovar_res_allinfo$LOC)
annovar_res_allinfo$LOC <- gsub("upstream", "Upstream/Downstream", annovar_res_allinfo$LOC)
annovar_res_allinfo$LOC <- gsub("downstream", "Upstream/Downstream", annovar_res_allinfo$LOC)
annovar_res_allinfo$LOC <- gsub("intergenic", "Intergenes", annovar_res_allinfo$LOC)
annovar_res_allinfo$LOC <- gsub("intronic", "Introns", annovar_res_allinfo$LOC)
annovar_res_allinfo$LOC <- gsub("exonic", "Exons", annovar_res_allinfo$LOC)
annovar_res_allinfo$LOC <- gsub("UTR5", "5'UTR", annovar_res_allinfo$LOC)
annovar_res_allinfo$LOC <- gsub("UTR3", "3'UTR", annovar_res_allinfo$LOC)
annovar_res_allinfo$LOC <- gsub("ncRNA_intronic", "ncRNA_Introns", annovar_res_allinfo$LOC)
annovar_res_allinfo$LOC <- gsub("ncRNA_exonic", "ncRNA_Exons", annovar_res_allinfo$LOC)

annovar_res_allinfo$POS <- annovar_res_allinfo$Start - 1

## add info from annovar_res_allinfo that were not before

TE_freq_anno_selscan_rec_master_filtered <- TE_freq_anno_selscan_rec_master_filtered %>%
  left_join(
    annovar_res_allinfo %>% dplyr::select(CHR, Start, Gene.refGene, GeneDetail.refGene),
    by = c("CHR" = "CHR", "POS" = "Start"))

sum(is.na(TE_freq_anno_selscan_rec_master_filtered$LOC)) # 66
sum(is.na(TE_freq_anno_selscan_rec_master_filtered$GeneDetail.refGene)) # 66
sum(is.na(TE_freq_anno_selscan_rec_master_filtered$Gene.refGene)) # 66

TE_freq_anno_selscan_rec_master_filtered_mod <- TE_freq_anno_selscan_rec_master_filtered %>%
  left_join(
    annovar_res_allinfo %>% dplyr::select(CHR, POS, LOC, Gene.refGene, GeneDetail.refGene),
    by = c("CHR" = "CHR", "POS" = "POS"),
    suffix = c("", "_AnnoVar")) %>%
  mutate(
    LOC = ifelse(is.na(LOC), LOC_AnnoVar, LOC),
    Gene.refGene = ifelse(is.na(Gene.refGene), Gene.refGene_AnnoVar, Gene.refGene),
    GeneDetail.refGene = ifelse(is.na(GeneDetail.refGene), GeneDetail.refGene_AnnoVar, GeneDetail.refGene)
    ) %>%
  dplyr::select(-LOC_AnnoVar,-Gene.refGene_AnnoVar,-GeneDetail.refGene_AnnoVar)

sum(is.na(TE_freq_anno_selscan_rec_master$LOC)) # 743
sum(is.na(TE_freq_anno_selscan_rec_master_filtered_mod$LOC)) # 1
TE_freq_anno_selscan_rec_master_filtered_mod %>% filter(is.na(LOC))

## fix the chr11:112208430 manually which is "intronic	BCO2"

colnames(TE_freq_anno_selscan_rec_master_filtered_mod)

TE_freq_anno_selscan_rec_master_filtered_mod <- TE_freq_anno_selscan_rec_master_filtered_mod %>%
  mutate(
    LOC = case_when(
      CHR == "chr11" & POS == 112208430   ~ "Introns",
      TRUE ~ LOC),
    Gene.refGene = case_when(
      CHR == "chr11" & POS == 112208430   ~ "BCO2",
      TRUE ~ Gene.refGene))

sum(is.na(TE_freq_anno_selscan_rec_master_filtered_mod$LOC)) # 0
nrow(TE_freq_anno_selscan_rec_master_filtered_mod) # 183

colnames(TE_freq_anno_selscan_rec_master_filtered_mod)

# TE_freq_anno_selscan_rec_master_filtered_mod <- TE_freq_anno_selscan_rec_master_filtered_mod %>%
#   dplyr::select(-c("TE_allele", "FREQ_REF.AFR", "FREQ_ALT.AFR", "FREQ_REF.OOA", "FREQ_ALT.OOA"))

TE_freq_anno_selscan_rec_master_filtered_mod <- TE_freq_anno_selscan_rec_master_filtered_mod[, c(
    "CHR", "POS", "ANNO", "new_class",
    "FREQ_REF.AFR", "FREQ_ALT.AFR", "FREQ_REF.OOA", "FREQ_ALT.OOA",
    "TE_allele", "TE_FREQ", "TE_FREQ.AFR","TE_FREQ.OOA", "TE_size",
    "distTSS", "LOC", "Gene.refGene", "GeneDetail.refGene",
    "mean_recrate", "recomb_category",
    "selscans", "n_selscans",
    "fst", "p_fst", "q_fst",
    "ihs_AFR", "LOGPVALUE_AFR",
    "ihs_OOA", "LOGPVALUE_OOA",
    "xpehh", "LOGPVALUE",
    "lrt", "p_ohana", "q_ohana",
    "selection_res", "selection_res_comb")]

## is there a correlation between sel scan and rGREAT analysis?

genes_rGREAT_vec
TE_freq_anno_selscan_rec_master_filtered_mod$Gene.refGene

genes_selscan <- unique(unlist(strsplit(TE_freq_anno_selscan_rec_master_filtered_mod$Gene.refGene, ";")))
genes_selscan <- genes_selscan[genes_selscan != "NONE"]

inter_selscan_rGREAT <- lapply(genes_rGREAT_vec, function(x) {
  intersect(genes_selscan, x)
})

length(unlist(unique(genes_rGREAT_vec)))
# 2089
length(genes_selscan)
# 252
length(unlist(unique(inter_selscan_rGREAT)))
# 78

# selscan_rGREAT_genes <- TE_freq_anno_selscan_rec_master_filtered_mod %>%
#   mutate(Gene_refGene_split = strsplit(Gene.refGene, ";")) %>%
#   tidyr::unnest(Gene_refGene_split) %>%
#   mutate(
#     rGREAT_class = case_when(
#       Gene_refGene_split %in% inter_selscan_rGREAT$Rare &
#         Gene_refGene_split %in% inter_selscan_rGREAT$Polymorphic ~ "Rich_Rare_Polymorphic",
#       Gene_refGene_split %in% inter_selscan_rGREAT$Rare ~ "Rich_Rare",
#       Gene_refGene_split %in% inter_selscan_rGREAT$Polymorphic ~ "Rich_Polymorphic",
#       TRUE ~ NA_character_))

TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod <- TE_freq_anno_selscan_rec_master_filtered_mod %>%
  mutate(
    rGREAT_class = case_when(
      sapply(strsplit(Gene.refGene, ";"), function(x)
        any(x %in% inter_selscan_rGREAT$Rare) &
          any(x %in% inter_selscan_rGREAT$Polymorphic)) ~ "Rich_Rare_Polymorphic",
      
      sapply(strsplit(Gene.refGene, ";"), function(x)
        any(x %in% inter_selscan_rGREAT$Rare)) ~ "Rich_Rare",
      
      sapply(strsplit(Gene.refGene, ";"), function(x)
        any(x %in% inter_selscan_rGREAT$Polymorphic)) ~ "Rich_Polymorphic",
      
      TRUE ~ NA_character_))

## ***TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod*** is your final master TE df under +ve sel

## summary statistics table

TE_freq_anno_selscan_rec_rGREAT_master_154 <- TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod %>%
  filter(n_selscans > 0)

TE_master_pos_sel_3meth_18 <- TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod %>%
  filter(selection_res_comb == "positive_selection_3meth")

nrow(TE_freq_anno_selscan_rec_rGREAT_master_154)
nrow(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod)
nrow(TE_master_pos_sel_3meth_18)

## TE family
table(TE_freq_anno_selscan_rec_rGREAT_master_154$ANNO)
round(prop.table(table(TE_freq_anno_selscan_rec_rGREAT_master_154$ANNO)) * 100, 1)

## allele freq
table(TE_freq_anno_selscan_rec_rGREAT_master_154$new_class)
round(prop.table(table(TE_freq_anno_selscan_rec_rGREAT_master_154$new_class)) * 100, 1)

## rec rate
sum(is.na(TE_freq_anno_selscan_rec_rGREAT_master_154$mean_recrate)) # NA: 73 out of 154; 81 are not NA
table(TE_freq_anno_selscan_rec_rGREAT_master_154$recomb_category)
round(prop.table(table(TE_freq_anno_selscan_rec_rGREAT_master_154$recomb_category)) * 100, 1)

## distance to genes
sum(is.na(TE_freq_anno_selscan_rec_rGREAT_master_154$distTSS)) # NA: 4 out of 154
summary(TE_freq_anno_selscan_rec_rGREAT_master_154$distTSS)
TE_freq_anno_selscan_rec_rGREAT_master_154 %>%
  group_by(recomb_category) %>%
  summarise(mean_distTSS = mean(distTSS, na.rm = TRUE), .groups = "drop") %>%
  print(digits = Inf)

## genomic location
table(TE_freq_anno_selscan_rec_rGREAT_master_154$LOC)
round(prop.table(table(TE_freq_anno_selscan_rec_rGREAT_master_154$LOC)) * 100, 1)
TE_freq_anno_selscan_rec_rGREAT_master_154 %>%
  filter(LOC == "Intergenes") %>%
  summarise(mean_distTSS = mean(distTSS, na.rm = TRUE))

table(TE_master_pos_sel_3meth_18$ANNO)
table(TE_master_pos_sel_3meth_18$new_class)
table(TE_master_pos_sel_3meth_18$LOC)

table(TE_freq_anno_selscan_rec_rGREAT_master_154$rGREAT_class)

## func enrich of sel scan genes

## extract genes associated with intronic TE var under +ve sel 

TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod %>%
  filter(TE_FREQ == 0) %>%
  pull(selscans) %>%
  unique()

TE_freq_anno_selscan_rec_rGREAT_master_154 %>%
  filter(TE_FREQ == 0) %>%
  nrow

super_freq_anno_v2026 %>%
  filter(TE_FREQ == 0) %>%
  nrow

# TE_genes_pos <- TE_freq_anno_selscan_rec_rGREAT_master_154 %>%
TE_genes_pos <- TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod %>%
  # filter(!grepl(";", Gene.refGene)) %>%
  # filter(LOC == "Introns") %>%
  filter(LOC == "Introns" | LOC == "ncRNA_Introns") %>%
  pull(Gene.refGene) %>%
  unique()

## symbols --> Entrez

TE_genes_pos_map <- bitr(
  TE_genes_pos,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db)

TE_genes_pos_entrez <- unique(TE_genes_pos_map$ENTREZID)

length(unique(TE_genes_pos))
length(unique(TE_genes_pos_entrez))
setdiff(TE_genes_pos, TE_genes_pos_map$SYMBOL)

## GO BP enrichment

ORA_GO_BP_TE_genes_pos <- enrichGO(
  gene          = TE_genes_pos_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

ORA_GO_BP_TE_genes_pos_df <- as.data.frame(ORA_GO_BP_TE_genes_pos)

ORA_KEGG_TE_genes_pos <- enrichKEGG(
  gene = unique(TE_genes_pos_entrez),
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05)

ORA_KEGG_TE_genes_pos_df <- as.data.frame(ORA_KEGG_TE_genes_pos)

dotplot(ORA_GO_BP_TE_genes_pos, showCategory = 20)
dotplot(ORA_KEGG_TE_genes_pos,  showCategory = 20)

ggplot(ORA_KEGG_TE_genes_pos,
       aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment),
           fill = -log10(p.adjust))) +
  geom_col(width = 0.7) +
  scale_fill_viridis_c(name = expression(-log[10](q))) +
  labs(x = "Fold enrichment", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.grid.major.y = element_blank())

## extract genes associated with intergenic TE var under +ve sel 

# TE_genes_pos_intergenic <- TE_freq_anno_selscan_rec_rGREAT_master_154 %>%
TE_genes_pos_intergenic <- TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod %>%
  filter(LOC == "Intergenes") %>%
  mutate(genes = strsplit(Gene.refGene, ";"),
         dists = strsplit(GeneDetail.refGene, ";")) %>%
  mutate(gene_closest = map2_chr(genes, dists, ~{
      d <- as.numeric(gsub("dist=|NONE", "", .y))
      .x[which.min(ifelse(is.na(d), Inf, d))]})) %>%
  pull(gene_closest) %>%
  unique()

## symbols --> Entrez

TE_genes_pos_intergenic_map <- bitr(
  TE_genes_pos_intergenic,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db)

TE_genes_pos_intergenic_entrez <- unique(TE_genes_pos_intergenic_map$ENTREZID)

length(unique(TE_genes_pos_intergenic))
length(unique(TE_genes_pos_intergenic_entrez))
setdiff(TE_genes_pos_intergenic, TE_genes_pos_intergenic_map$SYMBOL)

## GO BP enrichment

ORA_GO_BP_TE_intergene_pos <- enrichGO(
  gene          = TE_genes_pos_intergenic_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

ORA_GO_BP_TE_intergene_pos_df <- as.data.frame(ORA_GO_BP_TE_intergene_pos)

ORA_KEGG_TE_intergene_pos <- enrichKEGG(
  gene = unique(TE_genes_pos_intergenic_entrez),
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05)

ORA_KEGG_TE_intergene_pos_df <- as.data.frame(ORA_KEGG_TE_intergene_pos)

dotplot(ORA_GO_BP_TE_intergene_pos, showCategory = 20)
dotplot(ORA_KEGG_TE_intergene_pos,  showCategory = 20)

ggplot(ORA_KEGG_TE_intergene_pos,
       aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment),
           fill = -log10(p.adjust))) +
  geom_col(width = 0.7) +
  scale_fill_viridis_c(name = expression(-log[10](q))) +
  labs(x = "Fold enrichment", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.grid.major.y = element_blank())

## lets separate high freq AFR and high freq OOA and then
## extract genes associated with TE var under +ve sel 

# for TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod

TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod <- TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod %>%
  mutate(
    freq_diff = TE_FREQ.OOA - TE_FREQ.AFR,
    # high_freq_diff = abs(freq_diff) >= quantile(abs(freq_diff), 0.75, na.rm = TRUE),
    high_freq_diff_AFR = freq_diff <= quantile(freq_diff, 0.10, na.rm = TRUE),
    high_freq_diff_OOA = freq_diff >= quantile(freq_diff, 0.90, na.rm = TRUE),
    high_freq_diff = high_freq_diff_AFR == T | high_freq_diff_OOA == T)

plot(density(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod$TE_FREQ.AFR))
plot(density(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod$TE_FREQ.OOA))
plot(density(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod$freq_diff))

ggplot(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod, aes(x = TE_FREQ.AFR, y = TE_FREQ.OOA, color = high_freq_diff)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#333333", "firebrick"), labels = c("Low difference", "High difference"), name = NULL) +
  labs(x = "AFR TE frequency", y = "Non-AFR TE frequency") +
  theme_minimal()

ggplot(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod, aes(x = TE_FREQ.AFR, y = TE_FREQ.OOA, color = high_freq_diff)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2.5) +
  geom_text_repel(
    data = subset(TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod, high_freq_diff == T & LOC == "Introns"),
    aes(label = Gene.refGene), color = "black", size = 3, box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf) +
  scale_color_manual(values = c("#333333", "firebrick"), labels = c("Low difference", "High difference"), name = NULL) +
  labs(x = "AFR TE frequency", y = "Non-AFR TE frequency") +
  theme_minimal()

# for TE_freq_anno_selscan_rec_rGREAT_master_154

TE_freq_anno_selscan_rec_rGREAT_master_154 <- TE_freq_anno_selscan_rec_rGREAT_master_154 %>%
  mutate(
    freq_diff = TE_FREQ.OOA - TE_FREQ.AFR,
    # high_freq_diff = abs(freq_diff) >= quantile(abs(freq_diff), 0.75, na.rm = TRUE),
    high_freq_diff_AFR = freq_diff <= quantile(freq_diff, 0.10, na.rm = TRUE),
    high_freq_diff_OOA = freq_diff >= quantile(freq_diff, 0.90, na.rm = TRUE),
    high_freq_diff = high_freq_diff_AFR == T | high_freq_diff_OOA == T)

plot(density(TE_freq_anno_selscan_rec_rGREAT_master_154$TE_FREQ.AFR))
plot(density(TE_freq_anno_selscan_rec_rGREAT_master_154$TE_FREQ.OOA))
plot(density(TE_freq_anno_selscan_rec_rGREAT_master_154$freq_diff))

figsupp_selscan_freqdiff <- ggplot(TE_freq_anno_selscan_rec_rGREAT_master_154, aes(x = TE_FREQ.AFR, y = TE_FREQ.OOA, color = high_freq_diff)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("#333333", "firebrick"), labels = c("Low difference", "High difference"), name = NULL) +
  labs(x = "AFR TE frequency", y = "Non-AFR TE frequency") +
  theme_minimal()

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_selscan_freqdiff.png"),
       plot = figsupp_selscan_freqdiff, width = 8, height = 6, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

ggplot(TE_freq_anno_selscan_rec_rGREAT_master_154, aes(x = TE_FREQ.AFR, y = TE_FREQ.OOA, color = high_freq_diff)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2.5) +
  geom_text_repel(
    data = subset(TE_freq_anno_selscan_rec_rGREAT_master_154, high_freq_diff == T & LOC == "Introns"),
    aes(label = Gene.refGene), color = "black", size = 3, box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf) +
  scale_color_manual(values = c("#333333", "firebrick"), labels = c("Low difference", "High difference"), name = NULL) +
  labs(x = "AFR TE frequency", y = "Non-AFR TE frequency") +
  theme_minimal()

## 4 AFR

TE_high_AFR_genes_pos <- TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod %>%
  filter(high_freq_diff_AFR == T) %>%
  filter(!grepl(";", Gene.refGene)) %>%
  pull(Gene.refGene) %>%
  unique()
length(unique(TE_high_AFR_genes_pos))

## symbols --> Entrez

TE_high_AFR_genes_pos_map <- bitr(
  TE_high_AFR_genes_pos,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db)
TE_high_AFR_genes_pos_entrez <- unique(TE_high_AFR_genes_pos_map$ENTREZID)
length(unique(TE_high_AFR_genes_pos))
length(unique(TE_high_AFR_genes_pos_entrez))
setdiff(TE_high_AFR_genes_pos, TE_high_AFR_genes_pos_map$SYMBOL)

## GO BP enrichment

ORA_GO_BP_TE_high_AFR_genes_pos <- enrichGO(
  gene          = TE_high_AFR_genes_pos_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
ORA_GO_BP_TE_high_AFR_genes_pos_df <- as.data.frame(ORA_GO_BP_TE_high_AFR_genes_pos)
ORA_KEGG_TE_high_AFR_genes_pos <- enrichKEGG(
  gene = unique(TE_high_AFR_genes_pos_entrez),
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05)
ORA_KEGG_TE_high_AFR_genes_pos_df <- as.data.frame(ORA_KEGG_TE_high_AFR_genes_pos)
dotplot(ORA_GO_BP_TE_high_AFR_genes_pos, showCategory = 10)
dotplot(ORA_KEGG_TE_high_AFR_genes_pos,  showCategory = 10)

## 4 OOA

TE_high_OOA_genes_pos <- TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod %>%
  filter(high_freq_diff_OOA == T) %>%
  filter(!grepl(";", Gene.refGene)) %>%
  pull(Gene.refGene) %>%
  unique()
length(unique(TE_high_OOA_genes_pos))

## symbols --> Entrez

TE_high_OOA_genes_pos_map <- bitr(
  TE_high_OOA_genes_pos,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db)
TE_high_OOA_genes_pos_entrez <- unique(TE_high_OOA_genes_pos_map$ENTREZID)
length(unique(TE_high_OOA_genes_pos))
length(unique(TE_high_OOA_genes_pos_entrez))
setdiff(TE_high_OOA_genes_pos, TE_high_OOA_genes_pos_map$SYMBOL)

## GO BP enrichment

ORA_GO_BP_TE_high_OOA_genes_pos <- enrichGO(
  gene          = TE_high_OOA_genes_pos_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
ORA_GO_BP_TE_high_OOA_genes_pos_df <- as.data.frame(ORA_GO_BP_TE_high_OOA_genes_pos)
ORA_KEGG_TE_high_OOA_genes_pos <- enrichKEGG(
  gene = unique(TE_high_OOA_genes_pos_entrez),
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05)
ORA_KEGG_TE_high_OOA_genes_pos_df <- as.data.frame(ORA_KEGG_TE_high_OOA_genes_pos)
dotplot(ORA_GO_BP_TE_high_OOA_genes_pos, showCategory = 10)
dotplot(ORA_KEGG_TE_high_OOA_genes_pos,  showCategory = 10)

ggplot(ORA_GO_BP_TE_high_OOA_genes_pos,
       aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment),
           fill = -log10(p.adjust))) +
  geom_col(width = 0.7) +
  scale_fill_viridis_c(name = expression(-log[10](q))) +
  labs(x = "Fold enrichment", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    legend.position = "right",
    panel.grid.major.y = element_blank())

## fun enrich of positive_selection_3meth_genes

positive_selection_3meth_genes <- TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod %>%
  filter(selection_res_comb == "positive_selection_3meth" &
           (LOC == "Introns" | LOC == "ncRNA_Introns")) %>%
  pull(Gene.refGene) %>%
  unique()

positive_selection_3meth_genes

# positive_selection_3meth_genes <- c("INCENP", "SESTD1", "PLCB1", "SULF1", "ZNF529")

positive_selection_3meth_genes_map <- bitr(
  positive_selection_3meth_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db)

positive_selection_3meth_genes_entrez <- unique(positive_selection_3meth_genes_map$ENTREZID)
length(unique(positive_selection_3meth_genes))
length(unique(positive_selection_3meth_genes_entrez))

## GO BP enrichment

ORA_GO_BP_positive_selection_3meth_genes <- enrichGO(
  gene          = positive_selection_3meth_genes_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)
ORA_GO_BP_positive_selection_3meth_genes_df <- as.data.frame(ORA_GO_BP_positive_selection_3meth_genes)
ORA_KEGG_positive_selection_3meth_genes <- enrichKEGG(
  gene = unique(positive_selection_3meth_genes_entrez),
  organism = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05)
ORA_KEGG_positive_selection_3meth_genes_df <- as.data.frame(ORA_KEGG_positive_selection_3meth_genes)
dotplot(ORA_GO_BP_positive_selection_3meth_genes, showCategory = 10) +
dotplot(ORA_KEGG_positive_selection_3meth_genes,  showCategory = 10)

ORA_GO_BP_positive_selection_3meth_genes_df_sig <- ORA_GO_BP_positive_selection_3meth_genes_df %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10)

ORA_KEGG_positive_selection_3meth_genes_df_sig <- ORA_KEGG_positive_selection_3meth_genes_df %>%
  filter(p.adjust < 0.05) %>%
  arrange(p.adjust) %>%
  slice_head(n = 10)

ORA_GO_BP_positive_selection_3meth_genes_plot <- ggplot(
  ORA_GO_BP_positive_selection_3meth_genes_df_sig,
  aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment), fill = -log10(p.adjust))) +
  geom_col(width = 0.7) +
  # scale_fill_viridis_c(name = expression(-log[10](p.adjust))) +
  scale_fill_gradient(name = expression(-log[10](p.adjust)), low = "firebrick2", high = "firebrick4") +
  labs(x = "Fold enrichment", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    # legend.title = element_text(size = 13),
    # legend.text = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.major.y = element_blank())
ORA_KEGG_positive_selection_3meth_genes_plot <- ggplot(
  ORA_KEGG_positive_selection_3meth_genes_df_sig,
  aes(x = FoldEnrichment, y = fct_reorder(Description, FoldEnrichment), fill = -log10(p.adjust))) +
  geom_col(width = 0.7) +
  # scale_fill_viridis_c(name = expression(-log[10](p.adjust))) +
  scale_fill_gradient(name = expression(-log[10](p.adjust)), low = "firebrick2", high = "firebrick4") +
  labs(x = "Fold enrichment", y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    # legend.title = element_text(size = 13),
    # legend.text = element_text(size = 12),
    legend.position = "bottom",
    panel.grid.major.y = element_blank())

figsupp_ORA_positive_selection_3meth_genes <-
  ORA_GO_BP_positive_selection_3meth_genes_plot + ORA_KEGG_positive_selection_3meth_genes_plot

ggsave(filename = file.path(sys_dir, "PanTE_human_2025_v2/paper_fig/HumGenom_Final", "figsupp_ORA_positive_selection_3meth_genes.png"),
       plot = figsupp_ORA_positive_selection_3meth_genes, width = 20, height = 12, units = "in", dpi = 600, bg = "white", limitsize = FALSE)

## I am not sure if this is a specific or random effect of genes
## create a random set of 5 genes five times and detected fun enrich
## sometimes yes and sometimes no; it's better not to do enrichment unless there is an underlying sig gene-set

random_genes <- annovar_res_allinfo %>%
  filter(LOC == "Introns") %>%
  filter(!is.na(Gene.refGene), Gene.refGene != "") %>%
  distinct(Gene.refGene) %>%
  pull(Gene.refGene)

GO_random <- list()
KEGG_random <- list()
genes_random <- list()

for (i in 1:5) {
  
  test <- sample(random_genes, 5)
  genes_random[[paste0("random_", i)]] <- test
  
  test_v2 <- bitr(
    test,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db)
  
  test_v2_entrez <- unique(test_v2$ENTREZID)
  
  ORA_GO_BP_TE_random <- enrichGO(
    gene          = test_v2_entrez,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    readable      = TRUE)
  
  GO_random[[paste0("random_", i)]] <- as.data.frame(ORA_GO_BP_TE_random) %>%
    mutate(random_set = paste0("random_", i))
  
  ORA_KEGG_TE_random <- enrichKEGG(
    gene          = test_v2_entrez,
    organism      = "hsa",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05)
  
  KEGG_random[[paste0("random_", i)]] <- as.data.frame(ORA_KEGG_TE_random) %>%
    mutate(random_set = paste0("random_", i))
}

GO_random_df <- bind_rows(GO_random)
KEGG_random_df <- bind_rows(KEGG_random)

## GO and KEGG plots

GO_random_plot <- GO_random_df %>%
  mutate(
    GeneRatio = sapply(strsplit(GeneRatio, "/"),
                       function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
  group_by(random_set) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  ggplot(aes(x = GeneRatio, y = reorder(Description, GeneRatio), size = Count, color = p.adjust)) +
  geom_point() +
  facet_wrap(~random_set, scales = "free_y") +
  labs(x = "Gene ratio", y = NULL, size = "Count", color = "Adjusted p-value") +
  theme_bw()

KEGG_random_plot <- KEGG_random_df %>%
  mutate(
    GeneRatio = sapply(strsplit(GeneRatio, "/"),
                       function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
  group_by(random_set) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  ggplot(aes(x = GeneRatio, y = reorder(Description, GeneRatio),size = Count, color = p.adjust)) +
  geom_point() +
  facet_wrap(~random_set, scales = "free_y") +
  labs(x = "Gene ratio", y = NULL, size = "Count", color = "Adjusted p-value") +
  theme_bw()

# ## should we supply the annovar missed POS with the proper "Start-1"
# ## select a random set of 10 var and test in the Genome Browser
# ## all PASSED
# 
# df_NA <- TE_freq_anno_selscan_rec_master_filtered[is.na(TE_freq_anno_selscan_rec_master_filtered$LOC), ]
# set.seed(2009)
# df_random <- df_NA[sample(nrow(df_NA), 10), ] 
# df_random <- df_random %>% arrange(CHR, POS)
# paste0(df_random$CHR, ":", df_random$POS - 2, "-", df_random$POS + 3)
# 
# annovar_res %>% filter(POS %in% c((df_random$POS)+1)) %>% arrange(CHR, POS)
# # CHR       POS        LOC
# # 1   chr1  78061629    Introns     OK
# # 2  chr11 110523823 Intergenes     OK
# # 3  chr18  14822957    Introns     OK
# # 4  chr19  46452845 Intergenes     OK
# # 5   chr2 183110365 Intergenes     OK - might be a gene?
# # 6  chr20  32311425    Introns     OK     
# # 7  chr21  36216094    Introns     OK  
# # 8   chr4  66711938 Intergenes     OK
# # 9   chr5 123236805 Intergenes     OK - might be a gene?
# # 10  chr9  34159868 Intergenes     OK
# 
# sum(is.na(TE_freq_anno_selscan_rec_master_filtered$LOC)) # 59
# table(TE_freq_anno_selscan_rec_master_filtered$LOC)
# 
# table(TE_freq_anno_selscan_rec_master_filtered$selection_res)
# table(TE_freq_anno_selscan_rec_master_filtered$comb_sel_scan_res)
# 
# # TE_freq_anno_selscan_rec_master_filtered <- TE_freq_anno_selscan_rec_master %>%
# #   # filter(n_selscans > 0) %>%
# #   mutate(selection_res = "positive_selection") %>%
# #   filter(POS %in% c(merged_ohana_fst_xpehh_filtered$pos))
# 
# missing_candidates <- merged_ohana_fst_xpehh_filtered %>%
#   mutate(CHR = paste0("chr", chr), POS = as.numeric(pos)) %>%
#   distinct(CHR, POS, .keep_all = TRUE) %>%
#   anti_join(TE_freq_anno_selscan_rec_master_filtered, by = c("CHR", "POS"))
# 
# table(TE_freq_anno_selscan_rec_master_filtered$n_selscans)
# 
# nrow(TE_freq_anno_selscan_rec_master_filtered)
# length(unique(TE_freq_anno_selscan_rec_master_filtered$POS)) # 13
# length(unique(merged_ohana_fst_xpehh_filtered$pos)) # 18
# 
# length(unique(candi_all$pos)) # 154
# length(unique(TE_freq_anno_selscan_rec_master_filtered$POS)) # 154
# 
# all(
#   paste(candi_all$chr, candi_all$pos) %in%
#     paste(TE_freq_anno_selscan_rec_master_filtered$CHR, TE_freq_anno_selscan_rec_master_filtered$POS)
#   )
# 
# # [this analysis is a bit old as the updated version already adds the missing POS anno]
# 
# ## for these sites, grep the original anno file (anno_g.hg38_multianno.txt) to check their location
# # chr19:36555290-36555291     OK
# # chr2:179200614-179200615	  OK
# # chr8:69626129-69626130      OK
# # chr4:66711937-66711938      OK
# # chr7:113714286-113714287    OK
# # chr8:116387572-116387573    OK
# # chr5:123236804-123236805	  OK
# ## checked genome browser and all are valid location in coordinate X and coordinates X+1
# 
# # grep -E '66711937|113714286|116387572|123236804|69626129|179200614|36555290' anno_g.hg38_multianno.txt 
# # ## NA
# # grep -E '66711938|113714287|116387573|123236805|69626130|179200615|36555291' anno_g.hg38_multianno.txt 
# # ## not NA
# # grep -E '62143188|8835848' anno_g.hg38_multianno.txt 
# # ## not NA
# 
# ## the gene introns are 
# # chr19	36555291	intronic	ZNF529
# # chr2	179200615	intronic	SESTD1
# # chr8	69626130	intronic	SULF1
# # chr11	62143188	intronic	INCENP
# # chr20	8835848	intronic	PLCB1

# Supp Tables - GB -------------------------------------------------------------

# Rdataframes_CT <- data.frame(
#     Sheet = c("Sheet1", "Sheet2", "Sheet3", "Sheet4", "Sheet5", "Sheet6", 
#               "Sheet7", "Sheet8", "Sheet9", "Sheet10", "Sheet11", "Sheet12", 
#               "Sheet13", "Sheet14", "Sheet15", "Sheet16", "Sheet17", "Sheet18"),
#     DataFrameName = c("global_freq", "global_freq_annotation", "genomic_location", 
#                       "TE_dataframe", "TE_divergence_fst", "TE_FULL_INFO", 
#                       "TE_size_chromosome_length_counted", "TE_size_chromosome_length", "TE_counts_per_chromosome", 
#                       "background_recrate", "major_singletons_shared_per_individual", "TE_counts_per_individual", 
#                       "ihs", "xpehh", "ohana", "candidate_TE_ihs", "candidate_TE_xpehh", "candidate_TE_all"))
# 
# write_xlsx(list(
#                 ContentTable = Rdataframes_CT,
#                 Sheet1 = super_freq_7,
#                 Sheet2 = super_freq_anno_counted,
#                 Sheet3 = TE_loc_c,
#                 Sheet4 = TE_dataframe,
#                 Sheet5 = TE_div_Fst_filtered,
#                 Sheet6 = TE_FULL_INFO,
#                 Sheet7 = size_length_counted,
#                 Sheet8 = size_length,
#                 Sheet9 = biTE_count_2,
#                 Sheet10 = background_rec,
#                 Sheet11 = maj_sin_shar_tidy,
#                 Sheet12 = TE_indv,
#                 Sheet13 = ihs_AFR_outofAFR,
#                 Sheet14 = xpehh,
#                 Sheet15 = ohana_postions,
#                 Sheet16 = ihs_AFR_outofAFR_candidate_TE,
#                 Sheet17 = xpehh_per_0199,
#                 Sheet18 = candi_all),
#            file.path(sys_dir,"PanTE_human/scripts/supp_data_Rdataframes.xlsx"))

# Supp Tables - NAR -------------------------------------------------------------

# Rdataframes_CT <- data.frame(
#   Sheet = c("Sheet1", "Sheet2", "Sheet3", "Sheet4", "Sheet5", "Sheet6", 
#             "Sheet7", "Sheet8", "Sheet9", "Sheet10", "Sheet11", "Sheet12", 
#             "Sheet13", "Sheet14", "Sheet15", "Sheet16", "Sheet17", "Sheet18", "Sheet19"),
#   
#   DataFrameName = c(
#     "global_freq",
#     "global_freq_annotation",
#     "genomic_location",
#     "TE_dataframe",
#     "TE_divergence_fst",
#     "TE_FULL_INFO",
#     "TE_size_chromosome_length_counted",
#     "TE_size_chromosome_length",
#     "TE_counts_per_chromosome",
#     "background_recrate",
#     "major_singletons_shared_per_individual",
#     "TE_counts_per_individual",
#     "ihs",
#     "xpehh",
#     "ohana",
#     "candidate_TE_ihs",
#     "candidate_TE_xpehh",
#     "candidate_TE_all",
#     "integrated_candidate_TE"
#   ),
#   
#   DataFrameNote = c(
#     "Genome-wide allele frequencies of all 6,407 biallelic TE variants across the HPRC v.1 individuals",
#     "TE allele frequencies annotated with repeat TE types from RepeatMasker",
#     "Genomic context of each TE insertion including intergenic, intronic, UTRs, and exonic regions",
#     "Master TE dataset containing TE coordinates, sizes, divergence and annotations per genomic window of 10 kb",
#     "TE divergence estimates integrated with Fst values",
#     "Complete table of TE features including frequency class, location, recombination rate and population metrics",
#     "TE counts per chromosome",
#     "TE counts per chromosome normalized to chromosome/TE size",
#     "TE counts per chromosome across the pangenome dataset",
#     "Random background genomic windows (10 kb) with recombination rate estimates",
#     "Counts of singleton, major alleles (frequency > 0.5), and alleles present in all individuals but absent from GRCh38, grouped by population",
#     "Total TE counts per individual across the HPRC dataset grouped by superpopulation",
#     "Haplotype-based selection scan summary table - iHS scores for TE loci to detect selection signatures",
#     "Haplotype-based selection scan summary table - xp-EHH scores comparing African vs non-African haplotypes",
#     "Ohana admixture-aware selection scan summary table - Ohana likelihood ratio statistics (p-value: X2 distribution, one degree of freedom; q-value: Benjamini–Hochberg method)",
#     "Candidate TE loci under positive selection identified by top 1% iHS scores",
#     "Candidate TE loci under positive selection identified by top 1% xp-EHH scores",
#     "Combined list of candidate TE loci from all selection scans (top 1% in Fst, iHS, xp-EHH, and Ohana LRT)",
#     "Integrated list of candidate TE loci significant for three selection scans (Fst > 0.25; xp-EHH and Ohana likelihood ratio p-value < 0.05)"
#   )
# )
# 
# openxlsx::write.xlsx(list(
#   ContentTable = Rdataframes_CT,
#   Sheet1 = super_freq_7,
#   Sheet2 = super_freq_anno_counted,
#   Sheet3 = TE_loc_c,
#   Sheet4 = TE_dataframe,
#   Sheet5 = TE_div_Fst_filtered,
#   Sheet6 = TE_FULL_INFO,
#   Sheet7 = size_length_counted,
#   Sheet8 = size_length,
#   Sheet9 = biTE_count_2,
#   Sheet10 = background_rec,
#   Sheet11 = maj_sin_shar_tidy,
#   Sheet12 = TE_indv,
#   Sheet13 = ihs_AFR_outofAFR,
#   Sheet14 = xpehh,
#   Sheet15 = ohana_postions_stats, # instead of ohana_postions
#   Sheet16 = ihs_AFR_outofAFR_candidate_TE,
#   Sheet17 = xpehh_per_0199,
#   Sheet18 = candi_all,
#   Sheet19 = merged_ohana_fst_xpehh_filtered
# ),
# file.path(sys_dir,"PanTE_human_2025_v2/scripts/supp_data_Rdataframes_v2.xlsx"))

# Supp Tables - HG -------------------------------------------------------------

## gather selection scans res together 2 report

selection_combined_2report <- Reduce(function(a, b) full_join(a, b, by = c("chr", "pos")), list(fst_tab, ihs_tab, xpehh_tab, ohana_tab)) %>% as.data.table()
selection_combined_2report <- selection_combined_2report[order(chr, pos)]
selection_combined_2report <- selection_combined_2report %>%
  dplyr::select(-c(TE))

identical(selection_combined$pos, selection_combined_2report$pos)

## gather AFR and non AFR freq res together 2 report

AFR_nonAFR_freq_3_v2026_2report <- merge(
  AFR_freq_3_v2026,
  outofAFR_freq_3_v2026,
  by = c("CHR", "POS"))

nrow(AFR_nonAFR_freq_3_v2026_2report)
identical(AFR_nonAFR_freq_3_v2026_2report$REF, AFR_nonAFR_freq_3_v2026_2report$REF.o)
identical(AFR_nonAFR_freq_3_v2026_2report$ALT, AFR_nonAFR_freq_3_v2026_2report$ALT.o)

## gather GO:BP and KEGG res together 2 report

ORA_GO_BP_positive_selection_3meth_genes_df_sig$category <- "GO Biological Process"
ORA_GO_BP_positive_selection_3meth_genes_df_sig$subcategory <- NA

ORA_GO_BP_KEGG_positive_selection_3meth_genes_df_sig_2report <- bind_rows(
  ORA_GO_BP_positive_selection_3meth_genes_df_sig,
  ORA_KEGG_positive_selection_3meth_genes_df_sig)

## run for the final SUPP TABLES

Rdataframes_CT <- data.frame(
  Sheet = c("Sheet1", "Sheet2", "Sheet3", "Sheet4", "Sheet5",
            "Sheet6", "Sheet7", "Sheet8", "Sheet9", "Sheet10",
            "Sheet11", "Sheet12", "Sheet13", "Sheet14", "Sheet15",
            "Sheet16", "Sheet17", "Sheet18", "Sheet19", "Sheet20"),
  
  DataFrameName = c(
    "super_freq_anno_v2026_chr_length",
    "AFR_nonAFR_freq_3_v2026_2report",
    "annovar_res_allinfo",
    "TE_newloc",
    
    "TE_candi_all_v2026",
    "TE_dataframe",
    "TE_div_Fst_filtered",
    "TE_FULL_INFO",
    
    "size_length_counted",
    "size_length",
    "biTE_count_2",
    "background_rec",
    "maj_sin_shar_tidy",
    "TE_indv_sin_tidy",
    
    "tb_allfreq_master_top10",
    
    "candidate_TE_div_Fst",
    "selection_combined_2report",
    "TE_freq_anno_selscan_rec_master",
    "TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod",
    "ORA_GO_BP_KEGG_positive_selection_3meth_genes_df_sig_2report"
  ),
  
  DataFrameNote = c(
    "Genome-wide TE allele frequencies of all 6,407 biallelic TE variants across the HPRC v.1 individuals with chromosome length information",
    "Population-specific TE allele frequencies in African and non-African populations",
    "ANNOVAR annotation of TE variants including genomic location, overlapping genes, and annotations based on GRCh38",
    "Counts of TE variants by genomic location, frequency class, and TE type",
    "TE annotations with allele frequencies, genomic location, divergence, age estimates, and candidate loci under positive selection",
    "TE genomic window dataset containing TE coordinates, annotations, sizes, divergence, and recombination rate per 10 kb genomic window",
    "TE divergence estimates integrated with Fst values",
    "Master TE genomic window dataset containing TE features including frequency class, genomic location, divergence, recombination rate, and population metrics",
    "TE counts per chromosome",
    "TE counts per chromosome normalized to chromosome/TE size",
    "TE counts per chromosome across the pangenome dataset annotated with TE type",
    "Random background genomic windows of 10 kb with recombination rate estimates",
    "Counts of singleton, major alleles (frequency > 0.5), and alleles present in all individuals but absent from GRCh38, grouped by population",
    "Total TE counts and singleton counts per individual across the HPRC dataset grouped by superpopulation",
    "rGREAT results with top 10 most enriched genomic regions for each TE frequency class (binomial or hypergeometric adjusted p-value ≤ 0.05)",
    "Candidate TE loci with high Fst values and young age estimates",
    "Combined selection scan results including Fst, iHS, xp-EHH and Ohana summary statistics for all TE loci included in the analyses",
    "TE dataset integrating TE frequency, size, genomic location, recombination rate, and selection scan results",
    "Master TE dataset for loci under positive selection (top 1% in Fst, iHS, xp-EHH, and Ohana LRT) integrating all annotations (e.g., frequencies, genomic location, recombination rate, selection scan statistics, and rGREAT results)",
    "Significantly enriched GO biological process terms and KEGG pathways for genes associated with intronic TE loci under positive selection"
  )
)

openxlsx::write.xlsx(list(
  ContentTable = Rdataframes_CT,
  Sheet1 = super_freq_anno_v2026_chr_length,
  Sheet2 = AFR_nonAFR_freq_3_v2026_2report,
  Sheet3 = annovar_res_allinfo,
  Sheet4 = TE_newloc,
  Sheet5 = TE_candi_all_v2026,
  Sheet6 = TE_dataframe,
  Sheet7 = TE_div_Fst_filtered,
  Sheet8 = TE_FULL_INFO,
  Sheet9 = size_length_counted,
  Sheet10 = size_length,
  Sheet11 = biTE_count_2,
  Sheet12 = background_rec,
  Sheet13 = maj_sin_shar_tidy,
  Sheet14 = TE_indv_sin_tidy,
  Sheet15 = tb_allfreq_master_top10,
  Sheet16 = candidate_TE_div_Fst,
  Sheet17 = selection_combined_2report,
  Sheet18 = TE_freq_anno_selscan_rec_master,
  Sheet19 = TE_freq_anno_selscan_rec_rGREAT_master_filtered_mod,
  Sheet20 = ORA_GO_BP_KEGG_positive_selection_3meth_genes_df_sig_2report
),
file.path(sys_dir,"PanTE_human_2025_v2/scripts/supp_data_Rdataframes_HumGenom.xlsx"))


