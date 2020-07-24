library(tidyverse) 
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
chr=args[1]
left=as.numeric(args[2])
right=as.numeric(args[3])

### NIC526 - p1 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/NIC526_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

NIC526 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("NIC526") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### JU1400 - p2 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/JU1400_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

JU1400 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("JU1400") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### JU310 - p3 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/JU310_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

JU310 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("JU310") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### JU2600 - p4 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/JU2600_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

JU2600 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("JU2600") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### MY2147 - p4 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/MY2147_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

MY2147 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("MY2147") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### NIC2 - p6 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/NIC2_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

NIC2 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("NIC2") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### MY2693 - p7 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/MY2693_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

MY2693 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("MY2693") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)


### EG4725 - p8 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/EG4725_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

EG4725 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("EG4725") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### JU2526 - p9 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/JU2526_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

JU2526 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("JU2526") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### QX1794 - p10 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/QX1794_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

QX1794 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("QX1794") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### DL238 - p11 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/DL238_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

DL238 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("DL238") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### CB4856 -p12 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/CB4856_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

CB4856 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("CB4856") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### ECA36 -p13 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/ECA36_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

ECA36 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("ECA36") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### ECA396 -p14 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/ECA396_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

ECA396 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.4) +
  geom_rect() +
  ggtitle("ECA396") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

### XZ1516 -p15 ####
df_pacbio <- read_tsv(file="../whole_genome_alignments/XZ1516_vs_N2.coords") %>%
  dplyr::rename(N2_start = '[S1]', N2_end = '[E1]', CB4856_start = '[S2]', CB4856_end = '[E2]', 
                N2_length = '[LEN 1]', CB4856_length = '[LEN 2]', Identity = "[% IDY]", 
                N2_chrom_length = '[LEN R]', CB4856_chrom_length='[LEN Q]', CHROM='[TAGS]')

XZ1516 <- df_pacbio %>%
  dplyr::filter(CHROM == chr, N2_length > 1000) %>%
  ggplot(.) +
  aes(xmin=N2_start/1e6, xmax=N2_end/1e6, ymax=Identity+0.5, ymin = Identity-0.5) +
  geom_rect() +
  ggtitle("XZ1516") + 
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_text(size=10, color='black'), 
        axis.title=element_text(size=11, color='black'), plot.margin = unit(c(0,2,0,2), "pt")) +
  labs(y="Identity (%)", fill = "") +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim=c(left/1e6, right/1e6)) +
  ylim(90,100.5)

# plot
pdf(file=paste(args[1], "_", args[2], "-", args[3], "_long_read_alignment.pdf", sep = ""), width=10.11, height=13.85)
grid.arrange(NIC2, ECA36, ECA396, QX1794, JU1400,  MY2693,  JU2526, JU2600, MY2147,  JU310,  NIC526, CB4856, EG4725, XZ1516, DL238, nrow = 15)
dev.off()

