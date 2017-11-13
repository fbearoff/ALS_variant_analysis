library(vcfR)
library(reshape2)
library(ggplot2)

vcf <- read.vcfR("BB.Chr17.recode.vcf")
dna <- ape::read.dna("/home/frank/GRCm38_construct/Mus_musculus.GRCm38.dna.chromosome.17.fa", format = "fasta")
gff <- read.table("/home/frank/GRCm38_construct/Mus_musculus.GRCm38.89.gff3.gz", sep="\t", quote="")
gff2 <- gff[grep("17", gff[,1]),]

#depth violin plot
dp <- extract.gt(vcf, element = 'DP', as.numeric = TRUE)
dpf <- melt(dp, varnames=c('Index', 'Sample'), value.name="Depth", na.rm=TRUE)
dpf <- dpf[ dpf$Depth > 0,]
p <- ggplot(dpf, aes(x=Sample, y=Depth)) +
  geom_violin(fill='#C0C0C0', adjust=1.0,
            scale = 'count', trim=TRUE)
p <- p + theme_bw()
p <- p + ylab('Read Depth (DP)')
p <- p + theme(axis.title.x = element_blank(),
               axis.text.x = element_text(angle = 60, hjust = 1))
p <- p + stat_summary(fun.data=mean_sdl,
                      geom='pointrange', color='black')
p <- p + scale_y_continuous(trans=scales::log2_trans(),
                            breaks=c(1, 10, 100, 1000))
p

heatmap.bp(dp[501:1500,], rlabels = FALSE)

chrom <- create.chromR(name="BB", vcf = vcf, seq = dna, ann = gff2)
chrom <- proc.chromR(chrom)
chromoqc(chrom, dp.alpha=20) 
#chromoqc(chrom, xlim=c(34956436, 34959238), dp.alpha=20) #xlim determines interval graphed
plot(chrom)

#setup some baseline filters
chrom_mask <- masker(chrom, min_QUAL=30, min_DP=3, max_DP=100, min_MQ=59, max_MQ=61)
plot(chrom_mask)
