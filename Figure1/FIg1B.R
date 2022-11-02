## Fig1B 
##  SNP and disease variant located in uORF
##  Only upload Freqtbl file

# library
library(ggplot2)
library(ggpubr) 
library(ggsci)

# Load Data
## intersect data
col.nm <- paste(rep(c("v", "u"), each = 6), rep(c("chromosome", "start", "stop", "name", "score", "strand"), 2), sep = "_")

all_bed <- unique(read.delim("All_5u_uORF.bed", col.names = col.nm))
snp_bed <- unique(read.delim("SNP151-uORF_hg38.bed", col.names = col.nm))

# Analysis
## Manipulate Table
invisible(lapply(c("all_", "snp_"), FUN = function(x){
  tmp <- get(paste0(x, "bed"))
  tmp$Same_Strand <- ifelse(tmp$v_strand==tmp$u_strand, "TRUE", "FALSE")
  tmp$v_score <- NULL -> tmp$u_score
  assign(paste0(x, "bed"), tmp, envir = .GlobalEnv)
}))

## Calculate Position
invisible(lapply(c("all_", "snp_"), FUN = function(x){
  tmp <- get(paste0(x, "bed"))
  tmp <- transform(tmp, 
                   Region = ifelse(u_strand=="+", 
                                   ifelse((v_stop-u_start)<4, "Start", ifelse((u_stop-v_start)<4, "Stop", "uCDS")), 
                                   ifelse((v_stop-u_start)<4, "Stop", ifelse((u_stop-v_start)<4, "Start", "uCDS"))),
                   Distance = ifelse(u_strand=="+", 
                                     round((v_stop-u_start)/(u_stop-u_start), digits = 3), 
                                     round((u_stop-v_start)/(u_stop-u_start), digits = 3))
  )
  assign(paste0(x, "bed"), tmp, envir = .GlobalEnv)
}))

# Plot
## Freqplot table
Freqtbl <- rbind(data.frame(Dist = all_bed$Distance, Grp = "All"), data.frame(Dist = snp_bed$Distance, Grp = "SNP151"))

## plot
snpUORF <- 
  ggplot(Freqtbl)+
  geom_freqpoly(aes(x = Dist, y = ..density.., group = Grp, color = Grp), bins = 100, size = 1.2, alpha = 0.7)+
  labs(x = "Relative Distance\n(Dist to uORF Start / Full Length uORF)", y = "Density (%)", title = "Variants located on uORF")+
  scale_color_manual(name = "Group", labels = c("All Collected 5' Variant (n=465)", "SNP151 (n=560474)"), values = pal_d3(palette = "category10")(4)[c(3,4)])+
  theme_pubclean()+
  theme(legend.position = "bottom")
