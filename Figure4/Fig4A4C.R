## Fig4A (5'UTR) & Fig4C (3'UTR)
##  TA content correlate with half-life, SH-SY5Y, Both 5/3'UTR
##  upload SHdiNT


#  library
library(DBI)
library(RSQLite)
library(ggplot2)
library(ggpubr)
library(ggthemes)

#  function
sequence_kmers <- function(sequence, k){
  k_mers <- lapply(sequence,function(x){
    seq_loop_size <- length(Biostrings::RNAString(x))-k+1
    
    kmers <- sapply(1:seq_loop_size, function(z){
      y <- z + k -1
      kmer <- substr(x=x, start=z, stop=y)
      return(kmer)
    })
    return(kmers)
  })
  
  uniq <- unique(unlist(k_mers))
  ind <- t(sapply(k_mers, function(x){
    tabulate(match(x, uniq), length(uniq))
  }))
  colnames(ind) <- uniq
  
  return(ind)
}

# prepare database
ftLink <- dbConnect(SQLite(), "LibFeature_220331.db")
stLink <- dbConnect(SQLite(), "SHSY5Yori_wlm_P0nI.db")

# sequence data
SHGCFE <- sapply(c("U5", "U3"), simplify = FALSE, FUN = function(utr){
  ## collect  
  qt <- quantile(unlist(dbGetQuery(stLink, sprintf("SELECT MSE FROM 'QC_%s_HS'", utr))), ifelse(utr=="U5", 0.95, ifelse(utr=="U3", 0.8, NA)))
  utrSeq <- dbGetQuery(stLink, sprintf("SELECT QC_%s_HS.* FROM QC_%s_HS INNER JOIN cpmQC_HS ON QC_%s_HS.seqName = cpmQC_HS.seqName WHERE t05 > 0 AND MSE < %s", utr, utr, utr, qt))
  utrSeq <- subset(utrSeq, subset = log(t05)<quantile(log(t05), 0.975))
  utrSeq$logT05 <- log(utrSeq$t05)
  ## merge sequence data
  utrSeq <- merge(x = utrSeq, y = dbGetQuery(ftLink, sprintf("SELECT * FROM seq_FE_GC WHERE seqName IN ('%s')", paste0(utrSeq$seqName, collapse = "','"))), by = "seqName")
  ## remove CDS sequence (because of primer)
  utrSeq$sequence <- sub(ifelse(utr=="U5", "AUGGUGA$", "^GGACGAGCUGUACAAGUAA"), "",utrSeq$sequence)
  
  return(utrSeq)
})

# calculate TA content
SHdiNT <- lapply(SHGCFE, FUN = function(x){
  tmp <- subset(x, select = c("seqName", "sequence", "GCcontent", "t05", "logT05"))
  tmp <- transform(tmp, halfLife = t05, logHalfLife = logT05, t05 = NULL, logT05 = NULL)
  
  diNtTbl <- sequence_kmers(tmp$sequence, 2)
  diNTSum <- rowSums(diNtTbl)
  
  diNt <- 
    data.frame(UA = diNtTbl[,"UA"], AU = diNtTbl[,"AU"], UAAU = rowSums(diNtTbl[,c("UA", "AU")]))
  
  diNt <- 
    transform(diNt,
              UAcontent = UA/diNTSum, 
              AUcontent = AU/diNTSum, 
              UAAUcontent = UAAU/diNTSum, NTSum = diNTSum)
  
  return(cbind(tmp, diNt))
})

dbDisconnect(ftLink)
dbDisconnect(stLink)

# corplot
## 5'UTR (Fig4A)
ggplot(SHdiNT$U5, aes(x = UAcontent, y = logHalfLife))+
  geom_jitter(width = 0.005, color = "#1f77b4", alpha = 0.7)+
  geom_smooth(method = "lm", formula = y~x, color = "#ff7f0e")+
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, 
           label = sprintf("Pearson's corr = %.3f\n\tpval = %.3e", cor.test(SHdiNT$U5$UAcontent, SHdiNT$U5$logHalfLife)$estimate, cor.test(SHdiNT$U5$UAcontent, SHdiNT$U5$logHalfLife)$p.value), size = 5)+
  theme_few()+
  theme(panel.border = element_rect(size = 1.5))+
  labs(x = "UA content", y = "half life (log)", subtitle = "SH-SY5Y 5'UTR")


## 3'UTR (Fig4C)
ggplot(SHdiNT$U3, aes(x = UAcontent, y = logHalfLife))+
  geom_jitter(width = 0.005, color = "#1f77b4", alpha = 0.7)+
  geom_smooth(method = "lm", formula = y~x, color = "#ff7f0e")+
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, 
           label = sprintf("Pearson's corr = %.3f\n\tpval = %.3e", cor.test(SHdiNT$U3$UAcontent, SHdiNT$U3$logHalfLife)$estimate, cor.test(SHdiNT$U3$UAcontent, SHdiNT$U3$logHalfLife)$p.value), size = 5)+
  theme_few()+
  theme(panel.border = element_rect(size = 1.5))+
  labs(x = "UA content", y = "half life (log)", subtitle = "SH-SY5Y 3'UTR")

