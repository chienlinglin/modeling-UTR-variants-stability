## Figure 2A 
##  UTR library half-life w./wo ARE motif in SH-SY5Y cell line

# library
library(DBI)
library(RSQLite)
library(coin)
library(ggplot2)
library(ggpubr) 
library(ggsci)  

# prepare data
## SH-SY5Y
stLink <- dbConnect(SQLite(), "SHSY5Yori_wlm_P0nI.db")
ftLink <- dbConnect(SQLite(), "LibFeature_220629.db")

SHAREYN <- 
  sapply(c("U5", "U3"), simplify = FALSE, FUN = function(utr){
    # load seqName/sequence/stability
    qt <- quantile(unlist(dbGetQuery(stLink, sprintf("SELECT MSE FROM 'QC_%s_HS'", utr))), ifelse(utr=="U5", 0.95, ifelse(utr=="U3", 0.8, NA)))
    utrSeq <- dbGetQuery(stLink, sprintf("SELECT QC_%s_HS.seqName,QC_%s_HS.t05 FROM QC_%s_HS INNER JOIN cpmQC_HS ON QC_%s_HS.seqName = cpmQC_HS.seqName WHERE t05 > 0 AND MSE < %s", utr, utr, utr, utr, qt))
    utrSeq <- subset(utrSeq, subset = log(t05)<quantile(log(t05), 0.975))
    
    # extract ARE info
    AREtbl <- 
      data.frame(
        dbGetQuery(ftLink, sprintf("SELECT * FROM ARE WHERE seqName IN ('%s')", paste0(utrSeq$seqName, collapse = "','"))), row.names = "seqName")
    
    AREtbl <- apply(AREtbl, MARGIN = 1, FUN = function(x){ifelse(all(x=="="),0,1)})
    AREtbl <- data.frame(seqName = names(AREtbl), grp = as.integer(AREtbl))
    
    # merge
    AREtbl <- merge(utrSeq, AREtbl, by = "seqName")
    
    # test
    AREtest <- oneway_test(AREtbl$t05~as.factor(AREtbl$grp))
    
    return(
      list(
        AREtbl = AREtbl, 
        AREtest = AREtest
      ))
  })


# plot
## SH-SY5Y
par(mfrow = c(2,1), mar = c(4,4,2,0), oma = c(1,1,0.5,0.5))

for(u in c("U5", "U3")){
  boxplot(log(SHAREYN[[u]]$AREtbl$t05)~SHAREYN[[u]]$AREtbl$grp, axes = FALSE, col = NA, boxcol = c("#1f77b4","#ff7f0e"), boxlwd = 2, boxwex = 0.4, medcol = c("#1f77b4","#ff7f0e"), whisklty = 1, whisklwd = 1.5,  outpch = 16,horizontal = TRUE, ylim = c(1,7), xlab = NA, ylab = NA)
  box(bty = "L", lwd = 2)
  axis(1, lwd = 0, lwd.ticks = 1.5)
  axis(2, lwd = 0, lwd.ticks = 1.5, at = c(1,2), labels = c("w/o\nARE", "w.\nARE"), las = 2)
  mtext(sprintf("SH-SY5Y %s'UTR", sub("U", "", u)), line = 0.5, adj = 0, cex = 1.2, font = 2)
  if(pvalue(SHAREYN[[u]]$AREtest)<2.2e-16){
    mtext("permutation test p.val < 2.22e-16", line = 0.5, adj = 0.4)  
  }else{
    mtext(sprintf("permutation test p.val = %.2e", pvalue(SHAREYN[[u]]$AREtest)), line = 0.5, adj = 0.4)  
  }
}

title(xlab = "half-life (log)", ylab = "with/without ARE", line = 0, outer = TRUE)

