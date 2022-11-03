## Figure 6E & 6F
##  GO (BP & MF) that both 5/3'UTR have significance high/low UA-dimer content
##  upload BPboth (fig6E) and MFboth (fig6F)

# library
library(DBI)
library(RSQLite)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(ggpubr)

# prepare data
dblink <- new.env()
dblink$GOdb <- dbConnect(SQLite(), "GOdb.db")
dblink$GOGCsliding <- dbConnect(SQLite(), "GOUA(pure)SlidingTest.db")
dblink$GOGCtest <- dbConnect(SQLite(), "GOUA(pure)test.db")

extTbl <- sapply(c("BP", "CC", "MF"), simplify = FALSE, FUN = function(GO){
  comm <- sprintf("SELECT Code,Name,GO_u5_mean,GO_u5_SD,All_u5_mean,All_u5_SD, `adj.pval_cds`, `adj.pval_u5`,`adj.pval_u3`, `adj.pval_all` ,GO_u3_mean,GO_u3_SD,All_u3_mean,All_u3_SD FROM %s WHERE `adj.pval_cds` > 0.01 AND (`adj.pval_u5` < 0.01 OR `adj.pval_u3` < 0.01 )", GO)
  
  return( dbGetQuery(dblink$GOGCtest, comm) )
})

# fig6E
 BPboth <- subset(extTbl$BP, `adj.pval_u5`<0.01&`adj.pval_u3`<0.01)
 invisible(lapply(c(3:6,10:13), FUN = function(idx){BPboth[,idx] <<- as.numeric(BPboth[,idx])}))
 
 BPboth <- transform(BPboth, 
                     adj.pval_u5 = ifelse(adj.pval_u5==0, 10^-20, adj.pval_u5),
                     adj.pval_u3 = ifelse(adj.pval_u3==0, 10^-20, adj.pval_u3),
                     diffU5 = GO_u5_mean - All_u5_mean,
                     diffU3 = GO_u3_mean - All_u3_mean,
                     area = ifelse((GO_u5_mean - All_u5_mean)>0 & (GO_u3_mean - All_u3_mean)>0, 1,
                                   ifelse((GO_u5_mean - All_u5_mean)<0 & (GO_u3_mean - All_u3_mean)>0, 2,
                                          ifelse((GO_u5_mean - All_u5_mean)<0 & (GO_u3_mean - All_u3_mean)<0, 3, 
                                                 ifelse((GO_u5_mean - All_u5_mean)>0 & (GO_u3_mean - All_u3_mean)<0, 4, NA)))))
 
 BPboth$area <- factor(BPboth$area, levels = c(1,3,2,4))
 
 GOBPplot <- 
   ggplot(BPboth)+
   geom_hline(yintercept = unique(BPboth$All_u3_mean), color = "grey30", linetype = "dashed", size = 0.9, alpha = 0.5)+
   geom_vline(xintercept = unique(BPboth$All_u5_mean), color = "grey30", linetype = "dashed", size = 0.9, alpha = 0.5)+
   geom_segment(aes(x = GO_u5_mean-GO_u5_SD, xend = GO_u5_mean+GO_u5_SD, y = GO_u3_mean, yend = GO_u3_mean), size = 0.3, col = "grey30", alpha = 0.5)+
   geom_segment(aes(x = GO_u5_mean, xend = GO_u5_mean, y = GO_u3_mean-GO_u3_SD, yend = GO_u3_mean+GO_u3_SD), size = 0.3, col = "grey30", alpha = 0.5)+
   geom_point(aes(x = GO_u5_mean, y = GO_u3_mean, size = -log10(adj.pval_u5), fill = -log10(adj.pval_u3)), color = "#2b2b2b", shape = 21)+
   geom_label_repel(aes(x = GO_u5_mean, y = GO_u3_mean, label = Name), max.overlaps = 10, size = 5.5, min.segment.length = 0.6, box.padding = 1, label.r = 0, label.size = NA, fill = alpha("white", 0.8))+
   scale_x_continuous(breaks = c(0,0.05,0.1,0.15), name = "5'UTR TA-dimer ratio (Mean\u00B1SD)")+
   scale_y_continuous(breaks = c(0,0.05,0.1,0.15), name = "3'UTR TA-dimer ratio (Mean\u00B1SD)", limits = c(0, NA))+
   scale_size(name = "5'UTR FDR\n(-log10)", range = c(1,10), breaks = c(5,10,15,20), guide = guide_legend(order = 1))+
   scale_fill_gradient(name = "3'UTR FDR\n(-log10)",limit = c(1,20), high = "#ff7f0e", low = "#FBFBFB", guide = guide_colorbar(reverse = FALSE))+
   theme_few()+
   labs(title = "Biologocal Process (BP)")+
   theme(panel.border = element_blank(), legend.position = "bottom", legend.spacing = unit(2, units = "cm"), plot.title = element_text(face = "bold", hjust = 0, size = 18))

 # fig6F
 MFboth <- subset(extTbl$MF, `adj.pval_u5`<0.01&`adj.pval_u3`<0.01)
 invisible(lapply(c(3:6,10:13), FUN = function(idx){MFboth[,idx] <<- as.numeric(MFboth[,idx])}))
 
 MFboth <- transform(MFboth, 
                     adj.pval_u5 = ifelse(adj.pval_u5==0, 10^-20, adj.pval_u5),
                     adj.pval_u3 = ifelse(adj.pval_u3==0, 10^-20, adj.pval_u3),
                     diffU5 = GO_u5_mean - All_u5_mean,
                     diffU3 = GO_u3_mean - All_u3_mean,
                     area = ifelse((GO_u5_mean - All_u5_mean)>0 & (GO_u3_mean - All_u3_mean)>0, 1,
                                   ifelse((GO_u5_mean - All_u5_mean)<0 & (GO_u3_mean - All_u3_mean)>0, 2,
                                          ifelse((GO_u5_mean - All_u5_mean)<0 & (GO_u3_mean - All_u3_mean)<0, 3, 
                                                 ifelse((GO_u5_mean - All_u5_mean)>0 & (GO_u3_mean - All_u3_mean)<0, 4, NA)))))
 
 MFboth$area <- factor(MFboth$area, levels = c(1,3,2,4))
 
 GOMFplot <- 
   ggplot(MFboth)+
   geom_hline(yintercept = unique(MFboth$All_u3_mean), color = "grey30", linetype = "dashed", size = 0.9, alpha = 0.5)+
   geom_vline(xintercept = unique(MFboth$All_u5_mean), color = "grey30", linetype = "dashed", size = 0.9, alpha = 0.5)+
   geom_segment(aes(x = GO_u5_mean-GO_u5_SD, xend = GO_u5_mean+GO_u5_SD, y = GO_u3_mean, yend = GO_u3_mean), size = 0.3, col = "grey30", alpha = 0.5)+
   geom_segment(aes(x = GO_u5_mean, xend = GO_u5_mean, y = GO_u3_mean-GO_u3_SD, yend = GO_u3_mean+GO_u3_SD), size = 0.3, col = "grey30", alpha = 0.5)+
   geom_point(aes(x = GO_u5_mean, y = GO_u3_mean, size = -log10(adj.pval_u5), fill = -log10(adj.pval_u3)), color = "#2b2b2b", shape = 21)+
   geom_label_repel(aes(x = GO_u5_mean, y = GO_u3_mean, label = Name), max.overlaps = 10, size = 5.5, min.segment.length = 0.6, box.padding = 1, label.r = 0, label.size = NA, fill = alpha("white", 0.8))+
   scale_x_continuous(breaks = c(0,0.05,0.1,0.15), name = "5'UTR TA-dimer ratio (Mean\u00B1SD)")+
   scale_y_continuous(breaks = c(0,0.05,0.1,0.15), name = "3'UTR TA-dimer ratio (Mean\u00B1SD)")+
   scale_size(name = "5'UTR FDR\n(-log10)", range = c(1,10), breaks = c(5,10,15,20), guide = guide_legend(order = 1))+
   scale_fill_gradient(name = "3'UTR FDR\n(-log10)",limit = c(1,20), high = "#ff7f0e", low = "#FBFBFB", guide = guide_colorbar(reverse = FALSE))+
   theme_few()+
   labs(title = "Molecular Function (MF)") +
   theme(panel.border = element_blank(), legend.position = "bottom", legend.spacing = unit(2, units = "cm"), plot.title = element_text(face = "bold", hjust = 0, size = 18))
 

