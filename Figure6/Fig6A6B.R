## Fig6A & 6B
##  Biological function have significance higher/lower UA dimer content in 5/3'UTR
##  upload BPU5tbl, BPU3tbl

# library
library(DBI)
library(RSQLite)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(ggpubr)

# prepare data from database
dblink <- new.env()
dblink$GOdb <- dbConnect(SQLite(), "GOdb.db")
dblink$GOGCsliding <- dbConnect(SQLite(), "GOUA(pure)SlidingTest.db")
dblink$GOGCtest <- dbConnect(SQLite(), "GOUA(pure)test.db")

# select CDS not significance & UTR significance GO
extTbl <- sapply(c("BP", "CC", "MF"), simplify = FALSE, FUN = function(GO){
  comm <- sprintf("SELECT Code,Name,GO_u5_mean,GO_u5_SD,All_u5_mean,All_u5_SD, `adj.pval_cds`, `adj.pval_u5`,`adj.pval_u3`, `adj.pval_all` ,GO_u3_mean,GO_u3_SD,All_u3_mean,All_u3_SD FROM %s WHERE `adj.pval_cds` > 0.01 AND (`adj.pval_u5` < 0.01 OR `adj.pval_u3` < 0.01 )", GO)
  
  return( dbGetQuery(dblink$GOGCtest, comm) )
})

# fig6A (5'UTR)

BPU5tbl <- head(extTbl$BP[order(extTbl$BP$adj.pval_u5),], 10)[,c("Code","Name","GO_u5_mean","GO_u5_SD","All_u5_mean","All_u5_SD","adj.pval_cds","adj.pval_u5","adj.pval_u3")]
BPU5tbl <- BPU5tbl[order(BPU5tbl$GO_u5_mean, decreasing = TRUE),]
invisible(lapply(c(3:6), FUN = function(idx){BPU5tbl[,idx] <<- as.numeric(BPU5tbl[,idx])}))

BPU5tbl$Name <- factor(BPU5tbl$Name, levels = BPU5tbl$Name[order(BPU5tbl$GO_u5_mean-BPU5tbl$All_u5_mean)])

fig6A <- 
  ggplot(BPU5tbl)+
  geom_vline(xintercept = unique(BPU5tbl$All_u5_mean), linetype = 2, color = "grey30")+
  geom_linerange(aes(x = GO_u5_mean, y = Name, xmin = GO_u5_mean-GO_u5_SD, xmax = GO_u5_mean+GO_u5_SD), color = "grey50", size = 1.1)+
  geom_point(aes(x = GO_u5_mean, y = Name, color = (GO_u5_mean-All_u5_mean)>0), size = 3.5, show.legend = FALSE)+
  geom_segment(aes(x = -0.20, xend = 0.3, y = 0, yend = 0), size = 1)+
  geom_text(aes(x = -0.55, y = Name, label = Name), hjust = 0, size = 5.3)+
  coord_cartesian(clip = "off")+
  scale_x_continuous(limits = c(-0.55, 0.35), breaks = seq(-0.2,0.3,0.1))+
  scale_color_manual(values = c("#e377c2", "#17becf"))+
  labs(x = "5'UTR TA-dimer ratio (Mean\u00B1SD)", y = NULL)+
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_line(size = 0.8),
    axis.ticks.length.x = unit(5, "points"),
    axis.text.x = element_text(size = 11),
    axis.ticks.length.y = unit(0, "points"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 15, vjust = -2, hjust = 0.75)  
  )


# fig6B (3'UTR)

BPU3tbl <- head(extTbl$BP[order(extTbl$BP$adj.pval_u3), ], 10)[,c("Code","Name","GO_u3_mean","GO_u3_SD","All_u3_mean","All_u3_SD","adj.pval_cds","adj.pval_u5","adj.pval_u3")]
BPU3tbl <- BPU3tbl[order(BPU3tbl$GO_u3_mean, decreasing = TRUE),]
invisible(lapply(c(3:6), FUN = function(idx){BPU3tbl[,idx] <<- as.numeric(BPU3tbl[,idx])}))

BPU3tbl$Name <- factor(BPU3tbl$Name, levels = BPU3tbl$Name[order(BPU3tbl$GO_u3_mean-BPU3tbl$All_u3_mean)])

fig6B <- 
  ggplot(BPU3tbl)+
  geom_vline(xintercept = unique(BPU3tbl$All_u3_mean), linetype = 2, color = "grey30")+
  geom_linerange(aes(x = GO_u3_mean, y = Name, xmin = GO_u3_mean-GO_u3_SD, xmax = GO_u3_mean+GO_u3_SD), color = "grey50", size = 1.1)+
  geom_point(aes(x = GO_u3_mean, y = Name, color = (GO_u3_mean-All_u3_mean)>0), size = 3.5, show.legend = FALSE)+
  geom_segment(aes(x = -0.10, xend = 0.3, y = 0, yend = 0), size = 1)+
  geom_text(aes(x = -0.45, y = Name, label = Name), hjust = 0, size = 5.3)+
  coord_cartesian(clip = "off")+
  scale_x_continuous(limits = c(-0.45, 0.35), breaks = seq(-0.1,0.3,0.1))+
  scale_color_manual(values = c("#17becf", "#e377c2"))+
  labs(x = "3'UTR TA-dimer ratio (Mean\u00B1SD)", y = NULL)+
  theme(
    panel.background = element_blank(),
    axis.ticks.x = element_line(size = 0.8),
    axis.ticks.length.x = unit(5, "points"),
    axis.ticks.length.y = unit(0, "points"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 15, vjust = -2, hjust = 0.8)  
  )
