## figure 1B 
##  volcano plot show statistical significance group
##  upload SHtbl, HEKgroup

# library
library(DBI)
library(RSQLite)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(qvalue)

# plot
## SH-SY5Y
stLink <- dbConnect(SQLite(), "SHSY5Yori_wlm_P0nI.db")
SHtbl <- dbGetQuery(stLink, "SELECT DISTINCT Group_No,UTR,Allele,Time_vary,`adj.Time_vary`,coef FROM M2_U3 UNION SELECT DISTINCT Group_No,UTR,Allele,Time_vary,`adj.Time_vary`,coef FROM M2_U5")
SHtbl <- 
  do.call(rbind, lapply(split(SHtbl, SHtbl$Group_No), FUN = function(x){
    hl <- setNames(x$coef, x$Allele)[c("Ref", "Mut")]
    hl <- log(2)/-(hl)
    return(data.frame(grpNo = unique(x$Group_No),UTR = unique(x$UTR), aovPval =  unique(x$Time_vary), aovPvalFDR = unique(x$adj.Time_vary), quotHL = hl["Mut"]/hl["Ref"], wtHL = hl["Ref"], mtHL = hl["Mut"]))
  }))

SHtbl <- 
  do.call(rbind, lapply(split(SHtbl, SHtbl$UTR), FUN = function(x){
    x$aovQval <- qvalue(x$aovPval)$qvalues
    return(x)
  }))

SHtbl$Qcol <- factor(ifelse(SHtbl$aovQval>0.2, NA, ifelse(SHtbl$quotHL<2&SHtbl$quotHL>0.5, NA, ifelse(SHtbl$UTR=="U3UTR", "3'UTR", "5'UTR"))), levels = c("5'UTR", "3'UTR"))
SHtbl$UTR <- sub("UTR$", "'UTR", sub("^U", "", SHtbl$UTR))

SHtblp <- split(SHtbl, is.na(SHtbl$Qcol))

SHtblnum <- table(SHtbl$Qcol)

SHqplot <- 
  ggplot(SHtblp$`FALSE`)+
  geom_point(aes(x = log2(quotHL), y = -log10(aovQval), color = Qcol), alpha = 0.8, size = 2)+
  geom_point(data = SHtblp$`TRUE`, mapping = aes(x = log2(quotHL), y = -log10(aovQval)), alpha = 0.8, color = "grey50", show.legend = FALSE, size = 2)+
  geom_hline(yintercept = -log10(0.2), color = "grey70", linetype = 2)+
  geom_vline(xintercept = log2(0.5), color = "grey70", linetype = 2)+
  geom_vline(xintercept = log2(2), color = "grey70", linetype = 2)+
  coord_cartesian(xlim = c(-10,10))+
  ggsci::scale_color_d3(name = sprintf("UTR (n=%d)", nrow(SHtbl)), na.value = "grey50", label = c(sprintf("5'UTR (n=%d)", as.integer(SHtblnum["5'UTR"])), sprintf("3'UTR (n=%d)", as.integer(SHtblnum["3'UTR"]))))+
  theme_pubclean()+
  labs(x = "log2 halflife quotient (mt/WT)", y = "-log10 qvalue", subtitle = "SH-SY5Y")+
  theme(legend.direction = "vertical", legend.position = c(0.85,0.85), legend.title.align = 0.5, panel.grid.major.y = element_blank())

## HEK293
HEKgroup <- 
  do.call(rbind,
          lapply(c(1,2), FUN = function(pg){
            readxl::read_xlsx("df_test_t05.xlsx", sheet = pg)
          }))

HEKgroup <- 
  transform(HEKgroup,
            quotHL = T05_mt/T05_wt,
            sig = NULL,
            T05_mt = NULL,
            T05_wt = NULL)

HEKgroup$colq <- factor(ifelse(HEKgroup$qvalue>0.2, NA, ifelse(HEKgroup$quotHL>0.5&HEKgroup$quotHL<2, NA, ifelse(grepl("^3[MP]", HEKgroup$WT), "3'UTR", "5'UTR"))), levels = c("5'UTR", "3'UTR"))
HEKgroup$WT <- NULL    # too big, too long, not essential
HEKgroup$mt <- NULL

HEKqTbl <- split(HEKgroup, is.na(HEKgroup$colq))

HEKtblnum <- table(HEKqTbl$`FALSE`$colq)

HEKqPlot <- 
  ggplot(HEKqTbl$`FALSE`)+
  geom_point(aes(x = log2(quotHL), y = -log10(qvalue), color = colq), alpha = 0.8, size = 2)+
  geom_point(data = HEKqTbl$`TRUE`, mapping = aes(x = log2(quotHL), y = -log10(qvalue)), alpha = 0.8, color = "grey50", show.legend = FALSE, size = 2)+
  geom_hline(yintercept = -log10(0.2), color = "grey70", linetype = 2)+
  geom_vline(xintercept = log2(2), color = "grey70", linetype = 2)+
  geom_vline(xintercept = log2(0.5), color = "grey70", linetype = 2)+
  coord_cartesian(xlim = c(-10,10))+
  ggsci::scale_color_d3(name = sprintf("UTR (n=%d)", nrow(HEKgroup)), na.value = "grey50", label = c(sprintf("5'UTR (n=%d)", as.integer(HEKtblnum["5'UTR"])), sprintf("3'UTR (n=%d)", as.integer(HEKtblnum["3'UTR"]))))+
  theme_pubclean()+
  labs(x = "log2 halflife quotient (mt/WT)", y = "-log10 qvalue", subtitle = "HEK293")+
  theme(legend.direction = "vertical", legend.position = c(0.85,0.85), legend.title.align = 0.5, panel.grid.major.y = element_blank())
