## Figure 7A & 7B
##  GDC data show potential UA-dimer disfunction
##  upload UCEC_DPP4_RNA, UCEC_DPP4_RPPA, STAD_CASP7_RNA, STAD_CASP7_RPPA

# library
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(coin)

# prepare data
UCEC <- as.environment(list(RNA = new.env(), RPPA = new.env()))
STAD <- as.environment(list(RNA = new.env(), RPPA = new.env()))
load("Project_DataSet_UCEC.rda", envir = UCEC$RNA)
load("Project-UCEC-20190214.rda", envir = UCEC$RPPA)
load("Project_DataSet_STAD.rda", envir = STAD$RNA)
load("Project-STAD-20190214.rda", envir = STAD$RPPA)

## DPP4
UCEC_DPP4 <- 
  list(
    RNA = subset(UCEC$RNA$Proj.Table$UCEC$RNA$`chr2-161992536-161992536T--`$Exp, Condition!="Control"),
    RPPA = UCEC$RNA$Proj.Table$UCEC$Protein$DPP4$Exp$`CD26.R.V$chr2-161992536-161992536T--`
  )

UCEC_DPP4$RNA$Condition <- factor(sub("Donater", "Donor", UCEC_DPP4$RNA$Condition), levels = c("Donor_wo", "Donor_w"))
UCEC_DPP4$RPPA$Condition <- factor(sub("Donater", "Donor", UCEC_DPP4$RPPA$Condition), levels = c("Donor_wo", "Donor_w"))

UCEC_DPP4_test <- 
  list(
    RNA = oneway_test(count~Condition, UCEC_DPP4$RNA),
    RPPA = oneway_test(Exp~Condition, UCEC_DPP4$RPPA)
  )

UCEC_DPP4_test <- 
  list(
    RNA_tbl = 
      data.frame(
        group1 = "Donor_wo", group2 = "Donor_w", p = pvalue(UCEC_DPP4_test$RNA), y.position = ceiling(max(log10(UCEC_DPP4$RNA$count+1)))*1.1), 
    RPPA_tbl = 
      data.frame(
        group1 = "Donor_wo", group2 = "Donor_w", p = pvalue(UCEC_DPP4_test$RPPA), y.position = ceiling(max(UCEC_DPP4$RPPA$Exp))*1.1)
  )


## CASP7
STAD_CASP7 <- 
  list(
    RNA = subset(STAD$RNA$Proj.Table$STAD$RNA$`chr10-113729598-113729598T--`$Exp, Condition!="Control"),
    RPPA = STAD$RNA$Proj.Table$STAD$Protein$CASP7$Exp$`Caspase.7_cleavedD198.R.C$chr10-113729598-113729598T--`
  )

STAD_CASP7$RNA$Condition <- factor(sub("Donater", "Donor", STAD_CASP7$RNA$Condition), levels = c("Donor_wo", "Donor_w"))
STAD_CASP7$RPPA$Condition <- factor(sub("Donater", "Donor", STAD_CASP7$RPPA$Condition), levels = c("Donor_wo", "Donor_w"))

STAD_CASP7_test <- 
    list(
      RNA = oneway_test(count~Condition, STAD_CASP7$RNA),
      RPPA = oneway_test(Exp~Condition, STAD_CASP7$RPPA)
    )
  
STAD_CASP7_test <- 
  list(
    RNA_tbl = 
      data.frame(
        group1 = "Donor_wo", group2 = "Donor_w", p = pvalue(STAD_CASP7_test$RNA), y.position = ceiling(max(log10(STAD_CASP7$RNA$count+1)))*1.1), 
    RPPA_tbl = 
      data.frame(
        group1 = "Donor_wo", group2 = "Donor_w", p = pvalue(STAD_CASP7_test$RPPA), y.position = ceiling(max(STAD_CASP7$RPPA$Exp))*1.1)
  )


# plot
## DPP4
UCEC_DPP4_plot <- 
  list(
    RNA = 
      ggplot(UCEC_DPP4$RNA)+
      geom_dotplot(aes(x = Condition, y = log10(count+1), fill = Condition, color = Condition), binaxis = "y", stackdir = "center", dotsize = 1, binwidth = 0.05, show.legend = FALSE)+
      geom_boxplot(aes(x = Condition, y = log10(count+1)), fill = NA, width = 0.125, size = 0.6,show.legend = FALSE, )+
      scale_x_discrete(label = c("A", "A>-"))+
      scale_y_continuous(limits = c(0, (ceiling(max(log10(UCEC_DPP4$RNA$count+1)))+1)))+
      stat_pvalue_manual(UCEC_DPP4_test$RNA_tbl, tip.length = 0, label = "pval = {sprintf('%.1e',p)}", bracket.size = 0.6)+
      scale_color_manual(values = c("#ff7f0e", "#1f77b4"))+
      scale_fill_manual(values = c("#ff7f0e", "#1f77b4"))+
      labs(x = "genotype of UCEC", y = "RNA Expression (log10)", subtitle = "DPP4[chr2:161992535-161992536(-)]")+
      theme_clean()
      ,
    RPPA = 
      ggplot(UCEC_DPP4$RPPA)+
      geom_dotplot(aes(x = Condition, y = Exp, fill = Condition, color = Condition), binaxis = "y", stackdir = "center", dotsize = 1, binwidth = 0.05, show.legend = FALSE)+
      geom_boxplot(aes(x = Condition, y = Exp), fill = NA, width = 0.125, size = 0.6,show.legend = FALSE)+
      scale_x_discrete(label = c("A", "A>-"))+
      scale_y_continuous(limits = c(-1.5, (ceiling(max(UCEC_DPP4$RPPA$Exp))+0.5)))+
      stat_pvalue_manual(UCEC_DPP4_test$RPPA_tbl, tip.length = 0, label = "pval = {sprintf('%.1e',p)}", bracket.size = 0.6)+
      scale_color_manual(values =c("#ff7f0e", "#1f77b4"))+
      scale_fill_manual(values =c("#ff7f0e", "#1f77b4"))+
      labs(x = "genotype of UCEC", y = "Protein (RPPA) Expression", subtitle = "DPP4[chr2:161992535-161992536(-)]")+
      theme_clean()
  )

## CASP7
STAD_CASP7_plot <-
  list(
    RNA = 
      ggplot(STAD_CASP7$RNA)+
      geom_dotplot(aes(x = Condition, y = log10(count+1), fill = Condition, color = Condition), binaxis = "y", stackdir = "center", dotsize = 0.7, binwidth = 0.05, show.legend = FALSE)+
      geom_boxplot(aes(x = Condition, y = log10(count+1)), fill = NA, width = 0.125, size = 0.6,show.legend = FALSE)+
      scale_x_discrete(label = c("T", "T>-"))+
      scale_y_continuous(limits = c(2, (ceiling(max(log10(STAD_CASP7$RNA$count+1)))+1)))+
      stat_pvalue_manual(STAD_CASP7_test$RNA_tbl, tip.length = 0, label = "pval = {sprintf('%.1e',p)}", bracket.size = 0.6)+
      scale_color_manual(values =c("#ff7f0e", "#1f77b4"))+
      scale_fill_manual(values =c("#ff7f0e", "#1f77b4"))+
      labs(x = "genotype of STAD", y = "RNA Expression (log10)", subtitle = "CASP7[chr10:113729597-113729598(+)]")+
      theme_clean()
      ,
    RPPA =
      ggplot(STAD_CASP7$RPPA)+
      geom_dotplot(aes(x = Condition, y = Exp, fill = Condition, color = Condition), binaxis = "y", stackdir = "center", dotsize = 1, binwidth = 0.05, show.legend = FALSE)+
      geom_boxplot(aes(x = Condition, y = Exp), fill = NA, width = 0.125, size = 0.6,show.legend = FALSE)+
      scale_x_discrete(label = c("T", "T>-"))+
      scale_y_continuous(limits = c(-2, (ceiling(max(STAD_CASP7$RPPA$Exp))+0.5)))+
      stat_pvalue_manual(STAD_CASP7_test$RPPA_tbl, tip.length = 0, label = "pval = {sprintf('%.1e',p)}", bracket.size = 0.6)+
      scale_color_manual(values =c("#ff7f0e", "#1f77b4"))+
      scale_fill_manual(values =c("#ff7f0e", "#1f77b4"))+
      labs(x = "genotype of STAD", y = "Protein (RPPA) Expression", subtitle = "CASP7[chr10:113729597-113729598(+)]")+
      theme_clean()
  )
