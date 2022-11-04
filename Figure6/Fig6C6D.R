## Fig 6C & 6D
##  several example for meta-analysis UA-dimer content 
##  upload plot6C_DF, plot6D_DF, combined dataframe for all GO

# library
library(DBI)
library(RSQLite)
library(ggplot2)
library(ggthemes)
library(ggrepel)
library(ggpubr)

# prepare data for meta-analysis
dblink <- new.env()
dblink$GOdb <- dbConnect(SQLite(), "GOdb.db")
dblink$GOGCsliding <- dbConnect(SQLite(), "GOUA(pure)SlidingTest.db")
dblink$GOGCtest <- dbConnect(SQLite(), "GOUA(pure)test.db")

load("20220829_UA(only)slidingNorm.rds")  

allTSsliding <- 
  list(
    CDS = data.frame(freg = c(1:100), do.call(rbind, apply(allUAslidingNorm$cds, 2, simplify = FALSE, FUN = function(x){return(c(allMean = mean(x), allSD = sd(x)))}))),
    U3 = data.frame(freg = c(1:100), do.call(rbind, apply(allUAslidingNorm$U3, 2, simplify = FALSE, FUN = function(x){return(c(allMean = mean(x), allSD = sd(x)))}))),
    U5 = data.frame(freg = c(1:100), do.call(rbind, apply(allUAslidingNorm$U5, 2, simplify = FALSE, FUN = function(x){return(c(allMean = mean(x), allSD = sd(x)))})))
  )


# fig6C (5'UTR)
plot6C <- 
  sapply(paste0("GO:",c("0016255", "0006941", "0099536", "0008284")), simplify = FALSE, FUN = function(GO){
    utr <- "u5"
    
    #goName <- unlist(dbGetQuery(dblink$GOdb, sprintf("SELECT Name FROM ListBP WHERE Code=='%s'", GO)))
    goName <- spaceWarp(unlist(dbGetQuery(dblink$GOdb, sprintf("SELECT Name FROM ListBP WHERE Code=='%s'", GO))), Equal = FALSE)
    goMean <- unlist(dbGetQuery(dblink$GOGCsliding, sprintf("SELECT * FROM '%s' WHERE GO=='%s'", paste0("BP_", utr, "_Mean"), GO))[,-1])
    goSD <- unlist(dbGetQuery(dblink$GOGCsliding, sprintf("SELECT * FROM '%s' WHERE GO=='%s'", paste0("BP_", utr, "_SD"), GO))[,-1])
    goNum <- unlist(dbGetQuery(dblink$GOdb, sprintf("SELECT %s FROM TSamt WHERE GOid=='%s'", utr, GO)))
    
    plotDF <- 
      data.frame(
        group = rep(c("All", "GO"), each = 100), 
        freg = rep(c(1:100), 2),
        mean = c(allTSsliding[[toupper(utr)]]$allMean, goMean),
        sd = c(allTSsliding[[toupper(utr)]]$allSD, goSD)
      )
    
    plotDF$group <- factor(plotDF$group, levels = c("GO", "All"))
    
    plotDF$freg <- rep(c(-100:-1), 2)
    
    goPlot <- 
      ggplot(plotDF)+
      geom_line(aes(x = freg, y = mean, color =group), size = 1.8)+
      scale_color_manual(name = NULL, values = c("#ee7800", "#00a1e9"), labels = c("GO transcript", "All transcript"))+
      annotate("text",x = -99, y = Inf, label = goName, vjust = 1.25, hjust = 0, size = 4.8)+
      annotate("text",x = -15, y = Inf, label = sprintf("n=%s", goNum), vjust = 1.25, hjust = 0, size = 4.8)+
      labs(title = GO, x = NULL, y = NULL)+  
      theme_clean()+
      guides(color = guide_legend(byrow = TRUE))+
      theme(
        legend.position = "bottom", 
        legend.background = element_rect(colour = NA), 
        legend.text = element_text(size = 18),
        legend.spacing.x = unit(2, "char"),
        plot.title = element_text(size = 15),
        plot.background = element_rect(color = NA))
    
    return(goPlot)
  })


fig6C <- 
  ggarrange(NULL,
            ggarrange(plotlist = plot6C, align = "v", common.legend = TRUE, legend = "none", vjust = 0), nrow = 2, heights = c(0.5, 10))
  
fig6C <- 
  annotate_figure(fig6C, 
                  left = text_grob("TA-dimer ratio (Mean)", size = 15, rot = 90), 
                  bottom = text_grob("position", size = 15), fig.lab = "5'UTR", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")


# fig6D (3'UTR)
plot6D <- 
  sapply(paste0("GO:",c("0006614", "0034121", "0019083", "0002181")), simplify = FALSE, FUN = function(GO){
    utr <- "u3"
    
    goName <- spaceWarp(unlist(dbGetQuery(dblink$GOdb, sprintf("SELECT Name FROM ListBP WHERE Code=='%s'", GO))), Equal = FALSE)
    goMean <- unlist(dbGetQuery(dblink$GOGCsliding, sprintf("SELECT * FROM '%s' WHERE GO=='%s'", paste0("BP_", utr, "_Mean"), GO))[,-1])
    goSD <- unlist(dbGetQuery(dblink$GOGCsliding, sprintf("SELECT * FROM '%s' WHERE GO=='%s'", paste0("BP_", utr, "_SD"), GO))[,-1])
    goNum <- unlist(dbGetQuery(dblink$GOdb, sprintf("SELECT %s FROM TSamt WHERE GOid=='%s'", utr, GO)))

    plotDF <- 
      data.frame(
        group = rep(c("All", "GO"), each = 100), 
        freg = rep(c(1:100), 2),
        mean = c(allTSsliding[[toupper(utr)]]$allMean, goMean),
        sd = c(allTSsliding[[toupper(utr)]]$allSD, goSD)
      )
    
    plotDF$group <- factor(plotDF$group, levels = c("GO", "All"))
    
    
    goPlot <- 
      ggplot(plotDF)+
      geom_line(aes(x = freg, y = mean, color =group), size = 1.8)+
      scale_color_manual(name = NULL, values = c("#ee7800", "#00a1e9"), labels = c("GO transcript", "All transcript"))+
      annotate("text",x = 1, y = Inf, label = goName, vjust = 1.25, hjust = 0, size = 4.8)+
      annotate("text",x = 85, y = Inf, label = sprintf("n=%s", goNum), vjust = 1.25, hjust = 0, size = 4.8)+
      labs(title = GO, x = NULL, y = NULL)+  
      theme_clean()+
      theme(
        legend.position = "bottom", 
        legend.background = element_rect(colour = NA), 
        legend.text = element_text(size = 18),
        plot.title = element_text(size = 15),
        plot.background = element_rect(color = NA))
    
    return(goPlot)
  })


fig6D <- 
  ggarrange(NULL,
            ggarrange(plotlist = plot6D, align = "v", common.legend = TRUE, legend = "none", vjust = 0), nrow = 2, heights = c(0.5, 10))

fig6D <- 
  annotate_figure(fig6D, 
                  left = text_grob("TA-dimer ratio (Mean)", size = 15, rot = 90), 
                  bottom = text_grob("position", size = 15), fig.lab = "3'UTR", fig.lab.pos = "top.left", fig.lab.size = 18, fig.lab.face = "bold")
