## Figure 1A 
##  SNP / disease variant locate in mRNA 
##  Due to file size, only upload summerized file (SNP151, ClinNew)


# library
library(ggsci)
library(ggplot2)

# data prepare
## SNP151
SNP151pos <- 
  rbind(
    read.delim("snp151Ref3U", header = FALSE),
    read.delim("snp151Ref5U", header = FALSE))

SNP151pos$utr <- regmatches(SNP151pos$V14, regexpr("utr[35]", SNP151pos$V14))  # extract UTR information

SNP151tbl <-                                                                   # split by strand
    do.call(rbind, lapply(
      split(SNP151pos, SNP151pos$V16), FUN = function(x){
        sig <- unique(x$V16)
        
        x$allDist <- x$V13-x$V12
        if(sig=="+"){
          x$toStr <- x$V2-x$V12
        }else if(sig=="-"){
          x$toStr <- x$V13-x$V3
        }
      
      x$perTo5 <- x$toStr/x$allDist
      
      return(x) 
    }
  ))[,c("utr", "perTo5")]

SNP151tbl <- subset(SNP151tbl, perTo5>=0)                                      # clean negative data

SNP151tbl <- split(SNP151tbl, SNP151tbl$utr)

SNP151 <-                                                                      # metagene position calculation
  do.call(rbind, lapply(SNP151tbl, FUN = function(x){
    utr <- unique(x$utr)
    freqTbl <- hist(x$perTo5, breaks = 101, plot = FALSE)
    
    tbl <- data.frame(pos = tail(freqTbl$breaks, 100), cnt = freqTbl$counts, density = freqTbl$density, utr = utr, data = "SNP151")
  
    tbl$smoothD <- predict(smooth.spline(tbl$pos, tbl$density, spar = 0.3))$y
  
    if(utr=="utr5"){
      tbl$pos <- tbl$pos*100
    }else if(utr=="utr3"){
      tbl$pos <- (tbl$pos*100)+150
    }
    
    return(tbl)
  }))

## disease variant
loadBED <- read.delim("utrInterSected.BED", header = FALSE)
loadBED$varTx <- sub("\\.[0-9]+", "", sapply(strsplit(loadBED$V6, split = "\\("), "[[", 1))
loadBED$utrTx <- sub("\\.[0-9]+", "", sub("_utr[35]+.+$", "", loadBED$V13))    

loadBED <- loadBED[loadBED$varTx==loadBED$utrTx,]
loadBED$utr <- regmatches(loadBED$V13, regexpr("utr[35]", loadBED$V13))

calcData <- 
  do.call(rbind, lapply(
    split(loadBED, loadBED$V15), FUN = function(x){
      sig <- unique(x$V15)
      
      x$allDist <- x$V12-x$V11
      if(sig=="+"){
        x$toStr <- x$V2-x$V11
      }else if(sig=="-"){
        x$toStr <- x$V12-x$V3
      }
      
      x$perTo5 <- x$toStr/x$allDist
      
      return(x) 
    }
  ))
  
calcData <- subset(calcData, perTo5>=0)
calcData <- split(calcData, calcData$utr)

ClinNew <- 
  do.call(rbind, lapply(calcData, FUN = function(x){
    utr <- unique(x$utr)
    freqTbl <- hist(x$perTo5, breaks = 101, plot = FALSE)
    
    tbl <- data.frame(pos = tail(freqTbl$breaks, 100), cnt = freqTbl$counts, density = freqTbl$density, utr = utr, data = "ClinVar(20220901)")
    
    tbl$smoothD <- predict(smooth.spline(tbl$pos, tbl$density, spar = 0.3))$y
    
    if(utr=="utr5"){
      tbl$pos <- tbl$pos*100
    }else if(utr=="utr3"){
      tbl$pos <- (tbl$pos*100)+150
    }
    
    return(tbl)
  }))


## CDS  ## create empty/spacing dataframe
CDS <- 
  data.frame(pos = c(101:150), cnt = NA_integer_, density = NA_integer_, utr = NA_character_, data = NA_character_, smoothD = NA_integer_)


## combine
combNew <- rbind(SNP151, CDS, ClinNew)
combNew <- transform(combNew, utr = factor(utr, levels = c("utr5", "utr3")), data = factor(data, levels = c("SNP151", "ClinVar(20220901)")))
combNew <- combNew[order(combNew$pos),]


# plot
combNewPlot <- 
  ggplot(combNew)+
  geom_line(data = subset(combNew, utr=="utr5"), aes(x = pos, y = smoothD, color = data, group = data), size = 1.1, alpha = 0.75)+
  geom_line(data = subset(combNew, utr=="utr3"), aes(x = pos, y = smoothD, color = data, group = data), size = 1.1, alpha = 0.75)+
  scale_y_continuous(name = "density")+
  scale_x_continuous(name = "region", breaks = c(0,100,150,250), labels = c(-100,0,0,100))+
  scale_color_d3(labels = c(sprintf("SNP151\n5'UTR n=%s\n3'UTR n=%s", nrow(SNP151tbl$utr5), nrow(SNP151tbl$utr3)), sprintf("ClinVar(20220901)\n5'UTR n=%s\n3'UTR n=%s", nrow(calcData$utr5), nrow(calcData$utr3))))+
  annotate("text", x = c(50, 125, 200), y = -Inf, label = c("5'UTR", "CDS", "3'UTR"), vjust = -0.01, size = 5)+
  ggpubr::theme_pubclean()+
  theme(legend.position = "bottom")
