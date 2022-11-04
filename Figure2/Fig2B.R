## Figure 2G
##  LASSO estimate unique ARE pattern effect
##  upload AREtbl_U5, AREtbl_U3

# library
library(DBI)
library(RSQLite)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggsci)
library(glmnet)
library(parallel)


# extract data
## database link
dbLink <- dbConnect(SQLite(), "LibFeature_211214.db")
shData <- dbConnect(SQLite(), "SHSY5Yori_wlm_P0nI.db")

## SHSY5Y half-life data
invisible(lapply(c("U5", "U3"), FUN = function(utr){
  qt <- quantile(unlist(dbGetQuery(shData, sprintf("SELECT MSE FROM 'QC_%s_HS'", utr))), ifelse(utr=="U5", 0.95, ifelse(utr=="U3", 0.8, NA)))
  utrSeq <- dbGetQuery(shData, sprintf("SELECT QC_%s_HS.* FROM QC_%s_HS INNER JOIN cpmQC_HS ON QC_%s_HS.seqName = cpmQC_HS.seqName WHERE t05 > 0 AND MSE < %s", utr, utr, utr, qt))
  utrSeq <- subset(utrSeq, subset = log(t05)<quantile(log(t05), 0.975))
  utrSeq$logT05 <- log(utrSeq$t05)
  assign(paste0("sh", Biostrings::reverse(utr)), utrSeq, envir = .GlobalEnv)
}))

## ARE pattern amount in each UTR sequence
AREdetail <- dbReadTable(dbLink, "AREdetail")
row.names(AREdetail) <- AREdetail$seqName
AREdetail$seqName <- NULL

ARE3U <- AREdetail[sh3U$seqName,] 
ARE5U <- AREdetail[sh5U$seqName,] 
rm(AREdetail)

## bootstrap sampling + LASSO estimate motif coeff
invisible(lapply(c("ARE5U", "ARE3U"), FUN = function(x){                  ## table cleaning
  tmp <- get(x)
  tmp <- sapply(tmp, FUN = function(dat){
    dat[dat=="="] <- 0
    
    tmp.dat <- gregexpr(";", dat[dat!="0"])
    
    dat[dat!="0"] <- unlist(lapply(tmp.dat, length))
    
    return(as.integer(dat))
  })
  
  tmp <- tmp[,apply(tmp, MARGIN = 2, function(x)sum(x!=0))>20]
  
  assign(x, tmp, envir = .GlobalEnv)
}))

invisible(lapply(c("5U", "3U"), FUN = function(utr){                      ## bootstrap sampling + LASSO
  all.x <- get(paste0("ARE", utr))
  all.y <- log(get(paste0("sh", utr))$t05)
  
  bchData <- 
    do.call(rbind,
            mclapply(c(1:1000), mc.cores = 4, FUN = function(x){
              idx <- sample(length(all.y), size = length(all.y)*0.15, replace = TRUE)
              
              Tr_x <- all.x[-idx,]
              Tr_y <- all.y[-idx]
              
              Tt_x <- all.x[idx,]
              Tt_y <- all.y[idx]
              
              mdl <- glmnet(x = Tr_x, y = Tr_y, alpha = 1)
              cv.mdl <- cv.glmnet(x = Tr_x, y = Tr_y, alpha = 1)
              
              pred_y <- predict(mdl, s = cv.mdl$lambda.min, Tt_x)
              
              temp <- coef.glmnet(mdl, s = cv.mdl$lambda.min)
              coeff <- as.numeric(temp)
              
              #return(c(sprintf("bch%03d", x), caret::R2(pred_y, Tt_y, form = "traditional"), caret::RMSE(pred_y, Tt_y), coeff))
              return(c(caret::R2(pred_y, Tt_y, form = "traditional"), caret::RMSE(pred_y, Tt_y), coeff))
            }))  
  
  assign(paste0("mdlCoef", utr), bchData, envir = .GlobalEnv)
}))

## summerize data for visualization
AREplotTbl <- list()

### 3'UTR
coef3U <- mdlCoef3U[,-c(1:3)]

colnames(coef3U) <- colnames(ARE3U)

coef3U <- coef3U[,!apply(coef3U, 2, FUN = function(x)all(x==0))]

coef3UIdx <- unlist(mclapply(seq_len(ncol(coef3U)), mc.cores = 2,FUN = function(N){
  tmpDat <- c(as.numeric(coef3U[,N]), rep(0,nrow(coef3U)))
  tmpFac <- rep(factor(c("data", "zero"), c("data", "zero")), each = nrow(coef3U))
  
  tmpTest <- coin::oneway_test(tmpDat~tmpFac)
  return(coin::pvalue(tmpTest))
}))
coef3UIdx <- setNames(p.adjust(coef3UIdx, method = "bonferroni"), colnames(coef3U))
coef3U <- coef3U[,names(head(sort(coef3UIdx[coef3UIdx<0.05]), n=60))]

coef3UMean <- colMeans(coef3U) |> sort()

AREplotTbl$U3 <- 
  transform(
    data.frame(
      utr = "3'UTR",
      ptrn = chartr("T", "U", head(names(sort(abs(coef3UMean), decreasing = TRUE)), n = 20)), 
      beta = coef3UMean[head(names(sort(abs(coef3UMean), decreasing = TRUE)), n = 20)]),
    effect = ifelse(beta>0,"pos","neg"))

AREplotTbl$U3 <- AREplotTbl$U3[order(AREplotTbl$U3$beta),]
AREplotTbl$U3$ptrn <- factor(AREplotTbl$U3$ptrn, levels = AREplotTbl$U3$ptrn)

### 5'UTR
coef5U <- mdlCoef5U[,-c(1:3)]

colnames(coef5U) <- colnames(ARE5U)

coef5U <- coef5U[,!apply(coef5U, 2, FUN = function(x)all(x==0))]

coef5UIdx <- unlist(mclapply(seq_len(ncol(coef5U)), mc.cores = 2,FUN = function(N){
  tmpDat <- c(as.numeric(coef5U[,N]), rep(0,nrow(coef5U)))
  tmpFac <- rep(factor(c("data", "zero"), c("data", "zero")), each = nrow(coef5U))
  
  tmpTest <- coin::oneway_test(tmpDat~tmpFac)
  return(coin::pvalue(tmpTest))
}))

coef5UIdx <- setNames(p.adjust(coef5UIdx, method = "bonferroni"), colnames(coef5U))
coef5U <- coef5U[,names(head(sort(coef5UIdx[coef5UIdx<0.05]), n=60))]

coef5UMean <- colMeans(coef5U) |> sort()

AREplotTbl$U5 <- 
  transform(
    data.frame(
      utr = "5'UTR",
      ptrn = chartr("T", "U", head(names(sort(abs(coef5UMean), decreasing = TRUE)), n = 20)), 
      beta = coef5UMean[head(names(sort(abs(coef5UMean), decreasing = TRUE)), n = 20)]),
    effect = ifelse(beta>0,"pos","neg"))

AREplotTbl$U5 <- AREplotTbl$U5[order(AREplotTbl$U5$beta),]
AREplotTbl$U5$ptrn <- factor(AREplotTbl$U5$ptrn, levels = AREplotTbl$U5$ptrn)
                                                                                        
AREtbl <-                                                                                      
  sapply(AREplotTbl, simplify = FALSE, FUN = function(utrTbl){
    utrTbl <- transform(utrTbl, ptrn = as.character(ptrn), abs = abs(beta))
    utrTbl <- head(utrTbl[order(utrTbl$abs, decreasing = TRUE),], 10)
    utrTbl <- utrTbl[order(utrTbl$beta, decreasing = TRUE),]
    utrTbl$ptrn <- factor(utrTbl$ptrn, levels = rev(utrTbl$ptrn))
    
    return(utrTbl)
  })

## plot
AREplot <- list(U5 = 
                  ggplot(AREtbl$U5)+
                  geom_col(aes(y = beta, x = as.numeric(ptrn)+0.1, fill = effect), show.legend = FALSE, width = 0.4)+
                  geom_text(aes(y = 0, x = as.numeric(ptrn)-0.4, label = as.character(ptrn)), size = 4.2)+
                  scale_y_continuous(limits = c(-0.3,0.3), breaks = c(-0.2,-0.1,0.1,0.2))+
                  scale_x_continuous(breaks = NULL)+
                  scale_fill_d3()+
                  theme_pubclean()+
                  coord_flip()+  
                  theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 7), axis.ticks.length.y = unit(0,"points"), plot.background = element_rect(fill = "transparent"), panel.background = element_rect(fill="transparent"), panel.grid.major.x = element_line(colour = "grey75", linetype = 2, size = 0.4))+
                  labs(y = NULL, x = NULL, subtitle = "5'UTR"),
                U3 = 
                  ggplot(AREtbl$U3)+
                  geom_col(aes(y = beta, x = as.numeric(ptrn)+0.1, fill = effect), show.legend = FALSE, width = 0.4)+
                  geom_text(aes(y = 0, x = as.numeric(ptrn)-0.4, label = as.character(ptrn)), size = 4.2)+
                  scale_y_continuous(limits = c(-0.3,0.3), breaks = c(-0.2,-0.1,0.1,0.2))+
                  scale_x_continuous(breaks = NULL)+   
                  scale_fill_d3()+
                  theme_pubclean()+
                  coord_flip()+
                  theme(axis.text.y = element_blank(), axis.text.x = element_text(size = 7), axis.ticks.length.y = unit(0,"points"), plot.background = element_rect(fill = "transparent"), panel.background = element_rect(fill="transparent"), panel.grid.major.x = element_line(colour = "grey75", linetype = 2, size = 0.4))+
                  labs(y = NULL, x = NULL, subtitle = "3'UTR")
)

AREarrange <-  ggarrange(plotlist = AREplot, ncol = 2, align = "hv")
AREarrange <- annotate_figure(AREarrange, bottom = "coefficient")
