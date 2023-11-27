require(ggrepel)
require(ggpubr)
require(tidyverse)

source("fig_data_code_GR/utils.R")

list_df_fig5J <- readRDS("fig_data_code_GR/fig_5J.rds")


list_df_fig5J %>%
  pmap(function(allDf, names_lst, filteredDf){
    ggplot(allDf, aes(x = Coef, y = -log10(pvalue), label = Gene)) +
      geom_point() + 
      geom_label_repel(data         = filteredDf,
                       size          = 4,
                       box.padding   = 0.5,
                       point.padding = 0.5,
                       force         = 10,
                       segment.size  = 0.2,
                       segment.color = "grey50",
                       max.overlaps = 15) +
      ggtitle(names_lst %>% str_replace_all("\\.", " - ")) +
      theme_Publication() +
      xlab("Effect on half-life") +
      theme_Publication()
  }) %>%
  ggarrange(plotlist = ., nrow = 2, ncol = 2)
