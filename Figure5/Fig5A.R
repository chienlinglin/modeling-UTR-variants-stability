require(rio)
require(SUmisc)
require(ggpubr)
require(tidyverse)

df_GC_TA_SH <- import_list(file.path(getwd(), "fig_data_code", "fig_5A.xlsx"))

draw_scatter_TA <- 
  function(x_axis, grp_by, top_title,
           label.x.npc = 0.76,
           label.y.npc = 0.99, df){
    
    df %>%
      map(function(df_plot){
        (ggscatter(df_plot, x = x_axis, y = "logT05",
                   add = "reg.line",
                   add.params = list(color = "blue", fill = "lightgray"),
                   conf.int = TRUE 
        ) + stat_cor(method = "pearson",
                     label.x.npc = label.x.npc,
                     label.y.npc = label.y.npc)) %>%
          facet(facet.by = grp_by, scales = 'free') +
          SUmisc:::theme_Publication() +
          theme(text = element_text(size = 18),
                axis.title.y = element_text(face = "bold"),
                plot.title = element_text(size = 20,
                                          face = "bold"))+ 
          ylab(bquote(bold(Log(t["1/2"]))))
      }) %>%
      map2(names(df), ~(.x + ggtitle(paste(top_title, .y))))
    
  }


draw_scatter_TA(x_axis = "kmer_TA_ratio",
                grp_by = "GC_grp",
                top_title = "SH-SY5Y",
                label.x.npc = 0.7,
                label.y.npc = 0.95,
                df = df_GC_TA_SH)
