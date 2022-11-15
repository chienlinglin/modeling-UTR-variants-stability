require(rio)
require(iheatmapr)
require(RColorBrewer)
require(nord)
require(plotly)
require(tidyverse)




df_heatmap <- 
  import_list(file.path(getwd(), "fig_data_code", "fig_4BD.xlsx")) %>%
  map(~.x %>% 
        column_to_rownames("...1"))

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
p <- 
  df_heatmap %>%
  map2(c("HEK 5' UTR",
         "HEK 3' UTR",
         "SH-SY5Y 5' UTR",
         "SH-SY5Y 3' UTR"),
       function(plot_df, p_title){
         main_heatmap(plot_df[,-(1:2)] %>% as.matrix() %>% t, zmid = 0, zmin = -3, zmax = 3, 
                      name = "Standardized level", colors =  myPalette(100) %>% replace(51, "white")) %>%
           add_row_labels(font = list(size = 14))  %>% 
           add_col_signal(signal = plot_df[,2], name = "Log(t<sub>1/2</sub>)", 
                          show_title = FALSE, colors = c(nord("snowstorm"), rev(nord("polarnight")))) %>%
           add_col_title(p_title, font = list(size = 20), side = "top", buffer = 0.05) %>%
           add_main_heatmap(plot_df[,-(1:2)] %>% as.matrix() %>% cor(method = "spearman"), 
                            name = "Spearman's correlation<br>between features", 
                            colors = c(low = "#0C6291", mid = "#FBFEF9", mid = "#FBFEF9", high = "#A63446"), 
                            zmid = 0, zmin = -1, zmax = 1) %>%
           modify_layout(list(font = list(size = 16),
                              margin = list(l = 100, r = 80))) 
       })
p
