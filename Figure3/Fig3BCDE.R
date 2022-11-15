require(rio)
require(ggplot2)
require(wesanderson)
require(ggpubr)
require(grid)
require(SUmisc)
require(tidyverse)



df_sig_coef_plot <- import_list(file.path(getwd(), "fig_data_code", "fig_3BCDE.xlsx"))

p_top_coef <- 
  list(x = df_sig_coef_plot,
       y = names(df_sig_coef_plot),
       z = list(list(bb = c(seq(-0.2, 0.1, 0.1) %>% round(1), 0.2), ll = c(-0.2, 0.2)),
                list(bb = c(seq(-0.03, -0.01, 0.01), 0, 0.01, 0.02, 0.03), ll = c(-0.03, 0.03)),
                list(bb = c(seq(-0.3, 0.1, 0.1) %>% round(1), 0.2), ll = c(-0.3, 0.2)),
                list(bb = c(seq(-0.3, 0.1, 0.1) %>% round(1), 0.2), ll = c(-0.3, 0.1)))) %>%
  pmap(function(x, y, z){
      plot <- 
      ggplot(x, aes(x = Effect_size, 
                    y = reorder(Coefficient, abs(Effect_size)), 
                    color = Effect)) + 
      geom_pointrange(aes(xmin = Lower_CI, xmax = Upper_CI)) +
      geom_errorbar(aes(xmin = Lower_CI, xmax = Upper_CI), 
                    width = 0.25, cex = 0.5) +
      ggtitle(y) +
      scale_color_manual(values= c(wes_palette("GrandBudapest1", n = 4)[2],
                                   wes_palette("Zissou1", n = 5)[2])) +
      geom_vline(xintercept = 0, linetype="dashed") +
      SUmisc:::theme_Publication() +
      rremove("legend") +
      rremove("xlab") +
      rremove("ylab") +
      scale_x_continuous(breaks = z$bb, limits = z$ll) 
    
    return(plot)
  })

ggarrange(plotlist = p_top_coef, ncol = 2, nrow = 2, common.legend = TRUE, align = "v") %>%
  annotate_figure(bottom = textGrob("Beta coefficient", gp = gpar(cex = 1.3)),
                  left = textGrob("Feature", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
