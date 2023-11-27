require(ggsci)
require(ggpubr)
require(nord)
require(ghibli)
require(readxl)
require(tidyverse)

source("fig_data_code_GR/utils.R")

df_5H <- read_xlsx(file.path(getwd(), "fig_data_code_GR", "fig_5HI.xlsx"), sheet = "fig_H")
df_5I <- read_xlsx(file.path(getwd(), "fig_data_code_GR", "fig_5HI.xlsx"), sheet = "fig_I")



p_5H <- 
  ggplot(df_5H, aes(x = log_t05_lib, 
                    y = log_t05)) + 
  geom_point(aes(size = weighted_sum, color = GeneSymbol)) +
  theme_Publication() +
  scale_color_manual(values = ghibli_palette("PonyoMedium")[c(2:4, 6)], guide = "none") +
  geom_label(aes(size = NULL, color = GeneSymbol, label = GeneSymbol), 
             nudge_y = 0.2, nudge_x = 0.01, show_guide  = FALSE,  size = 5) +
  scale_size_continuous(range = c(3.5, 12), name = "#RBP per TA-dinucleotide") +
  ylab(bquote(Log(t[1/2])~"-RNA stability assay")) +
  xlab(bquote(Log(t[1/2])~"-MPRA")) +
  theme(legend.position = "top") +
  stat_cor(aes(label = ..r.label..), method = "spearman", 
           label.x = 2.4, label.y = 7.1, cor.coef.name = "rho", size = 6) +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.4)) 

p_5H

p_5I <- 
  ggplot(df_5I, 
         aes(x = Time, y = MEAN, color = GeneSymbol)) + 
  geom_line(linewidth = 0.8) + 
  geom_point() + 
  scale_x_continuous(breaks = c(0, 2, 4)) +
  scale_color_manual(values = ghibli_palette("PonyoMedium")[c(2:4, 6)]) +
  theme_Publication() +
  ylab("EGFP (± SE)") + 
  xlab("ActD treatment (hr)") +
  geom_errorbar(aes(ymin = MEAN - (SD/sqrt(3)), ymax = MEAN + (SD/sqrt(3))), width=.2,
                position = position_dodge(0.08), alpha = 0.75) +
  theme(legend.title = element_blank()) +
  ggtitle("RNA stability") +
  theme(legend.position = "top")

p_5I
