require(readxl)
require(wesanderson)
require(ggpubr)
require(EnvStats)
require(tidyverse)


df_sh_5utr_feat_plot_5B <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_5B.xlsx"), sheet = "UTR5")
df_sh_3utr_feat_plot_5B <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_5B.xlsx"), sheet = "UTR3")

p_TAdelta_GC_5utr <- 
  ggplot(df_sh_5utr_feat_plot_5B, aes(x = kmer_TA_delta, y = logT05_delta)) + 
  geom_boxplot(aes(fill = kmer_TA_delta), width = .5, outlier.size = 0.8) +
  stat_summary(fun = median, geom = "line", aes(group = 1), linetype = "dashed") +
  facet_wrap(~GC_grp, strip.position = "bottom") + 
  xlab("Level of GC-content in wild-type") + 
  ylab(expression(Delta*Log(t["1/2"]))) +
  ggtitle("SH-SY5Y 5' UTR") +
  SUmisc:::theme_Publication() +
  stat_compare_means(label.y = 4) +
  stat_n_text() +
  scale_fill_manual(name = expression(Delta*"TA-dimer"),
                    values = c(wes_palette("Darjeeling2")[2], 
                               wes_palette("GrandBudapest1")[2])) +
  theme(strip.placement = "outside") +
  rremove("x.text") +
  rremove("x.ticks")

p_TAdelta_GC_3utr <- 
  ggplot(df_sh_3utr_feat_plot_5B, aes(x = kmer_TA_delta, y = logT05_delta)) + 
  geom_boxplot(aes(fill = kmer_TA_delta), width = .5, outlier.size = 0.8) +
  stat_summary(fun = median, geom = "line", aes(group = 1), linetype = "dashed") +
  facet_wrap(~GC_grp, strip.position = "bottom") + 
  xlab("Level of GC-content in wild-type") + 
  ylab(expression(Delta*Log(t["1/2"]))) +
  ggtitle("SH-SY5Y 3' UTR") +
  SUmisc:::theme_Publication() +
  stat_compare_means(label.y = 7) +
  stat_n_text() +
  scale_fill_manual(name = expression(Delta*"TA-dimer"),
                    values = c(wes_palette("Darjeeling2")[2], 
                               wes_palette("GrandBudapest1")[2])) +
  theme(strip.placement = "outside") +
  rremove("x.text") +
  rremove("x.ticks")

p_TAdelta_GC_all <- ggarrange(p_TAdelta_GC_5utr, p_TAdelta_GC_3utr, nrow = 1, common.legend = TRUE)

p_TAdelta_GC_all
