require(readxl)
require(SUmisc)
require(ggpubr)
require(wesanderson)
require(tidyverse)


df_sh_5utr_feat_plot <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_4F.xlsx"), sheet = "SH_UTR5")
df_sh_3utr_feat_plot <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_4F.xlsx"), sheet = "SH_UTR3")
p5_sign_df <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_4F.xlsx"), sheet = "SH_UTR5_pvs")
p3_sign_df <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_4F.xlsx"), sheet = "SH_UTR3_pvs")


p5 <- 
  ggplot(df_sh_5utr_feat_plot,  aes(x = sig_dir)
  ) + geom_col(aes(x = sig_dir, y = n, fill = kmer_TA_delta), position = 'dodge', width = 0.7) + 
  ggtitle("SH-SY5Y 5' UTR") +
  SUmisc:::theme_Publication() + 
  rremove("xlab") +
  rremove("ylab") +
  geom_text(data = p5_sign_df, aes(x = x, y = y, label = label)) +
  scale_fill_manual(name = expression(Delta*"TA-dimer"),
                    values = c(wes_palette("Darjeeling2")[2], 
                               wes_palette("GrandBudapest1")[2]))


p3 <- 
  ggplot(df_sh_3utr_feat_plot,  aes(x = sig_dir)
  ) + geom_col(aes(x = sig_dir, y = n, fill = kmer_TA_delta), position = 'dodge', width = 0.7) + 
  ggtitle("SH-SY5Y 3' UTR") +
  SUmisc:::theme_Publication() + 
  rremove("xlab") +
  rremove("ylab") +
  geom_text(data = p3_sign_df, aes(x = x, y = y, label = label)) +
  scale_fill_manual(name = expression(Delta*"TA-dimer"),
                    values = c(wes_palette("Darjeeling2")[2], 
                               wes_palette("GrandBudapest1")[2]))


p_sign_all <- 
  ggarrange(p5, p3, nrow = 1, common.legend = TRUE) %>%
  annotate_figure(bottom = text_grob("Testing for stability (Mutant vs. Wild-type)", 
                                     face = "bold", size = 15),
                  left = text_grob("Percentage", rot = 90,
                                   face = "bold", size = 15)) 

p_sign_all
