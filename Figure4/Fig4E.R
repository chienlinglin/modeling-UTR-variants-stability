require(readxl)
require(SUmisc)
require(ggpubr)
require(wesanderson)
require(tidyverse)


df_sh_3utr_feat_plot <- read_xlsx(file.path(getwd(), "fig_data_code_GR", "fig_4E.xlsx"), sheet = "SH_UTR3")
p3_sign_df <- read_xlsx(file.path(getwd(), "fig_data_code_GR", "fig_4E.xlsx"), sheet = "SH_UTR3_pvs")



p3 <- 
  ggplot(df_sh_3utr_feat_plot,  aes(x = sig_dir)
  ) + geom_col(aes(x = sig_dir, y = n, fill = kmer_TA_delta), position = 'dodge', width = 0.7) + 
  ggtitle("SH-SY5Y 3' UTR") +
  SUmisc:::theme_Publication() + 
  rremove("xlab") +
  ylab("Percentage") +
  ggtitle("MPRA results (Mutant vs. Wild-type)") +
  geom_text(data = p3_sign_df, aes(x = x, y = y, label = label)) +
  scale_fill_manual(name = expression(Delta*"TA-dimer"),
                    values = c(wes_palette("Darjeeling2")[2], 
                               wes_palette("GrandBudapest1")[2]))


p3
