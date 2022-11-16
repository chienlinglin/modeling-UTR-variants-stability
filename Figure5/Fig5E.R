require(readxl)
require(wesanderson)
require(ggpubr)
require(EnvStats)
require(SUmisc)
require(tidyverse)


SH_ft_hl_tbl_utr5 <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_5E.xlsx"), sheet = "UTR5")
SH_ft_hl_tbl_utr3 <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_5E.xlsx"), sheet = "UTR3")


p5 <- ggplot(SH_ft_hl_tbl_utr5, 
             aes(x = sum_grp, y = logT05, fill = sum_grp)) +
  geom_boxplot(width = .5, alpha = 0.8,
               position = position_dodge(.4),
               show.legend = TRUE,
               outlier.size = 0.8) +
  stat_n_text(size = 7) +
  stat_compare_means(label.y = 6.3, size = 7) +
  SUmisc:::theme_Publication()  +
  ggtitle("SH-SY5Y 5' UTR") +
  rremove("xlab") +
  ylab(bquote(Log(t[1/2])-WT)) + 
  labs(fill="Number of RBP-binding per TA-dimer  ") +
  rremove("x.text") +
  rremove("x.ticks")  + 
  scale_fill_manual(values = wesanderson::wes_palette("GrandBudapest1")[2:1]) +
  theme(text = element_text(size = 22)) +
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size = 21), 
        legend.text = element_text(size = 18)) 

p3 <- ggplot(SH_ft_hl_tbl_utr3, 
             aes(x = sum_grp, y = logT05, fill = sum_grp)) +
  geom_boxplot(width = .5, alpha = 0.8,
               position = position_dodge(.4),
               show.legend = TRUE,
               outlier.size = 0.8) +
  stat_n_text(size = 7) +
  stat_compare_means(label.y = 6.3, size = 7) +
  SUmisc:::theme_Publication()  +
  ggtitle("SH-SY5Y 3' UTR") +
  rremove("xlab") +
  ylab(bquote(Log(t[1/2])-WT)) + 
  labs(fill="Number of RBP-binding per TA-dimer  ") +
  rremove("x.text") +
  rremove("x.ticks") + 
  scale_fill_manual(values = wesanderson::wes_palette("GrandBudapest1")[2:1]) +
  theme(text = element_text(size = 22)) +
  theme(legend.key.size = unit(1, 'cm'), 
        legend.key.height = unit(1, 'cm'), 
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size = 21), 
        legend.text = element_text(size = 18)) 


p_RBP_all <- 
  ggarrange(p5, p3, nrow = 1, 
            common.legend = TRUE, legend = "bottom") %>%
  annotate_figure(top = text_grob("eCLIP", face = "bold", size = 28))

p_RBP_all
