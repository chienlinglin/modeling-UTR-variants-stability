require(SUmisc)
require(readxl)
require(ggplot2)

df_plot_utr5 <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_1E.xlsx"), sheet = "UTR5")
df_plot_utr3 <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_1E.xlsx"), sheet = "UTR3")


trend_plot_utr3 <- 
  ggplot(df_plot_utr3, aes(x = Time, y = mean_log_vt, color = Type, group = Type)) + 
  geom_line() +
  geom_errorbar(aes(ymin = mean_log_vt-log_vt_sd, ymax = mean_log_vt+log_vt_sd), width=.2,
                position=position_dodge(0.05))  +
  scale_x_continuous(limits = c(0, 95),
                     breaks = seq(0, 90, 45),
                     labels = seq(0, 90, 45)) +
  ylab("Log(Normalized Counts)") +
  facet_wrap(~Name) + 
  ggtitle("3' UTR") + 
  SUmisc:::theme_Publication() +
  xlab("minutes post-transfection") +
  scale_color_manual(values = c("black", "#990000"), labels = c("Wild-type", "Mutant")) +
  ylab("Log(Normalized Counts)") 

trend_plot_utr3

trend_plot_utr5 <- 
  ggplot(df_plot_utr5, aes(x = Time, y = mean_log_vt, color = Type, group = Type)) + 
  geom_line() +
  geom_errorbar(aes(ymin = mean_log_vt-log_vt_sd, ymax = mean_log_vt+log_vt_sd), width=.2,
                position=position_dodge(0.05)) +
  scale_x_continuous(limits = c(0, 95),
                     breaks = seq(0, 90, 45),
                     labels = seq(0, 90, 45)) +
  facet_wrap(~Name) + 
  ggtitle("5' UTR") + 
  SUmisc:::theme_Publication() +
  xlab("minutes post-transfection") +
  scale_color_manual(values = c("black", "#990000"), labels = c("Wild-type", "Mutant")) +
  ylab("Log(Normalized Counts)") 

trend_plot_utr5
