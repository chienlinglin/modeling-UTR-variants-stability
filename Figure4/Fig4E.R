require(readxl)
require(ComplexUpset)
require(ggplot2)
require(tidyverse)

sig_grps <- read_xlsx(file.path(getwd(), "fig_data_code", "fig_4E.xlsx"))


upset(data = sig_grps %>% as.data.frame(), 
      intersect = c("HEK 5' UTR",
                    "HEK 3' UTR",
                    "SH-SY5Y 5' UTR",
                    "SH-SY5Y 3' UTR") %>% rev, width_ratio = 0.4,
      name = "",
      sort_sets = FALSE,
      set_sizes=(
        upset_set_size(geom = geom_bar(width = 0.6, fill = "#2E3440"))
        + ylab("Number of selected features\n(ungrouped)")
        + geom_text(stat='count', aes(label=..count..), hjust = -1, color = "white")),
      base_annotations=list(
        'Composition ratio' =list(
          aes=aes(x=intersection, fill = Type),
          geom=list(
            geom_bar(stat='count', position='fill'), 
            scale_fill_manual(values = nord::nord("aurora")[1:4] %>% 
                                set_names(sig_grps$Type %>% unique),
                              name = ""),
            scale_y_continuous(labels=scales::percent_format())
          ) 
        ),
        'Size of the intersected\nfeatures (ungrouped)'=intersection_size(
          counts=FALSE,
          mapping=aes(fill=grp_TA)
        ) + scale_fill_manual(values = c('TA-correlated features' = nord::nord("lumina")[4],
                                         'Others' = "#2E3440"),
                              name = "")
      ))
