library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(MetBrewer)
library(ragg)
library(ggpubr)

# get date
date()

setwd(file.path("C:/Users/User/OneDrive/Documenten/Documenten-UU/Master/MCLS/", 
                "2022-2023/Writing_assignment/Data_visualization"))
Data_path <- file.path("C:/Users/User/OneDrive/Documenten/Documenten-UU/Master", 
                  "MCLS/2022-2023/Writing_assignment/Data_visualization/Data")
Figure_path <- file.path("C:/Users/User/OneDrive/Documenten/Documenten-UU/Master", 
                       "MCLS/2022-2023/Writing_assignment/Data_visualization/Figures")

## load critical appraisal sheet
appraisal <- read_delim(file.path(Data_path, "Critical_appraisal_sheet.tsv"), 
                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)
## load EWAS_list
EWAS_list <- read_delim(file.path(Data_path, "EWAS_list.tsv"), delim = "\t", 
                        escape_double = FALSE, trim_ws = TRUE)
EWAS_list %<>%
  as_tibble() %>%
  filter(!is.na(Exposure)) %>%
  mutate(Exposure = ifelse(Exposure == "Chlorinated persistent pollutants (POPs)", "POPs", 
                           Exposure))

# join 'Exposure' column to appraisal tibble via 'key' column
appraisal <- EWAS_list %>%
  select(key, Exposure, EWAS) %>%
  .[!(duplicated(EWAS_list$key)),] %>%
  full_join(appraisal, by = join_by(key)) %>%
  # manually add the categories of the three studies not performing EWAS and
  # only predicting environmental exposures
  mutate(Exposure = ifelse(key %in% c("rayyan-1000417471", "rayyan-1000417517",
                                      "rayyan-1000417624"), "Smoking", Exposure)) %>%
  # remove 'Confounders adjusted' column
  select(!`Confounders adjusted`)

# get critical appraisal category counts
category_counts <- appraisal %>%
  as_tibble() %>%
  select(!(c(key, Exposure, EWAS, Citation, title, `Sample size (discovery)`, Notes))) %>%
  mutate(`classifier prediction tested` = ifelse(is.na(`classifier prediction tested`), "-", `classifier prediction tested`)) %>%
  mutate(across(everything(), as.factor)) %>%
  pivot_longer(cols = everything()) %>%
  group_by(name, value) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = value, values_from = n) %>%
  mutate(`+/-` = ifelse(is.na(`+/-`), 0, `+/-`)) %>%
  rename(`Critical appraisal` = name)

# prepare critical appraisal df for plotting
factor_counts <- appraisal %>%
  as_tibble() %>%
  select(!(c(key, Exposure, EWAS, Citation, title, `Sample size (discovery)`, Notes))) %>%
  mutate(`classifier prediction tested` = ifelse(is.na(`classifier prediction tested`), "-", `classifier prediction tested`)) %>%
  mutate(across(everything(), as.factor)) %>%
  pivot_longer(cols = everything())

## create plots
# critical appraisal barplots
p1 <- ggplot(data = factor_counts, aes(x = value)) + 
  geom_bar(aes(fill = value)) + 
  facet_wrap(~factor(name, levels = c("Demographics", "Sample size calc", 
                                      "Quality control", "Confounder correction",
                                      "blood cell correction", "DNAm measuring technique",
                                      "environment collection", "replication attempt",
                                      "classifier prediction tested", "Meta-analysis"))) +
  scale_fill_manual(values = met.brewer(name = "Manet", n = 4, type = "continuous")) +
  theme(axis.text = element_text(size = 15, face = "bold"), axis.title.x = element_blank(),
        legend.position = c(0.93, 0.12), 
        legend.text = element_text(size = 15, face = "bold"),
        axis.title.y = element_text(size = 15),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 13)) +
  guides(fill = guide_legend(title = "Appraisal\ncategories")) +
  ylab("Number of articles")
p1

# combine p1 and category counts
table_1 <- ggtexttable(category_counts, rows = NULL, 
                       theme = ttheme("lBlack", 
                                      colnames.style = colnames_style(size = 12),
                                      tbody.style = tbody_style()))
table_1 <- table_1 %>% table_cell_font(row = 2:tab_nrow(table_1), 
                                       column = 1:tab_ncol(table_1), size = 12)
table_1

p1b <- ggarrange(p1, table_1, nrow = 2, ncol = 1, labels = c("A", "B"),
          heights = c(2,1))

# sample size histogram
p2 <- ggplot(data = appraisal, aes(x = `Sample size (discovery)`)) + 
  geom_histogram(binwidth = 500) +
  theme(axis.text = element_text(size = 12, face = "bold")) +
  xlab("Study sample size") +
  ylab("Frequency")
p2

p3 <- ggplot(data = appraisal, aes(x = `Sample size (discovery)`)) + 
  geom_histogram(aes(fill = Exposure), binwidth = 500, ) +
  scale_fill_manual(values = met.brewer(name = "Manet", 
                                        n = length(unique(appraisal$Exposure)), 
                                        type = "continuous")) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 11, face = "italic")) +
  labs(title = "binwidth = 500") +
  xlab("Study sample size") +
  ylab("Frequency")
p3

p3b <- ggplot(data = appraisal, aes(x = `Sample size (discovery)`)) + 
  geom_histogram(aes(fill = Exposure), binwidth = 30, ) +
  scale_fill_manual(values = met.brewer(name = "Manet", 
                                        n = length(unique(appraisal$Exposure)), 
                                        type = "continuous")) +
  theme(axis.text = element_text(size = 12, face = "bold"), 
        plot.title = element_text(size = 11, face = "italic")) +
  labs(title = "binwidth = 30") +
  xlab("Study sample size") +
  ylab("Frequency") +
  xlim(c(0,1000))
p3b

# combine histograms in one plot
p3c <- ggarrange(p3, p3b, ncol = 1, nrow = 2, align = "v", common.legend = TRUE, 
          legend = "right", labels = c("A", "B"))
p3c
# sample size boxplots
p4 <- ggplot(data = appraisal, aes(x = `Sample size (discovery)`, y = Exposure)) +
  geom_boxplot(aes(col = Exposure), show.legend = FALSE) + 
  scale_color_manual(values = met.brewer(name = "Manet", 
                                         n = length(unique(appraisal$Exposure)), 
                                         type = "continuous")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1), 
        axis.title.x = element_blank()) +
  xlab("Study sample size") +
  coord_flip()
p4

# save png plot
scaling_factor = 3
agg_png(filename = file.path(Figure_path, "Critical_appraisal_plot.png"), width = 700*scaling_factor, 
        height = 502*scaling_factor, scaling = scaling_factor)
p1
dev.off()

agg_png(filename = file.path(Figure_path, "Critical_appraisal_plot_table.png"), width = 800*scaling_factor, 
        height = 800*scaling_factor, scaling = scaling_factor)
p1b
dev.off()

agg_png(filename = file.path(Figure_path, "Sample_size_plot.png"), width = 501*scaling_factor,
        height = 501*scaling_factor, scaling = scaling_factor)
p2
dev.off()

agg_png(filename = file.path(Figure_path, "Sample_size_Exposure_plot.png"), width = 501*scaling_factor,
        height = 501*scaling_factor, scaling = scaling_factor)
p3
dev.off()

agg_png(filename = file.path(Figure_path, "Sample_size_Exposure_combined.png"), width = 800*scaling_factor,
        height = 574*scaling_factor, scaling = scaling_factor)
p3c
dev.off()

agg_png(filename = file.path(Figure_path, "Sample_size_Exposure_boxplot.png"), width = 501*scaling_factor,
        height = 501*scaling_factor, scaling = scaling_factor)
p4
dev.off()

# category counts
write_delim(category_counts, 
            file = file.path(Data_path, "Critical_appraisal_counts.tsv"), delim = "\t", 
            col_names = TRUE)
# EWAS_list + critical_appraisal_sheet joined together 
write_delim(appraisal, 
            file = file.path(Data_path, "Critical_appraisal_EWAS_list.tsv"), delim = "\t", 
            col_names = TRUE)

# get sessioninfo
date()
sessionInfo()
