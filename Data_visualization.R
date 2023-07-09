library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(MetBrewer)
library(plotly)
library(purrr)
library(ragg)
library(reticulate)

## get date
date()

## set working directory
setwd(file.path("C:/Users/User/OneDrive/Documenten/Documenten-UU/Master/MCLS/", 
                "2022-2023/Writing_assignment/Data_visualization"))
## load data
# give data and figure paths
Data_path <- file.path("C:/Users/User/OneDrive/Documenten/Documenten-UU/Master/", 
                  "MCLS/2022-2023/Writing_assignment/Data_visualization/Data")
Figure_path <- file.path("C:/Users/User/OneDrive/Documenten/Documenten-UU/Master", 
                         "MCLS/2022-2023/Writing_assignment/Data_visualization/Figures")

# load EWAS_list
EWAS_list <- read_delim(file.path(Data_path, "EWAS_list.tsv"), delim = "\t", 
                        escape_double = FALSE, trim_ws = TRUE)
# shorten long names
EWAS_list %<>%
  as_tibble() %>%
  filter(!is.na(Exposure)) %>%
  mutate(Exposure = ifelse(Exposure == "Chlorinated persistent pollutants (POPs)", "POPs", 
                           Exposure)) %>%
  mutate(Exposure = ifelse(Exposure == "phylloquinone (vitamin K1)", "phylloquinone", 
                           Exposure))
# load Prediction_list
prediction_df <- read_delim(file.path(Data_path, "Prediction_list.tsv"), delim = "\t", 
                            escape_double = FALSE, trim_ws = TRUE)

## create barplot of the number of EWAS per category

# create counts for reordering barplot
EWAS_counts <- count(EWAS_list, Exposure, sort = TRUE)
# get correct factor level order
cat_levels <- levels(reorder(EWAS_counts$Exposure, EWAS_counts$n, decreasing = TRUE))
# create plot
p1 <- ggplot(data = EWAS_list, aes(x = factor(Exposure, levels = cat_levels))) + 
  scale_fill_manual(values = met.brewer("Manet", n = 24, type = "continuous")) +
  geom_bar(aes(fill = factor(Exposure, levels = cat_levels)), show.legend = FALSE) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1), 
        axis.title.x = element_blank()) +
  ylab("Number of EWAS")
p1
## Make a sunburst diagram of the performed EWAS

# get exposure category counts
Exp_values <- sapply(unique(EWAS_list$Exposure), function(i) {
  len = length(EWAS_list$Exposure[EWAS_list$Exposure == i])
  return(len)
})

# get EWAS category counts
EWAS_values <- sapply(unique(EWAS_list$EWAS), function(i) {
  len = length(EWAS_list$EWAS[EWAS_list$EWAS == i])
  return(len)
})

# get parents vector 
parents <- rep("", length(unique(EWAS_list$Exposure)))
parents_v <- c(parents, EWAS_list$Exposure[!duplicated(EWAS_list$EWAS)])

# get label vector
labels_v <- c(unique(EWAS_list$Exposure), unique(EWAS_list$EWAS))

# get values vector
values_v <- unname(c(Exp_values, EWAS_values))

# generate sunburst plot
p2 <- plot_ly()

p2 <- p2 %>% add_trace(
  type = "sunburst",
  parents = parents_v,
  labels = labels_v,
  values = values_v,
  textinfo = "label+value",
  branchvalues = "total"
)
  
p2

## Visualization of significant CpGs and DMRs

# summary table of the DMPs
summary_DMP <- EWAS_list %>%
  group_by(Exposure) %>%
  summarise(EWAS = n(),
            Mean = mean(`Reported CpG-sites (different significance threshholds)`),
            SD = sd(`Reported CpG-sites (different significance threshholds)`),
            Median = median(`Reported CpG-sites (different significance threshholds)`),
            IQR = IQR(`Reported CpG-sites (different significance threshholds)`),
            Max = max(`Reported CpG-sites (different significance threshholds)`),
            Min = min(`Reported CpG-sites (different significance threshholds)`))
# summary table of the DMRs
summary_DMR <- EWAS_list %>%
  subset(!(is.na(DMRs))) %>%
  group_by(Exposure) %>%
  summarise(EWAS = n(),
            Mean = mean(DMRs),
            SD = sd(DMRs),
            Median = median(DMRs),
            IQR = IQR(DMRs),
            Max = max(DMRs),
            Min = min(DMRs))

# Scatterplot of the number of significant CpG sites (log10-transformed)
p3 <- ggplot(data = EWAS_list, aes(x = Exposure, 
                             y = log10(`Reported CpG-sites (different significance threshholds)`))) +
  geom_jitter(aes(col = Exposure, shape = Exposure), show.legend = TRUE) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank()) +
  scale_shape_manual(values = rep(c(15, 16, 17, 18, 19), 4)) +
  scale_color_manual(values = met.brewer("Manet", n = length(unique(EWAS_list$Exposure)), 
                                        type = "continuous")) +
  ylab("Log10(Number of significant DMPs)")
p3

# Scatterplot of the number of significant DMRs
p4 <- ggplot(data = subset(EWAS_list, !(is.na(DMRs))), aes(x = Exposure, 
                             y = log10(DMRs))) +
  geom_jitter(aes(col = Exposure, shape = Exposure), show.legend = TRUE) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank()) +
  scale_shape_manual(values = rep(c(15, 16, 17, 18, 19), 4)) +
  scale_color_manual(values = met.brewer("Manet", n = length(unique(EWAS_list$Exposure)), 
                                         type = "continuous")) +
  ylab("Log10(Number of significant DMRs)")
p4
# Prepare tibble for Barplot
EWAS_list <- EWAS_list %>%
  rename(Cpg_sites = `Reported CpG-sites (different significance threshholds)`) %>%
  # Get CpG sites and DMRs in one column
  pivot_longer(cols = Cpg_sites:DMRs, names_to = "Type", values_to = "Sign position") %>%
  # If not 0, log10 transform the data
  mutate(`Sign position` = ifelse(is.na(`Sign position`), 0, `Sign position`)) %>%
  mutate(`Sign position` = ifelse(`Sign position` != 0, log10(`Sign position`), `Sign position`))

# Barplot of bot the CpG sites and DMRs
p5 <- ggplot(data = EWAS_list, aes(x = Exposure, y = `Sign position`)) +
  geom_col(aes(fill = Type), position = position_dodge()) +
  ylab("Number of Significant positions (log10)") +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank()) +
  scale_fill_manual(values = c(met.brewer("Manet")[3], met.brewer("Manet")[7]),
                    labels = c("DMPs", "DMRs"),
                    name = "Site type")
p5


## create plots for predictions

p8 <- ggplot(prediction_df, aes(x = `AUC (replication)`, y = Exposure)) +
  geom_point(aes(col = Exposure, size = `Prediction CpGs`)) +
  scale_color_manual(values = met.brewer("Manet", n = 4, type = "continuous")) +
  scale_size_continuous(name = "Number of CpG-sites\nused in prediction")+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12)) +
  ylab("Environmental exposure") +
  xlab("AUC (validation set)")
p8

## create ggpubr plots

# for DMPs
table_1 <- ggtexttable(summary_DMP, rows = NULL, 
                       theme = ttheme("lBlack", 
                                      colnames.style = colnames_style(size = 8)))
table_1 <- table_1 %>% table_cell_font(row = 2:tab_nrow(table_1), column = 1:tab_ncol(table_1),
                  size = 8)

p6 <- ggarrange(p3 + theme(legend.position = c(1.13, 0.38),
                           plot.margin = unit(c(0,8,0,1), 'lines')), 
          table_1 + theme(plot.margin = unit(c(13.5,0,0,0), 'lines')), 
          nrow = 2, ncol = 1, labels = c("A", "B"),
          heights = c(2,1),
          vjust = c(1.5, 0.5)) +
  theme(plot.margin = unit(c(0,0,14,0), 'lines'))

# for DMRs

table_2 <- ggtexttable(summary_DMR, rows = NULL, 
                       theme = ttheme("lBlack", 
                                      colnames.style = colnames_style(size = 8)))
table_2 <- table_2 %>% table_cell_font(row = 2:tab_nrow(table_2), column = 1:tab_ncol(table_2),
                                       size = 8)

p7 <- ggarrange(p4 + theme(legend.position = c(1.11, 0.5),
                           plot.margin = unit(c(0,7,0,1), 'lines')), 
                table_2, 
                nrow = 2, ncol = 1, labels = c("A", "B"),
                heights = c(1,1))

# for predictions and DMP/DMR barplot

p9 <- ggarrange(p5,
                p8 + theme(plot.margin = unit(c(3,0,0,0), 'lines')),
          nrow = 2, ncol = 1, 
          labels = c("A", "B"),
          heights = c(2,1)) +
  theme(plot.margin = unit(c(0,0,1.5,0), 'lines'))
p9

# save png plots
scaling_factor = 3
agg_png(filename = file.path(Figure_path, "EWAS_number_barplot.png"), 
        width = 800*scaling_factor, 
        height = 574*scaling_factor, scaling = scaling_factor)
p1
dev.off()

agg_png(filename = file.path(Figure_path, "EWAS_CpGnumber_scatter.png"), 
        width = 850*scaling_factor, 
        height = 574*scaling_factor, scaling = scaling_factor)
p3
dev.off()

agg_png(filename = file.path(Figure_path, "EWAS_DMRnumber_scatter.png"), 
        width = 501*scaling_factor,
        height = 501*scaling_factor, scaling = scaling_factor)
p4
dev.off()

agg_png(filename = file.path(Figure_path, "EWAS_CpGandDMRnumber_barplot.png"), 
        width = 501*scaling_factor,
        height = 501*scaling_factor, scaling = scaling_factor)
p5
dev.off()

agg_png(filename = file.path(Figure_path, "EWAS_DMPs_table.png"), 
        width = 600*scaling_factor, 
        height = 800*scaling_factor, scaling = scaling_factor)
p6
dev.off()

agg_png(filename = file.path(Figure_path, "EWAS_DMRs_table.png"), 
        width = 600*scaling_factor, 
        height = 800*scaling_factor, scaling = scaling_factor)
p7
dev.off()

agg_png(filename = file.path(Figure_path, "Prediction_AUC_plot.png"), 
        width = 623*scaling_factor, 
        height = 319*scaling_factor, scaling = scaling_factor)
p8
dev.off()

agg_png(filename = file.path(Figure_path, "Prediction_AUC_CpG_barplot_combined.png"), 
        width = 623*scaling_factor, 
        height = 594*scaling_factor, scaling = scaling_factor)
p9
dev.off()

# save DMP and DMR summary tables
write_delim(summary_DMP, 
            file = file.path(Data_path, "Summary_DMP.tsv"), delim = "\t", 
            col_names = TRUE)
write_delim(summary_DMR, 
            file = file.path(Data_path, "Summary_DMR.tsv"), delim = "\t", 
            col_names = TRUE)

## save sunburst diagram
# load the kaleido conda environment
use_condaenv(condaenv = "kaleido", conda = "C:/Users/User/anaconda3/_conda.exe")

# import sys and plotly modules
reticulate::py_run_string("import sys")
reticulate::py_run_string("import plotly")
# save plotly image
save_image(p2, file.path(Figure_path, "EWAS_number_sunburst.png"), width = 501, 
           height = 501, scale = 6)

## get sessioninfo
date()
sessionInfo()