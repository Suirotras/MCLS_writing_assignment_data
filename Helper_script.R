library(magrittr)

selected_studies <- appraisal %>%
  filter(`classifier prediction tested` %in% c("+", "+/-"))

selected_studies <- appraisal %>%
  filter(Exposure %in% c("Pesticides"))

#selected_studies <- appraisal %>%
#  filter(key %in% c("rayyan-1000417268", "rayyan-1000417286", "rayyan-1000417344",
#                    "rayyan-1000428229", "rayyan-1000428324", "rayyan-1000428518",
#                    "rayyan-1000428603"))

keys <- selected_studies %>%
  select(key)

selected_EWAS <- EWAS_list %>%
  filter(key %in% keys$key)

study_counts <- count(selected_studies, Exposure, sort = TRUE)
# get correct factor level order
cat_levels <- levels(reorder(study_counts$Exposure, study_counts$n, decreasing = TRUE))
# create plot
ggplot(data = selected_studies, aes(x = factor(Exposure, levels = cat_levels))) + 
  scale_fill_manual(values = met.brewer("Manet", n = length(cat_levels), type = "continuous")) +
  geom_bar(aes(fill = factor(Exposure, levels = cat_levels)), show.legend = FALSE) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1), 
        axis.title.x = element_blank()) +
  ylab("Number of EWAS")


# select the EWAS exposure of interest
filter_EWAS <- FALSE
EWAS_exp <- ""

# for DMPs
selected_EWAS %>%
  {if(filter_EWAS) filter(., EWAS == EWAS_exp) else .} %>%
  summarise(mean = mean(`Reported CpG-sites (different significance threshholds)`),
            sd = sd(`Reported CpG-sites (different significance threshholds)`),
            median = median(`Reported CpG-sites (different significance threshholds)`),
            IQR = IQR(`Reported CpG-sites (different significance threshholds)`)) %>%
  View()
# for DMRs
selected_EWAS %>%
  {if(filter_EWAS) filter(., EWAS == EWAS_exp) else .} %>%
  filter(!(is.na(DMRs))) %>%
  summarise(mean = mean(DMRs),
            sd = sd(DMRs),
            median = median(DMRs),
            IQR = IQR(DMRs)) %>%
  View()

appraisal %>%
  select(`Sample size (discovery)`) %>%
  flatten_dbl() %>%
  IQR()
  