#Stable isotope mixing model (SIMMr) analysis #

#Load packages
library(R2jags)
library(simmr)
library(dplyr)

###SIMMr figure
library(simmr)
library(ggplot2)
library(reshape2)
library(patchwork)

#Load the CSV files containing target data, source data, and trophic enrichment factors (TEFs).
#Targets: csv structure -> 	d13C d15N ID nail portion sex
targets <- read.csv("targets ID.csv", header=T, sep = ";", dec = ",", stringsAsFactors = T)

#Sources: csv structure -> source_ID Meand13C Meand15N Sd13C Sd15N
sources <- read.csv("sources.csv", header=T, sep = ";", dec = ",", stringsAsFactors = T)
sources$sources<-as.character(sources$sources)

#TEFs: csv structure -> 	source_TEF mean13C mean15N sd13C sd15N
#TEFs values here used: mean13C = 1,18; mean15N =	2,76; sd13C =	0,39; sd15N =	0,55
TEFs <- read.csv("TEFs.csv", header=T, sep = ";", dec = ",", stringsAsFactors = T)

#Load data into the SIMMr model object, specifying ID as the grouping factor.
simmr_in = simmr_load(mixtures = as.matrix(targets[, 1:2]),
                      source_names = sources$sources,
                      source_means = as.matrix(sources[,2:3]),
                      source_sds = as.matrix(sources[,4:5]),
                      correction_means = as.matrix(TEFs[,2:3]),
                      correction_sds = as.matrix(TEFs[,4:5]),
                      group=as.factor(targets$ID)) #specify ID as the grouping factor so the model is fitted per individual

#Plot the SIMMr input data
plot(simmr_in, 
     xlab = expression(paste(delta^13, "C (\u2030)",
                             sep = "")), 
     ylab = expression(paste(delta^15, "N (\u2030)",
                             sep = "")))

#Run the simmr model
simmr_out = simmr_mcmc(simmr_in)

#Assess model fit
summary(simmr_out, type = "diagnostics")  
post_pred = posterior_predictive(simmr_out)
print(post_pred)

#Combine E. antarctica and T. newnesi into a single source category named "Fish".
final_model<-combine_sources(
  simmr_out,
  to_combine = c("E.antarctica", "T.newnesi"),
  new_source_name = "Fish"
)

#Save model results
save(final_model, file = "simmr_final_output.rda")


# Visualization of model results:
# Note that in the original dataset IDs range from 15–30 and 36,
# but SIMMr indexes groups from 1–17.

#Mean and quartiles
summary(final_model, type = 'statistics', group = 1)
summary(final_model, type = 'quantiles', group = 1)

#Compare posterior prey proportions among individuals
compare_groups(final_model, source_name = "E.superba", groups = 1:17)
compare_groups(final_model, source_name = "E.superba juv", groups = 1:17)
compare_groups(final_model, source_name = "Fish", groups = 1:17)
compare_groups(final_model, source_name = "Jellyfish", groups = 1:17)


#Extract and combine posterior samples for all individuals
ids <- c(15:30, 36)

results <- do.call(rbind, lapply(ids, function(i) {
  mat <- final_model[["output"]][[as.character(i)]][["BUGSoutput"]][["sims.matrix"]]
  cbind(mat, ID = i)
}))

write.csv(results, "simmr_final_output.csv", row.names = FALSE)


#Individual- and population-level metrics derived from simmr_final_output.csv
#Read posterior samples
DF <- read.csv("simmr_final_output.csv", header=T, sep = ",", stringsAsFactors = T)
DF$ID <- as.factor(DF$ID)

#Convert the ID variable into ordinal numeric values.
DF$ID_ord <- as.integer(factor(DF$ID, levels = sort(unique(DF$ID))))

#Individual-level posterior mean and 95% credible interval
diet_summary_ID_95IC <- DF %>%
  group_by(ID_ord) %>%
  summarise(
    mean_E_superba = round(mean(E.superba, na.rm = TRUE), 2),
    l95_E_superba = round(quantile(E.superba, 0.025, na.rm = TRUE), 2),
    u95_E_superba = round(quantile(E.superba, 0.975, na.rm = TRUE), 2),

    mean_E_superba_juv = round(mean(E.superba.juv, na.rm = TRUE), 2),
    l95_E_superba_juv = round(quantile(E.superba.juv, 0.025, na.rm = TRUE), 2),
    u95_E_superba_juv = round(quantile(E.superba.juv, 0.975, na.rm = TRUE), 2),
    
    mean_Fish = round(mean(Fish, na.rm = TRUE), 2),
    l95_Fish = round(quantile(Fish, 0.025, na.rm = TRUE), 2),
    u95_Fish = round(quantile(Fish, 0.975, na.rm = TRUE), 2),
    
    mean_Jellyfish = round(mean(Jellyfish, na.rm = TRUE), 2),
    l95_Jellyfish = round(quantile(Jellyfish, 0.025, na.rm = TRUE), 2),
    u95_Jellyfish = round(quantile(Jellyfish, 0.975, na.rm = TRUE), 2),
    
    .groups = "drop"
  )

#Population-level metrics
#Mean
pop_mean <- DF %>%
  summarise(
    mean_E_superba = round(mean(E.superba, na.rm = TRUE), 2),
    mean_E_superba_juv = round(mean(E.superba.juv, na.rm = TRUE), 2),
    mean_Fish = round(mean(Fish, na.rm = TRUE), 2),
    mean_Jellyfish = round(mean(Jellyfish, na.rm = TRUE), 2)
  )

#Range (min–max) of individual posterior means
pop_range <- diet_summary_ID_95IC %>%
  summarise(
    min_E_superba = round(min(mean_E_superba), 2),
    max_E_superba = round(max(mean_E_superba), 2),
    
    min_E_superba_juv = round(min(mean_E_superba_juv), 2),
    max_E_superba_juv = round(max(mean_E_superba_juv), 2),
    
    min_Fish = round(min(mean_Fish), 2),
    max_Fish = round(max(mean_Fish), 2),
    
    min_Jellyfish = round(min(mean_Jellyfish), 2),
    max_Jellyfish = round(max(mean_Jellyfish), 2)
  )

#Export data
write.csv(diet_summary_ID_95IC, file = "output_summary_ID_95IC.csv", row.names = FALSE)

write.csv(pop_mean, file = "output_summary_pop_mean.csv", row.names = FALSE)

write.csv(pop_range, file = "output_summary_pop_range.csv", row.names = FALSE)


# VISUALIZATION of individual boxplots from SIMMr mixing models #

#Load simmr output
df <- read.csv("simmr_final_output.csv", header=T, sep = ",", stringsAsFactors = T)

#Create a new dataframe for the figure
simmr_figure<-as.data.frame(df)

#Convert the ID variable into ordinal numeric values.
simmr_figure$ID_ord <- as.integer(factor(simmr_figure$ID, levels = sort(unique(simmr_figure$ID))))

#Reshape the data from wide to long format for boxplot visualization.
simmr_boxplot <- melt(
  simmr_figure,
  id.vars = "ID_ord",
  measure.vars = c("E.superba", "E.superba.juv", "Fish")
)
#Jellyfish will not be included in the final figure.

summary(simmr_boxplot)

#Create a subset of simmr_boxplot, since the dghlm model excludes ID 13 (sex = NA).
simmr_boxplot2<-subset(simmr_boxplot, ID_ord!="13")
summary(simmr_boxplot2)

#Check it turned out ok.
levels(factor(simmr_boxplot2$ID_ord)) 

#Reorder individuals to match the ordering in the dhglm model output.
manual_order <- c(9,10,7,8,15,14,1,2,4,12,11,6,3,17,16,5)

simmr_boxplot2$ID <- factor(simmr_boxplot2$ID_ord, levels = manual_order)

#Check it turned out ok.
levels(factor(simmr_boxplot2$ID_ord, levels = manual_order))

#Create a new column Prey based on the variable column.
simmr_boxplot2$Prey<-simmr_boxplot2$variable

#Figure
#Note that when using coord_flip(), all theme elements defined for the x-axis apply to the y-axis, and vice versa.
simmr_plot=ggplot(simmr_boxplot2, 
                  aes(
                    x = factor(ID_ord, levels = manual_order),
                    y = value,
                    fill = Prey,
                    color = Prey
                  ))+
  geom_boxplot(
    outlier.shape = NA, show.legend = FALSE, color = "gray12",                 
    size = 0.4,
    width = 0.7,
    position = position_dodge(width = 0.7) 
  ) +
  theme_classic() +
  labs(x = "Individual (ID)", y = "Proportion in diet") +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, size = 26, colour = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size=28),
    axis.title.y = element_blank(),
    legend.title = element_text(size=20),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    axis.ticks.y = element_blank()
  )+
  scale_fill_manual(
    values = c(
      "E.superba"     = "grey30",
      "E.superba.juv" = "grey60",
      "Fish"         = "white"
    ),
    labels = c(
      expression("Adult "*italic("E. superba")),
      expression("Juvenile "*italic("E. superba")),
      "Fish"
    )
  ) +  
  scale_color_manual(
    values = c(
      "E.superba"     = "black",
      "E.superba.juv" = "black",
      "Fish"         = "black"
    ),
    labels = c(
      expression("Adult "*italic("E. superba")),
      expression("Juvenile "*italic("E. superba")),
      "Fish"
    )
  ) +
  coord_flip()

#Add margins for visualization purposes
simmr_plot=simmr_plot+theme(
  plot.margin = margin(5, 25, 5, 5))

simmr_plot

