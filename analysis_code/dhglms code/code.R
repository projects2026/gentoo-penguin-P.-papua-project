# Load required packages for bayesian analysis, data manipulation, and visualization

library(brms) # Bayesian regression models using Stan
library(purrr) # Functional programming tools
library(bayestestR) # Bayesian model testing and estimation
library(ellipse) # Ellipse calculations for plotting
library(tidyverse) # Data manipulation and visualization (dplyr, ggplot2, tidyr, etc.)
library(patchwork) # Combine multiple ggplot objects
library(ggridges) # Ridge density plots
source("functions.R") # Load custom helper functions (plot_blups, repeatt, summ_bayes, plots_ell)

setwd("C:/Users/Dell/Desktop/code paper pinguis")

# Load individual-level isotope data (d13C and d15N measurements)
DF <- read.csv("papua.csv", header = T, sep = ";", stringsAsFactors = T, dec = ",")

# Define DHGLM (Double Hierarchical Generalized Linear Model) formulas
# Model both the mean (location) and variance (dispersion) of d13C
# form1: Mean structure - d13C varies by sex with random intercepts by individual
# form2: Variance structure - sigma (within-individual variance) also varies by sex and individual
form1<- bf(d13C ~ sexo + (1|a|id)) + lf (sigma ~ sexo + (1|a|id))
form2<- bf(d15N ~ sexo + (1|a|id)) + lf (sigma ~ sexo + (1|a|id))

# Fit the Bayesian DHGLM or load it if already fitted
fit1 <- brm(formula = form1 + form2, data = DF, iter = 8000, chains = 4, thin = 10, cores = 4, rescor=F, 
            control = list(adapt_delta = 0.95))

#### Save the fitted model to disk for future use
saveRDS(fit1, file = "papua_model.rds")
# Load pre-fitted model to save time
model1 <- readRDS("papua_model_sexo.rds")
summary(model1)


# ========== VISUALIZATION: ========== #
# Plot individual-level random effects showing how each individual deviates from population mean #
# Creates scatter plot of individual estimates for mean and variance components #

# Function to integrate fixed effects
# (Intercept + sexo / 2)
# -------------------------------
fn_int <- function(obj) {
  if (length(obj) == 1) return(obj[1])
  obj[1] + obj[2] / 2
}

# -------------------------------
# Function to extract BLUPs and plot
# -------------------------------
plot_blups <- function(model1){
  
  # ---- Extract fixed effects (posterior means)
  fix_coef <- fixef(model1)[, 1]
  # ---- Integrated fixed effects (mu)
  int_d13C <- fn_int(fix_coef[c("d13C_Intercept", "d13C_sexom")])
  int_d15N <- fn_int(fix_coef[c("d15N_Intercept", "d15N_sexom")])
  # ---- Integrated fixed effects (sigma, log scale)
  sig_d13C <- fn_int(fix_coef[c("sigma_d13C_Intercept", "sigma_d13C_sexom")])
  sig_d15N <- fn_int(fix_coef[c("sigma_d15N_Intercept", "sigma_d15N_sexom")])
  # ---- Extract BLUPs (id level)
  blups <- ranef(model1)$id[, 1, ] |>
    data.frame() |>
    rename_with(~ gsub("_Intercept", "", .x))
  # ---- Combine fixed + random effects
  bl_fix <- blups |>
    mutate(
      id = rownames(blups),
      # individual means
      int_d13C = d13C + int_d13C,
      int_d15N = d15N + int_d15N,
      # individual SDs (sigma = exp(eta))
      sig_d13C = exp(sigma_d13C + sig_d13C),
      sig_d15N = exp(sigma_d15N + sig_d15N))
  # ---- Build plot ---- #
  p <- ggplot(bl_fix, aes(x = int_d13C, y = int_d15N)) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(xmin = int_d13C - sig_d13C,
          xmax = int_d13C + sig_d13C),
      color = "darkblue",
      width = 0,
      linewidth = 1.3
    ) +
    geom_errorbar(
      aes(ymin = int_d15N - sig_d15N,
          ymax = int_d15N + sig_d15N),
      color = "coral",
      width = 0) +
    theme_classic() +
    labs(
      x = expression({delta}^13*C~'\u2030'),
      y = expression({delta}^15*N~'\u2030')
    )+
    theme(
      axis.text  = element_text(size = 11),
      axis.title = element_text(size = 11)
    ) +
    coord_cartesian(
      xlim = c(-25.5, -23),
      ylim = c(6, 9.5)
    ) +
    scale_x_continuous(breaks = seq(-25.5, -23, by = 0.5)) +
    scale_y_continuous(breaks = seq(6, 9.5, by = 0.5))
  
  return(p)
}

# -------------------------------
# RUN FUNCTION
# -------------------------------
p <- plot_blups(model1)
print(p)

# save the plot
tiff("figura_blups.tiff", width = 6, height = 6, units = "in", res = 300)
print(p)  
dev.off()

# ========== EXTRACT INDIVIDUAL-LEVEL INTRA-INDIVIDUAL VARIANCE (rIIV) ==========

# Select only columns for sigma (within-individual variance) random effects
# These represent individual-specific variance components for d13C and d15N

comunidad1<-posterior_samples(model1)

df_filtrado <- comunidad1 %>%
  select(starts_with("r_id__sigma_d"))

#In order to interpret rIIV in biological terms we backtransform rIIV by taking
# In DHGLMs, sigma is modeled on log scale, so exp() returns to sd scalecommunity.sp<- exp(df_filtrado)
community.sp<- exp(df_filtrado)
Sp.C<- community.sp[, c(1:16)]  
Sp.N<- community.sp[, c(17:32)]

# Split by isotope: columns 1-16 are Carbon (d13C), columns 17-32 are Nitrogen (d15N)
Sp.C <- community.sp[, c(1:16)] # Individual-specific variance for Carbon
Sp.N <- community.sp[, c(17:32)] # Individual-specific variance for Nitrogen

# ========== RESHAPE DATA FOR VISUALIZATION ==========
# Reshape Carbon data from wide to long format (one row per posterior draw per individual)

# for Carbon
df_Carbon <- Sp.C %>%
  pivot_longer(cols = everything(), 
               names_to = "id", 
               values_to = "Intraindividual Variance") %>% # Posterior samples of within-individual variance
  mutate(id = gsub("r_id__sigma_d13c\\[|,Intercept\\]", "", id))%>%
  arrange(id)  # Sort by individual ID

df_Carbon <- Sp.C %>%
  pivot_longer(cols = everything(), 
               names_to = "id", #Individual ID
               values_to = "Intraindividual Variance") %>%
  mutate(id = sub(".*\\[([0-9]+),.*\\]", "\\1", id)) %>%
  arrange(as.numeric(id)) # Sort by individual ID


df_Carbon <- as.data.frame(df_Carbon)
df_Carbon$id<- as.factor(df_Carbon$id)
str(df_Carbon)

# for Nitrogen
df_Nitrogen <- Sp.N %>%
  pivot_longer(cols = everything(), 
               names_to = "id", #individual ID
               values_to = "Intraindividual Variance") %>% # Posterior samples of within-individual variance
  mutate(id = sub(".*\\[([0-9]+),.*\\]", "\\1", id)) %>%
  arrange(as.numeric(id)) # Sort by individual ID

df_N <- as.data.frame(df_Nitrogen)
df_N$id<- as.factor(df_N$id)

# Combine Carbon and Nitrogen data into single dataframe
alls<- rbind(df_N,df_Carbon)
Isotope <- factor(c(rep("Nitrogen", 25600), rep("Carbon", 25600)))

# Create isotope type identifier (Nitrogen first, then Carbon)
df_largo4<- data.frame(Isotope = Isotope)
DF_plot<-cbind(alls,df_largo4)
str(DF_plot)

# Order individuals by their mean intra-individual variance (low to high)
# This makes the ridge plot easier to interpret
DF_plot$id <- factor(DF_plot$id, levels = DF_plot %>%
                       group_by(id) %>%
                       summarise(Mean_IV = mean(`Intraindividual Variance`)) %>%
                       arrange(Mean_IV) %>%
                       pull(id))

# VISUALIZATION OF INTRAINDIVIDUAL VARIANCE WITH GGPLOT2 #
spp_plot <- ggplot(DF_plot, aes(x = `Intraindividual Variance`, y = id, fill = Isotope)) +
  geom_density_ridges(alpha = 0.5, scale = 1.5) +
  theme_classic() +
  scale_fill_manual(values = c("#E44840", "#4D82BC")) +  # Rojo para C, Azul para N
  scale_x_continuous(limits = c(0, 4)) +
  labs(x = "Intraindividual variance", y = "Individuals (ID)", title = "") +
  theme(text = element_text(family = "Times New Roman")) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
    
  )+
  theme(legend.position = "bottom")
spp_plot

# save the pot
ggsave("pinguis_ind_ level.tiff", spp_plot, units="px", width= 8000, height = 5000, device = tiff, dpi=350)

