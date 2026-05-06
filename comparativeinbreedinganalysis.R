##ANNA CARVAJAL##


##set working direcotry and load packages
setwd("~/Downloads/EEB331/PAPERRR/PLOTS")
getwd()

library("pedtools")
library("optiSel") 
library("dplyr")
library("ggplot2")

###pedigreeinfo###

# read files
ped_data <- read.csv("marmot_Parentage_Informed_Pedigree.csv")
LH <- read.csv("RADseq_marmot_metadata.csv")

# clean pedigree
ped_clean <- ped_data %>%
  select(id, dam, sire)

# replace blank parent IDs with NA
ped_clean$dam[ped_clean$dam == ""] <- NA
ped_clean$sire[ped_clean$sire == ""] <- NA

# clean sex data
LH <- LH %>%
  rename(id = file_name) %>%
  mutate(
    sex = case_when(
      sex %in% c("M", "male", "Male") ~ 1,
      sex %in% c("F", "female", "Female") ~ 2,
      TRUE ~ NA_real_
    )
  )

# add sex data to pedigree file
ped_final <- ped_clean %>%
  left_join(LH %>% select(id, sex), by = "id")

#for missing sex 
ped_final$sex[ped_final$id %in% ped_final$sire] <- 1
ped_final$sex[ped_final$id %in% ped_final$dam] <- 2

#for missing parent
one_parent <- xor(is.na(ped_final$sire), is.na(ped_final$dam))
ped_final$sire[one_parent] <- NA
ped_final$dam[one_parent] <- NA

#format pedigree
ped_opti <- ped_final
colnames(ped_opti) <- c("Indiv", "Dam", "Sire", "Sex")

#prepare ped and calculate F 
ped_opti <- prePed(ped_opti)
F <- pedInbreeding(ped_opti)

head(F)

###ROHinfooo###

# genome size for FROH
genome_size <- 2319465413  # in bp

# read and clean
lines <- readLines("Marmots_maf3_geno90_mind20_ROH.txt")
lines <- lines[!grepl("^#", lines)]
rg_lines <- lines[grepl("^RG", lines)]

# turn into table
roh_segments <- read.table(
  text = rg_lines,
  header = FALSE,
  stringsAsFactors = FALSE
)

# add column names
colnames(roh_segments) <- c(
  "type",
  "sample",
  "chr",
  "start",
  "end",
  "length_bp",
  "n_markers",
  "quality"
)

# summarize ROH by individual
roh_summary <- roh_segments %>%
  mutate(length_bp = as.numeric(end) - as.numeric(start)) %>%
  group_by(sample) %>%
  summarise(
    NSEG = n(),
    total_bp = sum(length_bp),
    mean_bp = mean(length_bp),
    mean_length = total_bp / NSEG
  ) %>%
  mutate(FROH = total_bp / genome_size)

head(roh_summary)

##theme## 
# set up plot aesthetics we can use for all figures
eeb_theme <- function(base_size = 9, base_family = "Arial") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      # Set all text to Arial
      text = element_text(family = "Arial"),
      # Remove grid lines and background
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      # Plot title
      plot.title = element_text(size = base_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = base_size, hjust = 0.5, margin = margin(b = 10)),
      # Axes
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      axis.title = element_text(size = 9),
      # Legend
      legend.background = element_rect(fill = "transparent", color = NA),
      legend.key = element_blank(),
      legend.text = element_text(size = 9),
      legend.title = element_text(size = 9),
    ) 
}

##MERGINGF&FROH###
head(F$Indiv)
head(roh_summary$sample)

#clean sample names 
roh_summary$sample <- gsub("_sorted_UID_.*$", "", roh_summary$sample)
roh_summary$sample <- gsub("_sorted_RMBL_.*$", "", roh_summary$sample)

#merge
plot_data <- merge(F, roh_summary, by.x = "Indiv", by.y = "sample")
head(plot_data)
dim(plot_data)

#remove indv missing ped or gen inbreeding
plot_data <- plot_data %>%
  filter(!is.na(Inbr), !is.na(FROH))
head(plot_data)
dim(plot_data)

### PLOT: PEDIGREE INBREEDING VS GENOMIC INBREEDING ###
f_vs_froh <- ggplot(plot_data, aes(x = Inbr, y = FROH)) +
  geom_jitter(color = "black", size = 2.2, alpha = 0.5,
              width = 0.003, height = 0) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80") +
  labs(
    title = "Pedigree vs Genomic Inbreeding",
    x = "Pedigree Inbreeding Coefficient (F)",
    y = "Genomic Inbreeding Coefficient (FROH)"
  ) +
  eeb_theme()

f_vs_froh

#saveplot
ggsave("f_vs_froh.png", plot = f_vs_froh,
       width = 7, height = 4, units = "in", dpi = 500, bg = "white")

##stats for plot one##
cor.test(plot_data$Inbr, plot_data$FROH)
model <- lm(FROH ~ Inbr, data = plot_data)
summary(model)

##FIGURETWOOOO, Distributions###

#packages
library(cowplot)
library(ggplot2)
library(dplyr)

#use merged data set
fig2_data <- plot_data %>%
  filter(!is.na(Inbr), !is.na(FROH))
head(fig2_data)

#panel A ped inbreeding
ped_f_hist <- ggplot(fig2_data, aes(x = Inbr)) +
  geom_histogram(binwidth = 0.02, fill = "gray70", color = "black") +
  labs(
    title = "Pedigree Inbreeding",
    x = "Pedigree Inbreeding Coefficient (F)",
    y = "Number of Individuals"
  ) +
  eeb_theme()

ped_f_hist

#panel B genomic inbreeding
froh_hist <- ggplot(fig2_data, aes(x = FROH)) +
  geom_histogram(binwidth = 0.02, fill = "gray70", color = "black") +
  labs(
    title = "Genomic Inbreeding",
    x = "Genomic Inbreeding Coefficient (FROH)",
    y = "Number of Individuals"
  ) +
  eeb_theme()

froh_hist

#combine panels
fig2 <- plot_grid(ped_f_hist, froh_hist, ncol = 2, labels = c("A", "B"))
fig2

#save fig
ggsave("figure2_inbreeding_distributions.png", plot = fig2,
       width = 8, height = 4, units = "in", dpi = 500, bg = "white")

#stats
#ped
summary(fig2_data$Inbr)
sd(fig2_data$Inbr)
#genomic
summary(fig2_data$FROH)
sd(fig2_data$FROH)
#ind F=0
mean(fig2_data$Inbr == 0)

##PLOTTHHREEEE, FROH vs NSEG ###

#Prepare ROH Data
roh_plot <- roh_summary %>%
  filter(!is.na(FROH), !is.na(NSEG)) #keeping indv w FROH and NSEG

head(roh_plot)
dim(roh_plot)

#plot, shows how genomic inbreeding relates to number ROH segs per indv
froh_nseg <- ggplot(roh_plot, aes(x = NSEG, y = FROH)) +
  geom_point(color = "black", size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", fill = "gray80") +
  labs(
    title = "Genomic Inbreeding vs ROH Structure",
    x = "Number of ROH Segments",
    y = "Genomic Inbreeding Coefficient (FROH)"
  ) +
  eeb_theme()

froh_nseg

#save
ggsave("figure3_froh_vs_nseg.png", plot = froh_nseg,
       width = 6.5, height = 4.5, units = "in", dpi = 500, bg = "white")

#stats
cor.test(roh_plot$NSEG, roh_plot$FROH)

model_nseg <- lm(FROH ~ NSEG, data = roh_plot)
summary(model_nseg)










