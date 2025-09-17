library(corrplot)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(Hmisc)
library(openxlsx)
library(reshape)
library(tidypaleo)
library(vegan)

theme_set(theme_paleo(8))

init_project <- function() {
    rm(list = ls())
    dir.create("Output", recursive = TRUE, showWarnings = FALSE)
}


read_raw_data <- function(folder_path, section) {
  raw_files <- list.files(path = folder_path)
  data <- data.frame(matrix(ncol = 1, nrow = 0, dimnames = list(NULL, "compound")))
  for (filename in raw_files) {
    fullfilename <- paste(folder_path, filename, sep = "/")
    sample <- read.table(fullfilename, sep = "\t", header = TRUE)
    data <- merge(data, sample, by = "compound", all = TRUE)
  }
  
  rownames(data) <- data$compound
  data <- data[, !names(data) %in% "compound"]
  data <- as.data.frame(t(data))
  data$sample_ID <- rownames(data)
  metadata <- read.xlsx("Raw_data/Metadata.xlsx", sheet = section) # import metadata
  log_height <- data.frame(cbind(sample_ID = metadata$sample_ID, log_height = metadata$log_height))
  data <- merge(data, log_height, by = "sample_ID")
  data <- data[, !names(data) %in% "sample_ID"]
  rownames(data) <- data$log_height
  data <- data[, !names(data) %in% "log_height"]
  data <- data[, c(ncol(data), 1:ncol(data)-1)]
  
  return(data)
}


CPI_calc <- function(a_data) {
  # calculates the carbon preference index from the integrated peaks of the compounds
  
  CPI <- ((a_data$`nC21` + a_data$`nC23` + a_data$`nC25` + a_data$`nC27` + a_data$`nC29` + a_data$`nC31` + a_data$`nC33`) +
    (a_data$`nC23` + a_data$`nC25` + a_data$`nC27` + a_data$`nC29` + a_data$`nC31` + a_data$`nC33` + a_data$`nC35`)) /
    (2*(a_data$`nC22` + a_data$`nC24` +a_data$`nC26` + a_data$`nC28` + a_data$`nC30` + a_data$`nC32` + a_data$`nC34`))
    
  CPI <- as.data.frame(cbind(a_data$log_height, CPI))
  names(CPI) <- c("log_height", "CPI")
  
  return(CPI)
}

TAR_calc <- function(a_data) {
  # calculates the terrigenous-aquatic ratio (TAR) from the integrated peaks of the compounds
  
  TAR <- (a_data$`nC27` + a_data$`nC29` + a_data$`nC31`)/(a_data$`nC15` + a_data$`nC17` + a_data$`nC19`)
  
  TAR <- as.data.frame(cbind(a_data$log_height, TAR))
  names(TAR) <- c("log_height", "TAR")
  
  return(TAR)
}

MPI_calc <- function(MS_data) {
  # calculates the methylphenanthrene index (MPI) from the integrated peaks of the compounds
  # calculation after Cassani et al. (1988):
  # MPI = 1.89(2-MP + 3-MP)/[P + 1.26(1-MP + 9-MP)].
  
  MPI <- (1.89*(MS_data$"2MP" + MS_data$"3MP")/(MS_data$phenanthrene + 1.26*(MS_data$"1MP" + MS_data$"9MP")))
  
  MPI <- as.data.frame(cbind(MS_data$log_height, MPI))
  names(MPI) <- c("log_height", "MPI")
  
  return(MPI)
}

show_metadata <- function(section) {
  # show a vector with column names of metadata file
  
  all_metadata <- read.xlsx("Raw_data/Metadata.xlsx", sheet = section)
  metadata <- colnames(all_metadata)
  
  return(metadata)
}

get_metadata <- function(metadata_variable, section) {
  # gets a variable from the metadata file
  # and creates a data frame with log height and the respective variable
  
  all_metadata <- read.xlsx("Raw_data/Metadata.xlsx", sheet = section)
  metadata <- select(all_metadata, c("log_height", all_of(metadata_variable)))
  
  return(metadata)
}

add_metadata <- function(data, metadata_variable, section) {
  # creates a data frame with the input data and merges metadata with it
  # based on log_height
  
  all_metadata <- read.xlsx("Raw_data/Metadata.xlsx", sheet = section)
  metadata <- select(all_metadata, c("log_height", all_of(metadata_variable)))
  
  merged <- merge(data, metadata, by = "log_height")
  
  return(merged)
}

long_df <- function(data) {
  # creates a data frame with the input data and merges metadata with it
  # based on log_height
  
  wide_to_long <- melt(data,  id.vars = 'log_height', variable.name = 'compound') # long format
  names(wide_to_long) <- c("log_height", "compound", "mass")
  
  return(wide_to_long)
}

# DEFINE A FUNCTION FOR CUSTOMIZING PLOTS:
#   variables:
#   scale_x, scale_y: expand scale limits to 0
#   y_name: label y axis
#   extinction horizon: horizontal line at extinction horizon
plot_common_parameters <- function(scale_x = c(0, 0), scale_y = c(0, 0), x_name, y_name, extinction_horizon = TRUE, axis_ticks = "black", rev_ax = TRUE, ticks_labels = "black", axis_limit = NULL) {
  theme <- theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 14, color = axis_ticks), # labels in middle of tick (vjust)
                 axis.title.x = element_text(vjust = 0.5, size = 14),
                 axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14, colour = ticks_labels),
                 axis.title.y = element_text(size = 14),
                 plot.title = element_text(size=14, face="bold", hjust = 0.5),
                 strip.text.x = element_text(angle = 70), # rotates facet labels
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = "black"),
                 panel.grid = element_blank(),
                 panel.border = element_rect(colour = "grey70", fill=NA),
                 panel.background = element_blank())
  
  extinction_horizon_gg <- geom_vline(xintercept = 0, colour = "brown1", linetype = "dashed", lwd = 0.7)
  
  rev_ax_gg <- coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on") # reverse axes
  
  parts <- list(theme)
  
  limits <- expand_limits(x = 0, y = 0)
  parts <- c(parts, limits)
  
  if(rev_ax == TRUE) {
    parts <- c(parts, rev_ax_gg)
  }

  if(class(scale_x) == "numeric") {
    scale_x <- scale_x_continuous(expand = scale_x, name = x_name, breaks = c(-25:25))
    parts <- c(parts, scale_x)
  }
  
  if(class(scale_y) == "numeric") {
     scale_y <- scale_y_continuous(expand = scale_y, name = y_name, limits = axis_limit)
   } else {
     scale_y <- scale_y_continuous(name = y_name)
   }
   parts <- c(parts, scale_y)
  
  if(extinction_horizon == TRUE) {
    parts <- c(parts, extinction_horizon_gg)
  }
  
  plot_common_parameters <- parts
}

# remove y axis ticks and labels
remove_axis <- function(plot_position) {
  theme <- plot_position + theme(axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank())
  remove_axis <- theme
}

# prepares data frame with area of pristane and phytane for plotting
# Pr/nC17 vs. Ph/nC218 plot
# input is data frame with areas from FID
prep_depositional_environment_plot <- function(a_data, section) {
  all_metadata <- read.xlsx("Raw_data/Metadata.xlsx", sheet = section)
  data <- merge(as.data.frame(cbind(log_height = a_data$log_height,
                                    Pristane = a_data$Pristane,
                                    Phytane = a_data$Phytane,
                                    nC17 = a_data$nC17,
                                    nC18 = a_data$nC18)), as.data.frame(cbind(log_height = all_metadata$log_height,
                                                                              member = all_metadata$member)), by = "log_height")
  data$Pr_nC17 <- data$Pristane / data$nC17 
  data$Ph_nC18 <- data$Phytane / data$nC18 
  data$member <- as.factor(data$member)
  data$section <- rep(section, nrow(data))
  
  return(data)
}

# function for Pr_nC17 vs. Ph_nC18 plot
# input: data frame that contains log_height, Pr_nC17 and Ph_nC18 columns
plot_common_parameters_PrPh <- function(Pr_Ph_data) {
  theme <- theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, size = 14, color = "black"),
                 axis.text.y = element_text(vjust = 0.5, hjust = 0.5, size = 14, colour = "black"),
                 axis.title.y = element_text(size = 14),
                 axis.title.x = element_text(size = 14),
                 axis.line.x = element_line(colour = "black"),
                 axis.line.y = element_line(colour = "black"),
                 
                 panel.grid = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_rect(colour = "black", fill=NA),
                 
                 legend.position = "bottom",
                 legend.direction = "vertical",
                 legend.text = element_text(size = 14),
                 legend.title = element_blank(),
                 legend.background = element_blank(),
                 legend.key = element_blank(),
                 legend.spacing.y = unit(0, "mm"))
  
  parts <- list(theme)
  
  x_scale_c <- scale_x_continuous(name = expression(paste("Ph/", italic("n"), "-C"[18])), # scales (x and) in log10
                                trans = "log10",
                                limits = c(0.1, 10), expand = c(0, 0))
  parts <- c(parts, x_scale_c)
  
  y_scale_c <- scale_y_continuous(name = expression(paste("Pr/", italic("n"), "-C"[17])),
                                  trans = "log10",
                                  limits = c(0.1, 10), expand = c(0, 0))
  parts <- c(parts, y_scale_c)

  line1 <- geom_abline(intercept = 1, linetype = 2)
  line2 <- geom_abline(intercept = 0.5, linetype = 2)
  line3 <- geom_abline(intercept = 0.1, linetype = 2)
  line4 <- geom_abline(intercept = -0.5, linetype = 2)
  parts <- c(parts, line1, line2, line3, line4)
  
  col <- scale_color_manual(values = colors)
  parts <- c(parts, col)
  
  plot_common_parameters_PrPh <- parts
}