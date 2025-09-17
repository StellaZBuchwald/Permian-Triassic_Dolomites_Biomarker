# This script performs all calculations and reproduces all figures from the manuscript:
# Buchwald SZ, Birgel D, Kustatscher E, Prinoth H, Galasso F, Karapunar B, Frank AB,
# Gómez Correa MA, Lahajnar N, Peckmann J, Foster WJ
# "Molecular fossils record shallow marine ecosystem changes prior to and across the
# Permian–Triassic mass extinction in the Dolomites (Italy)"

# This script and the script "PTr_Dolomites_functions.R" need to be in the working
# directory to be properly loaded. Create a sub-folder titled "Raw_data" in the
# working directory, which should contain all raw and metadata.

# When initiating the project, a folder "Output" is created in the working directory,
# that will contain all data produced.

setwd("") # set working directory
source('PTr_Dolomites_functions.R')
init_project()

################################################################################
###                            FID data Siusi                                ###
################################################################################

Siusi_a_FID <- read_raw_data("Raw_data/FID/Siusi", "Siusi") # integrated peak areas
Siusi_a_FID$log_height <- as.numeric(rownames(Siusi_a_FID))

Siusi_CPI <- CPI_calc(Siusi_a_FID) # calculate carbon preference index (CPI) after Marzie et al. (1993)
Siusi_TAR <- TAR_calc(Siusi_a_FID) ## calculate terrigenous-aquatic ratio (TAR)

Siusi_ticks <- c("black", NA, NA, NA,NA)

# plot TOC
show_metadata("Siusi")
Siusi_TOC <- get_metadata("TOC", "Siusi")

Siusi_plot_TOC <- ggplot(Siusi_TOC, aes(log_height, TOC)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Siusi_TOC, aes(log_height, TOC, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "height (m)",
                         y_name = "TOC (wt%)",
                         ticks_labels = Siusi_ticks,
                         axis_limit = c(0, 4.6))

# plot TAR
Siusi_plot_TAR <- ggplot(Siusi_TAR, aes(log_height, TAR)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Siusi_TAR, aes(log_height, TAR, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = "TAR",
                         ticks_labels = Siusi_ticks,
                         axis_limit = c(0, 1.9))

# pristane, phytane, nC17, nC18
Siusi_Pr_Ph_FID <- Siusi_a_FID[c("Pristane", "Phytane", "nC17", "nC18", "log_height")]
Siusi_Pr_Ph_FID$Pr_Ph <- Siusi_Pr_Ph_FID$Pristane / Siusi_Pr_Ph_FID$Phytane

# plot Pr/Ph
Siusi_plot_Pr_Ph_FID <- ggplot(Siusi_Pr_Ph_FID, aes(log_height, Pr_Ph)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Siusi_Pr_Ph_FID, aes(log_height, Pr_Ph, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  geom_hline(yintercept = 1.0, colour = "grey50", linetype = "dashed", lwd = 0.6) + # vertical line at Pr/Ph = 1
  scale_color_manual(values =  "grey13") +
  theme(axis.ticks.y = element_blank()) +
  plot_common_parameters(x_name = "",
                         y_name = "Pr/Ph",
                         ticks_labels = Siusi_ticks,
                         axis_limit = c(0, 2.7))

# plot (Pr+Ph)/(nC17+nC18)
Siusi_Pr_Ph_FID$PrPh_nC17nC18 <- (Siusi_Pr_Ph_FID$Pristane + Siusi_Pr_Ph_FID$Phytane) / (Siusi_Pr_Ph_FID$nC17 + Siusi_Pr_Ph_FID$nC18)
Siusi_plot_photoautotrophs <- ggplot(Siusi_Pr_Ph_FID, aes(log_height, PrPh_nC17nC18)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Siusi_Pr_Ph_FID, aes(log_height, PrPh_nC17nC18, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  theme(axis.ticks.y = element_blank()) +
  plot_common_parameters(x_name = "",
                         y_name = expression(paste("(Pr+Ph)/(", italic("n"), "-C"[17], "+", italic("n"), "-C"[18], ")")),
                         ticks_labels = Siusi_ticks,
                         axis_limit = c(0, 2.8))

################################################################################
###                              MS data Siusi                               ###
################################################################################

Siusi_ID <- get_metadata("sample_ID", "Siusi")

Siusi_MS_data_F1 <- read.xlsx("Raw_data/MS_data.xlsx", sheet = "Siusi_F1")
Siusi_MS_data_F1 <- merge(Siusi_MS_data_F1, Siusi_ID, by = "sample_ID")

Siusi_MS_data_F2 <- read.xlsx("Raw_data/MS_data.xlsx", sheet = "Siusi_F2")
Siusi_MS_data_F2 <- merge(Siusi_MS_data_F2, Siusi_ID, by = "sample_ID")

# Homohopanes: m/z = 191 (F1)
Siusi_MS_data_F1$HHI31 <- Siusi_MS_data_F1$`abS-C31-HH` / (Siusi_MS_data_F1$`abS-C31-HH` + Siusi_MS_data_F1$`abR-C31-HH`)
Siusi_MS_data_F1$HHI32 <- Siusi_MS_data_F1$`abS-C32-HH` / (Siusi_MS_data_F1$`abS-C32-HH` + Siusi_MS_data_F1$`abR-C32-HH`)

# 2-Methylhopane Index = 2a-methyl-17a,21b-hopane/(2a-methyl-17a,21b-hopane + 17a,21b-hopane)
# 2a-methyl-17a,21b-hopane: m/z = 205 (F1)
# 17a,21b-hopane: m/z = 191 (F1)
Siusi_MS_data_F1$MeHopaneIndex <- (Siusi_MS_data_F1$`2a-methyl-17a21b-hopane`/(Siusi_MS_data_F1$`2a-methyl-17a21b-hopane` + Siusi_MS_data_F1$`17a21b-hopane`))*100
Siusi_MS_data_F1_2MHI_line <- Siusi_MS_data_F1[!is.na(Siusi_MS_data_F1$MeHopaneIndex),]
Siusi_plot_MeHopaneIndex <- ggplot(Siusi_MS_data_F1, aes(log_height, MeHopaneIndex)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Siusi_MS_data_F1_2MHI_line, aes(log_height, MeHopaneIndex, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = expression(paste("2", alpha, "-MHI (%)")),
                         ticks_labels = Siusi_ticks)

# Steranes: m/z = 217 (F1)
Siusi_Steranes <- as.data.frame(cbind(Cholestane = Siusi_MS_data_F1$`5a-Cholestane`,
                                      Ergosterane = Siusi_MS_data_F1$`5a-Ergosterane`,
                                      Stigmastane = Siusi_MS_data_F1$`5a-Stigmastane`))
Siusi_Steranes[is.na(Siusi_Steranes)] <- 0 # replace NAs by 0
Siusi_Steranes$sum <- rowSums(Siusi_Steranes)
Siusi_Steranes$C27_rel <- (Siusi_Steranes$Cholestane/Siusi_Steranes$sum)*100 # relative abundance of Cholestane
Siusi_Steranes$C28_rel <- (Siusi_Steranes$Ergosterane/Siusi_Steranes$sum)*100 # relative abundance of Ergosterane
Siusi_Steranes$C29_rel <- (Siusi_Steranes$Stigmastane/Siusi_Steranes$sum)*100 # relative abundance of Stigmastane
Siusi_Steranes_long <- long_df(as.data.frame(cbind(log_height = Siusi_MS_data_F1$log_height,
                                                   C27 = Siusi_Steranes$C27_rel,
                                                   C28 = Siusi_Steranes$C28_rel,
                                                   C29 = Siusi_Steranes$C29_rel)))
Siusi_plot_steranes_rel <- ggplot(Siusi_Steranes_long, aes(log_height, mass)) +
  geom_point(aes(colour = compound, shape = compound), size = 3, pch = 21, stroke = 1, show.legend = FALSE) +
  geom_line(data = Siusi_Steranes_long, aes(log_height, mass, colour = compound), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(name="",values = c("#99CC00", "#4775FF", "#990080"))+
  plot_common_parameters(y_name = "regular steranes (%)",
                         x_name = "",
                         ticks_labels = Siusi_ticks)

# Pyrene/Phenanthrene
# Pyrene: m/z = 202 (F2)
# Phenanthrene: m/z = 178 (F2)
Siusi_MS_data_F2$pyr_phe <- Siusi_MS_data_F2$pyrene / Siusi_MS_data_F2$phenanthrene

Siusi_plot_pyrene <- ggplot(Siusi_MS_data_F2, aes(log_height, pyr_phe)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Siusi_MS_data_F2, aes(log_height, pyr_phe, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = "pyrene/phenanthrene",
                         ticks_labels = Siusi_ticks)

# Coronene/Phenanthrene
# Coronene: m/z = 300 (F2)
# Phenanthrene: m/z = 178 (F2)
Siusi_MS_data_F2$cor_phe <- Siusi_MS_data_F2$coronene / Siusi_MS_data_F2$phenanthrene

Siusi_plot_coronene <- ggplot(Siusi_MS_data_F2, aes(log_height, cor_phe)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Siusi_MS_data_F2, aes(log_height, cor_phe, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = "coronene/phenanthrene",
                         ticks_labels = Siusi_ticks)

# Perylene/Phenanthrene
# Perylene: m/z = 252 (F2)
# Phenanthrene: m/z = 178 (F2)
Siusi_MS_data_F2$per_phe <- Siusi_MS_data_F2$perylene / Siusi_MS_data_F2$phenanthrene
Siusi_MS_data_F2$per_phe[is.na(Siusi_MS_data_F2$per_phe)] <- 0 # replace NAs by 0
Siusi_plot_perylene <- ggplot(Siusi_MS_data_F2, aes(log_height, per_phe)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Siusi_MS_data_F2, aes(log_height, per_phe, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = "perylene/phenanthrene",
                         ticks_labels = Siusi_ticks)

# correlations
Siusi_Pr_Ph_FID <- merge(Siusi_TOC, Siusi_Pr_Ph_FID, by = "log_height")
cor.test(Siusi_Pr_Ph_FID$TOC,Siusi_Pr_Ph_FID$Pr_Ph) # TOC ~ Pr/Ph -> r = -0.21; p = = 0.43
cor.test(Siusi_Pr_Ph_FID$TOC,Siusi_Pr_Ph_FID$PrPh_nC17nC18) # TOC ~ PrPh_nC17nC18 -> r = 0.58; p = = 0.02

Siusi_Steranes_TAR <- Siusi_Steranes
Siusi_Steranes_TAR$log_height <- Siusi_MS_data_F1$log_height
Siusi_Steranes_TAR <- merge(Siusi_TAR, Siusi_Steranes_TAR, by = "log_height")
cor.test(Siusi_Steranes_TAR$C29_rel, Siusi_Steranes_TAR$TAR) #  TAR ~  %C29 sterane -> r = -0.26; p = 0.31
Siusi_Steranes_TAR_pre <-Siusi_Steranes_TAR[-which(Siusi_Steranes_TAR$log_height > 0),] # remove log heights > 0
cor.test(Siusi_Steranes_TAR_pre$C29_rel, Siusi_Steranes_TAR_pre$TAR) #  TAR ~ %C29 sterane where steranes detected -> r = -0.26; p = 0.31

################################################################################
###                              FID data Seres                              ###
################################################################################

Seres_a_FID <- read_raw_data("Raw_data/FID/Seres", "Seres") # integrated peak areas
Seres_a_FID$log_height <- as.numeric(rownames(Seres_a_FID))

Seres_CPI <- CPI_calc(Seres_a_FID) # calculate carbon preference index (CPI) after Marzie et al. (1993)
Seres_TAR <- TAR_calc(Seres_a_FID) ## calculate terrigenous-aquatic ratio (TAR)

show_metadata("Seres")
Seres_TOC <- get_metadata("TOC", "Seres")

Seres_ticks <- c(NA, NA, "black", NA, NA)

# plot TOC
Seres_plot_TOC <- ggplot(Seres_TOC, aes(log_height, TOC)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Seres_TOC, aes(log_height, TOC, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "height (m)",
                         y_name = "TOC (wt%)",
                         ticks_labels = Seres_ticks,
                         axis_limit = c(0, 4.6))

# plot TAR
Seres_plot_TAR <- ggplot(Seres_TAR, aes(log_height, TAR)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Seres_TAR, aes(log_height, TAR, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = "TAR",
                         ticks_labels = Seres_ticks,
                         axis_limit = c(0, 1.9))

# pristane, phytane, nC17, nC18
Seres_Pr_Ph_FID <- Seres_a_FID[c("Pristane", "Phytane", "nC17", "nC18", "log_height")]
Seres_Pr_Ph_FID$Pr_Ph <- Seres_Pr_Ph_FID$Pristane / Seres_Pr_Ph_FID$Phytane

# plot Pr/Ph
Seres_plot_Pr_Ph_FID <- ggplot(Seres_Pr_Ph_FID, aes(log_height, Pr_Ph)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Seres_Pr_Ph_FID, aes(log_height, Pr_Ph, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  geom_hline(yintercept = 1.0, colour = "grey50", linetype = "dashed", lwd = 0.6) + # vertical line at TAR = 1
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = "Pr/Ph",
                         ticks_labels = Seres_ticks,
                         axis_limit = c(0, 2.7))

# plot (Pr+Ph)/(nC17+nC18)
Seres_Pr_Ph_FID$PrPh_nC17nC18 <- (Seres_Pr_Ph_FID$Pristane + Seres_Pr_Ph_FID$Phytane) / (Seres_Pr_Ph_FID$nC17 + Seres_Pr_Ph_FID$nC18)
Seres_plot_photoautotrophs <- ggplot(Seres_Pr_Ph_FID, aes(log_height, PrPh_nC17nC18)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Seres_Pr_Ph_FID, aes(log_height, PrPh_nC17nC18, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = expression(paste("(Pr+Ph)/(", italic("n"), "-C"[17], "+", italic("n"), "-C"[18], ")")),
                         ticks_labels = Seres_ticks,
                         axis_limit = c(0, 2.8))

################################################################################
###                              MS data Seres                               ###
################################################################################

Seres_ID <- get_metadata("sample_ID", "Seres")

Seres_MS_data_F1 <- read.xlsx("Raw_data/MS_data.xlsx", sheet = "Seres_F1")
Seres_MS_data_F1 <- merge(Seres_MS_data_F1, Seres_ID, by = "sample_ID")

Seres_MS_data_F2 <- read.xlsx("Raw_data/MS_data.xlsx", sheet = "Seres_F2")
Seres_MS_data_F2 <- merge(Seres_MS_data_F2, Seres_ID, by = "sample_ID")

# Homohopanes: m/z = 191 (F1)
Seres_MS_data_F1$HHI31 <- Seres_MS_data_F1$`abS-C31-HH` / (Seres_MS_data_F1$`abS-C31-HH` + Seres_MS_data_F1$`abR-C31-HH`)
Seres_MS_data_F1$HHI32 <- Seres_MS_data_F1$`abS-C32-HH` / (Seres_MS_data_F1$`abS-C32-HH` + Seres_MS_data_F1$`abR-C32-HH`)

# Pyrene/Phenanthrene
# Pyrene: m/z = 202 (F2)
# Phenanthrene: m/z = 178 (F2)
Seres_MS_data_F2$pyr_phe <- Seres_MS_data_F2$pyrene / Seres_MS_data_F2$phenanthrene

Seres_plot_pyrene <- ggplot(Seres_MS_data_F2, aes(log_height, pyr_phe)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Seres_MS_data_F2, aes(log_height, pyr_phe, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = "pyrene/phenanthrene",
                         ticks_labels = Seres_ticks)

# Coronene/Phenanthrene
# Coronene: m/z = 300 (F2)
# Phenanthrene: m/z = 178 (F2)
Seres_MS_data_F2$cor_phe <- Seres_MS_data_F2$coronene / Seres_MS_data_F2$phenanthrene

Seres_plot_coronene <- ggplot(Seres_MS_data_F2, aes(log_height, cor_phe)) +
  geom_point(pch = 21, size = 3, stroke = 1, show.legend = FALSE) +
  geom_line(data = Seres_MS_data_F2, aes(log_height, cor_phe, colour = "grey13"), linetype = "dotted", show.legend = FALSE) +
  scale_color_manual(values = "grey13") +
  plot_common_parameters(x_name = "",
                         y_name = "coronene/phenanthrene",
                         ticks_labels = Seres_ticks)

# correlations
Seres_Pr_Ph_FID <- merge(Seres_TOC, Seres_Pr_Ph_FID, by = "log_height")
cor.test(Seres_Pr_Ph_FID$TOC,Seres_Pr_Ph_FID$Pr_Ph) # TOC ~ Pr/Ph -> r = -0.30; p = = 0.40
cor.test(Seres_Pr_Ph_FID$TOC,Seres_Pr_Ph_FID$PrPh_nC17nC18) # TOC ~ PrPh_nC17nC18 -> r = -0.28; p = = 0.43

################################################################################
###                              Combined plot                               ###
################################################################################

# list the plots in the same order as in ggarrange
plots_sat <- list(Siusi_plot_TOC, Siusi_plot_TAR,  Siusi_plot_Pr_Ph_FID,  Siusi_plot_photoautotrophs,  Siusi_plot_MeHopaneIndex,
                  Siusi_plot_steranes_rel,
                  Seres_plot_TOC, Seres_plot_TAR,  Seres_plot_Pr_Ph_FID,  Seres_plot_photoautotrophs)
# remove y axis on all plots except those in the first column (plots 1 and 6)
plots_sat[[2]] <- remove_axis(plots_sat[[2]])
plots_sat[[3]] <- remove_axis(plots_sat[[3]])
plots_sat[[4]] <- remove_axis(plots_sat[[4]])
plots_sat[[5]] <- remove_axis(plots_sat[[5]])
plots_sat[[6]] <- remove_axis(plots_sat[[6]])
plots_sat[[8]] <- remove_axis(plots_sat[[8]])
plots_sat[[9]] <- remove_axis(plots_sat[[9]])
plots_sat[[10]] <- remove_axis(plots_sat[[10]])
# combine plots
combi_plot_sat <- ggarrange(plots_sat[[1]], plots_sat[[2]], plots_sat[[3]], plots_sat[[4]], plots_sat[[5]],
                            plots_sat[[6]], plots_sat[[7]], plots_sat[[8]], plots_sat[[9]], plots_sat[[10]], 
  nrow = 2, ncol = 6, align = "hv", labels = c("A","", "", "" ,"", "", "B"))
ggsave(combi_plot_sat, file = "Output/Combined_plot_sat.pdf", width = 36, height = 25, units = "cm")

# list plots in the same order as in ggarrange
plots_arom <- list(Siusi_plot_pyrene, Siusi_plot_coronene,
                  Siusi_plot_perylene,
                  Seres_plot_pyrene, Seres_plot_coronene)
# remove y axis on all plots except those in the first column (plots 1 and 6)
plots_arom[[2]] <- remove_axis(plots_arom[[2]])
plots_arom[[3]] <- remove_axis(plots_arom[[3]])
plots_arom[[5]] <- remove_axis(plots_arom[[5]])

# combine plots
combi_plot_arom <- ggarrange(plots_arom[[1]], plots_arom[[2]], plots_arom[[3]],  plots_arom[[4]],
                             plots_arom[[5]],
                            nrow = 2, ncol = 3, align = "hv", labels = c("A","", "", "B"))
ggsave(combi_plot_arom, file = "Output/Combined_plots_arom.pdf", width = 18, height = 22, units = "cm")

################################################################################
###                           Pr_nC17 vs. Ph_nC18                            ###
################################################################################

Siusi_dep_env_data <- prep_depositional_environment_plot(Siusi_a_FID, "Siusi") # for Pr/nC17 vs. Ph/nC18
Seres_dep_env_data <- prep_depositional_environment_plot(Seres_a_FID, "Seres") # for Pr/nC17 vs. Ph/nC18
dep_env_data <- rbind(Siusi_dep_env_data, Seres_dep_env_data)

colors <- c("#9C9EDE", "#8CA252", "#AD494A", "#E7BA52")

plot_enviproxy_FID <- ggplot(dep_env_data, aes(x = Ph_nC18, y = Pr_nC17)) +
  geom_point(aes(shape = section, colour = member), size = 4) +
  plot_common_parameters_PrPh(dep_env_data)
ggsave(plot_enviproxy_FID, file = "Output/PrnC17_PhnC18_FID.pdf", width = 20, height = 15, units = "cm")

################################################################################
###                               MPI Boxplot                                ###
################################################################################

# Methylphenanthrene index (MPI):
# methylphenanthrenes: m/z = 192 (F1)
# Phenanthrene: m/z = 178 (F1)
Siusi_MPI <- MPI_calc(Siusi_MS_data_F2)
Seres_MPI <- MPI_calc(Seres_MS_data_F2)

# all maturity proxies
Siusi_boxplot_prep <- merge(Siusi_MPI, Siusi_CPI, by = "log_height")
Siusi_boxplot_prep <- merge(Siusi_boxplot_prep, as.data.frame(cbind(HHI31 = Siusi_MS_data_F1$HHI31,
                                                                    HHI32 = Siusi_MS_data_F1$HHI32,
                                                                    log_height = Siusi_MS_data_F1$log_height)), by = "log_height")
Siusi_boxplot_prep$section <- rep("Siusi", nrow(Siusi_boxplot_prep))

Seres_boxplot_prep <- merge(Seres_MPI, Seres_CPI, by = "log_height")
Seres_boxplot_prep <- merge(Seres_boxplot_prep, as.data.frame(cbind(HHI31 = Seres_MS_data_F1$HHI31,
                                                                    HHI32 = Seres_MS_data_F1$HHI32,
                                                                    log_height = Seres_MS_data_F1$log_height)), by = "log_height")
Seres_boxplot_prep$section <- rep("Seres", nrow(Seres_boxplot_prep))

maturity <- rbind(Siusi_boxplot_prep, Seres_boxplot_prep)
maturity$section <- with(maturity, factor(section, levels = unique(section))) # to have order as in data frame later
maturity <- maturity[, !names(maturity) %in% c("log_height")] # remove column
maturity_long <- melt(maturity,  id.vars = 'section', variable.name = 'compound') # long format
maturity_long <- na.omit(maturity_long)

maturity_plot<- ggplot(maturity_long, aes(x = variable, y = value, colour = section)) +
  geom_boxplot(lwd = 0.3) +
  scale_y_continuous(name = "") +
  xlab("") +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title.y = element_text(size = 14, color = "black"),
        axis.line = element_line(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = c(0.88, 0.9)) +
scale_color_manual(name="",values = c("coral2", "dodgerblue3")) +
scale_x_discrete(labels=c("CPI1" = "CPI", "HHI31" = expression(paste(HHI ["C31"])),
                          "HHI32" = expression(paste(HHI ["C32"]))))
ggsave(maturity_plot, file = "Output/Maturity_boxplot.pdf", width = 10, height = 12, units = "cm")

################################################################################
###                                    NMDS                                  ###
################################################################################

#prepare Siusi dataset for nmds plot:
Siusi_nmds_prep <- data.frame(cbind(log_height = Siusi_Pr_Ph_FID$log_height, Pr_Ph = Siusi_Pr_Ph_FID$Pr_Ph, photoautotrophs = Siusi_Pr_Ph_FID$PrPh_nC17nC18))
Siusi_nmds_prep <- merge(Siusi_nmds_prep, data.frame(cbind(log_height = Siusi_TAR$log_height, TAR =  Siusi_TAR$TAR)), by = "log_height") # add TAR
Siusi_nmds_prep <- merge(Siusi_nmds_prep, data.frame(cbind(log_height = Siusi_CPI$log_height, CPI =  Siusi_CPI$CPI)), by = "log_height") # add CPI
Siusi_nmds_prep <- merge(Siusi_nmds_prep, data.frame(cbind(log_height = Siusi_MS_data_F2$log_height,
                                                           coronene = Siusi_MS_data_F2$cor_phe, # add coronene/phenanthrene
                                                           pyrene = Siusi_MS_data_F2$pyr_phe)), by = "log_height") # add pyrene/phenanthrene
Siusi_nmds_prep <- merge(Siusi_nmds_prep, Siusi_MPI, by = "log_height")

# prepare Seres dataset for nmds plot:
Seres_nmds_prep <- data.frame(cbind(log_height = Seres_Pr_Ph_FID$log_height, Pr_Ph = Seres_Pr_Ph_FID$Pr_Ph, photoautotrophs = Seres_Pr_Ph_FID$PrPh_nC17nC18))
Seres_nmds_prep <- merge(Seres_nmds_prep, data.frame(cbind(log_height = Seres_TAR$log_height, TAR =  Seres_TAR$TAR)), by = "log_height") # add TAR
Seres_nmds_prep <- merge(Seres_nmds_prep, data.frame(cbind(log_height = Seres_CPI$log_height, CPI =  Seres_CPI$CPI)), by = "log_height") # add CPI
Seres_nmds_prep <- merge(Seres_nmds_prep, data.frame(cbind(log_height = Seres_MS_data_F2$log_height,
                                                           coronene = Seres_MS_data_F2$cor_phe, # add coronene/phenanthrene
                                                          pyrene = Seres_MS_data_F2$pyr_phe)), by = "log_height") # add pyrene/phenanthrene
Seres_nmds_prep <- merge(Seres_nmds_prep, Seres_MPI, by = "log_height")

# combine to one data frame
nmds_prep <- rbind(Siusi_nmds_prep, Seres_nmds_prep)

rownames(nmds_prep) <-nmds_prep$log_height
nmds_prep <- nmds_prep[, !names(nmds_prep) %in% c("log_height")] # remove column
nmds_prep[is.na(nmds_prep)] <- 0 # replace NAs by 0

# chose normalization method
data_normalized <- as.data.frame(decostand(nmds_prep, method = "normalize"))

data_nmds <- metaMDS(data_normalized, # sample-by-compound matrix
                     k=2, # The number of reduced dimensions
                     trymax=100, # number of iterations
                     distance = "bray") # by default Bray-Curtis distance
# plot
stressplot(data_nmds)
plot(data_nmds) # open circles = sites; crosses = compounds

# customize plot
nmds_samples <- as.data.frame(data_nmds$points)
nmds_samples$section <- c(rep("Siusi", nrow(Siusi_nmds_prep)), rep("Seres", nrow(Seres_nmds_prep)))
nmds_samples$log_height <- as.numeric(rownames(nmds_samples))

Siusi_members <- get_metadata("member", "Siusi")
Seres_members <- get_metadata("member", "Seres")
members <- rbind(Siusi_members, Seres_members)
nmds_samples <- merge(nmds_samples, members, by = "log_height")

for(rownumber in 1:nrow(nmds_samples)) {
  if (nmds_samples[rownumber, "log_height"] > 0) {
    nmds_samples$pre_post[rownumber] <- "post-extinction"
  } else {
      nmds_samples$pre_post[rownumber] <- "pre-extinction"}
}

nmds_variables <- as.data.frame(data_nmds$species)
stress <- formatC(data_nmds$stress, digits =  1)

plot_NMDS <- ggplot(nmds_variables, aes(x = MDS1, y =MDS2)) +
  geom_point(data = nmds_variables,  size = 3, pch = 3, color = "grey56") + # compound data
  geom_point(data = nmds_samples, aes(x = MDS1, y = MDS2, shape = section, colour = member), size = 4, alpha = 0.8) +
  geom_text_repel(data = nmds_variables, # label compounds
                  aes(label=rownames(nmds_variables),
                      point.size = 10)) +
  theme(axis.text = element_text(size = 15, colour = "black"),
        axis.title = element_text(size = 15, colour = "black"),
        legend.title = element_blank(),
        legend.position = c(0.83, 0.88),
        legend.text = element_text(size = 14),
        panel.grid = element_blank()) +
  scale_color_manual(values = colors)+
  scale_shape_manual(name = "", values = c(16, 17)) +
  scale_size_manual(name = "", breaks=c("pre-extinction", "post-extinction")) +
  annotate("text", x = -0.2,  y = 0.5,
           label = paste("2D stress = ", stress),
           fontface = "bold", size = 5,
           colour = "gray66")
ggsave(plot_NMDS, file = "Output/nmds.pdf", width = 12, height = 12, units = "cm")

