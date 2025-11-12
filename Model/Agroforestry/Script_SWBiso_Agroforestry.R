### Conceptual plot-scale isotope SoilWaterBalance model, now named isoEcoPlot after Stevenson et al. (submitted) ###
## Christian Birkel, OACG-UCR ##

graphics.off()# clear all histories
rm(list=ls(all=TRUE))
cat("\014") 

library(parallel)
library(data.table)
require(FME)         # package used for calibration and uncertainty estimations
library(tidyverse)   # package used for data manipulation
require(gridExtra)   # package used for graphics with ggplot2    installed
require(corrplot)
require(gplots)
require(RColorBrewer)
require(factoextra)
library(ggplot2)
require(ggpubr)
require(ggsci)
library(scales)
library (lubridate)
#library(patchwork)
library(cowplot)
library(hydroGOF)

# single file for Tempisque plot SW2 model input and evaluation #
inp <- read.csv("Demnitz_inpmod_Ag.csv",na.strings="", sep = ",")     		    # read input into dataframe from .csv file with , as decimal separator
names(inp) <- c("Date","P","NR","AT","RH","LAI","PET","SW_10","SW_30","SW_100", "TF","Tobs","P_D","SW_10_D","SW_30_D","SW_100_D", "TempS_10", "TempS_30", "TempS_100")			            # define column names for direct use in model code
attach(inp,warn.conflicts = F)		
head(inp)
dim(inp)

source("CAL_SWBiso_Forest.R")


# parameter ranges used for calibration
parRange <- data.frame(min=c(-0.6, 0.1, 40, 40,  1,  1,  1,  50, 250, 1.0, 1.0, 1.0, 0.1, 0.5, 1,  3,  10,  0.25, 0.25, 0),
                       max=c(-0.1, 2.0, 60, 60, 20, 20, 20, 100, 450, 5.0, 5.0, 5.0, 0.9, 1.0, 20, 40, 100, 0.90, 0.75, 0))     # specify parameter ranges from min to max for alpha and beta parameters of the gamma df
rownames(parRange) <- c("rE", "alpha", "Smax", "Ic", "ks1", "ks2","ks3", "GWmax","Lmax", "g1", "g2","g3","PF_Scale", "INTp", "stoSp","gwSp","lowSp","k","x","beta")  # parameter names, needs to be consistent with your model and names

write.csv(data.frame(parRange), file = "parRange_forest.csv") 

LH<-Latinhyper(parRange,100000)     # Latin Hypercube MC sampling
#pairs(Latinhyper(parRange, 10000), main = "Latin hypercube")   # check this in a pairwise matrix plot

# Ensure ks1 > ks2 > ks3
ks_indices <- which(rownames(parRange) %in% c("ks1", "ks2", "ks3"))
LH[, ks_indices] <- t(apply(LH[, ks_indices], 1, sort, decreasing = TRUE))

# Ensure Smax/100 > GWmax/200 > Lmax/700
smax_idx <- which(rownames(parRange) == "Smax")
gwmax_idx <- which(rownames(parRange) == "GWmax")
lmax_idx <- which(rownames(parRange) == "Lmax")

# Scale the values and enforce constraints
scaled_values <- cbind(LH[, smax_idx] / 100, LH[, gwmax_idx] / 200, LH[, lmax_idx] / 700)
sorted_scaled_values <- t(apply(scaled_values, 1, sort, decreasing = TRUE))

# Re-scale back to original ranges
LH[, smax_idx] <- sorted_scaled_values[, 1] * 100
LH[, gwmax_idx] <- sorted_scaled_values[, 2] * 200
LH[, lmax_idx] <- sorted_scaled_values[, 3] * 700


## Cong Jiang, 2024.11

# Define the log file path
log_file <- "Ecoplot.log"

# Check if the log file exists and delete it
if (file.exists(log_file)) {
  file.remove(log_file)  # Delete the previous log file
}


Q_DE <- function(pa,inp){model_output <- SWMc (rE = pa[1], alpha=pa[2], Smax=pa[3], Ic=pa[4], ks1=pa[5], ks2=pa[6],ks3=pa[7],
                          GWmax=pa[8],Lmax=pa[9], g1=pa[10], g2=pa[11],g3=pa[12],PF_Scale=pa[13],INTp=pa[14], stoSp=pa[15], gwSp=pa[16],lowSp=pa[17], k=pa[18], x=pa[19], beta=pa[20],inp=inp)

 # Combine parameters and model output
  combined_output <- c(pa, model_output)

  # Return the combined output
  return(combined_output)

}

# Define the number of CPUs from SLURM environment variable
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

# Define a custom sensRange function for parallel processing
sensRange_parallel <- function(func, parInput, inp, num_cores) {
  cl <- makeCluster(num_cores) # Create a cluster
  clusterExport(cl, varlist = c("Q_DE", "inp", "SWMc")) # Export necessary variables to workers

  # Parallel processing: Apply the function to each parameter set
  results <- parLapply(cl, seq_len(nrow(parInput)), function(i) {
    pa <- parInput[i, ] # Extract parameter set
    func(pa, inp)       # Apply the model function
  })

  stopCluster(cl) # Stop the cluster

  # Combine results into a data frame
  results_df <- do.call(rbind, results) # Combine all results into a matrix/data.frame

  # Define column names
  param_names <- colnames(parInput)      # Parameter names

  output_names <-  c(
    "KGE1", "Q_Rsqu", "MAE1", "KGE2", "Q2_Rsqu", "MAE2", "Q3_Rsqu", "KGE3", "MAE3",
    "Q4_Rsqu", "KGE4", "MAE4", "Q5_Rsqu", "KGE5", "MAE5", "Q6_Rsqu", "KGE6", "MAE6",
    "Q7_Rsqu", "KGE7", "MAE7"
  )

  colnames(results_df) <- c(param_names, output_names) # Combine names

  return(results_df)

}

# Run the Modelsim with parallel sensRange
Modelsim <- sensRange_parallel(Q_DE, LH, inp, num_cores)

#write.csv(data.frame(Modelsim), file = "isoEcoPlot_LHres.csv")     # write model output

fwrite(data.frame(Modelsim), file = "isoEcoPlot_LHres.csv")


# analyze results
# inpR <- read.csv("isoEcoPlot_LHres.csv",na.strings="")     		    # read model output into dataframe

# Read the CSV file, set the first row as column names, and treat "NA", "na", and empty strings as missing values (NA)
inpR <- read.csv("isoEcoPlot_LHres.csv", header = TRUE, na.strings = c("NA"))

# Change the first column name to "Id"
# when use the fwrite, there are no default column name
#colnames(inpR)[1] <- "Id"

head(inpR)   # view first 5 rows of data frame

#names(inpR) <- c("Id", "rE","alpha", "Smax", "Ic", "ks1", "ks2","ks3", "GWmax","Lmax", "g1", "g2", "g3","PF_Scale","INTp","stoSp","gwSp", "lowSp","k","x","beta",
#                 "KGE1","Q_Rsqu","MAE1","KGE2","Q2_Rsqu","MAE2","KGE3","Q3_Rsqu3","MAE3",
#                 "KGE4","Q2_Rsqu4","MAE4","KGE5","Q5_Rsqu","MAE5","Q6_Rsqu", "KGE6", "MAE6", "Q7_Rsqu", "KGE7", "MAE7")			            # define column names for direct use in model code

attach(inpR,warn.conflicts = F)	

summary(inpR)  # generate summary statistics

# Check to see if we have any NA values for objective functions that could be indicative of a failed model run using set parameter ranges
##sapply(inpR[11:16], function(x) sum(is.na(x)))

# Remove rows with any NA values in the inpR data frame
inpR <- inpR %>% filter(complete.cases(.)) 

# Filter out runs according to objective function criteria
Threshold_KGE1<- quantile(inpR$KGE1, c(0.6))
Threshold_KGE1<- as.numeric(Threshold_KGE1)

Threshold_KGE2<- quantile(inpR$KGE2, c(0.6))
Threshold_KGE2<- as.numeric(Threshold_KGE2)

Threshold_KGE3<- quantile(inpR$KGE3, c(0.6))
Threshold_KGE3<- as.numeric(Threshold_KGE3)

Threshold_KGE5<- quantile(inpR$KGE5, c(0.6))
Threshold_KGE5<- as.numeric(Threshold_KGE5)

Threshold_KGE6<- quantile(inpR$KGE6, c(0.6))
Threshold_KGE6<- as.numeric(Threshold_KGE6)

Threshold_KGE7<- quantile(inpR$KGE7, c(0.6))
Threshold_KGE7<- as.numeric(Threshold_KGE7)

inpRsel <- inpR %>% filter(KGE1 > Threshold_KGE1 & KGE2 > Threshold_KGE2 & KGE3 > Threshold_KGE3
                           & KGE5 > Threshold_KGE5 & KGE6 > Threshold_KGE6 & KGE7 > Threshold_KGE7)

summary(inpRsel) 
head(inpRsel)
dim(inpRsel)

#write.csv(inpRsel,file="isoEcoPlot_LHres_sel.csv")

fwrite(inpRsel, file = "isoEcoPlot_LHres_sel.csv")

# Calculate KGE_mean and select the top 100 rows with the highest KGE_mean
inpRsel100 <- inpRsel %>%
  mutate(KGE100_mean = (KGE1 + KGE2 + KGE3 + KGE5 + KGE6 + KGE7) / 6) %>%  # Create a new column for KGE_mean
  arrange(desc(KGE100_mean)) %>%                      # Arrange rows by KGE_mean in descending order
  slice(1:100)                                        # Select the top 100 rows

print(summary(inpRsel100))
print(head(inpRsel100))
print(dim(inpRsel100))

# Write the filtered data frame to a CSV file
#write.csv(inpRsel100, file = paste0("isoEcoPlot_LHres_sel_100sims.csv"))

fwrite(inpRsel100, file = "isoEcoPlot_LHres_sel_100sims.csv")

# Plot parameter values against objective functions
## KGE1
Lmax_KGE1<- ggplot(inpR, aes(x = Lmax, y = KGE1)) + 
  geom_point() + theme_classic(base_size = 14) +
  geom_smooth(method='lm', colour = "red", se = FALSE)

ks3_KGE1<- ggplot(inpR, aes(x = ks3, y = KGE1)) + 
  geom_point() + theme_classic(base_size = 14)  +
  geom_smooth(method='lm', colour = "red", se = FALSE)

g3_KGE1<- ggplot(inpR, aes(x = g3, y = KGE1)) + 
  geom_point() + theme_classic(base_size = 14)  +
  geom_smooth(method='lm', colour = "red", se = FALSE)

rE_KGE1<- ggplot(inpR, aes(x = rE, y = KGE1)) + 
  geom_point() + theme_classic(base_size = 14)  +
  geom_smooth(method='lm', colour = "red", se = FALSE)

alpha_KGE1<- ggplot(inpR, aes(x = alpha, y = KGE1)) + 
  geom_point() + theme_classic(base_size = 14)  +
  geom_smooth(method='lm', colour = "red", se = FALSE)

PF_Scale_KGE1<- ggplot(inpR, aes(x = PF_Scale, y = KGE1)) + 
  geom_point() + theme_classic(base_size = 14)  +
  geom_smooth(method='lm', colour = "red", se = FALSE)

lowSp_KGE1<- ggplot(inpR, aes(x = lowSp, y = KGE1)) + 
  geom_point() + theme_classic(base_size = 14)  +
  geom_smooth(method='lm', colour = "red", se = FALSE)

k_KGE1<- ggplot(inpR, aes(x = k, y = KGE1)) + 
  geom_point() + theme_classic(base_size = 14)  +
  geom_smooth(method='lm', colour = "red", se = FALSE)

x_KGE1<- ggplot(inpR, aes(x = x, y = KGE1)) + 
  geom_point() + theme_classic(base_size = 14)  +
  geom_smooth(method='lm', colour = "red", se = FALSE)


KGE1_All<- plot_grid(rE_KGE1, alpha_KGE1, Lmax_KGE1, ks3_KGE1, g3_KGE1, 
                     PF_Scale_KGE1, lowSp_KGE1, k_KGE1, x_KGE1, nrow = 3, ncol = 3)
KGE1_All

# Save if required
png("ForestA_Params_KGE1_Correlation3.png", width = 1250, height = 800)
KGE1_All
dev.off()

#### Simulate output with retained parameters ####

source("out_SWBiso_Forest.R")        # load model script in working directory for simulation

# multi model output
para <- read.csv("isoEcoPlot_LHres_sel.csv",na.strings="") 				# read accepted parameter sets as input into dataframe
para <- para[1:20]
names(para)<- c("rE","alpha", "Smax", "Ic", "ks1", "ks2","k3", "GWmax","Lmax", "g1", "g2", "g3","PF_Scale","INTp","stoSp","gwSp", "lowSp","k","x","beta")
attach(para,warn.conflicts = F)					# make columns of input table directly available

# apply parameter sets to simulate series
Q_sim <- function(pa){SWMsim (rE = pa[1], alpha=pa[2], Smax=pa[3], Ic=pa[4], ks1=pa[5], ks2=pa[6],ks3=pa[7],
                              GWmax=pa[8],Lmax=pa[9], g1=pa[10], g2=pa[11],g3=pa[12],PF_Scale=pa[13],
                              INTp=pa[14], stoSp=pa[15], gwSp=pa[16],lowSp=pa[17], k=pa[18], x=pa[19], beta=pa[20])}   # call and execute gamma model

SWMout <- sensRange(func = Q_sim,parInput = as.matrix(para), map=NULL, num=nrow(para))
#write.csv(data.frame(SWMout), file = "isoEcoPlot_sim.csv")

fwrite(data.frame(SWMout), file = "isoEcoPlot_sim.csv")

#SWMout2 <- as.numeric(unlist(SWMout))
#write.csv(as.array(summary(SWMout2)), file = "summaryisoEcoPlot_sim.csv")

#plot(as.array(summary(SWMout2)),quant=TRUE)

if (FALSE) {

#### Simulate output with retained parameters ####

# single file for Tempisque plot SW2 model input and evaluation #
inp <- read.csv("Demnitz_inpmod_Fo_scen.csv",na.strings="", sep = ",")                   # read input into dataframe from .csv file with , as decimal separator
names(inp) <- c("Date","P","NR","AT","RH","LAI","PET","SW_10","SW_30","SW_100", "TF","Tobs","P_D","SW_10_D","SW_30_D","SW_100_D", "TempS_10", "TempS_30", "TempS_100")                              # define column names for direct use in model code
attach(inp,warn.conflicts = F)
head(inp)
dim(inp)

source("out_SWBiso_Forest.R")        # load model script in working directory for simulation

# multi model output
para <- read.csv("isoEcoPlot_LHres_sel.csv",na.strings="")                              # read accepted parameter sets as input into dataframe
para <- para[3:22]
names(para)<- c("rE","alpha", "Smax", "Ic", "ks1", "ks2","k3", "GWmax","Lmax", "g1", "g2", "g3","PF_Scale","INTp","stoSp","gwSp", "lowSp","k","x","beta")
attach(para,warn.conflicts = F)                                 # make columns of input table directly available

# Define the original and new ranges for scaling
rE_original_range <- c(-0.6, -0.1)   # Original range for rE
rE_new_range <- c(-0.6, -0.1)        # New range for rE

alpha_original_range <- c(0.1, 0.5) # Original range for alpha
alpha_new_range <- c(0.1, 0.5)      # New range for alpha

# Rescale function
rescale <- function(x, old_min, old_max, new_min, new_max) {
  ((x - old_min) / (old_max - old_min)) * (new_max - new_min) + new_min
}

# Scale `rE` and `alpha` to the new ranges
para$rE <- rescale(para$rE, rE_original_range[1], rE_original_range[2], rE_new_range[1], rE_new_range[2])
para$alpha <- rescale(para$alpha, alpha_original_range[1], alpha_original_range[2], alpha_new_range[1], alpha_new_range[2])


# apply parameter sets to simulate series
Q_sim <- function(pa){SWMsim (rE = pa[1], alpha=pa[2], Smax=pa[3], Ic=pa[4], ks1=pa[5], ks2=pa[6],ks3=pa[7],
                              GWmax=pa[8],Lmax=pa[9], g1=pa[10], g2=pa[11],g3=pa[12],PF_Scale=pa[13],
                              INTp=pa[14], stoSp=pa[15], gwSp=pa[16],lowSp=pa[17], k=pa[18], x=pa[19], beta=pa[20])}   # call and execute gamma model

SWMout <- sensRange(func = Q_sim,parInput = as.matrix(para), map=NULL, num=nrow(para))
write.csv(data.frame(SWMout), file = "isoEcoPlot_sim_scen.csv")


}

########################added plots
#PLOT UP OBSERVED AND MODELLED VALUES

# First clear the global environment
#rm(list=ls())

# Clear all objects except SWMout
rm(list = setdiff(ls(), "SWMout"))

# Call Required datasets
observed<- read.csv("Demnitz_inpmod_Ag.csv", sep = ",")
#modelled<- read.csv("isoEcoPlot_sim.csv")

modelled <- SWMout

# Drop paramater values
modelled<- modelled[,-c(1:20)]


ndays = nrow(observed) - 366
#select simulations with KGE potential (soil moisture and isotopes)
modelled_STO <- as.data.frame(t(modelled[1:nrow(modelled),(16*ndays+1):(17*ndays)]))
modelled_STO$STO_m <- observed[367:nrow(observed), 8]
modelled_GW <- as.data.frame(t(modelled[1:nrow(modelled),(17*ndays+1):(18*ndays)]))
modelled_GW$GW_m <- observed[367:nrow(observed), 9]
modelled_LOW <- as.data.frame(t(modelled[1:nrow(modelled),(18*ndays+1):(19*ndays)]))
modelled_LOW$LOW_m <- observed[367:nrow(observed), 10]
modelled_fupCSTO_D <- as.data.frame(t(modelled[1:nrow(modelled),(24*ndays+1):(25*ndays)]))
modelled_fupCSTO_D$SW10D_m <- observed[367:nrow(observed), 14]
modelled_gwCQ <- as.data.frame(t(modelled[1:nrow(modelled),(25*ndays+1):(26*ndays)]))
modelled_gwCQ$GW30D_m <- observed[367:nrow(observed), 15]
modelled_lowCQ <- as.data.frame(t(modelled[1:nrow(modelled),(26*ndays+1):(27*ndays)]))
modelled_lowCQ$LoW100D_m <- observed[367:nrow(observed), 16]

# create empty data frame for KGEs
m <- modelled_STO[1,-c(ncol(modelled_STO))]

# function to estimated KGEs of the simulation
KGEsim <- function(z){
  for(i in 1:(ncol(z)-1))
  {m[,i] <- (KGE(z[,i],z[,ncol(z)], method=c("2012"), na.rm=TRUE))}
  return(m)
}

# use function to create dataframe containing all KGEs for simulated results
mt <- KGEsim(modelled_STO)
mt[2,] <- KGEsim(modelled_GW)
mt[3,] <- KGEsim(modelled_LOW)
mt[4,] <- KGEsim(modelled_fupCSTO_D)
mt[5,] <- KGEsim(modelled_gwCQ)
mt[6,] <- KGEsim(modelled_lowCQ)

# create mean KGE result of soil moisture and isotope KGEs
mt[8,] <- (mt[1,] + mt[2,] + mt[3,] + mt[4,] + mt[5,] + mt[6,])/6

# transpose and rename
mt2<-as.data.frame(t(mt))
names(mt2) <- c("KGE1", "KGE2", "KGE3", "KGE5", "KGE6", "KGE7", "KGEms", "KGEmsi")

#select 100 best simulations
mt_sel <- dplyr::top_n(mt2, 100, KGEmsi)
#get row names, equals each simulation
mt_sel$names <- rownames(mt_sel)
modelled$names <- rownames(modelled)

#Join selected and simulated data
sim100 <- dplyr::left_join(mt_sel, modelled, by = "names")

#test
sim100t <- sim100[,-c(1:9)]
#transpose to fit following script
sim100 <- as.data.frame(t(sim100))

#write.csv(sim100, "KGEs_simulations_100.csv")

fwrite(sim100, file = "KGEs_simulations_100.csv")


if (FALSE) {

#####

#modelled_scen<- read.csv("isoEcoPlot_sim_scen.csv", sep = "," row.names = 1)  # Use the first column as row names


modelled_scen<-SWMout

modelled_scen<- modelled_scen[,-c(1:20)]

head(modelled_scen)

observed_scen<- read.csv("Demnitz_inpmod_Fo_scen.csv", sep = ",")

# Ensure row names exist


#print(dim(mt_sel))          # Dimensions of mt_sel
#print(dim(modelled_scen))   # Dimensions of modelled_scen
#print(head(mt_sel$names))   # Names column in mt_sel
#print(rownames(modelled_scen))  # Row names in modelled_scen



modelled_scen$names <- rownames(modelled_scen)



sim100_scen <- dplyr::left_join(mt_sel, modelled_scen, by = "names")

#transpose to fit following script
sim100_scen <- as.data.frame(t(sim100_scen))

write.csv(sim100_scen, "KGEs_simulations_100_scen.csv")

}

####

N_Variables<- 27 # Number of variables outputted by SWMc_Output function

SCF_TimeSeries<- sim100[((10) : ((nrow(sim100) / N_Variables)+9)),]
names(SCF_TimeSeries) <- c(paste("SCF", 1:ncol(sim100), sep = "_"))

D_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables)) + 10) : ((((nrow(sim100)-9) / N_Variables) * 2)+9),]
names(D_TimeSeries) <- c(paste("D", 1:ncol(sim100), sep = "_"))

Qs_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 2) + 10) : ((((nrow(sim100)-9) / N_Variables) * 3)+9),]
names(Qs_TimeSeries) <- c(paste("Qs", 1:ncol(sim100), sep = "_"))

Th_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 3) + 10) : ((((nrow(sim100)-9) / N_Variables) * 4)+9),]
names(Th_TimeSeries) <- c(paste("Th", 1:ncol(sim100), sep = "_"))

PN_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 4) + 10) : ((((nrow(sim100)-9) / N_Variables) * 5)+9),]
names(PN_TimeSeries) <- c(paste("PN", 1:ncol(sim100), sep = "_"))

I_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 5) + 10) : ((((nrow(sim100)-9) / N_Variables) * 6)+9),]
names(I_TimeSeries) <- c(paste("I", 1:ncol(sim100), sep = "_"))

Ei_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 6) + 10) : ((((nrow(sim100)-9) / N_Variables) * 7)+9),]
names(Ei_TimeSeries) <- c(paste("Ei", 1:ncol(sim100), sep = "_"))

Es_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 7) + 10) : ((((nrow(sim100)-9) / N_Variables) * 8)+9),]
names(Es_TimeSeries) <- c(paste("Es", 1:ncol(sim100), sep = "_"))

Tr_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 8) + 10) : ((((nrow(sim100)-9) / N_Variables) * 9)+9),]
names(Tr_TimeSeries) <- c(paste("Tr", 1:ncol(sim100), sep = "_"))

Tr_Upper_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 9) + 10) : ((((nrow(sim100)-9) / N_Variables) * 10)+9),]
names(Tr_Upper_TimeSeries) <- c(paste("Tr_Upper", 1:ncol(sim100), sep = "_"))

Tr_Lower_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 10) + 10) : ((((nrow(sim100)-9) / N_Variables) * 11)+9),]
names(Tr_Lower_TimeSeries) <- c(paste("Tr_Lower", 1:ncol(sim100), sep = "_"))

Tr_Deep_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 11) + 10) : ((((nrow(sim100)-9) / N_Variables) * 12)+9),]
names(Tr_Deep_TimeSeries) <- c(paste("Tr_Deep", 1:ncol(sim100), sep = "_"))

Perc_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 12) + 10) : ((((nrow(sim100)-9) / N_Variables) * 13)+9),]
names(Perc_TimeSeries) <- c(paste("Perc", 1:ncol(sim100), sep = "_"))

Pref_Flow_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 13) + 10) : ((((nrow(sim100)-9) / N_Variables) * 14)+9),]
names(Pref_Flow_TimeSeries) <- c(paste("Pref_Flow", 1:ncol(sim100), sep = "_"))

Recharge_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 14) + 10) : ((((nrow(sim100)-9) / N_Variables) * 15)+9),]
names(Recharge_TimeSeries) <- c(paste("Recharge", 1:ncol(sim100), sep = "_"))

GWflow_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 15) + 10) : ((((nrow(sim100)-9) / N_Variables) * 16)+9),]
names(GWflow_TimeSeries) <- c(paste("GWflow", 1:ncol(sim100), sep = "_"))

STO_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 16) + 10) : ((((nrow(sim100)-9) / N_Variables) * 17)+9),]
names(STO_TimeSeries) <- c(paste("STO", 1:ncol(sim100), sep = "_"))

GW_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 17) + 10) : ((((nrow(sim100)-9) / N_Variables) * 18)+9),]
names(GW_TimeSeries) <- c(paste("GW", 1:ncol(sim100), sep = "_"))

Sdeep_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 18) + 10) : ((((nrow(sim100)-9) / N_Variables) * 19)+9),]
names(Sdeep_TimeSeries) <- c(paste("Sdeep", 1:ncol(sim100), sep = "_"))

Int_CD_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 19) + 10) : ((((nrow(sim100)-9) / N_Variables) * 20)+9),]
names(Int_CD_TimeSeries) <- c(paste("Int_CD", 1:ncol(sim100), sep = "_"))

Idl_H_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 20) + 10) : ((((nrow(sim100)-9) / N_Variables) * 21)+9),]
names(Idl_H_TimeSeries) <- c(paste("Idl_H", 1:ncol(sim100), sep = "_"))

fInt_CD_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 21) + 10) : ((((nrow(sim100)-9) / N_Variables) * 22)+9),]
names(fInt_CD_TimeSeries) <- c(paste("fInt_CD", 1:ncol(sim100), sep = "_"))

upCSTO_D_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 22) + 10) : ((((nrow(sim100)-9) / N_Variables) * 23)+9),]
names(upCSTO_D_TimeSeries) <- c(paste("upCSTO_D", 1:ncol(sim100), sep = "_"))

Sdl_H_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 23) + 10) : ((((nrow(sim100)-9) / N_Variables) * 24)+9),]
names(Sdl_H_TimeSeries) <- c(paste("Sdl_H", 1:ncol(sim100), sep = "_"))

fupCSTO_D_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 24) + 10) : ((((nrow(sim100)-9) / N_Variables) * 25)+9),]
names(fupCSTO_D_TimeSeries) <- c(paste("fupCSTO_D", 1:ncol(sim100), sep = "_"))

gwCQ_D_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 25) + 10) : ((((nrow(sim100)-9) / N_Variables) * 26)+9),]
names(gwCQ_D_TimeSeries) <- c(paste("gwCQ_D", 1:ncol(sim100), sep = "_"))

lowCQ_D_TimeSeries<- sim100[((((nrow(sim100)-9) / N_Variables) * 26) + 10) : ((((nrow(sim100)-9) / N_Variables) * 27)+9),]
names(lowCQ_D_TimeSeries) <- c(paste("lowCQ_D", 1:ncol(sim100), sep = "_"))


# Merge these columns into one for ease
full_timeseries<- cbind(
  observed[367:nrow(observed),], 
  PN_TimeSeries, D_TimeSeries, I_TimeSeries, SCF_TimeSeries, Ei_TimeSeries,
  Tr_TimeSeries, Es_TimeSeries, STO_TimeSeries, Qs_TimeSeries,
  GW_TimeSeries, Sdeep_TimeSeries, GWflow_TimeSeries, Recharge_TimeSeries, 
  Th_TimeSeries, Perc_TimeSeries, Pref_Flow_TimeSeries,
  Int_CD_TimeSeries, fInt_CD_TimeSeries, upCSTO_D_TimeSeries, Sdl_H_TimeSeries,
  fupCSTO_D_TimeSeries, Idl_H_TimeSeries, gwCQ_D_TimeSeries, lowCQ_D_TimeSeries,
  Tr_Upper_TimeSeries, Tr_Lower_TimeSeries, Tr_Deep_TimeSeries
)


#set Date as Date value in R
full_timeseries$Date<-as.Date(strptime(full_timeseries[,1], format="%d.%m.%Y",tz=""))
for (i in 2:(ncol(full_timeseries)-1)){
  full_timeseries[ ,i]<-as.numeric(as.character(full_timeseries[ ,i]))
}

# Export this timeseries
# # deltel fist
#write.csv(full_timeseries, "All_Sim_Obs_ForestA_Time_Series.csv", row.names = F)

fwrite(full_timeseries, file = "All_Sim_Obs_ForestA_Time_Series.csv", row.names = FALSE)


if (FALSE) {

####

N_Variables<- 27 # Number of variables outputted by SWMc_Output function

SCF_TimeSeries_scen<- sim100_scen[((10) : ((nrow(sim100_scen) / N_Variables)+9)),]
names(SCF_TimeSeries_scen) <- c(paste("SCF", 1:ncol(sim100_scen), sep = "_"))

D_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables)) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 2)+9),]
names(D_TimeSeries_scen) <- c(paste("D", 1:ncol(sim100_scen), sep = "_"))

Qs_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 2) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 3)+9),]
names(Qs_TimeSeries_scen) <- c(paste("Qs", 1:ncol(sim100_scen), sep = "_"))

Th_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 3) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 4)+9),]
names(Th_TimeSeries_scen) <- c(paste("Th", 1:ncol(sim100_scen), sep = "_"))

PN_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 4) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 5)+9),]
names(PN_TimeSeries_scen) <- c(paste("PN", 1:ncol(sim100_scen), sep = "_"))

I_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 5) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 6)+9),]
names(I_TimeSeries_scen) <- c(paste("I", 1:ncol(sim100_scen), sep = "_"))

Ei_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 6) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 7)+9),]
names(Ei_TimeSeries_scen) <- c(paste("Ei", 1:ncol(sim100_scen), sep = "_"))

Es_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 7) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 8)+9),]
names(Es_TimeSeries_scen) <- c(paste("Es", 1:ncol(sim100_scen), sep = "_"))

Tr_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 8) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 9)+9),]
names(Tr_TimeSeries_scen) <- c(paste("Tr", 1:ncol(sim100_scen), sep = "_"))

Tr_Upper_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 9) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 10)+9),]
names(Tr_Upper_TimeSeries_scen) <- c(paste("Tr_Upper", 1:ncol(sim100_scen), sep = "_"))

Tr_Lower_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 10) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 11)+9),]
names(Tr_Lower_TimeSeries_scen) <- c(paste("Tr_Lower", 1:ncol(sim100_scen), sep = "_"))

Tr_Deep_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 11) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 12)+9),]
names(Tr_Deep_TimeSeries_scen) <- c(paste("Tr_Deep", 1:ncol(sim100_scen), sep = "_"))

Perc_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 12) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 13)+9),]
names(Perc_TimeSeries_scen) <- c(paste("Perc", 1:ncol(sim100_scen), sep = "_"))

Pref_Flow_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 13) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 14)+9),]
names(Pref_Flow_TimeSeries_scen) <- c(paste("Pref_Flow", 1:ncol(sim100_scen), sep = "_"))

Recharge_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 14) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 15)+9),]
names(Recharge_TimeSeries_scen) <- c(paste("Recharge", 1:ncol(sim100_scen), sep = "_"))

GWflow_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 15) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 16)+9),]
names(GWflow_TimeSeries_scen) <- c(paste("GWflow", 1:ncol(sim100_scen), sep = "_"))

STO_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 16) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 17)+9),]
names(STO_TimeSeries_scen) <- c(paste("STO", 1:ncol(sim100_scen), sep = "_"))

GW_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 17) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 18)+9),]
names(GW_TimeSeries_scen) <- c(paste("GW", 1:ncol(sim100_scen), sep = "_"))

Sdeep_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 18) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 19)+9),]
names(Sdeep_TimeSeries_scen) <- c(paste("Sdeep", 1:ncol(sim100), sep = "_"))

Int_CD_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 19) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 20)+9),]
names(Int_CD_TimeSeries_scen) <- c(paste("Int_CD", 1:ncol(sim100_scen), sep = "_"))

Idl_H_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 20) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 21)+9),]
names(Idl_H_TimeSeries_scen) <- c(paste("Idl_H", 1:ncol(sim100_scen), sep = "_"))

fInt_CD_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 21) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 22)+9),]
names(fInt_CD_TimeSeries_scen) <- c(paste("fInt_CD", 1:ncol(sim100_scen), sep = "_"))

upCSTO_D_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 22) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 23)+9),]
names(upCSTO_D_TimeSeries_scen) <- c(paste("upCSTO_D", 1:ncol(sim100_scen), sep = "_"))

Sdl_H_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 23) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 24)+9),]
names(Sdl_H_TimeSeries_scen) <- c(paste("Sdl_H", 1:ncol(sim100), sep = "_"))

fupCSTO_D_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 24) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 25)+9),]
names(fupCSTO_D_TimeSeries_scen) <- c(paste("fupCSTO_D", 1:ncol(sim100_scen), sep = "_"))

gwCQ_D_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 25) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 26)+9),]
names(gwCQ_D_TimeSeries_scen) <- c(paste("gwCQ_D", 1:ncol(sim100_scen), sep = "_"))

lowCQ_D_TimeSeries_scen<- sim100_scen[((((nrow(sim100_scen)-9) / N_Variables) * 26) + 10) : ((((nrow(sim100_scen)-9) / N_Variables) * 27)+9),]
names(lowCQ_D_TimeSeries_scen) <- c(paste("lowCQ_D", 1:ncol(sim100_scen), sep = "_"))



#print(dim(observed_scen[366:nrow (observed_scen),]))   # Dimensions of modelled_scen
#print(dim(PN_TimeSeries_scen))


# Merge these columns into one for ease
full_timeseries_scen<- cbind(
  observed_scen[367:nrow (observed_scen),], 
  PN_TimeSeries_scen, D_TimeSeries_scen, I_TimeSeries_scen, SCF_TimeSeries_scen, Ei_TimeSeries_scen,
  Tr_TimeSeries_scen, Es_TimeSeries_scen, STO_TimeSeries_scen, Qs_TimeSeries_scen,
  GW_TimeSeries_scen, Sdeep_TimeSeries_scen, GWflow_TimeSeries_scen, Recharge_TimeSeries_scen, 
  Th_TimeSeries_scen, Perc_TimeSeries_scen, Pref_Flow_TimeSeries_scen,
  Int_CD_TimeSeries_scen, fInt_CD_TimeSeries_scen, upCSTO_D_TimeSeries_scen, Sdl_H_TimeSeries_scen,
  fupCSTO_D_TimeSeries_scen, Idl_H_TimeSeries_scen, gwCQ_D_TimeSeries_scen, lowCQ_D_TimeSeries_scen,
  Tr_Upper_TimeSeries_scen, Tr_Lower_TimeSeries_scen, Tr_Deep_TimeSeries_scen
)


#set Date as Date value in R
full_timeseries_scen$Date<-as.Date(strptime(full_timeseries_scen[,1], format="%d.%m.%Y",tz=""))
for (i in 2:(ncol(full_timeseries_scen)-1)){
  full_timeseries_scen[ ,i]<-as.numeric(as.character(full_timeseries_scen[ ,i]))
}

# Export this timeseries
write.csv(full_timeseries_scen[-1, ], "All_Sim_Obs_ForestA_Time_Series_scen.csv", row.names = F)


}


# parameters for water balance estimation
water_balance_param <- select(full_timeseries, starts_with(c("P_mm", "LAI", "SCF", "PN", 
                                                             "Qs_", "I_", "Th_", "Ei_",
                                                             "Tr_", "Es_", "Recharge")))
# set parameters to numeric
for (i in 1:ncol(water_balance_param)){
  water_balance_param[ ,i]<-as.numeric(as.character(water_balance_param[ ,i]))
}
# transpose
water_balance <- t(apply(water_balance_param, 2, FUN = sum, na.rm = T))
# create new data frame
water_balance2 <- data.frame(matrix(nrow = 5, ncol = 1))
water_balance2[,1] <- c("max", "mean", "min", "sd", "quantile50")
c("P_mm", "LAI", "PN", "Qs", "I", "Th", "Ei",
  "Tr", "Tr_Upper", "Tr_Lower", "Tr_Deep", "Es_", "Recharge")
#  estimate values for water balance
water_balance2$PN[1] <- max(c(water_balance[,(((ncol(water_balance)-2)/12)+3):(((ncol(water_balance)-2)/12)*2+2)]))
water_balance2$PN[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)+3):(((ncol(water_balance)-2)/12)*2+2)]))
water_balance2$PN[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)+3):(((ncol(water_balance)-2)/12)*2+2)]))
water_balance2$PN[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)+3):(((ncol(water_balance)-2)/12)*2+2)]))
water_balance2$PN[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)+3):(((ncol(water_balance)-2)/12)*2+2)]),0.5)

water_balance2$Qs[1] <- max(c(water_balance[,((((ncol(water_balance)-2)/12)*2)+3):(((ncol(water_balance)-2)/12)*3+2)]))
water_balance2$Qs[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)*2+3):(((ncol(water_balance)-2)/12)*3+2)]))
water_balance2$Qs[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)*2+3):(((ncol(water_balance)-2)/12)*3+2)]))
water_balance2$Qs[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)*2+3):(((ncol(water_balance)-2)/12)*3+2)]))
water_balance2$Qs[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)*2+3):(((ncol(water_balance)-2)/12)*3+2)]),0.5)

water_balance2$I[1] <- max(c(water_balance[,(((ncol(water_balance)-2)/12)*3+3):(((ncol(water_balance)-2)/12)*4+2)]))
water_balance2$I[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)*3+3):(((ncol(water_balance)-2)/12)*4+2)]))
water_balance2$I[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)*3+3):(((ncol(water_balance)-2)/12)*4+2)]))
water_balance2$I[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)*3+3):(((ncol(water_balance)-2)/12)*4+2)]))
water_balance2$I[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)*3+3):(((ncol(water_balance)-2)/12)*4+2)]),0.5)

water_balance2$Th[1] <- max(c(water_balance[,(((ncol(water_balance)-2)/12)*4+3):(((ncol(water_balance)-2)/12)*5+2)]))
water_balance2$Th[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)*4+3):(((ncol(water_balance)-2)/12)*5+2)]))
water_balance2$Th[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)*4+3):(((ncol(water_balance)-2)/12)*5+2)]))
water_balance2$Th[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)*4+3):(((ncol(water_balance)-2)/12)*5+2)]))
water_balance2$Th[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)*4+3):(((ncol(water_balance)-2)/12)*5+2)]),0.5)

water_balance2$Ei[1] <- max(c(water_balance[,(((ncol(water_balance)-2)/12)*5+3):(((ncol(water_balance)-2)/12)*6+2)]))
water_balance2$Ei[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)*5+3):(((ncol(water_balance)-2)/12)*6+2)]))
water_balance2$Ei[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)*5+3):(((ncol(water_balance)-2)/12)*6+2)]))
water_balance2$Ei[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)*5+3):(((ncol(water_balance)-2)/12)*6+2)]))
water_balance2$Ei[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)*5+3):(((ncol(water_balance)-2)/12)*6+2)]),0.5)

water_balance2$Tr[1] <- max(c(water_balance[,(((ncol(water_balance)-2)/12)*6+3):(((ncol(water_balance)-2)/12)*7+2)]))
water_balance2$Tr[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)*6+3):(((ncol(water_balance)-2)/12)*7+2)]))
water_balance2$Tr[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)*6+3):(((ncol(water_balance)-2)/12)*7+2)]))
water_balance2$Tr[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)*6+3):(((ncol(water_balance)-2)/12)*7+2)]))
water_balance2$Tr[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)*6+3):(((ncol(water_balance)-2)/12)*7+2)]),0.5)

water_balance2$Tr_Upper[1] <- max(c(water_balance[,(((ncol(water_balance)-2)/12)*7+3):(((ncol(water_balance)-2)/12)*8+2)]))
water_balance2$Tr_Upper[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)*7+3):(((ncol(water_balance)-2)/12)*8+2)]))
water_balance2$Tr_Upper[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)*7+3):(((ncol(water_balance)-2)/12)*8+2)]))
water_balance2$Tr_Upper[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)*7+3):(((ncol(water_balance)-2)/12)*8+2)]))
water_balance2$Tr_Upper[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)*7+3):(((ncol(water_balance)-2)/12)*8+2)]),0.5)

water_balance2$Tr_Lower[1] <- max(c(water_balance[,(((ncol(water_balance)-2)/12)*8+3):(((ncol(water_balance)-2)/12)*9+2)]))
water_balance2$Tr_Lower[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)*8+3):(((ncol(water_balance)-2)/12)*9+2)]))
water_balance2$Tr_Lower[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)*8+3):(((ncol(water_balance)-2)/12)*9+2)]))
water_balance2$Tr_Lower[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)*8+3):(((ncol(water_balance)-2)/12)*9+2)]))
water_balance2$Tr_Lower[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)*8+3):(((ncol(water_balance)-2)/12)*9+2)]),0.5)

water_balance2$Tr_Deeper[1] <- max(c(water_balance[,(((ncol(water_balance)-2)/12)*9+3):(((ncol(water_balance)-2)/12)*10+2)]))
water_balance2$Tr_Deeper[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)*9+3):(((ncol(water_balance)-2)/12)*10+2)]))
water_balance2$Tr_Deeper[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)*9+3):(((ncol(water_balance)-2)/12)*10+2)]))
water_balance2$Tr_Deeper[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)*9+3):(((ncol(water_balance)-2)/12)*10+2)]))
water_balance2$Tr_Deeper[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)*9+3):(((ncol(water_balance)-2)/12)*10+2)]),0.5)

water_balance2$Es[1] <- max(c(water_balance[,(((ncol(water_balance)-2)/12)*10+3):(((ncol(water_balance)-2)/12)*11+2)]))
water_balance2$Es[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)*10+3):(((ncol(water_balance)-2)/12)*11+2)]))
water_balance2$Es[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)*10+3):(((ncol(water_balance)-2)/12)*11+2)]))
water_balance2$Es[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)*10+3):(((ncol(water_balance)-2)/12)*11+2)]))
water_balance2$Es[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)*10+3):(((ncol(water_balance)-2)/12)*11+2)]),0.5)

water_balance2$Recharge[1] <- max(c(water_balance[,(((ncol(water_balance)-2)/12)*11+3):(((ncol(water_balance)-2)/12)*12+2)]))
water_balance2$Recharge[2] <- mean(c(water_balance[,(((ncol(water_balance)-2)/12)*11+3):(((ncol(water_balance)-2)/12)*12+2)]))
water_balance2$Recharge[3] <- min(c(water_balance[,(((ncol(water_balance)-2)/12)*11+3):(((ncol(water_balance)-2)/12)*12+2)]))
water_balance2$Recharge[4] <- sd(c(water_balance[,(((ncol(water_balance)-2)/12)*11+3):(((ncol(water_balance)-2)/12)*12+2)]))
water_balance2$Recharge[5] <- quantile(c(water_balance[,(((ncol(water_balance)-2)/12)*11+3):(((ncol(water_balance)-2)/12)*12+2)]),0.5)

#estimate ET and Recharge+ET+Qs
water_balance2$ET <- water_balance2$Ei + water_balance2$Es + water_balance2$Tr
water_balance2$Re_ET_Qs <- water_balance2$ET + water_balance2$Recharge + water_balance2$Qs

# write water balance file
write.csv(water_balance2, "water_balance_Forest_sim100.csv", row.names=FALSE)


#create min and max for STO presentation
full_timeseries <- full_timeseries %>%
  mutate(minSTO = select(., starts_with("STO")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxSTO = select(., starts_with("STO")) %>% do.call(pmax, .))
full_timeseries$meanSTO <- rowMeans(select(full_timeseries, starts_with(c("STO"))), na.rm = TRUE)

#create min and max for GW presentation
full_timeseries <- full_timeseries %>%
  mutate(minGW = select(., starts_with("GW_")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxGW = select(., starts_with("GW_")) %>% do.call(pmax, .))
full_timeseries$meanGW <- rowMeans(select(full_timeseries, starts_with(c("GW_"))), na.rm = TRUE)

#create min and max for Low presentation
full_timeseries <- full_timeseries %>%
  mutate(minLow = select(., starts_with("Sdeep")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxLow = select(., starts_with("Sdeep")) %>% do.call(pmax, .))
full_timeseries$meanLow <- rowMeans(select(full_timeseries, starts_with(c("Sdeep"))), na.rm = TRUE)

#create min and max for I presentation
full_timeseries <- full_timeseries %>%
  mutate(minI = select(., starts_with("I")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxI = select(., starts_with("I")) %>% do.call(pmax, .))
full_timeseries$meanI <- rowMeans(select(full_timeseries, starts_with(c("I"))), na.rm = TRUE)

#create min and max for Ei presentation
full_timeseries <- full_timeseries %>%
  mutate(minEi = select(., starts_with("Ei")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxEi = select(., starts_with("Ei")) %>% do.call(pmax, .))
full_timeseries$meanEi <- rowMeans(select(full_timeseries, starts_with(c("Ei"))), na.rm = TRUE)

#create min and max for Tr presentation
full_timeseries <- full_timeseries %>%
  mutate(minTr = select(., matches("^Tr_[0-9]+$")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxTr = select(., matches("^Tr_[0-9]+$")) %>% do.call(pmax, .))
full_timeseries$meanTr <- rowMeans(select(full_timeseries, matches("^Tr_[0-9]+$")), na.rm = TRUE)

#create min and max for Es presentation
full_timeseries <- full_timeseries %>%
  mutate(minEs = select(., starts_with("Es")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxEs = select(., starts_with("Es")) %>% do.call(pmax, .))
full_timeseries$meanEs <- rowMeans(select(full_timeseries, starts_with(c("Es"))), na.rm = TRUE)

#create min and max for Recharge presentation
full_timeseries <- full_timeseries %>%
  mutate(minRecharge = select(., starts_with("Recharge")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxRecharge = select(., starts_with("Recharge")) %>% do.call(pmax, .))
full_timeseries$meanRecharge <- rowMeans(select(full_timeseries, starts_with(c("Recharge"))), na.rm = TRUE)

# create min and max for fupCSTO_D presentation
full_timeseries <- full_timeseries %>%
  mutate(minfupCSTO_D = select(., starts_with("fupCSTO_D")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxfupCSTO_D = select(., starts_with("fupCSTO_D")) %>% do.call(pmax, .))
full_timeseries$meanfupCSTO_D <- rowMeans(select(full_timeseries, starts_with(c("fupCSTO_D"))), na.rm = TRUE)

# create min and max for gwCQ_D presentation 
full_timeseries <- full_timeseries %>%
  mutate(mingwCQ_D = select(., starts_with("gwCQ_D")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxgwCQ_D = select(., starts_with("gwCQ_D")) %>% do.call(pmax, .))
full_timeseries$meangwCQ_D <- rowMeans(select(full_timeseries, starts_with(c("gwCQ_D"))), na.rm = TRUE)

# create min and max for lowCQ_D presentation
full_timeseries <- full_timeseries %>%
  mutate(minlowCQ_D = select(., starts_with("lowCQ_D")) %>% do.call(pmin, .))
full_timeseries <- full_timeseries %>%
  mutate(maxlowCQ_D = select(., starts_with("lowCQ_D")) %>% do.call(pmax, .))
full_timeseries$meanlowCQ_D <- rowMeans(select(full_timeseries, starts_with(c("lowCQ_D"))), na.rm = TRUE)


if (FALSE) {

#create min and max for STO presentation
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(minSTO_scen = select(., starts_with("STO")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxSTO_scen = select(., starts_with("STO")) %>% do.call(pmax, .))
full_timeseries_scen$meanSTO_scen <- rowMeans(select(full_timeseries_scen, starts_with(c("STO"))), na.rm = TRUE)

#create min and max for GW presentation
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(minGW_scen = select(., starts_with("GW_")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxGW_scen = select(., starts_with("GW_")) %>% do.call(pmax, .))
full_timeseries_scen$meanGW_scen <- rowMeans(select(full_timeseries_scen, starts_with(c("GW_"))), na.rm = TRUE)

#create min and max for Low presentation
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(minLow_scen = select(., starts_with("Sdeep")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxLow_scen = select(., starts_with("Sdeep")) %>% do.call(pmax, .))
full_timeseries_scen$meanLow_scen <- rowMeans(select(full_timeseries_scen, starts_with(c("Sdeep"))), na.rm = TRUE)

#create min and max for I presentation
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(minI_scen = select(., starts_with("I")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxI_scen = select(., starts_with("I")) %>% do.call(pmax, .))
full_timeseries_scen$meanI_scen <- rowMeans(select(full_timeseries_scen, starts_with(c("I"))), na.rm = TRUE)

#create min and max for Ei presentation
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(minEi_scen = select(., starts_with("Ei")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxEi_scen = select(., starts_with("Ei")) %>% do.call(pmax, .))
full_timeseries_scen$meanEi_scen <- rowMeans(select(full_timeseries_scen, starts_with(c("Ei"))), na.rm = TRUE)

#create min and max for Tr presentation
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(minTr_scen = select(., matches("^Tr_[0-9]+$")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxTr_scen = select(., matches("^Tr_[0-9]+$")) %>% do.call(pmax, .))
full_timeseries_scen$meanTr_scen <- rowMeans(select(full_timeseries_scen, matches("^Tr_[0-9]+$")), na.rm = TRUE)

#create min and max for Es presentation
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(minEs_scen = select(., starts_with("Es")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxEs_scen = select(., starts_with("Es")) %>% do.call(pmax, .))
full_timeseries_scen$meanEs_scen <- rowMeans(select(full_timeseries_scen, starts_with(c("Es"))), na.rm = TRUE)

#create min and max for Recharge presentation
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(minRecharge_scen = select(., starts_with("Recharge")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxRecharge_scen = select(., starts_with("Recharge")) %>% do.call(pmax, .))
full_timeseries_scen$meanRecharge_scen <- rowMeans(select(full_timeseries_scen, starts_with(c("Recharge"))), na.rm = TRUE)

# create min and max for fupCSTO_D presentation
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(minfupCSTO_D_scen = select(., starts_with("fupCSTO_D")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxfupCSTO_D_scen = select(., starts_with("fupCSTO_D")) %>% do.call(pmax, .))
full_timeseries_scen$meanfupCSTO_D_scen <- rowMeans(select(full_timeseries_scen, starts_with(c("fupCSTO_D"))), na.rm = TRUE)

# create min and max for gwCQ_D presentation 
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(mingwCQ_D_scen = select(., starts_with("gwCQ_D")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxgwCQ_D_scen = select(., starts_with("gwCQ_D")) %>% do.call(pmax, .))
full_timeseries_scen$meangwCQ_D_scen <- rowMeans(select(full_timeseries_scen, starts_with(c("gwCQ_D"))), na.rm = TRUE)

# create min and max for lowCQ_D presentation
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(minlowCQ_D_scen = select(., starts_with("lowCQ_D")) %>% do.call(pmin, .))
full_timeseries_scen <- full_timeseries_scen %>%
  mutate(maxlowCQ_D_scen = select(., starts_with("lowCQ_D")) %>% do.call(pmax, .))
full_timeseries_scen$meanlowCQ_D_scen <- rowMeans(select(full_timeseries_scen, starts_with(c("lowCQ_D"))), na.rm = TRUE)

#########################################################

}

#########################################################
#plotting of the results

# Upper Soil Store Observed versus Simulated
P1 <- ggplot(data = full_timeseries) +
  # geom_line(aes(x = Date, y = Tr_Upper_Q50/0.04), linewidth = 1.5, colour = "green") +  
  geom_line(aes(x = Date, y = meanSTO), linewidth = 1.5, colour = "black") +
#  geom_line(aes(x = Date, y = full_timeseries_scen$meanSTO_scen), linewidth = 1.5, colour = "yellow") +
  geom_ribbon(aes(x = Date, ymin= minSTO, ymax= maxSTO), linetype=2, alpha=0.4, fill="red") +
  geom_point(aes(x = Date, y = Upper_SM_mm), size = 5, colour = "blue") +
  # geom_errorbar(aes(x = Date, ymin = Moi_mm_min, ymax = Moi_mm_max),width = 10, color = "blue")+
  labs(x = "", y = "Upper soil store (mm)") +
  ylim(0,50)+
  ggtitle("Forest") +
  theme_classic(base_size = 40) +
  theme(plot.title = element_text(hjust = 0.5)) +
  # scale_y_continuous(limits = c(0, 50), sec.axis = sec_axis(~.*0.04, name=""))+
  scale_x_date(date_breaks = "3 month", date_labels = "%b",
               limits = as.Date(c('2021-01-01','2021-12-31')))

P1

# Lower Soil Store Observed versus Simulated
P2 <- ggplot(data = full_timeseries) +
  #geom_line(aes(x = Date, y = Tr_Lower_Q50/(1/40)), linewidth = 1.5, colour = "green") +    
  geom_line(aes(x = Date, y = meanGW), linewidth = 1.5, colour = "black") +
  geom_line(aes(x = Date, y = full_timeseries_scen$meanGW_scen), linewidth = 1.5, colour = "yellow") +
#  geom_ribbon(aes(x = Date, ymin= minGW, ymax= maxGW), linetype=2, alpha=0.4, fill="red") +
  geom_line(aes(x = Date, y = Medium_SM_mm),  linewidth = 1.5, colour = "blue") +
  ylim(0,60)+
  labs(x = "", y = "Lower soil store (mm)") +
  theme_classic(base_size = 40) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # scale_y_continuous(limits = c(0, 80), sec.axis = sec_axis(~.*(1/40), name=""))+
  scale_x_date(date_breaks = "3 month", date_labels = "%b",
               limits = as.Date(c('2021-01-01','2021-12-31')))
P2

# Deeper Soil Store Observed versus Simulated
P3 <- ggplot(data = full_timeseries) +
  # geom_line(aes(x = Date, y = Tr_Deep_Q50/(2/300)), linewidth = 1.5, colour = "green") +    
  geom_line(aes(x = Date, y = meanLow), linewidth = 1.5, colour = "black") +
#  geom_line(aes(x = Date, y = full_timeseries_scen$meanLow_scen), linewidth = 1.5, colour = "yellow") +
  geom_ribbon(aes(x = Date, ymin= minLow, ymax= maxLow), linetype=2, alpha=0.4, fill="red") +
  geom_line(aes(x = Date, y = Lower_SM_mm),  linewidth = 1.5, colour = "blue") +
  labs(x = "", y = "Deeper soil store (mm)") +
  ylim(0,200)+
  theme_classic(base_size = 40) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  # scale_y_continuous(limits = c(0, 300), sec.axis = sec_axis(~.*(2/300), name=""))+
  scale_x_date(date_breaks = "3 month", date_labels = "%b", 
               limits = as.Date(c('2021-01-01','2021-12-31')))
P3


Plots_All<- plot_grid(P1, P2, P3,nrow = 3, ncol = 1)
# Save if required
png("Forest_moi.png", width = 3200, height = 1600)
Plots_All
dev.off()

# Isotope dynamics plots
# Plot for Upper Soil Stable Water Isotope Dynamics
P4 <- ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = meanfupCSTO_D), linewidth = 1.5, colour = "black") +
  geom_ribbon(aes(x = Date, ymin = minfupCSTO_D, ymax = maxfupCSTO_D), linetype = 2, alpha = 0.4, fill = "blue") +
  geom_point(aes(x = Date, y = Upper_2H), size = 5, colour = "green") +
  labs(x = "", y = "Upper soil isotope (‰)") +
  ylim(-120,0)+
  ggtitle("Forest") +
  theme_classic(base_size = 40) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "3 month", date_labels = "%b",
               limits = as.Date(c('2021-01-01', '2021-12-31')))

# Plot for Lower Soil Stable Water Isotope Dynamics
P5 <- ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = meangwCQ_D), linewidth = 1.5, colour = "black") +
  geom_ribbon(aes(x = Date, ymin = mingwCQ_D, ymax = maxgwCQ_D), linetype = 2, alpha = 0.4, fill = "blue") +
  geom_point(aes(x = Date, y = Medium_2H), size = 5, colour = "green") +
  labs(x = "", y = "Lower soil isotope (‰)") +
  ylim(-120,0)+
  theme_classic(base_size = 40) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "3 month", date_labels = "%b",
               limits = as.Date(c('2021-01-01', '2021-12-31')))

# Plot for Deeper Soil Stable Water Isotope Dynamics
P6 <- ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = meanlowCQ_D), linewidth = 1.5, colour = "black") +
  geom_ribbon(aes(x = Date, ymin = minlowCQ_D, ymax = maxlowCQ_D), linetype = 2, alpha = 0.4, fill = "blue") +
  geom_point(aes(x = Date, y = Lower_2H), size = 5, colour = "green") +
  labs(x = "", y = "Deeper soil isotope (‰)") +
  ylim(-120,0)+
  theme_classic(base_size = 40) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "3 month", date_labels = "%b",
               limits = as.Date(c('2021-01-01', '2021-12-31')))

# Combine and save isotope dynamics plots
Isotope_Plots <- plot_grid(P4, P5, P6, nrow = 3, ncol = 1)
png("Forest_isotope_dynamics.png", width = 3200, height = 1600)
print(Isotope_Plots)
dev.off()

#######################
#careful I did not edit those to fit the current selection process
#maybe it's still helpful and I can help you to edit those if required

# Simulated Net P
P4<-ggplot(data = full_timeseries[!is.na(full_timeseries$P_mm),]) +
  geom_line(aes(x = Date, y = PN_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= PN_Q10, ymax= PN_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "PN - Net Precipitation (mm)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P4          

# Simulated Surface Cover Fraction
P5<-ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = SCF_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= SCF_Q10, ymax= SCF_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "SCF - Surface Cover Fraction (-)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters") +
  geom_hline(yintercept = 1, colour = "green", linewidth = 1) +
  geom_hline(yintercept = 0.80, colour = "red", linewidth = 1, linetype = 4)
P5

# Simulated Interception
P6<-ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = D_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= D_Q10, ymax= D_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "D - Intercepted Volume Per Timestep (mm)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P6

# Interception Store
P7<-ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = I_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= I_Q10, ymax= I_Q90), linetype=2, alpha=0.4, fill="red") +
  geom_line(aes(x = Date, y = ThroughFall_Threshold), linewidth = 1, colour = "green") +
  labs(x = "Time", y = "Interception Store (mm)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters
         Shading = 10th and 90th Percentiles of Retained Parameters")
P7

# Simulated Interception Evaporation
P8<-ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = Ei_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= Ei_Q10, ymax= Ei_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Interception Evaporation (mm)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P8

# Simulated Transpiration
P9<-ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = Tr_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= Tr_Q10, ymax= Tr_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Transpiration (mm)") +
  scale_y_continuous(limits = c(0, 4)) +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P9

# Simulated Soil Evaporation
P10<-ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = Es_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= Es_Q10, ymax= Es_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Soil Evaporation (mm)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P10


# Recharge from Upper to Lower Store
P11<-ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = Recharge_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= Recharge_Q10, ymax= Recharge_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Recharge Upper to Lower Store (mm)") +
  scale_y_continuous(limits = c(0, 3)) +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters
          Shading = 10th and 90th Percentiles of Retained Parameters")
P11

# Simulated Flow out of Lower Store
P12<-ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = GWflow_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= GWflow_Q10, ymax= GWflow_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Flow out of lower store (mm)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters
          Shading = 10th and 90th Percentiles of Retained Parameters")
P12

# # Throughfall
P13<- ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = Th_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= Th_Q10, ymax= Th_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Throughfall (mm)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters
           Shading = 10th and 90th Percentiles of Retained Parameters")

P13

# Prefential flow
P14<- ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = Pref_Flow_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= Pref_Flow_Q10, ymax= Pref_Flow_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Preferential Flow (mm)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters
           Shading = 10th and 90th Percentiles of Retained Parameters")

P14

# Transpiratoin from different compartments
P15<- ggplot(data = full_timeseries) +
  geom_line(aes(x = Date, y = Tr_Upper_Q50), linewidth = 1, colour = "green") +
  geom_line(aes(x = Date, y = Tr_Lower_Q50), linewidth = 1, colour = "brown") +
  geom_line(aes(x = Date, y = Tr_Deep_Q50), linewidth = 1, colour = "blue") +
    labs(x = "Time", y = "% Trans Contribution   Higher Box Tr = Green, Lower Tr = Brown, \n Deeper Tr = Blue") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b")

P15

# Upper Storage volume D concentration
P16<-ggplot(data = full_timeseries[!is.na(full_timeseries$Upper_2H),]) +
  geom_line(aes(x = Date, y = upCSTO_D_Q50), linewidth = 1, colour = "black") +  
  geom_line(aes(x = Date, y = Upper_2H), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(x = Date, ymin= upCSTO_D_Q10, ymax= upCSTO_D_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Upper storage D concentration") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters, Blue Line = Measured 2H    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P16 

# Lower Storage volume D concentration
P17<-ggplot(data = full_timeseries[!is.na(full_timeseries$Medium_2H),]) +
  geom_line(aes(x = Date, y = gwCQ_D_Q50), linewidth = 1, colour = "black") +
  geom_line(aes(x = Date, y = Medium_2H), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(x = Date, ymin= gwCQ_D_Q10, ymax= gwCQ_D_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Lower D discharge concentration in permil") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P17 

# Deeper Storage volume D concentration
P18<-ggplot(data = full_timeseries[!is.na(full_timeseries$Lower_2H),]) +
  geom_line(aes(x = Date, y = lowCQ_D_Q50), linewidth = 1, colour = "black") +
  geom_line(aes(x = Date, y = Lower_2H), linewidth = 1, colour = "blue") +
  geom_ribbon(aes(x = Date, ymin= lowCQ_D_Q10, ymax= lowCQ_D_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Deeper D discharge concentration in permill") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P18 


#Overland flow
P20<-ggplot(data = full_timeseries[!is.na(full_timeseries$Qs_Q10),]) +
  geom_line(aes(x = Date, y = Qs_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= Qs_Q10, ymax= Qs_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Overland flow in mm") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P20 

#Vertical flux from upper into lower Storage
P21<-ggplot(data = full_timeseries[!is.na(full_timeseries$Perc_Q10),]) +
  geom_line(aes(x = Date, y = Perc_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= Perc_Q10, ymax= Perc_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Vertical flux from upper into lower in mm") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P21 

#Interception D isotopic composition
P22<-ggplot(data = full_timeseries[!is.na(full_timeseries$Upper_2H),]) +
  geom_line(aes(x = Date, y = Int_CD_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= Int_CD_Q10, ymax= Int_CD_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Interception D isotopic composition in permille of the residual liquid (? Int_CD)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P22 

#Interception vapor D isotopic composition
P23<-ggplot(data = full_timeseries[!is.na(full_timeseries$Upper_2H),]) +
  geom_line(aes(x = Date, y = fInt_CD_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= fInt_CD_Q10, ymax= fInt_CD_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Interception vapor D isotopic composition in permille") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P23 

#Upper storage vapor D composition
P24<-ggplot(data = full_timeseries[!is.na(full_timeseries$Upper_2H),]) +
  geom_line(aes(x = Date, y = fupCSTO_D_Q50), linewidth = 1, colour = "black") +
  geom_ribbon(aes(x = Date, ymin= fupCSTO_D_Q10, ymax= fupCSTO_D_Q90), linetype=2, alpha=0.4, fill="red") +
  labs(x = "Time", y = "Upper storage vapor D isotopic composition in permill (? fupCSTO_D)") +
  theme_classic(base_size = 12) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  ggtitle("", subtitle = "Black Line = Median of Retained Parameters    
          Shading = 10th and 90th Percentiles of Retained Parameters")
P24 


# Export if required
pdf("Model_Output_Plots_ForestA.pdf")
P1
P2 
P3 
P4 
P5 
P6 
P7 
P8 
P9 
P10 
P11 
P12 
P13
P14
P15
P16
P17
P18
P20
P21
P22
P23
P24
dev.off()
