## Code to reproduce the analyzes of the ms:
## Complex Patterns of Diversification in the Gray Zone of Speciation:
## Model-Based Approaches Applied to Patagonian Liolaemid Lizards (Squamata: Liolaemus kingii clade)

## by: KI Sánchez



setwd('/your/working/directory/')

# load('environment.Rdata')
# save.image('environment.Rdata')

# font utilities
# library(showtext)
# font_add("latin modern",
#          regular = "/media/kevin/KEVIN_HDD/software/linux/fonts/LMSans/lmsans-regular-webfont.ttf",
#          italic = "/media/kevin/KEVIN_HDD/software/linux/fonts/LMSans/lmsans-oblique-webfont.ttf",
#          bold = "/media/kevin/KEVIN_HDD/software/linux/fonts/LMSans/lmsans-bold-webfont.ttf",
#          bolditalic = "/media/kevin/KEVIN_HDD/software/linux/fonts/LMSans/lmsans-boldoblique-webfont.ttf")
# showtext_auto()
# par(family = "latin modern")

# GDI ----
# plots PDF 8w x 6h

# set color palette
archeforus <- "#bebada"
baguali <- "#fdb462"
escarchadosi <- "#ffed6f"
kingiiSS <- "#e31a1c"
scolaroi_zullyae <- '#bc80bd'
tari <- "#b3de69"
tristis <- '#fccde5'
kingii1 <- '#8dd3c7'
kingii2 <- '#fb8072'
deseado <- '#80b1d3'

# Generate GDI

# full

# import MCMC file
mcmcfile <- read.table("bpp/A00_MSC_gdi/100loci_mcmc.txt", header = TRUE)
# check stationarity lnL
plot(mcmcfile$Gen[seq(1, nrow(mcmcfile), 50)], mcmcfile$lnL[seq(1, nrow(mcmcfile), 50)], type = "l", xlab = "Generation", ylab = "lnL")

dens_tari <- density(1 - exp(-2*mcmcfile$tau_16/mcmcfile$theta_10))
dens_escarch <- density(1 - exp(-2*mcmcfile$tau_16/mcmcfile$theta_4))
dens_arch <- density(1 - exp(-2*mcmcfile$tau_19/mcmcfile$theta_1))
dens_sz <- density(1 - exp(-2*mcmcfile$tau_19/mcmcfile$theta_8))
dens_kiSS <- density(1 - exp(-2*mcmcfile$tau_21/mcmcfile$theta_7))
dens_ki2 <- density(1 - exp(-2*mcmcfile$tau_21/mcmcfile$theta_6))

# colaps1

# import MCMC file
mcmcfile <- read.table("bpp/A00_MSC_gdi/colaps1_mcmc.txt", header = TRUE)
# check stationarity lnL
plot(mcmcfile$Gen[seq(1, nrow(mcmcfile), 50)], mcmcfile$lnL[seq(1, nrow(mcmcfile), 50)], type = "l", xlab = "Generation", ylab = "lnL")

dens_baguali <- density(1 - exp(-2*mcmcfile$tau_12bagualiescarch_tari/mcmcfile$theta_2baguali))
dens_escarch_tari <- density(1 - exp(-2*mcmcfile$tau_12bagualiescarch_tari/mcmcfile$theta_4escarch_tari))
dens_tristis <- density(1 - exp(-2*mcmcfile$tau_14tristisarcheforus_sz/mcmcfile$theta_8tristis))
dens_arch_sz <- density(1 - exp(-2*mcmcfile$tau_14tristisarcheforus_sz/mcmcfile$theta_1archeforus_sz))
dens_deseado <- density(1 - exp(-2*mcmcfile$tau_15deseadokingii/mcmcfile$theta_3deseado))
dens_kingii <- density(1 - exp(-2*mcmcfile$tau_15deseadokingii/mcmcfile$theta_5kingii))

# colaps2

# import MCMC file
mcmcfile <- read.table("bpp/A00_MSC_gdi/colaps2_mcmc.txt", header = TRUE)
# check stationarity lnL
plot(mcmcfile$Gen[seq(1, nrow(mcmcfile), 50)], mcmcfile$lnL[seq(1, nrow(mcmcfile), 50)], type = "l", xlab = "Generation", ylab = "lnL")

dens_west <- density(1 - exp(-2*mcmcfile$tau_9westkingii_deseado/mcmcfile$theta_5west))
dens_kingii_deseado <- density(1 - exp(-2*mcmcfile$tau_9westkingii_deseado/mcmcfile$theta_1kingii_deseado))

# colaps3

# import MCMC file
mcmcfile <- read.table("bpp/A00_MSC_gdi/colaps3_mcmc.txt", header = TRUE)
# check stationarity lnL
plot(mcmcfile$Gen[seq(1, nrow(mcmcfile), 50)], mcmcfile$lnL[seq(1, nrow(mcmcfile), 50)], type = "l", xlab = "Generation", ylab = "lnL")

dens_south <- density(1 - exp(-2*mcmcfile$tau_7southkingii_deseado_west/mcmcfile$theta_3south))
dens_kingii_deseado_west <- density(1 - exp(-2*mcmcfile$tau_7southkingii_deseado_west/mcmcfile$theta_4kingii_deseado_west))

# colaps4

# import MCMC file
mcmcfile <- read.table("bpp/A00_MSC_gdi/colaps4_mcmc.txt", header = TRUE)
# check stationarity lnL
plot(mcmcfile$Gen[seq(1, nrow(mcmcfile), 50)], mcmcfile$lnL[seq(1, nrow(mcmcfile), 50)], type = "l", xlab = "Generation", ylab = "lnL")

dens_kingii1 <- density(1 - exp(-2*mcmcfile$tau_5kingii1colaps/mcmcfile$theta_2kingii1))
dens_colaps <- density(1 - exp(-2*mcmcfile$tau_5kingii1colaps/mcmcfile$theta_3colaps))

# Plots

# full
plot((dens_ki2), xlim = c(0, 1), col = kingii2, lwd = 4, xlab = "gdi", main = NA, cex.lab = 1.2, yaxt = "n", xaxt = 'n')
axis(1, at = c(0, 0.2, 0.7, 1), las = 2) # mostrar solo los valores que deseo
lines(dens_tari, lw = 4, col = tari)
lines(dens_escarch, lw = 4, col = escarchadosi)
lines(dens_arch, lw = 4, col = archeforus)
lines(dens_sz, lw = 4, col = scolaroi_zullyae)
lines(dens_kiSS, lw = 4, col = kingiiSS)
abline(v = 0.2, col = "black", lty = 2)
abline(v = 0.7, col = "black", lty = 2)
legend("topright", 10, legend = c("kingii 2",
                                  "tari",
                                  "escarchadosi",
                                  "archeforus",
                                  "scolaroi-zullyae",
                                  "kingii sensu stricto"),
       col = c(kingii2, tari, escarchadosi, archeforus, scolaroi_zullyae, kingiiSS), pch = 15, pt.cex = 2, bty = "n")

# colaps1
plot((dens_kingii), xlim = c(0, 1), col = kingii2, lwd = 4, xlab = "gdi", main = NA, cex.lab = 1.2, yaxt = "n", xaxt = 'n')
axis(1, at = c(0, 0.2, 0.7, 1), las = 2) # mostrar solo los valores que deseo
lines(dens_baguali, lw = 4, col = baguali)
lines(dens_deseado, lw = 4, col = deseado)
lines(dens_escarch_tari, lw = 4, col = escarchadosi)
lines(dens_arch_sz, lw = 4, col = archeforus)
lines(dens_tristis, lw = 4, col = tristis)
abline(v = 0.2, col = "black", lty = 2)
abline(v = 0.7, col = "black", lty = 2)
legend("topright", 10, legend = c("kingii sensu stricto + kingii 2",
                                  "baguali",
                                  "Deseado clade",
                                  "escarchadosi + tari",
                                  "archeforus + scolaroi-zullyae",
                                  "tristis"),
       col = c(kingii2, baguali, deseado, escarchadosi, archeforus, tristis), pch = 15, pt.cex = 2, bty = "n")

# colaps2
plot((dens_kingii_deseado), xlim = c(0, 1), col = deseado, lwd = 4, xlab = "gdi", main = NA, cex.lab = 1.2, yaxt = "n", xaxt = 'n')
axis(1, at = c(0, 0.2, 0.7, 1), las = 2) # mostrar solo los valores que deseo
lines(dens_west, lw = 4, col = archeforus)
abline(v = 0.2, col = "black", lty = 2)
abline(v = 0.7, col = "black", lty = 2)
legend("topright", 10, legend = c("kingii sensu stricto + kingii 2 + Deseado clade",
                                  "escarchadosi + tari + baguali"),
       col = c(deseado, archeforus), pch = 15, pt.cex = 2, bty = "n")

# colaps3
plot((dens_kingii_deseado_west), xlim = c(0, 1), col = deseado, lwd = 4, xlab = "gdi", main = NA, cex.lab = 1.2, yaxt = "n", xaxt = 'n')
axis(1, at = c(0, 0.2, 0.7, 1), las = 2) # mostrar solo los valores que deseo
lines(dens_south, lw = 4, col = baguali)
abline(v = 0.2, col = "black", lty = 2)
abline(v = 0.7, col = "black", lty = 2)
legend("topright", 10, legend = c("kingii sensu stricto + kingii 2 + Deseado clade + west",
                                  "escarchadosi + tari + baguali"),
       col = c(deseado, baguali), pch = 15, pt.cex = 2, bty = "n")

# colaps4
plot((dens_colaps), xlim = c(0, 1), col = deseado, lwd = 4, xlab = "gdi", main = NA, cex.lab = 1.2, yaxt = "n", xaxt = 'n')
axis(1, at = c(0, 0.2, 0.7, 1), las = 2) # mostrar solo los valores que deseo
lines(dens_kingii1, lw = 4, col = kingii1)
abline(v = 0.2, col = "black", lty = 2)
abline(v = 0.7, col = "black", lty = 2)
legend("topright", 10, legend = c("kingii 1",
                                  "colapsado"),
       col = c(kingii1, deseado), pch = 15, pt.cex = 2, bty = "n")

# 
# # plot mean, median and mode of gdi for A
# # I did not use this, but is useful
# n <- length(densityA$y)         
# dx <- mean(diff(densityA$x))    
# y.unit <- sum(densityA$y) * dx  
# dx <- dx / y.unit
# 
# x.mean <- sum(densityA$y * densityA$x) * dx
# y.mean <- densityA$y[length(densityA$x[densityA$x < x.mean])]
# x.mode <- densityA$x[i.mode <- which.max(densityA$y)]
# y.mode <- densityA$y[i.mode]
# y.cs <- cumsum(densityA$y)
# x.med <- densityA$x[i.med <- length(y.cs[2*y.cs <= y.cs[n]])]
# y.med <- densityA$y[i.med]
# temp <- mapply(function(x,y,c) lines(c(x,x), c(0,y), lwd=2, col=c), 
#                c(x.mean, x.med, x.mode), 
#                c(y.mean, y.med, y.mode), 
#                c("black", "green", "blue"))
# 
# legend("topright", 10,legend=c("mean", "mode", "median"),
#        col=c("red","blue","green"), pch=15, pt.cex=2, bty="n")
# 
# cat("Mean GDI for A is (red line): ", x.mean)
# cat("Median GDI for A is (green line): ", x.med)
# cat("Mode GDI for A is (blue line): ", x.mode)



# Posterior Predictive distributions of summary statistics of MSC parameters ----
# compare empirical distributions of MSC parameters with their posterior predictive distributions
library(e1071)

# load empirical mcmc file from bpp
mcmc_emp <- read.table('bpp/postPredSim_MSC/mcmc.txt', header = TRUE)
params <- colnames(mcmc_emp)

# set the number of simulations and the number of samples of each simulation
nsim = 100
nsam = 10000

# load posterior predictive simulations (mcmc files) as a list of dataframes
# here I omit the first file, it was loaded above
mcmc_pp_files <- list.files(pattern = "*sub.txt", recursive = TRUE)[-1]
mcmc_pp <- lapply(mcmc_pp_files, read.table, header = TRUE)

# this is the first simulation
mcmc_pp[[1]]

# check the match between column names (i.e. parameters)
sum(colnames(mcmc_emp) == colnames(mcmc_pp[[1]])) # 24 columns

# define function to calculate summary statistics for columns in dataframes
my_summary <- function(x){
  c(mean = mean(x),
    median = median(x),
    sd = sd(x),
    q025 = quantile(x, probs = 0.025),
    q25 = quantile(x, probs = 0.25),
    q75 = quantile(x, probs = 0.75),
    q975 = quantile(x, probs = 0.975),
    skew = skewness(x, type = 2),
    kurt = kurtosis(x, type = 2) )
  }
ls_stats <- lapply(mcmc_pp, sapply, my_summary)
rownames(ls_stats[[1]])
colnames(ls_stats[[1]])

par(mfrow = c(7, 3), mar = c(1, 1, 1, 1)) # plot PDF 5 x 8 in.

plot(density(unlist(lapply(ls_stats, function(x) x[1, 2]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$theta_1archef_sz), col = "red", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 15]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$theta_14tristisarchef_sz), col = "darkgreen", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 22]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$tau_14tristisarchef_sz), col = "blue", lwd = 3, lty = 2)

plot(density(unlist(lapply(ls_stats, function(x) x[1, 3]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n', xlim = c(0, mean(mcmc_emp$theta_2baguali))) 
abline(v = mean(mcmc_emp$theta_2baguali), col = "red", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 16]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$theta_15deseadoKingii), col = "darkgreen", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 23]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$tau_15deseadoKingii), col = "blue", lwd = 3, lty = 2)

plot(density(unlist(lapply(ls_stats, function(x) x[1, 4]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n', xlim = c(0, mean(mcmc_emp$theta_3deseado))) 
abline(v = mean(mcmc_emp$theta_3deseado), col = "red", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 14]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$theta_13tristisarchef_szdeseadoKingii), col = "darkgreen", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 21]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$tau_13tristisarchef_szdeseadoKingii), col = "blue", lwd = 3, lty = 2)

plot(density(unlist(lapply(ls_stats, function(x) x[1, 5]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$theta_4escarch_tari), col = "red", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 13]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$theta_12bagualiescarch_tari), col = "darkgreen", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 20]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$tau_12bagualiescarch_tari), col = "blue", lwd = 3, lty = 2)

plot(density(unlist(lapply(ls_stats, function(x) x[1, 6]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n', xlim = c(0, mean(mcmc_emp$theta_5Kingii))) 
abline(v = mean(mcmc_emp$theta_5Kingii), col = "red", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 12]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$theta_11bagualiescarch_taritristisarchef_szdeseadoKingii), col = "darkgreen", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 19]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$tau_11bagualiescarch_taritristisarchef_szdeseadoKingii), col = "blue", lwd = 3, lty = 2)

plot(density(unlist(lapply(ls_stats, function(x) x[1, 7]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n', xlim = c(0, mean(mcmc_emp$theta_6kingii1))) 
abline(v = mean(mcmc_emp$theta_6kingii1), col = "red", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 11]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$theta_10kingii1bagualiescarch_taritristisarchef_szdeseadoKingii), col = "darkgreen", lwd = 3, lty = 2)
plot(density(unlist(lapply(ls_stats, function(x) x[1, 18]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n') 
abline(v = mean(mcmc_emp$tau_10kingii1bagualiescarch_taritristisarchef_szdeseadoKingii), col = "blue", lwd = 3, lty = 2)

plot(density(unlist(lapply(ls_stats, function(x) x[1, 9]))), lwd = 2, main = NA, xlab = NA, ylab = NA, yaxt = 'n', xlim = c(0, mean(mcmc_emp$theta_8tristis))) 
abline(v = mean(mcmc_emp$theta_8tristis), col = "red", lwd = 3, lty = 2)



# Concordance Factors table for network estimation ----
# here I use the SNPs2CF() function

# source the function
source("/path/to/SNPs2CF/functions.R");

# Opcional: Convert phased VCF to phy
library(pegas)
inputVCF <- "/media/kevin/KEVIN_HDD/academico/papers_caps/kingiiClade_SDL/data/RADseq/matrices_review/clust90min37_46_subset/clust90_min37_46_subset_outfiles/clust90_min37_46_subset.usnps.vcf"
vcf <- read.vcf(inputVCF); vcf
vcf2phylip(vcf.name = inputVCF, total.SNPs = 1755, output.name = 'test.phy');

inputPHY <- "clust90_min37_46_subset_unphased.usnps.phy" # PUT THE .phy MATRIX IN YOUR WORKING DIRECTORY

# Generate CF table
# for Network estimation: n.quartets = 100, bootstrap = TRUE
# for goodness of fit tests: n.quartets = 1, bootstrap = FALSE
output <- SNPs2CF(seqMatrix = inputPHY,
                  ImapName = "imapfile.txt",
                  between.sp.only = TRUE,
                  max.SNPs = NULL,
                  bootstrap = FALSE,
                  outputName = "btw_sp_1quart_indels.csv",
                  save.progress = FALSE,
                  n.quartets = 1,
                  # rm.outgroup = TRUE,
                  # outgroupSp = "silvanae",
                  indels.as.fifth.state = FALSE,
                  boots.rep = 100,
                  starting.sp.quartet = 1,
                  max.quartets = 100000,
                  cores = 4)

# Broken-stick heuristics of pseudo-log-Lik ----
# I use this as first approach to select the optimal network

library(ggplot2)
data <- data.frame(y = c(42049.0535425142,
                         39091.782999292,
                         40878.299818464,
                         46789.3171745189),
                   x = 0:3)

ggplot(data, aes(x, y, group = 1)) +
  geom_point() +
  geom_line(size = 1) +
  theme_minimal(base_size = 22, base_family = 'LM Sans 10') +
  xlab('hmax') + ylab('-log-pseudolikelihood')



# Time calibration of MSC and MSCi models ----

library(bppr)
mcmc_msc <- read.table('bpp/A00_MSC_gdi/colaps1_mcmc.txt', header = TRUE)
mcmc_msci <- read.table('bpp/A00_MSCi/mcmc.txt', header = TRUE)

# calibrate using a prior on molecular substitution rate:
# this converts tau into absolute divergence times and theta into effective population sizes (N)

# I use the nuclear genome-wide substitution rate for Squamata
# taken from Perry et al. 2018 (DOI: 10.1093/gbe/evy157):
# 7.6e-4 substitutions/site/Myears -> transform this to per-generation mutation rate (1 gen ≈ 2 years)
subst_rate_gen <- 7.6e-4*2/1e6
cal_msc <- msc2time.r(mcmc_msc, u.mean = subst_rate_gen, u.sd = 1e-10, g.mean = 2, g.sd = 0.25)
cal_msci <- msc2time.r(mcmc_msci, u.mean = subst_rate_gen, u.sd = 1e-10, g.mean = 2, g.sd = 0.25)

# posterior means and 95 % HPD
means_msc <- as.data.frame(apply(cal_msc, 2, mean))
means_msci <- as.data.frame(apply(cal_msci, 2, mean))
postHPD_msc <- as.data.frame(coda::HPDinterval(coda::as.mcmc(cal_msc), ))
postHPD_msci <- as.data.frame(coda::HPDinterval(coda::as.mcmc(cal_msci), ))
summary_msc <- data.frame(parameter = rownames(means_msc),
                      lower95HPD = postHPD_msc$lower,
                      mean = means_msc[, 1],
                      upper95HPD = postHPD_msc$upper)
summary_msci <- data.frame(parameter = rownames(means_msci),
                          lower95HPD = postHPD_msci$lower,
                          mean = means_msci[, 1],
                          upper95HPD = postHPD_msci$upper)
write.table(summary_msc, 'bpp/calibration/summary_msc.csv', quote = FALSE, sep = ',', row.names = FALSE)
write.table(summary_msci, 'bpp/calibration/summary_msci.csv', quote = FALSE, sep = ',', row.names = FALSE)




# Demographic model comparison and parameter estimation with PipeMaster ----

library(PipeMaster); library(doMC); library(caret)

# Functions

# Plot prior distributions
plot.priors <- function(model, nsamples = 1000, mu.rates = NULL){
  param = NULL
  pb = txtProgressBar(min = 1, max = nsamples, initial = 1)
  for(j in 1:nsamples){
    param.samples <- as.numeric(PipeMaster:::msABC.commander(model, use.alpha = FALSE, arg = 1)[[2]][2, ])
    if(!(is.null(mu.rates))){
      rates <- do.call(mu.rates[[1]], args = c(1, mu.rates[2:length(mu.rates)]))
      param <- rbind(param, c(param.samples, rates))
    } else {
      rates <- PipeMaster:::sample.mu.rates(model)[[2]]
      param <- rbind(param, c(param.samples, rates))
    }
    setTxtProgressBar(pb, j)
  }
  if(!(is.null(mu.rates))){
    colnames(param) <- c(PipeMaster:::msABC.commander(model, use.alpha = FALSE, arg = 1)[[2]][1, ], "mu")
  } else {
    colnames(param) <- c(PipeMaster:::msABC.commander(model, use.alpha = FALSE, arg = 1)[[2]][1, ], "mean.mu", "sd.mu")
  }
  par(mfrow = c(3, 3))
  for(i in 1:ncol(param)){
    plot(density(param[, i]), main = colnames(param)[i], col = 2)
  }
}

# Plot simulated against observed data
plot.sim.obs <- function(sim, obs)
{
  mylabels <- colnames(sim)
  for (i in 1:ncol(sim)) {
    hist(sim[, i], breaks = 20, xlab = mylabels[i], main = "",
         xlim = c(min(c(na.omit(sim[, i]), obs[i])), max(c(na.omit(sim[, i]), obs[i]))))
    abline(v = obs[i], col = 2)
  }
}

# use this convert split a .alleles file from ipyrad into separate fastas
# iPyrad.alleles.loci2fasta('/media/kevin/KEVIN_HDD/academico/papers_caps/kingiiClade_SDL/data/RADseq/matrices_review/clust90min15_19_west/clust90min15_19_west_outfiles/clust90min15_19_west.alleles', '/media/kevin/KEVIN_HDD/academico/papers_caps/kingiiClade_SDL/data/RADseq/matrices_review/clust90min15_19_west/clust90min15_19_west_outfiles/fastas/')

# Notes:
#   Time is measured in generations -> bppr: years/gen_time
#   Pop. sizes measured as number of individuals (Ne); time increases backwards (e.g. Ne1 is ancestral to Ne0) -> can be obtained from bppr
#   Per generation mutation rate: hyperparameter for mean and SD of rates (for genomic data)
#   usar 'u' obtenido de bppr
#   Migration: 4Nm_(ij) -> specified in the backwards sense, M1_2 = migrants that receive 1 from 2

# In all cases I will contrast three models:
# Isolation only (Is)
# Isolation + migration (IM)
# Isolation + secondary contact (IMsc)

#### kingii vs. Deseado ####
# 1: kingii (sensu stricto + kingii2), 2: "Deseado" clade

setwd('kingii_deseado')
fastas <- 'fastas/directory/'
length(list.files(fastas, pattern = '*.fas'))
# 1879 loci

Is_kingii_des <- main.menu() # isolation
IM_kingii_des <- main.menu(Is_kingii_des) # IM
IMsc_kingii_des <- main.menu(IM_kingii_des) # IMsc

IMsc_kingii_des <- main.menu(IMsc_kingii_des) # edit with menu

# plot models (PDF 5 x 5 in.)
PlotModel(model = Is_kingii_des, use.alpha = TRUE, average.of.priors = FALSE)
PlotModel(model = IM_kingii_des, use.alpha = FALSE, average.of.priors = FALSE)
PlotModel(model = IMsc_kingii_des, use.alpha = FALSE, average.of.priors = FALSE)

# plot priors (PDF lscape 10 x 6 in.)
plot.priors(IMsc_kingii_des, nsamples = 1000)

# Replicate empirical data structure

popmap_kingii_des <- read.table("popmap.txt", header = TRUE, sep = ',') # assignments

Is_kingii_des <- get.data.structure(model = Is_kingii_des, path.to.fasta = fastas, pop.assign = popmap_kingii_des, sanger = FALSE)
IM_kingii_des <- get.data.structure(model = IM_kingii_des, path.to.fasta = fastas, pop.assign = popmap_kingii_des, sanger = FALSE)
IMsc_kingii_des <- get.data.structure(model = IMsc_kingii_des, path.to.fasta = fastas, pop.assign = popmap_kingii_des, sanger = FALSE)

# save models
dput(Is_kingii_des, "Is_kingii_des.txt")
dput(IM_kingii_des, "IM_kingii_des.txt")
dput(IMsc_kingii_des, "IMsc_kingii_des.txt")
# load models
Is_kingii_des <- dget("Is_kingii_des.txt") 
IM_kingii_des <- dget("IM_kingii_des.txt")
IMsc_kingii_des <- dget("IMsc_kingii_des.txt")

# Summary statistics calculation

obs_kingii_des <- obs.sumstat.ngs(model = Is_kingii_des, path.to.fasta = fastas, pop.assign = popmap_kingii_des)
# 'model' is only specified to build the sumstats names correctly
# check msABC manual for descriptions of summary statistics
exclude <- c(grep("thomson", colnames(obs_kingii_des)),
             grep("pairwise_fst", colnames(obs_kingii_des)),
             grep("Fay", colnames(obs_kingii_des)),
             grep("fwh", colnames(obs_kingii_des)),
             grep("_dv", colnames(obs_kingii_des)),
             grep("_s_", colnames(obs_kingii_des)),
             grep("_ZnS", colnames(obs_kingii_des))) # exclude some statistics
obs_kingii_des <- t(data.frame(obs_kingii_des[, -exclude]))
colnames(obs_kingii_des) # check retained statistics
write.table(obs_kingii_des, "obs_kingii_des.txt", quote = FALSE, col.names = TRUE, row.names = FALSE) # save to file

# Simulate data
# number of simulations = nsim.blocks x block.size x ncores -> 100K
# a block.size of 1000 will be good for most cases, only modify ncores and nsim.blocks

sim.msABC.sumstat(Is_kingii_des, nsim.blocks = 10, use.alpha = FALSE, output.name = "Is_kingii_des",
                  append.sims = FALSE, block.size = 1000, ncores = 10)
sim.msABC.sumstat(IM_kingii_des, nsim.blocks = 10, use.alpha = FALSE, output.name = "IM_kingii_des",
                  append.sims = FALSE, block.size = 1000, ncores = 10)
sim.msABC.sumstat(IMsc_kingii_des, nsim.blocks = 10, use.alpha = FALSE, output.name = "IMsc_kingii_des",
                  append.sims = FALSE, block.size = 1000, ncores = 10)

# Visualization

# plot simulations against empirical data
Is_kingii_des_sims <- read.table("SIMS_Is_kingii_des.txt", header = TRUE)
IM_kingii_des_sims <- read.table("SIMS_IM_kingii_des.txt", header = TRUE)
IMsc_kingii_des_sims <- read.table("SIMS_IMsc_kingii_des.txt", header = TRUE)

obs_kingii_des <- read.table("obs_kingii_des.txt", header = TRUE) # empirical data
# match the simulations sumstats to the observed
Is_kingii_des_sims <- Is_kingii_des_sims[, colnames(Is_kingii_des_sims) %in% colnames(obs_kingii_des)]
IM_kingii_des_sims <- IM_kingii_des_sims[, colnames(IM_kingii_des_sims) %in% colnames(obs_kingii_des)]
IMsc_kingii_des_sims <- IMsc_kingii_des_sims[, colnames(IMsc_kingii_des_sims) %in% colnames(obs_kingii_des)]
# plot
plot.sim.obs(IM_kingii_des_sims, obs_kingii_des)

# plot PCA of simulations against empirical data
models <- rbind(Is_kingii_des_sims, IM_kingii_des_sims, IMsc_kingii_des_sims)
data <- c(rep("Is", nrow(Is_kingii_des_sims)), rep("IM", nrow(IM_kingii_des_sims)), rep("IMsc", nrow(IMsc_kingii_des_sims)))
plotPCs(model = models, index = data, observed = obs_kingii_des, subsample = 1) # PNG 1000w x 800h in.

# Model classification -> SML

# set up number of cores for SML
registerDoMC(4)

# combine simulations and index
models <- cbind(models, data)

# setup the outcome (name of the models, categories)
outcomeName <- 'data'

# set up predictors (summary statistics)
predictorsNames <- names(models)[names(models) != outcomeName]

# split the data into training and testing sets (75 % for trainning, 25 % testing)
splitIndex <- createDataPartition(models[, outcomeName], p = 0.75, list = FALSE, times = 1)
train <- models[splitIndex, ]
test  <- models[-splitIndex, ]

# bootstraps and other controls
objControl <- trainControl(method = 'boot', number = 10, returnResamp = 'final', classProbs = TRUE)
# train the algorithm
nnetModel_select <- train(train[, predictorsNames], train[, outcomeName],
                          method = "nnet",
                          maxit = 5000,
                          trControl = objControl,
                          metric = "Accuracy",
                          preProc = c("center", "scale"))

# predict outcome for testing data, classification
predictions <- predict(object = nnetModel_select, test[, predictorsNames], type = 'raw')

# calculate accuracy in model classification
accu <- postResample(pred = predictions, obs = as.factor(test[, outcomeName]))

# predict probabilities of the models for the observed data
pred <- predict(object = nnetModel_select, observed, type = 'prob')

# visualize results
t(c(pred, accu))

# write results to file
write.table(c(pred, accu), "model_selection.txt")

# Model classification -> ABC
# you can use this to contrast with SML

# prob <- postpr(target = obs_kingii_des, sumstat = models, index = data, method = "rejection", tol = 0.1)
# summary(prob)
# 
# # cross-validation
# cv <- cv4postpr(sumstat = models, index = data, method = "rejection", tol = 0.1, nval = 20)
# summary(cv)
# plot(cv)
# 
# # overall accuracy
# acc <- summary(cv)
# sum(diag(acc$conf.matrix$tol0.1))/60

# Parameter estimation -> ABC + NNet

# read selected model
Is_kingii_des_sims <- read.table("SIMS_Is_kingii_des.txt", header = TRUE)

# separate summary statistics from parameters
sims <- Is_kingii_des_sims[, colnames(Is_kingii_des_sims) %in% colnames(obs_kingii_des)]
colnames(Is_kingii_des_sims) # Ne0.pop1, Ne0.pop2, join1_2, mean.rate, sd.rate
param <- Is_kingii_des_sims[, 1:5]

# estimate posterior distribution of parameters

post <- abc(target = obs_kingii_des,
            param = param,
            sumstat = sims,
            sizenet = 20,
            method = "neuralnet",
            MaxNWts = 5000,
            tol = 0.1) 
# adjust tolerance level according to the total number of simulations.

# write results to file
write.table(summary(post), "parameters_Is.txt")

# plot posterior probabilities against prior
pdf('plot_posteriors.pdf')
extrafont::loadfonts()
par(mfrow = c(2, 3), family = "latin modern")
for(i in 1:ncol(param)){
  plot(density(post$unadj.values[, i]), col = 2, main = colnames(param)[i])
  lines(density(param[, i]))
}
dev.off()

# cross-validation for parameter estimates
cv <- cv4abc(param = param,
             sumstat = sims,
             nval = 20,
             sizenet = 20,
             method = "neuralnet",
             MaxNWts = 5000,
             tol = 0.1)
plot(cv) # pdf 5 x 5 
write.table(summary(cv), "cv_error.txt")


#### Southern clade ####
# 1: baguali, 2: escarchadosi_tari

setwd('south/')
fastas <- 'fastas/directory/'
length(list.files(fastas, pattern = '*.fas'))
# 2107 loci (I removed one locus for not having variable sites)

# Create models

Is_south <- main.menu(Is_kingii_des) # isolation
IM_south <- main.menu(IM_kingii_des) # IM
IMsc_south <- main.menu(IMsc_kingii_des) # IMsc

# IMsc_south <- main.menu(IMsc_south) # edit with menu

# plot models (PDF 5 x 5 in.)
PlotModel(model = Is_south, use.alpha = TRUE, average.of.priors = FALSE)
PlotModel(model = IM_south, use.alpha = FALSE, average.of.priors = FALSE)
PlotModel(model = IMsc_south, use.alpha = FALSE, average.of.priors = FALSE)

# plot priors (PDF lscape 10 x 6 in.)
plot.priors(IMsc_south, nsamples = 1000)

# Replicate empirical data structure

popmap_south <- read.table("popmap.txt", header = TRUE, sep = ',') # assignments

Is_south <- get.data.structure(model = Is_south, path.to.fasta = fastas, pop.assign = popmap_south, sanger = FALSE)
IM_south <- get.data.structure(model = IM_south, path.to.fasta = fastas, pop.assign = popmap_south, sanger = FALSE)
IMsc_south <- get.data.structure(model = IMsc_south, path.to.fasta = fastas, pop.assign = popmap_south, sanger = FALSE)

# save models
dput(Is_south, "Is_south.txt")
dput(IM_south, "IM_south.txt")
dput(IMsc_south, "IMsc_south.txt")
# load models
Is_south <- dget("Is_south.txt")
IM_south <- dget("IM_south.txt")
IMsc_south <- dget("IMsc_south.txt")

# Summary statistics calculation

obs_south <- obs.sumstat.ngs(model = Is_south, path.to.fasta = fastas, pop.assign = popmap_south)
# 'model' is only specified to build the sumstats names correctly
# check msABC manual for descriptions of summary statistics
exclude <- c(grep("thomson", colnames(obs_south)),
             grep("pairwise_fst", colnames(obs_south)),
             grep("Fay", colnames(obs_south)),
             grep("fwh", colnames(obs_south)),
             grep("_dv", colnames(obs_south)),
             grep("_s_", colnames(obs_south)),
             grep("_ZnS", colnames(obs_south))) # exclude some statistics
obs_south <- t(data.frame(obs_south[, -exclude]))
colnames(obs_south) # check retained statistics
write.table(obs_south, "obs_south.txt", quote = FALSE, col.names = TRUE, row.names = FALSE) # save to file

# Simulate data
# number of simulations = nsim.blocks x block.size x ncores -> 100K
# a block.size of 1000 will be good for most cases, only modify ncores and nsim.blocks

sim.msABC.sumstat(Is_south, nsim.blocks = 10, use.alpha = FALSE, output.name = "Is_south",
                  append.sims = FALSE, block.size = 1000, ncores = 10)
sim.msABC.sumstat(IM_south, nsim.blocks = 10, use.alpha = FALSE, output.name = "IM_south",
                  append.sims = FALSE, block.size = 1000, ncores = 10)
sim.msABC.sumstat(IMsc_south, nsim.blocks = 10, use.alpha = FALSE, output.name = "IMsc_south",
                  append.sims = FALSE, block.size = 1000, ncores = 10)

# Visualization

# plot simulations against empirical data
Is_south_sims <- read.table("SIMS_Is_south.txt", header = TRUE)
IM_south_sims <- read.table("SIMS_IM_south.txt", header = TRUE)
IMsc_south_sims <- read.table("SIMS_IMsc_south.txt", header = TRUE)

obs_south <- read.table("obs_south.txt", header = TRUE) # empirical data
# match the simulations sumstats to the observed
Is_south_sims <- Is_south_sims[, colnames(Is_south_sims) %in% colnames(obs_south)]
IM_south_sims <- IM_south_sims[, colnames(IM_south_sims) %in% colnames(obs_south)]
IMsc_south_sims <- IMsc_south_sims[, colnames(IMsc_south_sims) %in% colnames(obs_south)]
# plot
plot.sim.obs(Is_south_sims, obs_south)

# plot PCA of simulations against empirical data
models <- rbind(Is_south_sims, IM_south_sims, IMsc_south_sims)
data <- c(rep("Is", nrow(Is_south_sims)), rep("IM", nrow(IM_south_sims)), rep("IMsc", nrow(IMsc_south_sims)))
plotPCs(model = models, index = data, observed = obs_south, subsample = 1) # PNG 1000w x 800h in.

# Model classification -> SML

# set up number of cores for SML
registerDoMC(2)

# combine simulations and index
models <- cbind(models, data)

# setup the outcome (name of the models, categories)
outcomeName <- 'data'

# set up predictors (summary statistics)
predictorsNames <- names(models)[names(models) != outcomeName]

# split the data into training and testing sets (75 % for trainning, 25 % testing)
splitIndex <- createDataPartition(models[, outcomeName], p = 0.75, list = FALSE, times = 1)
train <- models[splitIndex, ]
test  <- models[-splitIndex, ]

# bootstraps and other controls
objControl <- trainControl(method = 'boot', number = 10, returnResamp = 'final', classProbs = TRUE)
# train the algorithm
nnetModel_select <- train(train[, predictorsNames], train[, outcomeName],
                          method = "nnet",
                          maxit = 5000,
                          trControl = objControl,
                          metric = "Accuracy",
                          preProc = c("center", "scale"))

# predict outcome for testing data, classification
predictions <- predict(object = nnetModel_select, test[, predictorsNames], type = 'raw')

# calculate accuracy in model classification
accu <- postResample(pred = predictions, obs = as.factor(test[, outcomeName]))

# predict probabilities of the models for the observed data
pred <- predict(object = nnetModel_select, observed, type = 'prob')

# visualize results
t(c(pred, accu))

# write results to file
write.table(c(pred, accu), "model_selection.txt")

# Parameter estimation -> ABC + NNet

# read selected model
IMsc_south_sims <- read.table("SIMS_IMsc_south.txt", header = TRUE)

# separate summary statistics from parameters
sims <- IMsc_south_sims[, colnames(IMsc_south_sims) %in% colnames(obs_south)]
colnames(IMsc_south_sims) # Ne0.pop1, Ne0.pop2, join1_2,
                          # t.mig1.1_2, t.mig1.2_1, mig0.1_2, mig0.2_1, mig1.1_2, mig1.2_1,
                          # mean.rate, sd.rate
param <- IMsc_south_sims[, 1:11]

# estimate posterior distribution of parameters

post <- abc(target = obs_south,
            param = param,
            sumstat = sims,
            sizenet = 20,
            method = "neuralnet",
            MaxNWts = 5000,
            tol = 0.1) 
# adjust tolerance level according to the total number of simulations.

# write results to file
write.table(summary(post), "parameters_IMsc.txt")

# plot posterior probabilities against prior
pdf('plot_posteriors.pdf')
extrafont::loadfonts()
par(mfrow = c(3, 4), family = "latin modern")
for(i in 1:ncol(param)){
  plot(density(post$unadj.values[, i]), col = 2, main = colnames(param)[i])
  lines(density(param[, i]))
}
dev.off()

# cross-validation for parameter estimates
cv <- cv4abc(param = param,
             sumstat = sims,
             nval = 20,
             sizenet = 20,
             method = "neuralnet",
             MaxNWts = 5000,
             tol = 0.1)
plot(cv) # pdf 5 x 5 
write.table(summary(cv), "cv_error.txt")


#### Western clade ####
# 1: archeforus_sz, 2: tristis

setwd('west/')
fastas <- 'directory/fastas'
length(list.files(fastas, pattern = '*.fas'))
# 2522 loci

Is_west <- main.menu(Is_south) # IS
IM_west <- main.menu(IM_south) # IM
IMsc_west <- main.menu(IMsc_south) # IMsc

IMsc_west <- main.menu(IMsc_west) # edit with menu

# plot models (PDF 5 x 5 in.)
PlotModel(model = Is_west, use.alpha = TRUE, average.of.priors = FALSE)
PlotModel(model = IM_west, use.alpha = FALSE, average.of.priors = FALSE)
PlotModel(model = IMsc_west, use.alpha = FALSE, average.of.priors = FALSE)

# plot priors (PDF lscape 10 x 6 in.)
plot.priors(IMsc_west, nsamples = 1000)

# Replicate empirical data structure

popmap_west <- read.table("popmap.txt", header = TRUE, sep = ',') # assignments

Is_west <- get.data.structure(model = Is_west, path.to.fasta = fastas, pop.assign = popmap_west, sanger = FALSE)
IM_west <- get.data.structure(model = IM_west, path.to.fasta = fastas, pop.assign = popmap_west, sanger = FALSE)
IMsc_west <- get.data.structure(model = IMsc_west, path.to.fasta = fastas, pop.assign = popmap_west, sanger = FALSE)

# save models
dput(Is_west, "Is_west.txt")
dput(IM_west, "IM_west.txt")
dput(IMsc_west, "IMsc_west.txt")
# load models
Is_west <- dget("Is_west.txt")
IM_west <- dget("IM_west.txt")
IMsc_west <- dget("IMsc_west.txt")

# Summary statistics calculation

obs_west2 <- obs.sumstat.ngs(model = IMsc_west, path.to.fasta = fastas, pop.assign = popmap_west)
# 'model' is only specified to build the sumstats names correctly
# check msABC manual for descriptions of summary statistics
exclude <- c(grep("thomson", colnames(obs_west)),
             grep("pairwise_fst", colnames(obs_west)),
             grep("Fay", colnames(obs_west)),
             grep("fwh", colnames(obs_west)),
             grep("_dv", colnames(obs_west)),
             grep("_s_", colnames(obs_west)),
             grep("_ZnS", colnames(obs_west))) # exclude some statistics
obs_west <- t(data.frame(obs_west[, -exclude]))
colnames(obs_west) # check retained statistics
write.table(obs_west, "obs_west.txt", quote = FALSE, col.names = TRUE, row.names = FALSE) # save to file

# Simulate data
# number of simulations = nsim.blocks x block.size x ncores -> 100K
# a block.size of 1000 will be good for most cases, only modify ncores and nsim.blocks

sim.msABC.sumstat(Is_west, nsim.blocks = 10, use.alpha = FALSE, output.name = "Is_west",
                  append.sims = FALSE, block.size = 1000, ncores = 10)
sim.msABC.sumstat(IM_west, nsim.blocks = 10, use.alpha = FALSE, output.name = "IM_west",
                  append.sims = FALSE, block.size = 1000, ncores = 10)
sim.msABC.sumstat(IMsc_west, nsim.blocks = 10, use.alpha = FALSE, output.name = "IMsc_west",
                  append.sims = FALSE, block.size = 1000, ncores = 10)

# Visualization

# plot simulations against empirical data
Is_west_sims <- read.table("SIMS_Is_west.txt", header = TRUE)
IM_west_sims <- read.table("SIMS_IM_west.txt", header = TRUE)
IMsc_west_sims <- read.table("SIMS_IMsc_west.txt", header = TRUE)

obs_west <- read.table("obs_west.txt", header = TRUE) # empirical data
# match the simulations sumstats to the observed
Is_west_sims <- Is_west_sims[, colnames(Is_west_sims) %in% colnames(obs_west)]
IM_west_sims <- IM_west_sims[, colnames(IM_west_sims) %in% colnames(obs_west)]
IMsc_west_sims <- IMsc_west_sims[, colnames(IMsc_west_sims) %in% colnames(obs_west)]
# plot
plot.sim.obs(Is_west_sims, obs_west)

# plot PCA of simulations against empirical data
models <- rbind(Is_west_sims, IM_west_sims, IMsc_west_sims)
data <- c(rep("Is", nrow(Is_west_sims)), rep("IM", nrow(IM_west_sims)), rep("IMsc", nrow(IMsc_west_sims)))
plotPCs(model = models, index = data, observed = obs_west, subsample = 1) # PNG 1000w x 800h in.

# Model classification -> SML

# set up number of cores for SML
registerDoMC(1)

# combine simulations and index
models <- cbind(models, data)

# setup the outcome (name of the models, categories)
outcomeName <- 'data'

# set up predictors (summary statistics)
predictorsNames <- names(models)[names(models) != outcomeName]

# split the data into training and testing sets (75 % for trainning, 25 % testing)
splitIndex <- createDataPartition(models[, outcomeName], p = 0.75, list = FALSE, times = 1)
train <- models[splitIndex, ]
test <- models[-splitIndex, ]

# bootstraps and other controls
objControl <- trainControl(method = 'boot', number = 10, returnResamp = 'final', classProbs = TRUE)
# train the algorithm
nnetModel_select <- train(train[, predictorsNames], train[, outcomeName],
                          method = "nnet",
                          maxit = 5000,
                          trControl = objControl,
                          metric = "Accuracy",
                          preProc = c("center", "scale"))

# predict outcome for testing data, classification
predictions <- predict(object = nnetModel_select, test[, predictorsNames], type = 'raw')

# calculate accuracy in model classification
accu <- postResample(pred = predictions, obs = as.factor(test[, outcomeName]))

# predict probabilities of the models for the observed data
pred <- predict(object = nnetModel_select, observed, type = 'prob')

# visualize results
t(c(pred, accu))

# write results to file
write.table(c(pred, accu), "model_selection.txt", quote = FALSE)

# Parameter estimation -> ABC + NNet

# read selected model
IMsc_west_sims <- read.table("SIMS_IMsc_west.txt", header = TRUE)

# separate summary statistics from parameters
sims <- IMsc_south_sims[, colnames(IMsc_south_sims) %in% colnames(obs_south)]
colnames(IMsc_south_sims) # Ne0.pop1, Ne0.pop2, join1_2,
                          # t.mig1.1_2, t.mig1.2_1, mig0.1_2, mig0.2_1, mig1.1_2, mig1.2_1,
                          # mean.rate, sd.rate
param <- IMsc_south_sims[, 1:11]

# estimate posterior distribution of parameters

post <- abc(target = obs_west,
            param = param,
            sumstat = sims,
            sizenet = 20,
            method = "neuralnet",
            MaxNWts = 5000,
            tol = 0.1) 
# adjust tolerance level according to the total number of simulations.

# write results to file
write.table(summary(post), "parameters_IMsc.txt")

# plot posterior probabilities against prior
pdf('plot_posteriors.pdf')
extrafont::loadfonts()
par(mfrow = c(3, 4), family = "latin modern")
for(i in 1:ncol(param)){
  plot(density(post$unadj.values[, i]), col = 2, main = colnames(param)[i])
  lines(density(param[, i]))
}
dev.off()

# cross-validation for parameter estimates
cv <- cv4abc(param = param,
             sumstat = sims,
             nval = 20,
             sizenet = 20,
             method = "neuralnet",
             MaxNWts = 5000,
             tol = 0.1)
plot(cv) # pdf 5 x 5 
write.table(summary(cv), "cv_error.txt")




sessionInfo()

# R version 4.2.2 Patched (2022-11-10 r83330)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Linux Mint 20.1
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=es_AR.UTF-8   
# [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=es_AR.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=es_AR.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] caret_6.0-93     lattice_0.20-45  PipeMaster_0.0.9 msm_1.7          abc_2.2.1        locfit_1.5-9.6   MASS_7.3-58.2    quantreg_5.94   
# [9] SparseM_1.81     nnet_7.3-18      abc.data_1.0     PopGenome_2.7.5  ff_4.0.7         bit_4.0.4        phyclust_0.1-32  bppr_0.6.2      
# [17] ggplot2_3.4.0    pegas_1.1        ape_5.6-2        doMC_1.3.8       iterators_1.0.14 foreach_1.5.2    e1071_1.7-12    
# 
# loaded via a namespace (and not attached):
# [1] fastmatch_1.1-3         PBSmapping_2.73.2       RcppEigen_0.3.3.9.2     plyr_1.8.7              igraph_1.3.1           
# [6] sp_1.5-0                splines_4.2.2           listenv_0.8.0           usethis_2.1.6           digest_0.6.30          
# [11] fansi_1.0.3             magrittr_2.0.3          phytools_1.0-3          memoise_2.0.1           remotes_2.4.2          
# [16] recipes_1.0.1           globals_0.16.1          gower_1.0.0             extrafont_0.18          extrafontdb_1.0        
# [21] hardhat_1.2.0           prettyunits_1.1.1       colorspace_2.0-3        dplyr_1.0.9             callr_3.7.1            
# [26] crayon_1.5.1            survival_3.4-0          phangorn_2.8.1          glue_1.6.2              gtable_0.3.1           
# [31] ipred_0.9-13            MatrixModels_0.5-1      pkgbuild_1.3.1          Rttf2pt1_1.3.10         future.apply_1.9.1     
# [36] maps_3.4.0              scales_1.2.1            mvtnorm_1.1-3           DBI_1.1.2               Rcpp_1.0.9             
# [41] plotrix_3.8-2           tmvnsim_1.0-2           units_0.8-0             proxy_0.4-27            stats4_4.2.2           
# [46] lava_1.6.10             prodlim_2019.11.13      RColorBrewer_1.1-3      ellipsis_0.3.2          pkgconfig_2.0.3        
# [51] utf8_1.2.2              tidyselect_1.1.2        rlang_1.0.6             ggspatial_1.1.6         reshape2_1.4.4         
# [56] munsell_0.5.0           tools_4.2.2             cachem_1.0.6            cli_3.4.1               generics_0.1.2         
# [61] devtools_2.4.3          stringr_1.5.0           fastmap_1.1.0           ModelMetrics_1.2.2.2    processx_3.7.0         
# [66] fs_1.5.2                tess3r_1.1.0            purrr_0.3.4             future_1.28.0           nlme_3.1-162           
# [71] brio_1.1.3              compiler_4.2.2          rstudioapi_0.14         testthat_3.1.4          clusterGeneration_1.3.7
# [76] tibble_3.1.8            stringi_1.7.8           ps_1.7.1                desc_1.4.1              rgeos_0.5-9            
# [81] Matrix_1.5-1            classInt_0.4-7          vctrs_0.5.1             pillar_1.8.1            lifecycle_1.0.3        
# [86] combinat_0.0-8          data.table_1.14.4       raster_3.5-21           R6_2.5.1                KernSmooth_2.23-20     
# [91] parallelly_1.32.1       sessioninfo_1.2.2       codetools_0.2-19        assertthat_0.2.1        pkgload_1.2.4          
# [96] rprojroot_2.0.3         withr_2.5.0             mnormt_2.0.2            expm_0.999-6            terra_1.5-34           
# [101] quadprog_1.5-8          grid_4.2.2              rpart_4.1.19            timeDate_4021.106       coda_0.19-4            
# [106] class_7.3-21            ggnewscale_0.4.7        sf_1.0-8                pROC_1.18.0             numDeriv_2016.8-1.1    
# [111] scatterplot3d_0.3-41    rEEMSplots_0.0.1        lubridate_1.8.0