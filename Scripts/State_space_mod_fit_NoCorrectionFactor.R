# State space population model for harbor seal stock assessment - revised as per reviewer comments 
# Written by Staci Amburgey, WDFW Research Scientist
rm(list=ls())

# LOAD libraries
require(gtools)
require(ggplot2)
require(dplyr)
require(readxl)
require(bayesplot)
# require(cmdstanr)
require(rstan)
require(posterior)
require(parallel)
require(fitdistrplus)
require(loo)
require(devtools)
require(lubridate)
# cmdstanr::set_cmdstan_path()
# options(cmdstanr_write_stan_file_dir = "C:/rtemp")

## Set working directly.----
# setwd("S:/WP/Science/Westside Team/Science Datasets/Pinnipeds/Aerial surveys/Analysis")

## Specify which years and haulouts to select from each stock.----
source("Scripts/data_selection.r")

## Load data.----
#Group three of the four stocks together as they have the same correction factor
alldata <- pv.df.seldatesyears %>%
  filter(case_when(
    #Washington Coast partial stock, limited to specific years by region
    Stock == "Coastal" & Year %in% CEYears & Year %in% OCYears ~ T,
    #Northern Inland stock, limited to specific years by region
    Stock == "Northern Inland" & Year %in% EBYears & Year %in% SJFYears & Year %in% SJIYears ~ T,
    #South Puget Sound stock
    Stock == "Southern Puget Sound" ~ T,
    # #Hood Canal stock
    # Stock == "Hood Canal" ~ T,
    T ~ F
  )) %>%
  mutate(Stock = if_else(Sitecode == 10, "Southern Puget Sound", Stock)) %>%     #10 which is unknown location in SPS was put in NI presumably by mistake - this sets to SPS (should be fixed in data now but still run just in case)
  filter(Sitecode != 10.89) %>%     #remove counts on log booms being transported through the area as not clear where they came from
  mutate(tidetype = 0, MFAug15 = 0, hour = 0, MinFHi.Cat = 0)

#Hood Canal has a separate correction factor
HCdata=droplevels(pv.df.seldatesyears[pv.df.seldatesyears$Stock=="Hood Canal",])
#add hi,lo,other factor variable for sites 
HCdata$tidetype=rep("Other",nrow(HCdata))
HCdata$tidetype[HCdata$Sitecode%in%c(8.04,8.05,8.06,8.08,8.13,8.14,8.15)]="High"
HCdata$tidetype[HCdata$Sitecode%in%c(8.16,8.19)]="Low"
HCdata$tidetype=factor(HCdata$tidetype)
#This is what was used in London et al although not strictly correct since 227 in non leap and 228 in leap years
#it is the decimal months from 15 August
HCdata$MFAug15 <- (HCdata$Julian- 226)/31
HCdata$hour = factor(hour(HCdata$Survey.time),levels=0:23)
#For non-high tide haulout locations use the peak value for high tide;
#hopefully this will minimize any positive bias in abundance estimate due to differences.
#Comment out the following line to examine the effect of this assumption
HCdata$time_from_high_minutes[HCdata$tidetype!="High"]=65
#Create high tide categories defined in London et al
HCdata$MinFHi.Cat=cut(HCdata$time_from_high_minutes,seq(-448-16,482+16,32),labels=seq(-448,482,32))

## Calculate correction factor for Washington Coast, Northern Inlands, and Southern Puget Sound.----
#Based on data from Huber et al. 2001, Table 2 detailing p_hat, number of animals, and number of surveys used to calculate correction factor (1.53)
Hubdat <- as.data.frame(matrix(c(c(0.70,0.62,0.69,0.54,0.73,0.66),
                                 c(34,17,18,16,13,26),
                                 c(4,2,2,3,3,5)),
                               nrow=6, ncol=3, dimnames=list(1:6,c("p","animals","nsurv"))))
#Calculate the variance ((p*[1-p])/n) and standard deviation (sqrt(var)) of p_hat from the Huber data
var_p <- (Hubdat[,"p"]*(1-Hubdat[,"p"]))/Hubdat[,"animals"]
sd_p <- sqrt(var_p)
#Calculate a weighted (by the number of surveys) mean of the variance and the standard deviation
var_pweighted <- weighted.mean(var_p,Hubdat[,"nsurv"]/max(Hubdat[,"nsurv"]))
sd_pweighted <- weighted.mean(sd_p,Hubdat[,"nsurv"]/max(Hubdat[,"nsurv"]))
#Calculate sigma_p ([(p*(1-p))/var]-1) from p_hat
# Hubdat$sigma_p <- ((Hubdat$p*(1-Hubdat$p))/var_p)-1
#Could bootstrap using these var and sd values to generate a distribution of vars and sds
# var_boot <- matrix(rbeta(6000,Hubdat$p*Hubdat$phi,(1-Hubdat$p)*Hubdat$sigma_p),6)  #could repeat per each number of surveys to get a weighted bootstrap std
# sd_boot <- sd(colMeans(pboot))

## Calculate correction factor for Hood Canal.----
# Load results from London et al to create survey specific correction factors
load("Data/all.fit.1.rda")
# get fixed effect table
fe=all.fit.1$fixed.effects
# create design matrix for survey data and add weighting of 1/3 for each of three years to use the mean value
dm=cbind(rep(1,nrow(HCdata)),model.matrix(~-1+hour,data=HCdata),model.matrix(~-1+MinFHi.Cat,data=HCdata),
         model.matrix(~-1+MFAug15,data=HCdata),model.matrix(~-1+hour:MFAug15,data=HCdata))
dm=cbind(dm[,1:56],matrix(1/3,nrow=nrow(dm),ncol=3),dm[,57:ncol(dm)])
#estimate of expected proportions hauled out at the covariate values when counts occurred
p=as.vector(plogis(dm%*%fe$estimate))
#variance-covariance matrix of proportion estimates
vc_p=(dm[,!is.na(fe$std.err)]*p*(1-p))%*%all.fit.1$covb%*%t(dm[,!is.na(fe$std.err)]*p*(1-p))
#multiplicative correction factor for counts
cf=1/p
#variance-covariance matrix of correction factors
vc_cf=cf^2%*%t(cf^2)*vc_p
## Pull predicted values from p and vc_p (diagonal is sd)
HCdata2 <- cbind(HCdata[,c("Sitecode","Stock","Year","Count.total")],p=p,sdp=sqrt(diag(vc_p)))

## Combine all datasets together.----
yall <- rbind(cbind(alldata[,c("Sitecode","Stock","Year","Count.total")],p=weighted.mean(Hubdat[,"p"],Hubdat[,"nsurv"]/max(Hubdat[,"nsurv"])),sdp=sd_pweighted),HCdata2)

## Calculate sigma_p (precision) in order to derive alpha and beta (shape parameters) to describe the beta distribution.----
yall$sigma_p <- ((yall$p*(1-yall$p))/(yall$sdp)^2)-1
yall$a <- yall$p*yall$sigma_p       #alpha
yall$b <- (1-yall$p)*yall$sigma_p   #beta

## Go from at least one or multiple counts at every haulout every year -> alpha and beta values for each stock each year.----
#First, take the mean count, alpha, and beta value at each haulout each year
yhaul <- aggregate(cbind(Count.total,a,b) ~ Stock + Year + Sitecode, data=yall, FUN = mean)
#Second, sum up counts across all haulouts in a stock over each year
ystock <- merge(yhaul,aggregate(Count.total ~ Stock + Year, data=yhaul, FUN = sum),by=c("Stock","Year"))
#variation across haulout counts each year, could consider using as random effect
# vary <-aggregate(Count.total ~ Stock + Sitecode + Year, data=yall, FUN = function(x) c(mean(x),sd(x)))
#Number of surveys of each haulout each year
# vary$nhaul <- rowSums(table(vary$Stock,vary$Sitecode) > 0)[vary$Stock]
# plot(vary$Count.total.y[,1],vary$Count.total.y[,2])   #comparison of stock level var vs. haulout level
#Third, sum up counts across the entire stock for each year; also generate summed V1 = a, V2 = b
yfinal <- aggregate(cbind(a*(Count.total.x/Count.total.y),b*(Count.total.x/Count.total.y)) ~ Stock + Year + Count.total.y, data=ystock, FUN = sum)
colnames(yfinal)[3:5] <- c("Count.total","alpha","beta")
## Visualizing - proportion of count at haulout to total at stock level
# ggplot(data=vary, aes(x=Count.total[,1], y=Count.total[,2],group=Stock, fill=Stock)) +
#   geom_point(pch=21) + geom_abline(aes(intercept=0,slope=1))

## Prepare data for analysis.----
# Mean count each year for each haulout in each stock
y <- yfinal %>%
  group_by(Stock) %>%
  mutate(occ = dense_rank(Year)) %>%
  ungroup() %>%
  mutate(Stock = factor(Stock, levels = c("Coastal","Hood Canal","Northern Inland","Southern Puget Sound"))) %>%
  arrange(Stock,Year) #%>%
  # filter(Stock != "Hood Canal") %>%
  # mutate(Stock = factor(Stock, levels = c("Coastal","Northern Inland","Southern Puget Sound")))

# Year before earliest survey
Year0 <- min(y$Year) - 1
# Obs = vector with stock survey counts (not adjusted for proportion not hauled out), 
Obs <- as.integer(y$Count.total)
#  yr = integer vector with year of each count, 
yr <- as.integer(y$Year - Year0)
#  st = integer vector with stock ID number of each count,
st <- as.integer(y$Stock)
# Nstk = number stocks,
Nstk <- max(st)
# Nyrs = number of years of dynamics (including years with no surveys),
Nyrs <- max(y$Year) - Year0
# Nsrv = total number surveys at the stock level (e.g. if 5 stocks were each surveyed in 5 different years, Nsrv=25)
Nsrv <- dim(y)[1] #length(table(y$Year))
# Get values for informative priors; however, since no longer correcting counts prior to the model, we may need to adjust these slightly
dfsum = y %>% group_by(Stock) %>%
  summarise(N0_pri = min(Count.total),
            N_max = max(Count.total))
N_min = dfsum$N0_pri
N_max = dfsum$N_max

## Set up model run.----
stan.data = list(Nstk=Nstk,
                 Nyrs=Nyrs,
                 Nsrv=Nsrv,
                 N_min=N_min,
                 N_max=N_max,
                 a=y$alpha,
                 b=y$beta,
                 yr=yr,
                 st=st,
                 Obs=Obs) #,
                 # cf=y$alpha/(y$alpha+y$beta))

# Fit model-------------------------------------------------
fitmodel = "Scripts/State_space_example_noarea_revised.stan"

parms = c("rmax","z","sigma_K","sigma_r","phi","K","MNPL","N","OSP")

init_fun = function(){
  list(rmax = runif(Nstk,0.2,0.3),
       z = runif(1,1,3),
       sigma_K = runif(1,0.25,1),
       sigma_r = runif(1,0,0.25),
       log_N0 = runif(Nstk,0.9 * log(dfsum$N0_pri),1.1 * log(dfsum$N0_pri) ),
       phi = runif(1,10,30) ) }

## try constant cf of a/(a+b) > check if beta should be narrower, check if cf is messing it up
## have beta parameters add higher values to keep constrained and not let too high or too low values
## remove Hood Canal and see if error still happens, or try phi three stocks and separate HC phi
## try phi cauchy with 1 or less sigma to keep large values from happening
## could try adding one to N * cf to see


#Rstan - trying these to deal with permission issues
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
fit <- stan(file = fitmodel, data = stan.data, chains=6, warmup = 1000, iter = 5000, pars = parms,  #chains = 6
            control=list(adapt_delta=0.99))  #, max_treedepth=12
stan_trace(fit, pars = parms[1])
saveRDS(fit, file="Results/SSM_revision.rds")
save(fit, file="Results/SSM_revision.Rdata")

# #subset chains that are stuck
# keepchains <- c(1:2,4:6)


#Cmdstanr - issues with writing files (permissions)
# nburnin = 2000#7000                    # number warm-up (burn-in) samples
# nsamples = 8000#50000                 # desired total number samples
# cores = detectCores()
# ncore = min(20,cores-2)
# Niter = round(nsamples/ncore)
# mod <- cmdstanr::cmdstan_model(fitmodel)   # compiles model (if necessary)
# suppressMessages(                # Suppress messages/warnings (if desired)
#   suppressWarnings (
#     fit <- mod$sample(
#       data = stan.data,
#       seed = 123,
#       chains = ncore,
#       parallel_chains = ncore,
#       refresh = 100,
#       init=init_fun,
#       iter_warmup = nburnin,
#       iter_sampling = Niter, 
#       adapt_delta = 0.99
#     )
#   )
# )
# Note: if error, run: fit$output(1)

# Save output
# saveRDS(fit, file = "Results/fit_SSM.rds")  #this will betray you and not save the tmp files needed to access results from model runs
# save(fit, file="Results/fit_SSM.Rdata")     #USE THIS INSTEAD

load("Results/SSM_revision.Rdata")


## CMDSTANR
# Create table of summary stats for params of interest:
# source("Scripts/cmdstan_sumstats.r")
# 
# 
# mcmc_trace(fit$draws("rmax"))
# mcmc_trace(fit$draws("z"))
# mcmc_trace(fit$draws("sigma_r"))
# mcmc_trace(fit$draws("sigma_K"))
# mcmc_trace(fit$draws("phi"))
# mcmc_trace(fit$draws("K"))
# 
# mcmc_areas(fit$draws(variables = c("K")), 
#            area_method="equal height",
#            prob = 0.8) + 
#   ggtitle(paste0("Posterior distributions, K")) +
#   theme_classic()
# 
# ifelse(sumstats$rhat > 1.1, print("Failed to converge"), print("Less than 1.1"))

#convert all iteration output to data frame
# df_res <- as_draws_df(fit)

# Plot expected vs observed

# N_exp = sumstats[which(startsWith(vn,"N[")),6]
# N_exp_lo = sumstats[which(startsWith(vn,"N[")),4]
# N_exp_hi = sumstats[which(startsWith(vn,"N[")),8]
# df_Nfit = data.frame(Stock = rep(seq(1,Nstk),Nyrs),
#                      Year = rep(seq(1,Nyrs),each = Nstk) + Year0,
#                      N_exp = N_exp,
#                      N_exp_lo = N_exp_lo,
#                      N_exp_hi = N_exp_hi)
# ii = which((N_exp_hi - N_exp_lo)>10*N_exp)  #no instances where this is true so omit this step
# df_Nfit = df_Nfit[-ii,]


# plt_trnd <- ggplot() +
#   geom_ribbon(data=df_Nfit,aes(ymin=N_exp_lo,ymax=N_exp_hi),alpha = 0.3) +
#   geom_line(data=df_Nfit,aes(x=Year,y=N_exp)) +
#   geom_point(data=y, aes(y=Abundance, x=Year)) +
#   theme_classic() +
#   facet_grid(vars(Stock), scales = "free")
# print(plt_trnd)

## RSTAN
source("Scripts/cmdstan_sumstats.r")

plot(fit, pars = "rmax", show_density = TRUE)
plot(fit, pars = "z", show_density = TRUE)
plot(fit, pars = "sigma_r", show_density = TRUE)
plot(fit, pars = "sigma_K", show_density = TRUE)
plot(fit, pars = "phi", show_density = TRUE)
plot(fit, pars = "K", show_density = TRUE)

#manually plot subset
# parint <- as.data.frame(subset_draws(as_draws(fit),chain=keepchains)[,,"rmax[1]"])
# 
# plot(parint[,1], type="l", col="red")
# lines(parint[,2], type="l", col="blue")
# lines(parint[,3], type="l", col="yellow")
# lines(parint[,4], type="l", col="green")
# lines(parint[,5], type="l", col="purple")

ifelse(sumstats$rhat > 1.1, print("Failed to converge"), print("Less than 1.1"))

#write summary of model results
write.csv(sumstats, "Results/HarborSeals_SSM_summary.csv")

## Updated plots by S. Amburgey (RSTAN)
start <- min(y$Year)
end <- max(y$Year)


Cyears <- subset(y, Stock == "Coastal")$Year
CSSM <- cbind(sumstats[which(startsWith(sumstats$variable,"N[1")),c(1:2,6,9)],as.data.frame(start:end))

png(filename="Results/Washington_Coast_SSM_rev.png", height = 8, width = 8, units = "in", res = 300)

ggplot() +
  geom_ribbon(data=CSSM, aes(ymin = `2.5%`, ymax = `97.5%`, x = `start:end`), fill = "grey", alpha = 0.5) +
  geom_line(data=CSSM, aes(x = `start:end`, y = mean), col="black") +
  geom_point(data = subset(CSSM, `start:end` %in% Cyears), aes(x = `start:end`, y = mean), pch=23, cex=4, fill="#06c5c9") +
  scale_x_continuous(limits = c(start,end)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  ggtitle("Washington Coast") +
  xlab("Year") + ylab("Abundance")

dev.off()


HCyears <- subset(y, Stock == "Hood Canal")$Year
HCSSM <- cbind(sumstats[which(startsWith(sumstats$variable,"N[2")),c(1:2,6,9)],as.data.frame(start:end))


png(filename="Results/Hood_Canal_Stock_SSM_rev.png", height = 8, width = 8, units = "in", res = 300)

ggplot() +
  annotate("rect", xmin = 1977, xmax = 2023, ymin = as.numeric(sumstats[which(startsWith(sumstats$variable,"K[2")),6]), ymax = as.numeric(sumstats[which(startsWith(sumstats$variable,"K[2")),9]), alpha = .1,fill = "#4F9573") +
  annotate("rect", xmin = 1977, xmax = 2023, ymin = as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[2")),6]), ymax = as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[2")),9]), alpha = .1,fill = "#bd5702") +
  geom_ribbon(data=HCSSM, aes(ymin = `2.5%`, ymax = `97.5%`, x = `start:end`), fill = "grey", alpha = 0.5) +
  geom_line(data=HCSSM, aes(x = `start:end`, y = mean), col="black") +
  geom_point(data = subset(HCSSM, `start:end` %in% HCyears), aes(x = `start:end`, y = mean), pch=23, cex=4, fill="#e82547") +
  geom_hline(yintercept=as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[2")),2]), col="#bd5702", lwd=0.65) +
  annotate(geom="text", x=1980, y=as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[2")),2]+200), label="MNPL", color="#bd5702") +
  geom_hline(yintercept=as.numeric(sumstats[which(startsWith(sumstats$variable,"K[2")),2]),col="#4F9573", lwd=0.65) +
  annotate(geom="text", x=1980, y=as.numeric(sumstats[which(startsWith(sumstats$variable,"K[2")),2]+200), label="K", color="#4F9573") +
  scale_x_continuous(limits = c(start,end)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  ggtitle("Hood Canal Stock") +
  xlab("Year") + ylab("Abundance")

dev.off()



NIyears <- subset(y, Stock == "Northern Inland")$Year
NISSM <- cbind(sumstats[which(startsWith(sumstats$variable,"N[3")),c(1:2,6,9)],as.data.frame(start:end))

png(filename="Results/Northern_Inland_Stock_SSM_rev.png", height = 8, width = 8, units = "in", res = 300)

ggplot() +
  annotate("rect", xmin = 1977, xmax = 2023, ymin = as.numeric(sumstats[which(startsWith(sumstats$variable,"K[3")),6]), ymax = as.numeric(sumstats[which(startsWith(sumstats$variable,"K[3")),9]), alpha = .1,fill = "#4F9573") +
  annotate("rect", xmin = 1977, xmax = 2023, ymin = as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[3")),6]), ymax = as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[3")),9]), alpha = .1,fill = "#bd5702") +
  geom_ribbon(data=NISSM, aes(ymin = `2.5%`, ymax = `97.5%`, x = `start:end`), fill = "grey", alpha = 0.5) +
  geom_line(data=NISSM, aes(x = `start:end`, y = mean), col="black") +
  geom_point(data = subset(NISSM, `start:end` %in% NIyears), aes(x = `start:end`, y = mean), pch=23, cex=4, fill="#af2bb5") +
  geom_hline(yintercept=as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[3")),2]), col="#bd5702", lwd=0.65) +
  annotate(geom="text", x=1980, y=as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[3")),2]+500), label="MNPL", color="#bd5702") +
  geom_hline(yintercept=as.numeric(sumstats[which(startsWith(sumstats$variable,"K[3")),2]),col="#4F9573", lwd=0.65) +
  annotate(geom="text", x=1980, y=as.numeric(sumstats[which(startsWith(sumstats$variable,"K[3")),2]+500), label="K", color="#4F9573") +
  scale_x_continuous(limits = c(start,end)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  ggtitle("Northern Inland Stock") +
  xlab("Year") + ylab("Abundance")

dev.off()



SPSyears <- subset(y, Stock == "Southern Puget Sound")$Year
SPSSSM <- cbind(sumstats[which(startsWith(sumstats$variable,"N[4")),c(1:2,6,9)],as.data.frame(start:end))

png(filename="Results/Southern_Puget_Sound_Stock_SSM_rev.png", height = 8, width = 8, units = "in", res = 300)

ggplot() +
  annotate("rect", xmin = 1977, xmax = 2023, ymin = as.numeric(sumstats[which(startsWith(sumstats$variable,"K[4")),6]), ymax = as.numeric(sumstats[which(startsWith(sumstats$variable,"K[4")),9]), alpha = .1,fill = "#4F9573") +
  annotate("rect", xmin = 1977, xmax = 2023, ymin = as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[4")),6]), ymax = as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[4")),9]), alpha = .1,fill = "#bd5702") +
  geom_ribbon(data=SPSSSM, aes(ymin = `2.5%`, ymax = `97.5%`, x = `start:end`), fill = "grey", alpha = 0.5) +
  geom_line(data=SPSSSM, aes(x = `start:end`, y = mean), col="black") +
  geom_point(data = subset(SPSSSM, `start:end` %in% SPSyears), aes(x = `start:end`, y = mean), pch=23, cex=4, fill="#028dd6") +
  geom_hline(yintercept=as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[4")),2]), col="#bd5702", lwd=0.65) +
  annotate(geom="text", x=1980, y=as.numeric(sumstats[which(startsWith(sumstats$variable,"MNPL[4")),2]+100), label="MNPL", color="#bd5702") +
  geom_hline(yintercept=as.numeric(sumstats[which(startsWith(sumstats$variable,"K[4")),2]),col="#4F9573", lwd=0.65) +
  annotate(geom="text", x=1980, y=as.numeric(sumstats[which(startsWith(sumstats$variable,"K[4")),2]+100), label="K", color="#4F9573") +
  scale_x_continuous(limits = c(start,end)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text = element_text(size=12), axis.title = element_text(size=14)) +
  ggtitle("Southern Puget Sound Stock") +
  xlab("Year") + ylab("Abundance")

dev.off()


# ## Updated plots by S. Amburgey (CMDSTANR)
# start <- min(y$Year)
# end <- max(y$Year)
# 
# 
# Cresults$Year <- as.numeric(as.character(Cresults$Year))
# 
# png(filename="Results/Washington_Coast_SSM.png", height = 8, width = 8, units = "in", res = 300)
# 
# ggplot() +
#   # annotate("rect", xmin = 1977, xmax = 2023, ymin = sumstats[which(startsWith(vn,"K[1")),4], ymax = sumstats[which(startsWith(vn,"K[1")),8], alpha = .1,fill = "#4F9573") +
#   # annotate("rect", xmin = 1977, xmax = 2023, ymin = sumstats[which(startsWith(vn,"MNPL[1")),4], ymax = sumstats[which(startsWith(vn,"MNPL[1")),8], alpha = .1,fill = "#bd5702") +
#   geom_ribbon(data=cbind(sumstats[which(startsWith(vn,"N[1")),c(4,8)],as.data.frame(start:end)), aes(ymin = `2.5%`, ymax = `97.5%`, x = `start:end`), fill = "grey", alpha = 0.5) +
#   geom_line(data=cbind(sumstats[which(startsWith(vn,"N[1")),1],as.data.frame(start:end)), aes(x = `start:end`, y = `sumstats[which(startsWith(vn, "N[1")), 1]`), col="black") +
#   geom_point(data=Cresults, aes(x=Year, y=Abundance), pch=21, cex=2.5, fill="#06c5c9") +
#   geom_errorbar(data=Cresults, aes(x=Year, ymin = LCL, ymax = UCL), width = 0) +
#   # geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[1")),1], col="#bd5702", lwd=0.65) +
#   # annotate(geom="text", x=1980, y=(sumstats[which(startsWith(vn,"MNPL[1")),1]+800), label="MNPL", color="#bd5702") +
#   # geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[1")),4],lty=3, col="#bd5702", lwd=0.65) +
#   # geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[1")),8],lty=3, col="#bd5702", lwd=0.65) +
#   # geom_hline(yintercept=sumstats[which(startsWith(vn,"K[1")),1],col="#4F9573", lwd=0.65) +
#   # annotate(geom="text", x=1980, y=(sumstats[which(startsWith(vn,"K[1")),1]+800), label="K", color="#4F9573") +
#   scale_x_continuous(limits = c(start,end)) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         axis.text = element_text(size=12), axis.title = element_text(size=14)) +
#   ggtitle("Washington Coast") +
#   xlab("Year") + ylab("Abundance")
# 
# dev.off()
# 
# 
# 
# HCresults$Year <- as.numeric(as.character(HCresults$Year))
# 
# png(filename="Results/Hood_Canal_Stock_SSM.png", height = 8, width = 8, units = "in", res = 300)
# 
# ggplot() +
#   annotate("rect", xmin = 1977, xmax = 2023, ymin = sumstats[which(startsWith(vn,"K[2")),4], ymax = sumstats[which(startsWith(vn,"K[2")),8], alpha = .1,fill = "#4F9573") +
#   annotate("rect", xmin = 1977, xmax = 2023, ymin = sumstats[which(startsWith(vn,"MNPL[2")),4], ymax = sumstats[which(startsWith(vn,"MNPL[2")),8], alpha = .1,fill = "#bd5702") +
#   geom_ribbon(data=cbind(sumstats[which(startsWith(vn,"N[2")),c(4,8)],as.data.frame(start:end)), aes(ymin = `2.5%`, ymax = `97.5%`, x = `start:end`), fill = "grey", alpha = 0.5) +
#   geom_line(data=cbind(sumstats[which(startsWith(vn,"N[2")),1],as.data.frame(start:end)), aes(x = `start:end`, y = `sumstats[which(startsWith(vn, "N[2")), 1]`), col="black") +
#   geom_point(data=HCresults, aes(x=Year, y=Abundance), pch=21, cex=2.5, fill="#e82547") +
#   geom_errorbar(data=HCresults, aes(x=Year, ymin = LCL, ymax = UCL), width = 0) +
#   geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[2")),1], col="#bd5702", lwd=0.65) +
#   annotate(geom="text", x=1980, y=(sumstats[which(startsWith(vn,"MNPL[2")),1]+200), label="MNPL", color="#bd5702") +
#   # geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[2")),4],lty=3, col="#bd5702", lwd=0.65) +
#   # geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[2")),8],lty=3, col="#bd5702", lwd=0.65) +
#   geom_hline(yintercept=sumstats[which(startsWith(vn,"K[2")),1],col="#4F9573", lwd=0.65) +
#   annotate(geom="text", x=1980, y=(sumstats[which(startsWith(vn,"K[2")),1]+200), label="K", color="#4F9573") +
#   scale_x_continuous(limits = c(start,end)) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         axis.text = element_text(size=12), axis.title = element_text(size=14)) +
#   ggtitle("Hood Canal Stock") +
#   xlab("Year") + ylab("Abundance")
# 
# dev.off()
# 
# 
# 
# NIresults$Year <- as.numeric(as.character(NIresults$Year))
# 
# png(filename="Results/Northern_Inland_Stock_SSM.png", height = 8, width = 8, units = "in", res = 300)
# 
# ggplot() +
#   annotate("rect", xmin = 1977, xmax = 2023, ymin = sumstats[which(startsWith(vn,"K[3")),4], ymax = sumstats[which(startsWith(vn,"K[3")),8], alpha = .1,fill = "#4F9573") +
#   annotate("rect", xmin = 1977, xmax = 2023, ymin = sumstats[which(startsWith(vn,"MNPL[3")),4], ymax = sumstats[which(startsWith(vn,"MNPL[3")),8], alpha = .1,fill = "#bd5702") +
#   geom_ribbon(data=cbind(sumstats[which(startsWith(vn,"N[3")),c(4,8)],as.data.frame(start:end)), aes(ymin = `2.5%`, ymax = `97.5%`, x = `start:end`), fill = "grey", alpha = 0.5) +
#   geom_line(data=cbind(sumstats[which(startsWith(vn,"N[3")),1],as.data.frame(start:end)), aes(x = `start:end`, y = `sumstats[which(startsWith(vn, "N[3")), 1]`), col="black") +
#   geom_point(data=NIresults, aes(x=Year, y=Abundance), pch=21, cex=2.5, fill="#af2bb5") +
#   geom_errorbar(data=NIresults, aes(x=Year, ymin = LCL, ymax = UCL), width = 0) +
#   geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[3")),1], col="#bd5702", lwd=0.65) +
#   annotate(geom="text", x=1980, y=(sumstats[which(startsWith(vn,"MNPL[3")),1]+500), label="MNPL", color="#bd5702") +
#   # geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[3")),4],lty=3, col="#bd5702", lwd=0.65) +
#   # geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[3")),8],lty=3, col="#bd5702", lwd=0.65) +
#   geom_hline(yintercept=sumstats[which(startsWith(vn,"K[3")),1],col="#4F9573", lwd=0.65) +
#   annotate(geom="text", x=1980, y=(sumstats[which(startsWith(vn,"K[3")),1]+500), label="K", color="#4F9573") +
#   scale_x_continuous(limits = c(start,end)) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         axis.text = element_text(size=12), axis.title = element_text(size=14)) +
#   ggtitle("Northern Inland Stock") +
#   xlab("Year") + ylab("Abundance")
# 
# dev.off()
# 
# 
# 
# SPSresults$Year <- as.numeric(as.character(SPSresults$Year))
# 
# png(filename="Results/Southern_Puget_Sound_Stock_SSM.png", height = 8, width = 8, units = "in", res = 300)
# 
# ggplot() +
#   annotate("rect", xmin = 1977, xmax = 2023, ymin = sumstats[which(startsWith(vn,"K[4")),4], ymax = sumstats[which(startsWith(vn,"K[4")),8], alpha = .1,fill = "#4F9573") +
#   annotate("rect", xmin = 1977, xmax = 2023, ymin = sumstats[which(startsWith(vn,"MNPL[4")),4], ymax = sumstats[which(startsWith(vn,"MNPL[4")),8], alpha = .1,fill = "#bd5702") +
#   geom_ribbon(data=cbind(sumstats[which(startsWith(vn,"N[4")),c(4,8)],as.data.frame(start:end)), aes(ymin = `2.5%`, ymax = `97.5%`, x = `start:end`), fill = "grey", alpha = 0.5) +
#   geom_line(data=cbind(sumstats[which(startsWith(vn,"N[4")),1],as.data.frame(start:end)), aes(x = `start:end`, y = `sumstats[which(startsWith(vn, "N[4")), 1]`), col="black") +
#   geom_point(data=SPSresults, aes(x=Year, y=Abundance), pch=21, cex=2.5, fill="#028dd6") +
#   geom_errorbar(data=SPSresults, aes(x=Year, ymin = LCL, ymax = UCL), width = 0) +
#   geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[4")),1], col="#bd5702", lwd=0.65) +
#   annotate(geom="text", x=1980, y=(sumstats[which(startsWith(vn,"MNPL[4")),1]+100), label="MNPL", color="#bd5702") +
#   # geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[4")),4],lty=3, col="#bd5702", lwd=0.65) +
#   # geom_hline(yintercept=sumstats[which(startsWith(vn,"MNPL[4")),8],lty=3, col="#bd5702", lwd=0.65) +
#   geom_hline(yintercept=sumstats[which(startsWith(vn,"K[4")),1],col="#4F9573", lwd=0.65) +
#   annotate(geom="text", x=1980, y=(sumstats[which(startsWith(vn,"K[4")),1]+100), label="K", color="#4F9573") +
#   scale_x_continuous(limits = c(start,end)) +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5, face = "bold"),
#         axis.text = element_text(size=12), axis.title = element_text(size=14)) +
#   ggtitle("Southern Puget Sound Stock") +
#   xlab("Year") + ylab("Abundance")
# 
# dev.off()

