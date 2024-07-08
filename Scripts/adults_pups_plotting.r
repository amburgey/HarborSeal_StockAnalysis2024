#Selection of years and dates for surveys

rm(list=ls())
setwd("S:/WP/Science/Westside Team/Science Datasets/Pinnipeds/Aerial surveys/Analysis")
library(lubridate);library(tidyverse);library(ggstream);library(ggpubr);library(ggplot2)

## For figure without correction factor
source("Scripts/data_selection.r")

## Run previous data prep and adjustment scripts for data analysis
#Create combined dataset for all stocks
source("Scripts/HC Estimate.r")
allData <- HCdata[,c("Sitecode","Day","Count.total","Count.pups","Count.nonpup","Estimated.total","Estimated.pups","Region","Stock","year","AbundanceEstimate")]
cfall <- as.data.frame(matrix(NA, nrow = (nrow(HCdata)+3), ncol = 2, dimnames = list(c(1:(nrow(HCdata)+3)), c("Stock","cf"))))
cfall[1:length(cf),1] <- "Hood Canal"
cfall[1:length(cf),2] <- cf
source("Scripts/Coastal Estimate.r")
allData <- rbind(allData,cbind(Cdata[,c("Sitecode","Day","Count.total","Count.pups","Count.nonpup","Estimated.total","Estimated.pups","Region","Stock","year")],as.data.frame(matrix(rep("NA", times=nrow(Cdata)), ncol=1, nrow = nrow(Cdata), dimnames = list(1:nrow(Cdata),"AbundanceEstimate")))))
cfall[nrow(HCdata)+1,1] <- "Coastal"
cfall[nrow(HCdata)+1,2] <- cf
source("Scripts/SPS_NI Estimate.r")
allData <- rbind(allData,cbind(NIdata[,c("Sitecode","Day","Count.total","Count.pups","Count.nonpup","Estimated.total","Estimated.pups","Region","Stock","year")],as.data.frame(matrix(rep("NA", times=nrow(NIdata)), ncol=1, nrow = nrow(NIdata), dimnames = list(1:nrow(NIdata),"AbundanceEstimate")))))
cfall[nrow(HCdata)+2,1] <- "Northern Inlands"
cfall[nrow(HCdata)+2,2] <- cf
allData <- rbind(allData,cbind(SPSdata[,c("Sitecode","Day","Count.total","Count.pups","Count.nonpup","Estimated.total","Estimated.pups","Region","Stock","year")],as.data.frame(matrix(rep("NA", times=nrow(SPSdata)), ncol=1, nrow = nrow(SPSdata), dimnames = list(1:nrow(SPSdata),"AbundanceEstimate")))))
cfall[nrow(HCdata)+3,1] <- "Southern Puget Sound"
cfall[nrow(HCdata)+3,2] <- cf


## Take the mean raw adult count vs. mean raw pup count by every year at every haulout and sum for every stock
allmean <- allData %>%
  mutate(year = factor(year, levels = sort(unique(year(allData$Day))))) %>%
  mutate(Count.pups = if_else(Count.pups == 0 & Count.nonpup == 0 & Count.total > 0, NA_real_, Count.pups)) %>%   #there are years (mostly early) where count.total didn't distinguish between pups and adults
  mutate(Count.nonpup = if_else(is.na(Count.pups) & Count.nonpup == 0 & Count.total > 0, NA_real_, Count.nonpup)) %>%
  group_by(Stock,Sitecode,year) %>%
  filter(!is.na(Count.pups)) %>%   #can't calculate ratios or means of data when pups and adults unable to be separated
  summarize(mnCtpups = mean(Count.pups, na.rm = TRUE), mnCtadults = mean(Count.nonpup, na.rm = TRUE), mnCtttl = mean(Count.total, na.rm = TRUE),
            minpup = min(Count.pups, na.rm = TRUE), minadult = min(Count.nonpup, na.rm = TRUE), mintotal = min(Count.total, na.rm = TRUE),
            maxpup = max(Count.pups, na.rm = TRUE), maxadult = max(Count.nonpup, na.rm = TRUE), maxtotal = max(Count.total, na.rm = TRUE)) %>%
  group_by(Stock, year) %>%
  summarize(summnPup = sum(mnCtpups), summnAd = sum(mnCtadults), summnTtl = sum(mnCtttl),
    summinpup = sum(minpup), summinadult = sum(minadult), summintotal = sum(mintotal),
    summaxpup = sum(maxpup), summaxadult = sum(maxadult), summaxtotal = sum(maxtotal)) %>%
  mutate(pupadult = round(summnPup/summnAd, 2)) %>%
  mutate(puptotal = round(summnPup/summnTtl,2)) %>%
  mutate(Stock = as.factor(Stock)) %>%
  mutate(Stock = fct_relevel(Stock, c("Coastal","Northern Inland","Hood Canal","Southern Puget Sound")))

  
## Plot for the sum of mean haulout counts per stock in every year
#To make y-axes match up for two larger and two larger stocks, manually add max value for a single year
blankdat <- data.frame(Stock = rep(factor(c("Coastal","Northern Inland","Hood Canal","Southern Puget Sound"), levels = c("Coastal","Northern Inland","Hood Canal","Southern Puget Sound")),times=2), 
                       year = as.factor(c(rep(1977,times=4),rep(2022,times=4))),
                       summnAd = rep(c(19000,19000,1650,1650),times=2))

countplot1 <- ggplot(data=subset(allmean, Stock %in% c("Coastal","Northern Inland"))) +
  geom_ribbon(aes(x=year, ymin=summinadult, ymax=summaxadult, group=Stock), fill = "#2712c6", alpha=0.3) +
  geom_ribbon(aes(x=year, ymin=summinpup*2.2, ymax=summaxpup*2.2, group=Stock), fill = "#7f9c0b", alpha=0.3) +
  geom_point(aes(x=year, y=summnAd, group=Stock), col="#2712c6", cex=1.5) +
  geom_line(aes(x=year, y=summnAd, group=Stock), col="#2712c6", linewidth=1) +
  geom_point(aes(x=year, y=summnPup*2.2, group=Stock), col="#7f9c0b", cex=1.5) +
  geom_line(aes(x=year, y=summnPup*2.2, group=Stock), col="#7f9c0b", linewidth=1) +
  scale_y_continuous("Summed Mean Count of Adults", sec.axis = sec_axis(~ . / 2.2, name = "Summed Mean Count of Pups")) +
  xlab("Year") +
  facet_wrap(vars(Stock)) +
  theme(axis.text.x = element_text(angle=90, size=8), axis.text.y = element_text(size=12), 
        axis.title = element_text(size=14), strip.text.x = element_text(size=12),
        axis.title.y.right = element_text(margin = margin(t=0, r=0, b=0, l=10))) +
  geom_point(data = subset(blankdat, Stock %in% c("Coastal","Northern Inland")), aes(x = year, y = summnAd), color = "red", alpha = 0)  #blanked out point to get scales the same

countplot2 <- ggplot(data=subset(allmean, Stock %in% c("Hood Canal","Southern Puget Sound"))) +
  geom_ribbon(aes(x=year, ymin=summinadult, ymax=summaxadult, group=Stock), fill = "#2712c6", alpha=0.3) +
  geom_ribbon(aes(x=year, ymin=summinpup*2.2, ymax=summaxpup*2.2, group=Stock), fill = "#7f9c0b", alpha=0.3) +
  geom_point(aes(x=year, y=summnAd, group=Stock), col="#2712c6", cex=1.5) +
  geom_line(aes(x=year, y=summnAd, group=Stock), col="#2712c6", linewidth=1) +
  geom_point(aes(x=year, y=summnPup*2.2, group=Stock), col="#7f9c0b", cex=1.5) +
  geom_line(aes(x=year, y=summnPup*2.2, group=Stock), col="#7f9c0b", linewidth=1) +
  scale_y_continuous("Summed Mean Count of Adults", sec.axis = sec_axis(~ . / 2.2, name = "Summed Mean Count of Pups")) +
  xlab("Year") +
  facet_wrap(vars(Stock)) +
  theme(axis.text.x = element_text(angle=90, size=8), axis.text.y = element_text(size=12), 
        axis.title = element_text(size=14), strip.text.x = element_text(size=12),
        axis.title.y.right = element_text(margin = margin(t=0, r=0, b=0, l=10))) +
  geom_point(data = subset(blankdat, Stock %in% c("Hood Canal","Southern Puget Sound")), aes(x = year, y = summnAd), color = "red", alpha = 0)  #blanked out point to get scales the same

jpeg(file = "Results/Adult_Pup_Counts_Stock.jpg", width=10, height=6, units="in", res=600)
ggarrange(countplot1, countplot2, ncol = 1, nrow = 2)
dev.off()

## Plot for the sum of mean haulout counts per stock in every year showing total count line
countplot3 <- ggplot(data=subset(allmean, Stock %in% c("Coastal","Northern Inland"))) +
  geom_ribbon(aes(x=year, ymin=summintotal, ymax=summaxtotal, group=Stock), fill = "#5A5A5A", alpha=0.3) +
  geom_ribbon(aes(x=year, ymin=summinadult, ymax=summaxadult, group=Stock), fill = "#2712c6", alpha=0.3) +
  geom_ribbon(aes(x=year, ymin=summinpup, ymax=summaxpup*2.2, group=Stock), fill = "#7f9c0b", alpha=0.3) +
  geom_point(aes(x=year, y=summnAd, group=Stock), col="#2712c6", cex=1.5) +
  geom_line(aes(x=year, y=summnAd, group=Stock), col="#2712c6", linewidth=1) +
  geom_point(aes(x=year, y=summnPup, group=Stock), col="#7f9c0b", cex=1.5) +
  geom_line(aes(x=year, y=summnPup, group=Stock), col="#7f9c0b", linewidth=1) +
  facet_grid(cols = vars(Stock)) +
  xlab("Year") + ylab("Summed Mean Count") +
  theme(axis.text.x = element_text(angle=90, size=8), axis.text.y = element_text(size=12), 
        axis.title = element_text(size=14), strip.text.x = element_text(size=12)) +
  geom_point(data = subset(blankdat, Stock %in% c("Coastal","Northern Inland")), aes(x = year, y = summnAd), color = "red", alpha = 0)  #blanked out point to get scales the same

countplot4 <- ggplot(data=subset(allmean, Stock %in% c("Hood Canal","Southern Puget Sound"))) +
  geom_ribbon(aes(x=year, ymin=summintotal, ymax=summaxtotal, group=Stock), fill = "#5A5A5A", alpha=0.3) +   #, color = "#404040"
  geom_ribbon(aes(x=year, ymin=summinadult, ymax=summaxadult, group=Stock), fill = "#2712c6", alpha=0.3) +
  geom_ribbon(aes(x=year, ymin=summinpup, ymax=summaxpup, group=Stock), fill = "#7f9c0b", alpha=0.3) +
  geom_point(aes(x=year, y=summnAd, group=Stock), col="#2712c6", cex=1.5) +
  geom_line(aes(x=year, y=summnAd, group=Stock), col="#2712c6", linewidth=1) +
  geom_point(aes(x=year, y=summnPup, group=Stock), col="#7f9c0b", cex=1.5) +
  geom_line(aes(x=year, y=summnPup, group=Stock), col="#7f9c0b", linewidth=1) +
  facet_grid(cols = vars(Stock), scales = "free_x") +
  xlab("Year") + ylab("Summed Mean Count") +
  theme(axis.text.x = element_text(angle=90, size=8), axis.text.y = element_text(size=12), 
        axis.title = element_text(size=14), strip.text.x = element_text(size=12)) +
  geom_point(data = subset(blankdat, Stock %in% c("Hood Canal","Southern Puget Sound")), aes(x = year, y = summnAd), color = "red", alpha = 0)  #blanked out point to get scales the same

jpeg(file = "Results/Adult_Pup_Counts_Total_Stock.jpg", width=10, height=6, units="in", res=600)
ggarrange(countplot3, countplot4, ncol = 1, nrow = 2)
dev.off()


## Plot of proportion of pups to total animals
cols <- c("#13AFEF","#7f9c0b","#2712c6")
cols2 <- c("#0d76a1","#556907","#180b78")

df1 <- allmean %>%
  select(Stock,year,summnPup) %>%
  rename(MeanCount = summnPup) %>%
  mutate(Group = "Pups")
df2 <- allmean %>%
  select(Stock,year,summnAd) %>%
  rename(MeanCount = summnAd) %>%
  mutate(Group = "Adults")
df3 <- allmean %>%
  select(Stock,year,summnTtl) %>%
  rename(MeanCount = summnTtl) %>%
  mutate(Group = "Total")
allmean_prop <- rbind(df1,df2,df3)
allmean_prop <- allmean_prop[order(allmean_prop$Stock,allmean_prop$year,allmean_prop$Group),]
allmean_prop$Group <- factor(allmean_prop$Group, levels = c("Total","Pups","Adults"))

#Something is amiss with this. Y axis is way too big and Northern Inland values are not centered at 0
# ggplot(allmean_prop, aes(x = as.numeric(as.character(year)), y = MeanCount, fill = Group, color = Group)) +
#   facet_wrap(vars(Stock), scales = "free") +
#   geom_stream(type="mirror") +
#   scale_fill_manual(values = cols) +
#   scale_color_manual(values = cols2) +
#   theme(axis.text.x = element_text(angle=90))

blankdat2 <- data.frame(Stock = rep(factor(c("Coastal","Northern Inland","Hood Canal","Southern Puget Sound"), levels = c("Coastal","Northern Inland","Hood Canal","Southern Puget Sound")), times=46), 
                       year = rep(1977:2022, times = 4),
                       MeanCount = rep(c(40000,40000,3500,3500), times = 46))

countplot5 <- ggplot() +
  geom_area(data=allmean_prop, aes(x = as.numeric(as.character(year)), y = MeanCount, fill = Group, color = Group)) +
  facet_wrap(vars(Stock), scales = "free_y") +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols2) +
  theme_pubclean() +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  xlab("Year") + ylab("Summed Mean Count") +
  geom_area(data=blankdat2, aes(x = as.numeric(as.character(year)), y = MeanCount), alpha = 0)  #blanked out point to get scales the same

jpeg(file = "Results/AdultPuptoTotal_Stock.jpg", width=10, height=6, units="in", res=600)
countplot5
dev.off()


## Plot proportion values only
cols3 <- c("#06c5c9", "#af2bb5", "#e82547", "#028dd6")

prop1 <- ggplot(allmean, aes(x = year, y = puptotal, fill = Stock)) +
  geom_point(pch=22, cex=3) +
  facet_wrap(vars(Stock)) +
  xlab("Year") + ylab("Pup:Total Count Ratio") +
  theme_bw() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), 
        strip.text.x = element_text(size = 14, face = "bold"), axis.text.x = element_text(angle = 90), 
        legend.position = "none") +
  scale_fill_manual(values=cols3) +
  scale_x_discrete(breaks=seq(1977,2022,5), limits=factor(1977:2022))

jpeg(file = "Results/PropPupAdult_Stock.jpg", width=10, height=6, units="in", res=600)
prop1
dev.off()


## Plot total raw count (adults + pups) vs. estimated abundance (count * correction factor) > basically same content as Supplements
HCmean2 <- allData %>%
  filter(Stock == "Hood Canal") %>%
  mutate(year = factor(year, levels = sort(unique(year(allData$Day))))) %>%
  group_by(Sitecode,year) %>%
  summarize(Count = mean(Count.total), EstAbundance = mean(as.numeric(AbundanceEstimate)),
            minCount = min(Count.total), minAbund = min(as.numeric(AbundanceEstimate)),
            maxCount = max(Count.total), maxAbund = max(as.numeric(AbundanceEstimate))) %>%
  group_by(year) %>%
  summarize(summnTtl = sum(Count), summnAbund = sum(EstAbundance),
            summinTtl = sum(minCount), summinAbund = sum(minAbund),
            summaxTtl = sum(maxCount), summaxAbund = sum(maxAbund))

AnnualEstimates=as.vector(HCdata$AbundanceEstimate%*%X)
names(AnnualEstimates)=sort(unique(HCdata$year))

ggplot(HCmean2) +
  geom_ribbon(aes(x=year, ymin=LCL, ymax=UCL), fill = "#3A6B35", alpha=0.3, group=1) +  #summinAbund, summaxAbund
  geom_ribbon(aes(x=year, ymin=summinTtl, ymax=summaxTtl), fill = "#CBD18F", alpha=0.3, group=1) +
  geom_line(aes(x=year, y=summnTtl), col="#CBD18F", linewidth=1, group=1) +
  geom_line(aes(x=year, y=AnnualEstimates), col="#3A6B35", linewidth=1, group=1) +  #sumnAbund
  geom_point(aes(x=year, y=summnTtl), pch=21, col="#686b48", fill="#CBD18F", cex=3) +
  geom_point(aes(x=year, y=AnnualEstimates), pch=21, col="#1b3319", fill="#3A6B35", cex=3) +  #sumnAbund
  ylab("Value")

notHCmean2 <- allData %>%
  filter(Stock != "Hood Canal") %>%
  mutate(year = factor(year, levels = sort(unique(year(allData$Day))))) %>%
  group_by(Stock, Sitecode,year) %>%
  summarize(Count = mean(Count.total), minCount = min(Count.total), maxCount = max(Count.total)) %>%
  group_by(Stock, year) %>%
  summarize(summnTtl = sum(Count), summinTtl = sum(minCount), summaxTtl = sum(maxCount)) %>%
  bind_cols(., Cresults[,c("Abundance","LCL","UCL")])


## Plot correction factors (cf)

cfplot <- ggplot(cfall, aes(y=Stock, x=cf, fill=cf)) +
  geom_jitter(data=subset(cfall, Stock == "Hood Canal"), pch=21, cex=4) +
  geom_point(data=subset(cfall, Stock != "Hood Canal"), pch=21, cex=4) +
  xlab(label = "Correction Factor") +
  scale_fill_continuous(type="viridis") +
  theme(axis.text = element_text(size=12), axis.title = element_text(size=16), legend.position = "none")

jpeg(file = "Results/CorrectionFactor_Stock.jpg", width=10, height=6, units="in", res=600)
cfplot
dev.off()



