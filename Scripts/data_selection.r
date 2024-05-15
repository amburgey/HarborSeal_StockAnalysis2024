library(lubridate)

#Pre-analysis code cleaning; if surveys at a haulout were done on the same day, they have to be >50min apart to count as separate surveys
#See Jeff Laake's code from pvCode.R in Archived_preOuterCoastAdditions

#Selection of years and dates for surveys

# Read in pvdf.csv (issues reading date in due to non-ambiguous format so use lubridate as well to specify)
pv.df=read.csv(file="Data/pvdf_FINALIZED_2023.csv",colClasses=c("numeric","character","character","character",rep("integer",5),rep("factor",3),"numeric","numeric","factor",
rep("numeric",5))) %>%
  mutate(Day = mdy(Day)) %>%
  mutate(Survey.time = mdy_hm(Survey.time))

#assign use = TRUE if in date range and FALSE otherwise
nonhood.dates <- 198:245
hood.dates <- 227:314
coastal.dates <- 121:181
nonhood <- c("Northern Inland", "Southern Puget Sound")
#Piece out target stocks by target dates
pv.df$use=FALSE
pv.df$use[pv.df$Stock %in% nonhood & pv.df$Julian %in% nonhood.dates]=TRUE
pv.df$use[pv.df$Stock=="Hood Canal" & pv.df$Julian %in% hood.dates]=TRUE
pv.df$use[pv.df$Stock=="Coastal" & pv.df$Julian %in% coastal.dates]=TRUE

# Years to use by region
HCYears=c(1977,1978,1984,1988,1991,1992,1993,1996,1998,1999,2000,2002,2003,2004,2005,2010,2013,2019,2022,2023)
SPSYears=c(1985,1991,1992,1993,1994,1996,1998,2000,2004,2005,2008,2010,2013,2014,2019,2023)
EBYears=c(1983,1984,1985,1986,1987,1988,1989,1991,1992,1993,1994,1995,1996,1998,1999,2000,2001,2002,2004,2005,2006,2007,2008,2010,2013,2014,2019,2022)
SJIYears=c(1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1998,1999,2000,2001,2002,2004,2005,2006,2010,2013,2014,2019,2022)
SJFYears=c(1983,1984,1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,1995,1996,1998,1999,2000,2004,2013,2014,2019,2022)
CEYears=c(1975:1978,1980:1989,1991:1997,1999:2001,2004:2005,2007,2014,2022)
OCYears=c(1980:1981,1983,1986:1987,1989,1991:1997,1999:2001,2004:2005,2007,2014,2022)

#restrict to dates
pv.df.seldates=pv.df[pv.df$use,]

# restrict to dates and years used
pv.df.seldatesyears=droplevels(pv.df.seldates[(pv.df.seldates$Region=="Hood Canal" & pv.df.seldates$Year%in%HCYears) |  
                                                (pv.df.seldates$Region=="San Juan Islands" & pv.df.seldates$Year%in%SJIYears) | 
                                                (pv.df.seldates$Region=="Strait of Juan de Fuca" & pv.df.seldates$Year%in%SJFYears) | 
                                                (pv.df.seldates$Region=="Eastern Bays" & pv.df.seldates$Year%in%EBYears) |
                                                (pv.df.seldates$Region=="Puget Sound" & pv.df.seldates$Year%in%SPSYears)|
                                                (pv.df.seldates$Region%in%c("Grays Harbor","Willapa Bay")& pv.df.seldates$Year%in%CEYears)|
                                                (pv.df.seldates$Region%in%c("Olympic Coast N","Olympic Coast S")& pv.df.seldates$Year%in%OCYears)
                                              ,])


