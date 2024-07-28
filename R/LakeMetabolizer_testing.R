library(LakeMetabolizer)
library(tidyverse)

#### fake data Example from help package ----
datetime <- seq(as.POSIXct("2014-06-16 00:00:00", tz="GMT"),
                as.POSIXct("2014-06-17 23:55:00", tz="GMT"), length.out=288*2)
do.obs <- 2*sin(2*pi*(1/288)*(1:(288*2))+1.1*pi) + 8 + rnorm(288*2, 0, 0.5)
wtr <- 3*sin(2*pi*(1/288)*(1:(288*2))+pi) + 17 + rnorm(288*2, 0, 0.15)
do.sat <- LakeMetabolizer::o2.at.sat.base(wtr, 960)
irr <- (1500*sin(2*pi*(1/288)*(1:(288*2))+1.5*pi) +650 + rnorm(288*2, 0, 0.25)) *
  ifelse(is.day(datetime, 42.3), 1, 0)
k.gas <- 0.4
z.mix <- 1
# plot time series
plot(wtr, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
par(new=TRUE); plot(do.obs, type="l", col="blue", xaxt="n", yaxt="n", xlab="", ylab="")
par(new=TRUE); plot(irr, type="l", col="orange", xaxt="n", yaxt="n", xlab="", ylab="")
abline(v=144, lty="dotted")
abline(v=288)
legend("topleft", legend=c("wtr", "do.obs", "irr"), lty=1,
       col=c("black", "blue", "orange"), inset=c(0.08, 0.01))
# put data in a data.frame
data <- data.frame(datetime=datetime, do.obs=do.obs, do.sat=do.sat, k.gas=k.gas,
                   z.mix=z.mix, irr=irr, wtr=wtr)
# run each metabolism model
m.bk <- metab(data, "bookkeep", lake.lat=42.6)
m.bk <- metab(data, lake.lat=42.6) # no method d

#### metab.bookkeep example ----
Sys.setenv(TZ='GMT')
doobs = load.ts(system.file('extdata',
                            'sparkling.doobs', package="LakeMetabolizer"))
wtr = load.ts(system.file('extdata',
                          'sparkling.wtr', package="LakeMetabolizer"))
wnd = load.ts(system.file('extdata',
                          'sparkling.wnd', package="LakeMetabolizer"))

#Subset a day
mod.date = as.POSIXct('2009-07-08', 'GMT')
doobs = doobs[trunc(doobs$datetime, 'day') == mod.date, ]
wtr = wtr[trunc(wtr$datetime, 'day') == mod.date, ]
wnd = wnd[trunc(wnd$datetime, 'day') == mod.date, ]
k.gas = k600.2.kGAS.base(k.cole.base(wnd[,2]), wtr[,3], 'O2')
do.sat = o2.at.sat.base(wtr[,3], altitude=300)
# Must supply 1 for daytime timesteps and 0 for nighttime timesteps
irr = as.integer(is.day(doobs[,1], 45))
metab.bookkeep(doobs[,2], do.sat, k.gas, z.mix=1, irr, datetime=doobs$datetime)


#### trying w/ FCR ----

devtools::source_url("https://raw.githubusercontent.com/dwh77/FCR_Metab/main/Scripts/02_Metab_Model/model_functions/calcZMixDens.R")
devtools::source_url("https://raw.githubusercontent.com/dwh77/FCR_Metab/main/Scripts/02_Metab_Model/model_functions/fillHoles.R")

file_1 <- "https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-catwalk-data-qaqc/fcre-waterquality_L1.csv" #catwalk github
file_2 <- "https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-metstation-data-qaqc/FCRmet_L1.csv" #met station github

#get data together 
fcr_do_24 <- read.csv(file_1)
fcr_met_24 <- read.csv(file_2)

#get met and water Q
met_var <- fcr_met_24 |> 
  select(DateTime, WindSpeed_Average_m_s, PAR_umolm2s_Average, BP_Average_kPa) |> 
  rename(dateTime = DateTime) |> 
  # mutate(datetime = ymd_hms(datetime)) |> 
  mutate(BP_Average_millibar = BP_Average_kPa*10) |> 
  select(dateTime, WindSpeed_Average_m_s, PAR_umolm2s_Average, BP_Average_millibar)


water_var <- fcr_do_24 |> 
  # mutate(DateTime = ymd_hms(DateTime)) |> 
  select(DateTime, EXODO_mgL_1, EXOTemp_C_1) |> 
  rename(dateTime = DateTime, 
         do.obs = EXODO_mgL_1,
         wtr = EXOTemp_C_1)

# do <- water_var |>  select(dateTime, do.obs)
# 
# do <- fillHoles(dataIn = do, maxLength = 60, timeStep = 10)
# 
# water_var$do.obs <- fillHoles(water_var$do.obs, maxLength=60, timeStep=10)

#get z mix
dataTempProfile <- fcr_do_24 |> 
  select(DateTime, 4:13) |> 
   #mutate(DateTime = ymd_hms(DateTime)) |> 
  rename(dateTime = DateTime,
         temp0.1 = ThermistorTemp_C_surface,
         temp1.0 = ThermistorTemp_C_1,
         temp2.0 = ThermistorTemp_C_2,
         temp3.0 = ThermistorTemp_C_3,
         temp4.0 = ThermistorTemp_C_4,
         temp5.0 = ThermistorTemp_C_5,
         temp6.0 = ThermistorTemp_C_6,
         temp7.0 = ThermistorTemp_C_7,
         temp8.0 = ThermistorTemp_C_8,
         temp9.0 = ThermistorTemp_C_9)

dataDensProfile<-dataTempProfile

for (j in 2:ncol(dataDensProfile)){
  dataDensProfile[,j]=1000*(1 - (dataTempProfile[,j]+288.9414)/(508929.2*(dataTempProfile[,j]+68.12963))*(dataTempProfile[,j]-3.9863)^2)
}

dataZMix <- calcZMixDens(dataDensProfile, thresh = 0.1)


#bring data together and make final conversions
fcr_all_vars <- left_join(water_var, met_var, by = "dateTime")
fcr_all_vars <- left_join(fcr_all_vars, dataZMix, by = c("dateTime"))

fcr_vars <- fcr_all_vars |> 
  rename(datetime = dateTime) |> 
  mutate(datetime = as.POSIXct(datetime, tz = "GMT")) |>
  mutate(do.sat = (o2.at.sat.base(wtr, BP_Average_millibar, altitude = (1663 * 0.3048) ))) |> 
  mutate(k_600_cole = k.cole.base(WindSpeed_Average_m_s)) |> 
  mutate(k600GAS = k600.2.kGAS.base(k_600_cole, wtr, gas="O2")) |> 
  mutate(nighttime = is.night(datetime, lat = 37.30),
         day_night = ifelse(nighttime == T, 0, 1))  
#   # filter(!is.na(do.obs)) |> 
# 
#   filter(datetime >= ymd_hms("2024-04-24 00:00:00"),
#          datetime < ymd_hms("2024-04-25 00:00:00")) 
# 
# fcr_vars$do.obs <- zoo::na.approx(fcr_vars$do.obs)
# fcr_vars$do.sat <- zoo::na.approx(fcr_vars$do.sat)
# fcr_vars$k600GAS <- zoo::na.approx(fcr_vars$k600GAS)
# 
# 
# 
# ## try bookkeeping 
# bk <- metab.bookkeep(datetime = fcr_vars$datetime, do.obs = fcr_vars$do.obs, do.sat = fcr_vars$do.sat, 
#                      k.gas = fcr_vars$k600GAS, z.mix = fcr_vars$zMix, irr = fcr_vars$day_night)
# 
# bk
# 
# metab(fcr_vars, "bookkeep", lake.lat=37.30)

#### Trying to make loop that goes by date, and then within a date will fill each gap

dates <- seq(min(as.Date(fcr_vars$datetime)), max(as.Date(fcr_vars$datetime)), by = "day")

metab_outputs <- data.frame(datetime = dates, GPP = NA, R = NA, NEP = NA)

for (i in dates) {
  
  fcr_vars <- fcr_vars |> 
    filter(as.Date(datetime) == i)
  
  fcr_vars$do.obs <- zoo::na.approx(fcr_vars$do.obs)
  fcr_vars$do.sat <- zoo::na.approx(fcr_vars$do.sat)
  fcr_vars$k600GAS <- zoo::na.approx(fcr_vars$k600GAS)
  
  bk <- metab.bookkeep(datetime = fcr_vars$datetime, do.obs = fcr_vars$do.obs, do.sat = fcr_vars$do.sat, 
                      k.gas = fcr_vars$k600GAS, z.mix = fcr_vars$zMix, irr = fcr_vars$day_night)
  
  metab_outputs <- left_join(metab_outputs, bk, by = "datetime")
  
  return(metab_outputs)
  
}







