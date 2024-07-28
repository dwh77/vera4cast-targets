#Metabolism variables (NEP, GPP, R) from FCRE catwalk and met station data 
#DWH

##libraries needed 
library(tidyverse)
library(suncalc)
library(LakeMetabolizer)


# function for calculating stratified depth in main function 

devtools::source_url("https://raw.githubusercontent.com/dwh77/FCR_Metab/main/Scripts/02_Metab_Model/model_functions/calcZMixDens.R")
devtools::source_url("https://raw.githubusercontent.com/dwh77/FCR_Metab/main/Scripts/02_Metab_Model/model_functions/fillHoles.R")

## these four lines are to define files needed and then test that function works 
file_1 <- "https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-catwalk-data-qaqc/fcre-waterquality_L1.csv" #catwalk github
file_2 <- "https://raw.githubusercontent.com/FLARE-forecast/FCRE-data/fcre-metstation-data-qaqc/FCRmet_L1.csv" #met station github
file_3 <- "https://pasta.lternet.edu/package/data/eml/edi/271/8/fbb8c7a0230f4587f1c6e11417fe9dce"   #catwalk EDI
file_4 <- "https://pasta.lternet.edu/package/data/eml/edi/389/8/d4c74bbb3b86ea293e5c52136347fbb0"  #met station EDI



target_generation_metab_function <- function(file_1, file_2, file_3, file_4){
  
  ## read in current data files from github 
  catwalk_git <- read_csv(file_1)
  
  met_git <- read_csv(file_2)
  
  ## read in historical data files on EDI
  inUrl1  <- file_3
  infile1 <- tempfile()
  try(download.file(inUrl1,infile1,method="curl"))
  if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")
  catwalk_edi <-read_csv(infile1)
  
  inUrl2  <- file_4
  infile2 <- tempfile()
  try(download.file(inUrl2,infile2,method="curl"))
  if (is.na(file.size(infile2))) download.file(inUrl2,infile2,method="auto")
  met_edi <-read_csv(infile2)
  
  ## manipulate the data files to match each other 
  ## bind the two files using row.bind()
  
  catwalk <- bind_rows(catwalk_edi, catwalk_git)
  
  met <- bind_rows(met_edi, met_git)
  
  ## select variables needed and summarize to hourly
  
  catwalk_hourly <- catwalk %>% 
    select(3:13, EXOTemp_C_1, EXODO_mgL_1) %>% 
    mutate(Date = as.Date(DateTime),
           Hour = hour(DateTime)) %>% 
    group_by(Date, Hour) %>% 
    summarise(across(.cols = where(is.numeric), .fns = ~mean(.x, na.rm = TRUE)), .groups = "drop")  %>% 
    mutate(DateTime = ymd_h(paste(Date, Hour, sep = " "))) %>% 
    select(DateTime, 3:14)
  
  colnames(catwalk_hourly) <- c('DateTime', 'temp0.1', 'temp1.0', 'temp2.0', 'temp3.0', 'temp4.0', 'temp5.0', 'temp6.0', 'temp7.0', 'temp8.0', 'temp9.0', 'EXOTemp_C_1', 'EXODO_mgL_1') #rename so calcZmix function works below
  
  
  met_hourly <- met %>% 
    select(DateTime, PAR_umolm2s_Average, WindSpeed_Average_m_s)%>% 
    mutate(Date = as.Date(DateTime),
           Hour = hour(DateTime)) %>% 
    group_by(Date, Hour) %>% 
    summarise(across(.cols = where(is.numeric), .fns = ~mean(.x, na.rm = TRUE)), .groups = "drop")  %>% 
    mutate(DateTime = ymd_h(paste(Date, Hour, sep = " "))) %>% 
    select(DateTime, PAR_umolm2s_Average, WindSpeed_Average_m_s) %>% 
    rename(PAR = PAR_umolm2s_Average)
  
  met_hourly$PAR[met_hourly$PAR<2]<-0
  met_hourly$PAR <- met_hourly$PAR/1000

  
  hourly_data <- left_join(catwalk_hourly, met_hourly, by = "DateTime") %>%
    mutate(DateTime = ymd_hms(DateTime)) %>% 
    filter(DateTime >= ymd("2021-01-01"))
  

  #Get the sunrise and sunset times  
  daysVec <- seq(as.Date(min(hourly_data$DateTime)), as.Date(max(hourly_data$DateTime)), "1 day")

  SunriseSunsetTimes<-getSunlightTimes(date=as.Date(daysVec),keep=c("sunrise","sunset"),lat=37.30,lon=-79.84,tz="EST")  
  #Create data.frame with sunrise, sunset times for each day
  sun <- data.frame(day=daysVec, sunrise=SunriseSunsetTimes$sunrise, sunset=SunriseSunsetTimes$sunset)

  
  #get z mix depth 
  dataTempProfile <- hourly_data %>% 
    select(1:11)
  dataDensProfile <- dataTempProfile |> rename(dateTime = 1)
  
  ##Loop through the columns
  for (j in 2:ncol(dataDensProfile))
  {
    dataDensProfile[,j]=1000*(1 - (dataTempProfile[,j]+288.9414)/(508929.2*(dataTempProfile[,j]+68.12963))*(dataTempProfile[,j]-3.9863)^2)
    ##Calculation for converting temperature to density:
    ##rho = 1000(1 - (T+288.9414)/(508929.2*(T+68.12963))*(T-3.9863)^2)
    # rho_t_24 = 1000*(24 - (24+288.9414)/(508929.2*(24+68.12963))*(24-3.9863)^2)
    ##End of column for loop
  }
  
  
  dataZMix <- calcZMixDens(dataDensProfile)
  
  fluxDummy <- as.numeric(dataZMix$zMix>1.6) #if zmix depth is greater than depth EXO is deployed
  
  flux_zmix_dataframe <- cbind(dataZMix, fluxDummy)
  
  ##Metab calcs
  
  data1 <- left_join(hourly_data, dataZMix, by = c("DateTime" = "dateTime")) %>% 
    select(DateTime, EXODO_mgL_1, EXOTemp_C_1, PAR, WindSpeed_Average_m_s, zMix)
  
  
  #Constants
  Pb <- 101325        #static pressure, pascals
  Tb <- 288.15        #standard temp, K
  Lb <- -0.0065       #standard temp lapse rate, K m-1
  h <- 1663 * 0.3048           #elevation above sea level, m
  hb <- 0             #elevation at bottom of atmospheric layer 0, m (note layer 0 extends to 11000 masl)
  Rstar <-  8.31432   #universal gas constant, N m mol-1 K-1 (equiv to J K-1 mol-1)  SI: 8.314472
  g0 <- 9.80665       #acceleration of gravity, m s-1
  M <- 0.0289644      #molar mass of Earth's air, kg mol-1
  
  #Pressure, in Pa (pascals)
  P <- Pb * (Tb/(Tb+Lb*(h-hb)))^(g0*M/(Rstar*Lb))
  # In mmHg
  atmPres <- P*0.00750061683
  
  
  ##
  #Calculate DO saturation
  #Use eqn from Weiss 1970 Deep Sea Res. 17:721-735; simplified since salinity=0
  # ln DO = A1 + A2 100/T + A3 ln T/100 + A4 T/100
  data2 <- data1 %>% 
    rename(sensorTemp = EXOTemp_C_1)
  
  attach(data2)
  

  
  #Convert sensorTemp to Kelvin
  sensorTempK <- sensorTemp + 273.15
  
  #Weiss equation
  A1 <- -173.4292;  A2 <- 249.6339;  A3 <- 143.3483;  A4 <- -21.8492
  DOSat <- exp(((A1 + (A2*100/sensorTempK) + A3*log(sensorTempK/100) + A4*(sensorTempK/100))))
  
  #Correction for local average atmospheric pressure
  u <- 10^(8.10765 - (1750.286/(235+sensorTemp)))
  DOSat <- (DOSat*((atmPres-u)/(760-u)))   #ml/L
  DOSat <- DOSat/1000                      #L/L
  
  #Convert using standard temperature and pressure. 
  #Similar to calculating saturation DO at STP in ml/L, converting to mg?L (at STP),
  #and then doing the above temperature and pressure conversions.
  R <- 0.082057  #L atm deg-1 mol-1
  O2molWt <- 15.999*2
  convFactor <- O2molWt*(1/R)*(1/273.15)*(760/760) #g/L
  DOSat <- DOSat*convFactor*1000                   #mg/L
  
  ##
  #Calculate kO2
  wp <- 0.15                       #exponent of wind profile power relationship, Smith 1985 Plant, Cell & Environment 8:387-398
  #wind10 <- (10/windHeight)^wp * windSpeed
  wind10 <- WindSpeed_Average_m_s #changing this to just equal wind speed since wind Height was normalized in published met data
  
  k600 <- 2.07 + 0.215*wind10^1.7  #k600 in cm hr-1 per Cole and Caraco 1998;
  k600 <- k600*24/100              #k600 in m day-1
  schmidt <- 1800.6 - 120.1*sensorTemp + 3.7818*sensorTemp^2 - 0.047608*sensorTemp^3
  kO2 <- k600*(schmidt/600)^-0.5   #Jahne et al. 87. exp could also be -.67
  kO2 <- kO2*(60/1440)       #change kO2 to units of m/(timeStep*min)
  
  detach(data2)
  
  ###testing to make calcs in tidyverse/a sepearte function, circle back to later
  # A1 <- -173.4292;  A2 <- 249.6339;  A3 <- 143.3483;  A4 <- -21.8492
  # data3 <- data1 %>% 
  #   mutate(sensorTempK = EXOTemp_C_1 + 273.15,
  #          DOSat = exp(((A1 + (A2*100/sensorTempK) + A3*log(sensorTempK/100) + A4*(sensorTempK/100)))),
  #          u <- 10^(8.10765 - (1750.286/(235+sensorTemp))),
  #          DOSat <- (DOSat*((atmPres-u)/(760-u))),   #ml/L
  #          DOSat <- DOSat/1000,                      #L/L
  #                     )
  

  data3 <- data.frame(dateTime=data2$DateTime, DOObs=data2$EXODO_mgL_1, DOSat=DOSat, irr=data2$PAR, kO2=kO2, zMix=data2$zMix, fluxDummy=fluxDummy)
  
  attach(data3)
  
  # calculate gasflux for each datapoint
  gasflux <- kO2/60 *(DOSat - DOObs) / zMix  # grams O2 m^-3 minute^-1 , 60 is our minute timestep
  
  # Find delta O2 / delta time
  dDO<-diff(DOObs)
  dT=diff(dateTime)
  dDOdT <- dDO / as.numeric(dT,units="mins") # grams O2 m^-3 minute^-1
  
  # calcuate "metabolism" value for each timestep, where m = R in dark, NEP in light
  # since metabolism is calculated for a time interval (opposed to a point), use 
  #   mean of zMix, and gasflux from start and end of each interval
  m <- dDOdT - gasflux[1:length(gasflux)-1]+diff(gasflux)/2
  # grams O2 m^-3 minute^-1
  
  
  # For each night, find areal rate of respiration
  nDays <- dim(sun)[1] - 1  #Not sure if this indexing works appropriately for all lakes
  dateRange <- c(sun$day[1],sun$day[nDays])
  outDays <- seq(dateRange[1],dateRange[2],"1 day")
  
  R_nightly <- rep(NA,nDays)
  for (i in 1:nDays)
  {
    
    r <- which(dateTime>sun$sunset[i] & dateTime<sun$sunrise[i+1])
    R_nightly[i] <- mean(m[r],na.rm=T)*1440    # grams O2 m^-3 day^-1
    rm(r)
  }
  
  # For each daylight period, average R from night before and after
  Rmean<-rep(NA,length(R_nightly)-1)
  for (i in 1:length(R_nightly)-1)
  {
    Rmean[i]<-(R_nightly[i]+R_nightly[i+1])/2 # grams O2 m^-3 day^-1
  }
  
  # Fill in R values for first and last days of deployment based on single nights
  R <- c(R_nightly[1], Rmean)  # grams O2 m^-3 day^-1
  
  # Calculate daylight NEP (this is NOT 24hr NEP)
  Nd <- rep(NA,nDays)
  for (i in 1:length(Nd))
  {
    r <- which(dateTime>sun$sunrise[i] & dateTime<sun$sunset[i])
    Nd[i] <- mean(m[r],na.rm=T) * (as.double(sun$sunset[i]-sun$sunrise[i],units="mins"))
    # grams O2 m^-3 daylight period^-1
    # line above multiplies the rate per minute by the number of daylight minutes
    rm(r)
  }
  
  # Calculate GPP & true NEP (24hr NEP)  # grams O2 m^-3 day^-1    
  GPP <- Nd + (-R * as.double(sun$sunset[1:nDays]-sun$sunrise[1:nDays],units="days") )
  # line above multiplies R (per day) by fraction of day that is daylight.
  NEP <- GPP + R
  
  #Redefine R so that positive respiration rates are displayed as positive
  R <- -R
  
  #Detach data2
  detach(data3)
  
  #Write out results of bookkeeping estimates
  bookOut <- data.frame(outDays,GPP,R)
  
  bookOut_qaqc <- bookOut %>% 
    mutate(GPP_qaqc = ifelse(GPP < 0, NA, GPP),
           R_qaqc = ifelse(R > 0, NA, R),
           NEP_qaqc = GPP_qaqc + R_qaqc)
  

  vera_metab_target <- bookOut_qaqc %>% 
    select(outDays, GPP_qaqc, R_qaqc, NEP_qaqc) %>% 
    rename(datetime = outDays,
           GPP_mgO2_L_day = GPP_qaqc,
           R_mgO2_L_day = R_qaqc,
           NEP_mgO2_L_day = NEP_qaqc) %>% 
    pivot_longer(-1) %>% 
    mutate(site_id = "fcre",
           depth = 1.6) %>% 
    rename(variable = name, observation = value) %>% 
    select(datetime, site_id, depth, observation, variable)
  
  return(vera_metab_target)  

  } #end of target generation bracket 


#test function
test <- target_generation_metab_function(file_1, file_2, file_3, file_4)

test %>%
  # pivot_wider(names_from = variable, values_from = observation) |> 
  ggplot(aes(x = datetime, y = observation, color = variable))+
  geom_point()
  
  





### Looking at bookkeeping from metab chapter 

# bk2019 <- read.table("C:/Users/dwh18/OneDrive/Desktop/R_Projects/FCR_Metabolism_OLD/Data/Model_Output/2019/FCR2019 bookOut.txt")
# bk2020 <- read.table("C:/Users/dwh18/OneDrive/Desktop/R_Projects/FCR_Metabolism_OLD/Data/Model_Output/2020/FCR2020 bookOut.txt")
# bk2021 <- read.table("C:/Users/dwh18/OneDrive/Desktop/R_Projects/FCR_Metabolism_OLD/Data/Model_Output/2021/FCR2021 bookOut.txt")
# 
# bk19_21 <- rbind(bk2019, bk2020, bk2021)
# 
# bk19_21 %>%
#   mutate(solarDay = ymd(solarDay)) %>% 
#   ggplot()+
#   geom_point(aes(x = solarDay, y = GPP), col = "green")+
#   geom_point(aes(x = solarDay, y = R), col = "black")

















