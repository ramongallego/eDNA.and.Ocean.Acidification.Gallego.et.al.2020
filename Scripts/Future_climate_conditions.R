#library(tidyverse)
#library(propagate)
library (fitdistrplus)
library (tidyverse)
library (vegan)
library (tidybayes)
library (here)

now.env.data <- read_csv(here("Input","Combined_Biol_Env_Plankton.csv")) %>% 
  select(event, Area, mean.TA, mean.DIC = mean.DIC.new.sal,Salinity = Salinity.new ,Temperature, Omega.aragonite, pH = pH_new ) %>% 
  distinct() %>%  
  filter(pH > 7.4) 


#calculations based upon Khangaonkar, T., Nugraha, A., Xu, W., & Balaguru, K. (2019). Salish Sea response to global climate change, sea level rise, and future nutrient loads. Journal of Geophysical Research: Oceans, 124. https://doi.org/10.1029/ 2018JC014670
#and in particular, upon Figure 14 in that paper. 

#1. The paper represents a single model run focusing on the Salish Sea. Hence, there are no available error estimates for changes in mean parameters over time.
#2. But because the model shows a good deal of spatial variation within the Salish Sea for that single model run, we can do the following:
#3. Assume spatial variation of means over the year is greater than the likely variation among model runs.
#4. Use (large) spatial variation as a conservative estimate of within-space variation among model runs. This is the best we can do.

#interpreted from Khangaonkar, Figure 14:

# pH.2000.mean <- 7.81
# pH.2000.CI25 <- 7.73
# pH.2000.CI75 <- 7.9
# pH.2095.mean <- 7.63
# pH.2095.CI25 <- 7.54
# pH.2095.CI75 <- 7.78
# 
# T.2000.mean <- 10.01
# T.2000.CI25 <- 9.0
# T.2000.CI75 <- 11.8
# T.2095.mean <- 11.52
# T.2095.CI25 <- 10.6
# T.2095.CI75 <- 13.5
# 
# 
# #####TO COERCE TO NORMAL DISTRIBUTIONS
# 
#     #forcing these into symmetrical distributions by empirically sampling different normal distributions to match
#     # rnorm(1000, mean = pH.2000.mean, sd = 0.25) %>% 
#     #   quantile(c(.25, .75))
#     # rnorm(1000, mean = pH.2095.mean, sd = 0.3) %>% 
#     #   quantile(c(.25, .75))
#     # rnorm(1000, mean = T.2000.mean, sd = 2.4) %>% 
#     #   quantile(c(.25, .75))
#     # rnorm(1000, mean = T.2095.mean, sd = 2.7) %>% 
#     #   quantile(c(.25, .75))
#     
#     #or by approximation using the t-quantile method, in which
#     #SE = (CIupper - mean)/tquantile
#     pH.2000.sd <- (pH.2000.CI75 - pH.2000.mean)/.75
#     pH.2095.sd <- (pH.2095.CI75 - pH.2095.mean)/.75
#     T.2000.sd <- (T.2000.CI75 - T.2000.mean)/.75
#     T.2095.sd <- (T.2095.CI75 - T.2095.mean)/.75
# 

###PROPAGATE ERROR
    
# #Estimate the total change in parameter values between 2000 and 2095, using package "propagate"
#   pH.year2000 <- c(pH.2000.mean, pH.2000.sd)
#   pH.year2095 <- c(pH.2095.mean, pH.2095.sd)
#     pH.df <- cbind(pH.year2000,pH.year2095)
# pH.delta <- propagate(expr = expression(pH.year2095 - pH.year2000),
#          data = pH.df)
# 
#   T.year2000 <- c(T.2000.mean, T.2000.sd)
#   T.year2095 <- c(T.2095.mean, T.2095.sd)
#     T.df <- cbind(T.year2000,T.year2095)
# T.delta <- propagate(expr = expression(T.year2000 - T.year2095),
#                       data = T.df)

    
    
###BUT given that our observed data are from non-random corners of the Salish Sea, these observed data are already way off of what 
    #Khangaonkar would suggest (in particular, our observed SD is far greater than they'd predict even for 2095). So perhaps the best way forward is simply taking our observed data and moving the means according to the average forcasted in the paper:
    #so, between 2000 and 2095, deltaT = 1.51, deltapH = -0.18. If we assume linear change, that is 0.01589474 degrees C per year, and -0.001894737 pH units, so
    #in between our observed year, 2017, and 2095 (78 years) we should see 1.24 degrees C change and -0.148 pH units. This model will assume the variance of climate params is constant over the coming century.
    
    get_pH_distribution_simple <- function(year, area) { # for any year between 2000 and 2095
      # returns mean and sd for Salish Sea for any given year, by simply shifting the observed data appropriately.
      # bases results on two underlying normal approximations, one for each of the Areas in our observed dataset
      
      if (year < 2000 | year > 2095) {
        print("year outside range 2000 - 2095")
      } else {
        
        
        ph.normal.params <- now.env.data %>% 
          filter(Area == area) %>% 
          pull(pH) %>% 
          fitdist(distr = "norm") 
        
        pH.delta.mean.model <- lm(deltapH ~ Year, data = data.frame(
          deltapH = c(0, -0.18),
          Year = c(2000, 2095)
        )) 
        
        predict(pH.delta.mean.model, newdata = data.frame(Year = 2020:2080)) %>% fitdist(distr = "norm") -> increase.model.params
        return(data.frame(
          mean =  ph.normal.params$estimate[1] +  #observed
            predict(pH.delta.mean.model, newdata = data.frame(Year = year)), #predicted change
          sd = ph.normal.params$estimate[2] + increase.model.params$estimate[2]  #observed
        ))
      }
    }
    #get_pH_distribution_simple(2017, "SJI")  #example
   
    
    
    
  r_pH_year <- function(year, area, ndraws) { # generate n random draws from the pH distribution in the Salish Sea in a given year
    if (year < 2000 | year > 2095) {
      print("year outside range 2000 - 2095")
    } else {
      tmp.mean <- as.numeric(
        get_pH_distribution_simple(year, area)[1]
      )
      tmp.sd <- as.numeric(
        get_pH_distribution_simple(year, area)[2]
      )
      return(rnorm(ndraws, mean = tmp.mean, sd = tmp.sd))
    }
  }
#    r_pH_year(2020,"Hood Canal", 100)   #example
  
  r_pH_year_STD <- function(year, area, ndraws, REFyear) { # 
    #generate n random draws from the pH distribution in the Salish Sea in a given year
    #IN STANDARDIZED UNITS, standardized to REFyear
    if (year < 2000 | year > 2095) {
      print("year outside range 2000 - 2095")
    } else {
      
      distribution.realization <- get_pH_distribution_simple(year, area)
      distribution.reference <- get_pH_distribution_simple(REFyear, area)
      
      
      tmp.mean <- as.numeric(
        distribution.realization[1]
      )
      tmp.sd <- as.numeric(
        distribution.realization[2]
      )
      
      tmp.distr <- rnorm(ndraws, mean = tmp.mean, sd = tmp.sd)
      
      return(
      (tmp.distr - 
          as.numeric(distribution.reference[1])) / 
        as.numeric(distribution.reference[2])
      )
    }
  }
  #r_pH_year_STD(2095, "SJI", 100, 2017)  #example 
  
  
  
  get_Temperature_distribution_simple <- function(year, area) { # for any year between 2000 and 2095
    # returns mean and sd for Salish Sea for any given year, by simply shifting the observed data appropriately.
    # bases results on two underlying normal approximations, one for each of the Areas in our observed dataset
    
    if (year < 2000 | year > 2095) {
      print("year outside range 2000 - 2095")
    } else {
      
      
      Temperature.normal.params <- now.env.data %>% 
        filter(Area == area) %>% 
        pull(Temperature) %>% 
        fitdist(distr = "norm") 
      
      Temperature.delta.mean.model <- lm(deltaTemperature ~ Year, data = data.frame(
        deltaTemperature = c(0, 1.51),
        Year = c(2000, 2095)
      ))
      
      predict(Temperature.delta.mean.model, newdata = data.frame(Year = 2020:2080)) %>% fitdist(distr = "norm") -> increase.model.params
      return(data.frame(
        mean =  Temperature.normal.params$estimate[1] +  #observed
          predict(Temperature.delta.mean.model, newdata = data.frame(Year = year)), #predicted change
        sd = Temperature.normal.params$estimate[2]   #observed
      ))
    }
  }
  #get_Temperature_distribution_simple(2089, "Hood Canal")  #example
  
  
  
  
  r_Temperature_year <- function(year, area, ndraws) { # generate n random draws from the Temperature distribution in the Salish Sea in a given year
    if (year < 2000 | year > 2095) {
      print("year outside range 2000 - 2095")
    } else {
      tmp.mean <- as.numeric(
        get_Temperature_distribution_simple(year, area)[1]
      )
      tmp.sd <- as.numeric(
        get_Temperature_distribution_simple(year, area)[2]
      )
      return(rnorm(ndraws, mean = tmp.mean, sd = tmp.sd))
    }
  }
  #    r_Temperature_year(2020, "SJI", 100)   #example
  
  r_Temperature_year_STD <- function(year, area, ndraws, REFyear) { # 
    #generate n random draws from the Temperature distribution in the Salish Sea in a given year
    #IN STANDARDIZED UNITS, standardized to REFyear
    if (year < 2000 | year > 2095) {
      print("year outside range 2000 - 2095")
    } else {
      
      distribution.realization <- get_Temperature_distribution_simple(year, area)
      distribution.reference <- get_Temperature_distribution_simple(REFyear, area)
      
      
      tmp.mean <- as.numeric(
        distribution.realization[1]
      )
      tmp.sd <- as.numeric(
        distribution.realization[2]
      )
      
      tmp.distr <- rnorm(ndraws, mean = tmp.mean, sd = tmp.sd)
      
      return(
        (tmp.distr - 
           as.numeric(distribution.reference[1])) / 
          as.numeric(distribution.reference[2])
      )
    }
  }
  #r_Temperature_year_STD(2095, "SJI", 100, 2017)  #example 
  


  
  