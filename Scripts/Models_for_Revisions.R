library(tidyverse)
library(stats)
library(rstanarm)
library(loo)
library(bayesplot)
library(tidybayes)
library(modelr)
library(here)
library(styler)

#uncomment below, if sourcing at command line (because this apparently doesn't respect here())
#setwd("/Volumes/GoogleDrive/My Drive/Kelly_Lab/Projects/OA_eDNA/Analysis/Splitting.data.by.lifestyle")


####LOAD DATA
options(mc.cores = parallel::detectCores())
data.for.stan <- read_csv("Input/data.for.stan.csv")
now.env.data <- read_csv("Input/env.data.for.log.regression.csv") %>%
  filter(pH > 7.4) %>%
  mutate(
    Season = case_when(
      str_detect(sample, "1703|1710|1711|1801|1803") ~ "Winter",
      TRUE ~ "Summer"
    ),
    Area = case_when(
      str_detect(sample, "CP|LK|FH") ~ "SJI",
      TRUE ~ "Hood Canal"
    ))
  
######FUNCTIONS
### DEFINE linear models to switch back and forth between standardized and natural-scale data
value2std <- function(df){lm(std~value, data = df)}
std2value <- function(df){lm(value~std, data = df)}

ModelPlot <- function(df, variable) {
  require(tidyverse)
  temp.plot <-
    df %>%
    mutate(Posterior = list(apply(posterior_linpred(df$model[[1]], transform = T), 2, median))) %>%
    select(-model) %>%
    unnest() %>%
    ggplot(aes(x = eval(as.symbol(variable)), y = binary, color = Area)) +
    geom_point() +
    geom_smooth(aes(y = Posterior)) +
    facet_grid(Area ~ .) +
    ggtitle(label = df$taxa[1]) +
    xlab(variable) +
    ylab("Probability of Occurrence")
  
  return(temp.plot)
}
#
ProjectionPlot <- function(taxon, modelset = mod.proj.phylum) {
  temp.plot <- modelset %>%
    mutate(meanRichness = round(meanRichness, 1)) %>%
    group_by(phylum, Area) %>%
    mutate(ScaledMeanRichness = meanRichness / max(meanRichness)) %>%
    filter(phylum == taxon) %>%
    ggplot(aes(x = Temperature, y = pH, fill = ScaledMeanRichness)) +
    geom_raster(interpolate = T) +
    facet_grid(. ~ Area) +
    ggtitle(label = taxon)
  
  return(temp.plot)
}
######


######DEFINE MODEL FOR STAN

           safe.stan.glm <- possibly(stan_glmer, otherwise = "Error")
           t_prior <- student_t(df = 7, location = 0, scale = 2.5)
          
          #  ModelFun <- function(df) {
          #   safe.stan.glm(binary ~ (1 | Season) + (1 + Temperature + pH | Area),
          #     family = "binomial",
          #     data = df,
          #     chains = 3,
          #     iter = 2000,
          #     prior = t_prior,
          #     prior_intercept = t_prior,
          #     adapt_delta = .99
          #   )
          # }
          # 
          #  #seasonality, but no universal intercept
          #  ModelFun1 <- function(df) {
          #    safe.stan.glm(binary ~ 0 + (1 | Season) + (1 + Temperature + pH | Area),
          #                  family = "binomial",
          #                  data = df,
          #                  chains = 3,
          #                  iter = 2000,
          #                  prior = t_prior,
          #                  prior_intercept = t_prior,
          #                  adapt_delta = .99
          #    )
          #  }
           
          #no universal intercept or seasonality
           ModelFun2 <- function(df) {
             safe.stan.glm(binary ~ 0 + (1 + Temperature + pH | Area),
                           family = "binomial",
                           data = df,
                           chains = 3,
                           iter = 2000,
                           prior = t_prior,
                           prior_intercept = t_prior,
                           adapt_delta = .99
             )
           } 
           
           # #universal intercept but no seasonality
           # ModelFun3 <- function(df) {
           #   safe.stan.glm(binary ~ 1 + (1 + Temperature + pH | Area),
           #     family = "binomial",
           #     data = df,
           #     chains = 3,
           #     iter = 2000,
           #     prior = t_prior,
           #     prior_intercept = t_prior,
           #     adapt_delta = .99
           #   )
           # }
           
                     
          # # far faster modeling, but imprecise; do not use for real results, but use to get a sense of what will work
          # QuickModelFun <- function(df) {
          #   safe.stan.glm(binary ~ (1 | Season) + (1 + Temperature + pH | Area),
          #     family = "binomial",
          #     data = df,
          #     adapt_delta = .99,
          #     algorithm = "meanfield"
          #   )
          # }


###################################
###################################

#take observed data and standardize environmental variables to be centered at zero w sd = 1
        now.df <- now.env.data %>% 
          ungroup() %>% 
          # select(-sample) %>% 
          gather(key = "key", value = "value", -Season,-Area, -sample) %>% 
          group_by(key, Area) %>%   #calculating different scale for each area; standardizing to different means
          mutate(std = vegan::decostand(value, method = "standardize")) %>%  
          ungroup() %>% 
          mutate(Season = as.factor(Season),
                 Area = as.factor(Area))

        
        #to help conversions from raw to standardized data and back again; custom models for each data subset (e.g., by Area)
        Temperature.mod <- now.df %>% 
          filter(key == "Temperature")  %>% 
          nest(dataset = c(value, std, Season)) %>% 
          mutate(toStd = map(dataset, value2std),
                 toValue = map(dataset, std2value))
        pH.mod <- now.df %>% 
          filter(key == "pH")  %>% 
          nest(dataset = c(value, std, Season)) %>% 
          mutate(toStd = map(dataset, value2std),
                 toValue = map(dataset, std2value))
        
        now.df %>% 
          filter(key %in% c("Temperature", "pH"))  %>% 
          nest(dataset = c(sample, value, std, Season)) %>% 
          mutate(toStd = map(dataset, value2std),
                 toValue = map(dataset, std2value)) %>% write_rds("Input/std_to_real.rds")

        #set up paramater space in natural units
        paramSpace.SJI <- expand.grid(Temperature = seq(7, 15, 1),
                                 pH = seq(7.3, 8.3, 0.1),
                                 Area = "SJI")
        paramSpace.HC <- expand.grid(Temperature = seq(7, 25, 1),
                                  pH = seq(7.6, 8.6, 0.1),
                                  Area = "Hood Canal")
        
        
        #convert this parameter space to standardized units, based on current observations
        
          new <- rbind(paramSpace.SJI, paramSpace.HC) %>%
                  #group_by(Season, Area) %>%
                  mutate(Temperature = case_when(Area == "SJI" ~ 
                                                   predict(Temperature.mod$toStd[[1]], newdata = data.frame(value = .$Temperature)),
                                                 Area == "Hood Canal" ~ 
                                                   predict(Temperature.mod$toStd[[2]], newdata = data.frame(value = .$Temperature)))
                         ) %>% 
                mutate(pH = case_when(Area == "SJI" ~ 
                                        predict(pH.mod$toStd[[1]], newdata = data.frame(value = .$pH)),
                                      Area == "Hood Canal" ~ 
                                        predict(pH.mod$toStd[[2]], newdata = data.frame(value = .$pH)))
                        )
          new %>% 
            distinct()
          #visualize for sanity check
          
          now.df %>% filter (key %in% c("pH", "Temperature")) %>% 
            ggplot(aes(x = value, fill = Area)) + 
            geom_density(alpha=0.5) +
            facet_grid(~key, scales = "free")
            
          new %>% 
            pivot_longer(cols = c(Temperature, pH)) %>% 
            ggplot(aes(x = value, fill = Area)) + 
            geom_bar()+
          # geom_density(alpha=0.5) +
            facet_grid(~name, scales = "free")
          
          
          
          

###################################
###################################
      #FIT MODELS

      # to read in existing model output
      #a <- readRDS("ModelOutputs/full_model_output_rpk_TpH.standardized_2019-08-13_1567.RDS")
          data.for.stan %>% summarise(n_distinct(sample))
      # to run new model
      a <- data.for.stan %>%
        group_by(Area) %>% 
        mutate(Temperature = vegan::decostand(Temperature, method = "standardize"),
               pH = vegan::decostand(pH, method = "standardize")) %>%
        group_by(taxa) %>%
        mutate(Noccur = sum(binary)) %>%
        filter(Noccur > 7 & Noccur < 67) %>%
        group_by(taxa) %>%
        nest() %>%
        mutate(model = map(data, ModelFun2)) #model without seasonality
      
          # to save data 
          saveRDS(a, file = here("ModelOutputs", paste0("full_model_output_rpk_TpH.standardized_nonSeasonal",
                                   Sys.Date(),  #add date to model output
                                   "_",
                                   sample(1:10000,1), #add random integer to label, in case you run a bunch of models in a day
                                   ".RDS")))
          
          data.for.stan %>%
            group_by(Area) %>% 
            mutate(Temperature = vegan::decostand(Temperature, method = "standardize"),
                   pH = vegan::decostand(pH, method = "standardize")) %>%
            group_by(taxa) %>%
            mutate(Noccur = sum(binary)) %>% 
            filter(Noccur < 7 | Noccur > 67) %>% 
            group_by(taxa) %>%
            nest() %>%
            mutate(model = map(data, ModelFun2)) %>% 
            #model without seasonality
          
          # to save data 
          saveRDS( file = here("ModelOutputs", paste0("full_model_output_rpk_TpH.standardized_nonSeasonal_outliers",
                                                        Sys.Date(),  #add date to model output
                                                        "_",
                                                        sample(1:10000,1), #add random integer to label, in case you run a bunch of models in a day
                                                        ".RDS")))       
          
          
      a <- read_rds("Input/full_model_output_rpk_TpH.standardized_nonSeasonal2020-06-24_3262.RDS")   
          
      # to visualize models individually
      # i = 16
      # a$data[[i]] %>%
      #   #group_by(Area, Season) %>%
      #   #data_grid(T.resid = seq_range(T.resid, n = 20),  #tidy equiv of "newdata"
      #   #          pH.resid = seq_range(pH.resid, n = 20)) %>%
      #   add_fitted_draws(a$model[[i]]) %>%  #tidy equiv of posterior_linpred
      #   ggplot(aes(x = Temperature, y = binary, color = Area)) +
      #   stat_lineribbon(aes(y = .value)) +
      #   geom_point(data = a$data[[i]]) +
      #   scale_fill_brewer(palette = "Greys") +
      #   scale_color_brewer(palette = "Set2") +
      #   facet_grid(Area~.) +
      #   ggtitle(label = a$taxa[i])
      
      # if (length(which(sapply(a$model, is.character))) > 0){a <- a[-which(sapply(a$model, is.character)),]}
      #


###################################
###################################
          #project models into the future
          
          #add phylum info to dataset; there's a better way to deal w nesting, but whatever
          phyla <- a %>%
            unnest(data) %>%
            select(taxa, phylum) %>%
            unique() %>%
            pull(phylum)
          a$phylum <- as.factor(phyla)
          levels(a$phylum) <-  c(levels(a$phylum),  "Unid")
          a$phylum[is.na(a$phylum)] <- "Unid"
          
          
          tmp <- a %>%
            group_by(taxa, phylum) %>%
            mutate(projections = map(model,
                                     add_predicted_draws,
                                     newdata = new,
                                     n = 100)) %>% 
            select(taxa, phylum, projections) %>%
            unnest(cols = projections) %>%
            ungroup() %>% 
            select(taxa, phylum, Temperature, pH, Area, .draw, .prediction)
          
          
          scenarios <- unique(tmp[,c("Temperature", "pH", "Area")]) %>% 
            group_by(Area) %>% 
            nest() %>% 
            mutate(data = map (data, ~.x %>% rownames_to_column("Scenario"))) %>% 
            unnest(data) %>% 
            unite(Area, Scenario, col = "Scenario", sep = "_", remove = F) 
          
          mod.proj.overall <- left_join(tmp, scenarios) %>% 
            group_by(Temperature, pH, Area, .draw) %>%
            dplyr::summarize(Richness = sum(.prediction)) %>% 
            group_by(Temperature, pH, Area) %>%
            summarize(meanRichness = mean(Richness)) 
          
          mod.proj.phylum <- left_join(tmp, scenarios) %>% 
            group_by(phylum, Temperature, pH, Area, .draw) %>%
            dplyr::summarize(Richness = sum(.prediction)) %>% 
            group_by(phylum, Temperature, pH, Area) %>%
            summarize(meanRichness = mean(Richness)) 
          
          mod.proj.species <- left_join(tmp, scenarios) %>% 
            group_by(taxa, phylum, Temperature, pH, Area, .draw) %>%
            dplyr::summarize(Prob = sum(.prediction)) %>% 
            group_by(taxa, phylum, Temperature, pH, Area) %>%
            summarize(meanProb = mean(Prob)) 
          
          
          # mod.proj <- readRDS("model_projections_rpk.RDS")
          

          mod.proj.phylum %>% 
            filter(phylum == "Bacillariophyta") %>% 
            filter(Area == "Hood Canal") %>% 
            ggplot(aes(x = Temperature, y = pH, fill = meanRichness)) +
              geom_tile()
          
          mod.proj.species %>% 
            filter(phylum == "Bacillariophyta") %>% 
            filter(Area == "Hood Canal") %>% 
            #separate(taxa, into = c("Family", "Genus", "Species"))
            ggplot(aes(x = Temperature, y = pH, fill = meanProb)) +
              geom_tile() + 
              facet_wrap(~as.factor(taxa))
          
          
          mod.proj.overall %>% 
            mutate(meanRichness = round(meanRichness, 1)) %>%
            group_by(Area) %>%
            mutate(ScaledMeanRichness = meanRichness / max(meanRichness)) %>% 
            #filter(Area == "SJI", Season == "Winter") %>% 
            ggplot(aes(x = Temperature, y = pH, fill = ScaledMeanRichness)) +
              geom_raster(interpolate = T) +
              facet_grid(~ Area) + 
              ggtitle(label = "Overall Richness") 
          
          
          
# (projectionplot <- ProjectionPlot("Bacillariophyta", modelset = mod.proj.phylum))
# 
#           std.to.val.Temp <- scales::trans_new(name = "std.to.val.Temp",
#                   transform = function(x) predict(Temperature.mod$toValue[[1]], newdata = data.frame(value = x)),
#                   inverse = function(x) predict(Temperature.mod$toStd[[1]], newdata = bind_rows(stats::model.frame(x), data.frame(std = x))))    


pdf("ProjectionPlots/full.projectionPlot.standardized.pdf")
  lapply(names(head(sort(table(a$phylum), decreasing = T), 10)),
         ProjectionPlot)
dev.off()




### which taxa are most affected by change in pH?
#
# create a df w model coefficients
coef.df <- data.frame(matrix(NA, nrow = nrow(a), ncol = 6))
for (i in 1:nrow(a)) {
  coef.df[i, ] <- a$model[[i]]$coefficients
}
coef.df$taxa <- a$taxa
names(coef.df) <- c(names(a$model[[i]]$coefficients), "taxa")
names(coef.df) <- gsub(":", "_", names(coef.df))


# now, easy to see and plot most extreme cases
coef.df %>%
  select(taxa, `b[pH Area_Hood_Canal]`) %>%
  mutate(sign = ifelse(`b[pH Area_Hood_Canal]` < 0, "-", "+")) %>% 
  mutate(`b[pH Area_Hood_Canal]` = abs(`b[pH Area_Hood_Canal]`)) %>% 
  arrange(desc(`b[pH Area_Hood_Canal]`)) %>%
  head(15)

a %>%
  filter(taxa == "Noelaerhabdaceae|Emiliania|Emiliania huxleyi") %>%
  ModelPlot(variable = "Temperature")


# b <- data.for.stan %>%
#   filter(taxa == "Gonyaulacaceae|Alexandrium|NA")
# 
# summary(glm(binary ~ pH, data = b, family = "binomial"))
# 
# 
# b %>%
# ggplot(aes(x = pH, y = binary, color = Area)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = F) +
#   facet_grid(Area~Season)




