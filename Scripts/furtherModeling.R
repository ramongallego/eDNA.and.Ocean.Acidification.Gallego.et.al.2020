library (lme4)
library (here)
library (tidyverse)
library (tidybayes)
library (rstanarm)
library (rstan)

data.for.stan <- read_csv(here("Input", "data.for.stan.csv"))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(108)
taxVec <- unique(data.for.stan$taxa)
subsetTax <- sample(taxVec, 20)

stan.input <- data.for.stan %>% 
  filter(taxa %in% subsetTax) %>% 
  group_by(Area) %>% 
  mutate(Temperature = vegan::decostand(Temperature, method = "standardize"),
         pH = vegan::decostand(pH, method = "standardize"),
         Salinity = vegan::decostand(Salinity, method = "standardize")) %>%
  group_by(taxa) %>%
  mutate(Noccur = sum(binary)) %>%
  filter(Noccur > 7 & Noccur < 69) %>%  #must be present or absent in at least 10% of 76 samples in dataset
  dplyr::select(-c(sample, phylum, DIC, Noccur)) %>% 
  ungroup() %>% 
  unite(Area, taxa, col = "AreaTaxa", remove = F) %>% 
  mutate(sampleIndex = 1:nrow(.),
         taxonIndex = match(taxa, unique(taxa)),
         areaIndex = match(Area, unique(Area)),
         seasonIndex = match(Season, unique(Season)),
         taxonAreaIndex = match(AreaTaxa, unique(AreaTaxa))) %>%   #unique combinations of taxon and area
  dplyr::select(-AreaTaxa) 


##rstanarm
      t_prior <- student_t(df = 7, location = 0, scale = 2.5)
      
      #model fitting one intercept per unique area/taxon combination, plus
      #a different slope for the effect of Temperature and of pH on each taxon.
      mod1.arm <- stan_glmer(binary ~ 0 + (0 + Temperature + pH | taxa) + (1 | taxonAreaIndex), 
                             data = stan.input, 
                             family = "binomial",
                             prior = t_prior,
                             prior_intercept = t_prior, 
                             chains = 2,
                             adapt_delta = 0.99)
      
      #model fitting one intercept per taxon, plus
      #a different slope for the effect of Temperature and of pH and of area on each taxon.
      mod2.arm <- stan_glmer(binary ~ 0 + (0 + Temperature + pH + Area| taxa) + (1 | taxa), 
                             data = stan.input, 
                             family = "binomial",
                             prior = t_prior,
                             prior_intercept = t_prior, 
                             chains = 2,
                             adapt_delta = 0.99)
      
      #model fitting one intercept per taxon/Area combination, plus
      #a different slope for the effect of Temperature and of pH and of Area on each taxon.
      mod3.arm <- stan_glmer(binary ~ 0 + (0 + Temperature + pH + Area| taxa) + (1 | taxonAreaIndex), 
                             data = stan.input, 
                             family = "binomial",
                             prior = t_prior,
                             prior_intercept = t_prior, 
                             chains = 2,
                             adapt_delta = 0.99)
      
      
      loo_compare(waic(mod1.arm),
                  waic(mod2.arm),
                  waic(mod3.arm))
      
      
      #now try one without Area , but with T and pH and Salinity, standardized across the whole dataset (not within area)
      stan.input.wholeSalish <- data.for.stan %>% 
        filter(taxa %in% subsetTax) %>% 
        #group_by(Area) %>% 
        mutate(Temperature = vegan::decostand(Temperature, method = "standardize"),
               pH = vegan::decostand(pH, method = "standardize"),
               Salinity = vegan::decostand(Salinity, method = "standardize")) %>%
        group_by(taxa) %>%
        mutate(Noccur = sum(binary)) %>%
        filter(Noccur > 7 & Noccur < 69) %>%  #must be present or absent in at least 10% of 76 samples in dataset
        dplyr::select(-c(sample, phylum, DIC, Noccur)) %>% 
        ungroup() %>% 
        unite(Area, taxa, col = "AreaTaxa", remove = F) %>% 
        mutate(sampleIndex = 1:nrow(.),
               taxonIndex = match(taxa, unique(taxa)),
               areaIndex = match(Area, unique(Area)),
               seasonIndex = match(Season, unique(Season)),
               taxonAreaIndex = match(AreaTaxa, unique(AreaTaxa))) %>%   #unique combinations of taxon and area
        dplyr::select(-AreaTaxa) 
      
      mod4.arm <- stan_glmer(binary ~ 0 + (0 + Temperature + pH | taxa) + (1 | taxonAreaIndex), 
                             data = stan.input.wholeSalish, 
                             family = "binomial",
                             prior = t_prior,
                             prior_intercept = t_prior, 
                             chains = 2,
                             adapt_delta = 0.99)
      mod5.arm <- stan_glmer(binary ~ 0 + (1 + Temperature + pH + Salinity | taxa), 
                             data = stan.input.wholeSalish, 
                             family = "binomial",
                             prior = t_prior,
                             prior_intercept = t_prior, 
                             chains = 2,
                             adapt_delta = 0.99)
      
      #OK, cool -- so there's something about area that isn't captured by T/Sal/pH, and so we want to keep area. 
      loo_compare(waic(mod4.arm),
                  waic(mod5.arm))

      
#Now, run best-fit model and export results      
      stan.input <- data.for.stan %>% 
        group_by(Area) %>% 
        mutate(Temperature = vegan::decostand(Temperature, method = "standardize"),
               pH = vegan::decostand(pH, method = "standardize"),
               Salinity = vegan::decostand(Salinity, method = "standardize")) %>%
        group_by(taxa) %>%
        mutate(Noccur = sum(binary)) %>%
        filter(Noccur > 7 & Noccur < 69) %>%  #must be present or absent in at least 10% of 76 samples in dataset
        dplyr::select(-c(sample, phylum, DIC, Noccur)) %>% 
        ungroup() %>% 
        unite(Area, taxa, col = "AreaTaxa", remove = F) %>% 
        mutate(sampleIndex = 1:nrow(.),
               taxonIndex = match(taxa, unique(taxa)),
               areaIndex = match(Area, unique(Area)),
             #  seasonIndex = match(Season, unique(Season)),
               taxonAreaIndex = match(AreaTaxa, unique(AreaTaxa))) %>%   #unique combinations of taxon and area
        dplyr::select(-AreaTaxa) 
      
      finalModel <- stan_glmer(binary ~ 0 + (0 + Temperature + pH | taxa) + (1 | taxonAreaIndex), 
                 data = stan.input, 
                 family = "binomial",
                 prior = t_prior,
                 prior_intercept = t_prior, 
                 chains = 4,
                 adapt_delta = 0.99)
      
      saveRDS(finalModel, file = here("ModelOutputs","finalModel_20200807.RDS" ))
      finalModel <- read_rds(here("ModelOutputs","finalModel_20200807.RDS"))
     
## And plot for sanity check     
      
      summary(finalModel$data == stan.input)
      
      stan.input$postpredict <- posterior_predict(finalModel) %>% colMeans()
      
      stan.input %>% 
        filter(taxa %in% unique(taxa)[1:5]) %>% 
        ggplot(aes(x = Temperature, y = binary)) +
          geom_point() +
          geom_point(aes(x = Temperature, y = postpredict), color = "blue") +
          facet_grid(taxa ~ Area)
      
##Full Stan
#I'm including the below (raw) Stan models in case they become useful later, but the rstanarm models above are far easier to work with for posterior predictions

# #model fitting one intercept per unique area/taxon combination, plus
# #a different slope for the effect of Temperature and of pH on each taxon.
# mod1 <- stan_model(file = "temperature_pH_by_taxon_int_by_AreaTaxon.stan")
# mod1.fit <- sampling(mod1, 
#      data = list(
#        N = nrow(stan.input),
#        Temperature = stan.input$Temperature[,1],
#        pH = stan.input$pH[,1],
#        Taxon = stan.input$taxonIndex,
#        Ntaxa = length(unique(stan.input$taxonIndex)),
#        AreaTaxon = stan.input$taxonAreaIndex,
#        NAreaTaxa = length(unique(stan.input$taxonAreaIndex)),
#        y = stan.input$binary
#      ),
#      chains = 4,
#      iter = 1000
#      )
# 
# #model fitting one intercept per taxon, plus
# #a different slope for the effect of Temperature and of pH and of area on each taxon.
# mod2 <- stan_model(file = "temperature_pH_Area_by_taxon.stan")
# mod2.fit <- sampling(mod2, 
#                      data = list(
#                        N = nrow(stan.input),
#                        Temperature = stan.input$Temperature[,1],
#                        pH = stan.input$pH[,1],
#                        Taxon = stan.input$taxonIndex,
#                        Ntaxa = length(unique(stan.input$taxonIndex)),
#                        Area = stan.input$areaIndex,
#                        Narea = length(unique(stan.input$areaIndex)),
#                        y = stan.input$binary
#                      ),
#                      chains = 4,
#                      iter = 1000
# )
# 
# #model fitting one intercept per taxon/Area combination, plus
# #a different slope for the effect of Temperature and of pH and of Area on each taxon.
# mod3 <- stan_model(file = "temperature_pH_Area_by_taxon_interceptTaxonArea.stan")
# mod3.fit <- sampling(mod3, 
#                      data = list(
#                        N = nrow(stan.input),
#                        Temperature = stan.input$Temperature[,1],
#                        pH = stan.input$pH[,1],
#                        Taxon = stan.input$taxonIndex,
#                        Ntaxa = length(unique(stan.input$taxonIndex)),
#                        Area = stan.input$areaIndex,
#                        Narea = length(unique(stan.input$areaIndex)),
#                        AreaTaxon = stan.input$taxonAreaIndex,
#                        NAreaTaxa = length(unique(stan.input$taxonAreaIndex)),
#                        y = stan.input$binary
#                      ),
#                      chains = 4,
#                      iter = 1000
#                      
# )
# 





