#############################################################################################################################################################
# These analyses correspond to the analyses presented in  "A global meta-analysis of temperature effects on marine fishes' digestion across trophic groups" #         
# Code used to create specific tables or figures in the main text or Appendix 2 are referenced.                                                             #
# All figures/tables from the Appendix start with S, e.g. Table S3.2.  Main text figures/tables have no S, e.g. Figure 2.                                   #                                                                                     
#############################################################################################################################################################

#Load packages and set seed

library(brms)
library(ggplot2)
library(ggpubr)
library(bayesplot)
library(reshape2)
library(tidybayes)
library(ggpubr)
library(dplyr)
library(loo)
library(scales)
library(plotrix)

set.seed(438)  #Final model results are from models run with this set.seed()

######################
#3. Gut passage time #
######################

GP.dat <- read.csv("data_gut_passage.csv")


#Models are saved so that they don't need to be rerun for every session.
GP.0 <- readRDS("GP0.rds")
GP.1 <- readRDS("GP1.rds")
GP.2 <- readRDS("GP2.rds")


##############
# Table S3.2 #
##############

GP.priors <- c(prior(normal(5, 1), class = "Intercept"), #Corresponds to the log-transformed longest recorded gut passage time
               prior(normal(0, 1), class = "b"),         #For fixed effects
               prior(exponential(1), class = "sd"),      #For random effects
               prior(gamma(0.01, 0.01),class = "shape")) #brms default; appears to behave well   

##############
# Model GP.0 #
##############

#Prior predictive checks

GP.0_prior <-brm(passage.time.h ~ temperature.c + (1|species),
                 data = GP.dat,
                 family = Gamma(link = "log"),
                 prior = GP.priors,
                 sample_prior = "only",
                 cores = 4,
                 control = list(adapt_delta = 0.999))
GP.0_prior
#See here for an explanation on divergent priors in prior predictive checks.
#https://discourse.mc-stan.org/t/meaning-of-divergences-in-prior-predictive-checks/10759/3
#It appears that despite divergent transitions in the priors, the subsequent models should perform fine

#Note that outcomes of prior predictive checks when nsamples = 11 may vary widely depending on the set.seed()

pp_check(GP.0_prior, type = "hist", nsamples = 11) #Prior predictive checks come out narrow in plots but prior model fit shows plenty of variance around 
                                                   #predictors, doesn't appear to be a problem based on model outcomes and posterior predictive checks


#Build model 

GP.0 <- brm(passage.time.h ~ temperature.c + (1|species),
                       data = GP.dat,
                       family = Gamma(link = "log"),
                       prior = GP.priors,
                       chains = 4, cores = 4,
                       iter = 8000, warmup = 4000)
       ##############
GP.0   # Table S3.4 #
       ##############

#Run LOO-PSIS cross-validation

GP.0_PK <- loo(GP.0)
GP.0_PK

#High Pareto K values plus p_loo < p suggests that the model has difficulty predicting individual observations because 
#there are relatively few observations per species (https://rdrr.io/github/stan-dev/loo/man/loo-glossary.html)
#This doesn't necessarily mean the model is misspecified, particularly given the good performance of the posterior
#predictive checks (see below; https://discourse.mc-stan.org/t/pareto-k-diagnostics-and-kfold-model-comparison/7409).
#However, estimates of random/group-level effects from the model will likely be sensitive to future data.

GP.0 <- add_criterion(GP.0, "loo", reloo = TRUE) #Use reloo prior to model weighting

                     ##############
GP.0$criteria$loo    # Table S3.3 #
                     ##############

#saveRDS(GP.0, file = "GP0.rds")  

#Posterior predictive checks

pp_check(GP.0, type="scatter_avg", nsamples=100)            #looks good                                       
pp_check(GP.0, type = "loo_pit_overlay", nsamples = NULL)   #looks good
pp_check(GP.0, type = "dens_overlay", nsamples = 99)        #looks good
pp_check(GP.0, type = "hist", nsamples = 11)                #looks good

##############
# Model GP.1 #
##############

#Prior predictive checks

GP.1_prior <-brm(passage.time.h ~ temperature.c + diet + (1|species),
                 data = GP.dat,
                 family = Gamma(link = "log"),
                 prior = GP.priors,
                 sample_prior = "only",
                 cores = 4,
                 control = list(adapt_delta = 0.999))

GP.1_prior #See GP.0 for explanation of divergent effects in priors

#Note that outcomes of prior predictive checks when nsamples = 11 may vary widely depending on the set.seed()

pp_check(GP.1_prior, type = "hist", nsamples = 11) #Prior predictive checks come out narrow in plots but prior model fit shows plenty of variance around 
                                                   #predictors, doesn't appear to be a problem based on model outcomes and posterior predictive checks


#Build model

GP.1 <- brm(passage.time.h ~ temperature.c + diet + (1|species),
                       data = GP.dat,
                       family = Gamma(link = "log"),
                       prior = GP.priors,
                       chains = 4, cores = 4,
                       iter = 8000, warmup = 4000)

       ##############
GP.1   # Table S3.4 #
       ##############

#Run LOO-PSIS cross-validation


GP.1_PK <- loo(GP.1) #See GP.0 for explanation of high Pareto-K values
GP.1_PK

GP.1 <- add_criterion(GP.1, "loo", reloo = TRUE) 

                    ##############
GP.1$criteria$loo   # Table S3.3 #
                    ##############


#saveRDS(GP.1, file = "GP1.rds")  

#Posterior predictive checks                                    
pp_check(GP.1, type="scatter_avg", nsamples=100)             #looks good         
pp_check(GP.1, type = "loo_pit_overlay",nsamples = NULL)     #looks good                                
pp_check(GP.1, type = "dens_overlay",nsamples = 99)          #looks good
pp_check(GP.1, type = "hist", nsamples = 11, binwidth = 10)  #looks good


##############
# Model GP.2 #
##############

#Prior predictive checks

GP.2_prior <-brm(passage.time.h ~ temperature.c * diet + (1|species),
                 data = GP.dat,
                 family = Gamma(link = "log"),
                 prior = GP.priors,
                 sample_prior = "only",
                 cores = 4,
                 control = list(adapt_delta = 0.999))
GP.2_prior #See GP.0 for explanation of divergent effects in priors

#Note that outcomes of prior predictive checks when nsamples = 11 may vary widely depending on the set.seed()

pp_check(GP.2_prior, type = "hist", nsamples = 11) #Priors capable of generating long tails but shouldn't be a problem


GP.2 <- brm(passage.time.h ~ temperature.c * diet + (1|species),
                       data = GP.dat,
                       family = Gamma(link = "log"),
                       prior = GP.priors,
                       chains = 4, cores = 4, 
                       iter = 8000, warmup = 4000)

       ############
GP.2   #Table S3.4#
       ############

#Run LOO-PSIS cross-validation

GP.2_PK <- loo(GP.2) #See GP.0 for explanation of high Pareto-K values
GP.2_PK

GP.2 <- add_criterion(GP.2, "loo", reloo = TRUE)

                    ##############
GP.2$criteria$loo   # Table S3.3 #
                    ##############

#saveRDS(GP.2, file = "GP2.rds")  

#Posterior predictive checks
pp_check(GP.2, type="scatter_avg",nsamples=100)                #looks good                             
pp_check(GP.2, type = "loo_pit_overlay",nsamples = NULL)       #looks good
pp_check(GP.2, type = "dens_overlay",nsamples = 99)            #looks good
pp_check(GP.2, type = "hist", nsamples = 11, binwidth = 10)    #looks good

############################################
# Gut passage model selection (Table S3.3) #
############################################

#Note: Because models with relatively small sample sizes can produce variable model weights, we reran models and model selection 10 times with a random seed 
#to quantify the variability around model weight estimates.  This was done for gut passage time, absorption efficiency and gut length models.

model_weights(GP.0, GP.1, GP.2, weights = "loo") #Compare model performance using LOO-PSIS as criterion

model_weights(GP.0, GP.1, GP.2, weights = "stacking") #Stacked model weights

#####################
# Assemble Figure 2 #
#####################


#Create predictor data across full temperature range for each diet classification
nd <- 
  tibble(temperature.c = seq(from = min(GP.dat$temperature.c), to = max(GP.dat$temperature.c), length.out = 100) %>% 
           rep(., times = 6), diet = c(rep("carnivore",100),rep("zooplankton",100),rep("diatoms",100),
                                            rep("seagrass",100), rep("macroalgae",100), rep("turf algae",100)), 
         species = rep(0, 600))

#Generate predictions from the best selected model (GP.2) corresponding to this new predictor data

f <-
  fitted(GP.2, newdata = nd,allow_new_levels=TRUE) %>%
  as_tibble() %>%
  bind_cols(nd)

#ln-transform estimates for ease of presentation
f$ln.Estimate<-log(f$Estimate)
f$ln.Q2.5 <- log(f$Q2.5)
f$ln.Q97.5 <- log(f$Q97.5)
f$diet <- factor(f$diet, levels = c("carnivore", "zooplankton", "macroalgae", "phytoplankton", "seagrass", "turf algae"))
f.sub<- subset(f, diet == "carnivore" | diet == "macroalgae") #only want to show 95% CI for carnivores and macroalgivores because no credible difference 
                                                              #between carnivores and trophic groups other than macroalgivore
#set colours, fill, shapes
GP_colours <- c("carnivore" = "red2", "macroalgae" = "steelblue2", "macroalgae (non-fermenter)" = "steelblue2", "macroalgae (fermenter)" = "black", "diatoms" = "steelblue2", "seagrass" = "steelblue2",
                "turf algae" = "steelblue2", "zooplankton" = "red2")
GP_fill <- c("carnivore" = "red2", "macroalgae" = "steelblue2", "macroalgae (non-fermenter)" = "steelblue2", "macroalgae (fermenter)" = "steelblue2", "diatoms" = "steelblue2", "seagrass" = "steelblue2",
                "turf algae" = "steelblue2", "zooplankton" = "red2")
GP_shapes <- c("carnivore" = 16, "macroalgae" = 16,"macroalgae (non-fermenter)" = 16, "macroalgae (fermenter)" = 21, "diatoms" = 0, "seagrass" = 1, "turf algae" = 2,"zooplankton" = 0)

xlab <- expression("Temperature " ( degree*C))

fig_2a <- ggplot(GP.dat,aes(x = temperature.c)) + 
  geom_smooth(data = f.sub,
              aes(y = ln.Estimate, ymin = ln.Q2.5, ymax = ln.Q97.5, color = diet, fill = diet),
              stat = "identity", 
              alpha = .15, size = 2, show.legend = FALSE) +
  geom_point(data = GP.dat, aes(y = log(passage.time.h), color = diet2, shape = diet2, fill = diet2),
             size = 5,alpha=.8, stroke = 1.5) +
  ylab("ln(gut passage time)") +xlab(xlab) + 
  scale_color_manual(values = GP_colours) + scale_fill_manual(values = GP_fill) + scale_shape_manual(values = GP_shapes) +
  theme_bw(base_size = 18) + labs(colour= "") +labs(shape= "") + theme(plot.margin = unit(c(0.1,0.1,0.5,0.1), "cm")) +
  theme(legend.position = c(0.2, 0.25)) + theme(legend.background = element_blank()) + guides(fill = FALSE)

fig_2a  #identified non-fermenters manually in powerpoint (only two)

#Plot estimated diet and temperature*diet effect sizes using posterior draws

GP.2_post <- posterior_samples(GP.2)
GP.2_post <- GP.2_post[,c(3:12)]
colnames(GP.2_post) <- c("Diatoms", "Macroalgae", 
                         "Seagrass", "Turf algae", "Zooplankton", "Diatoms-temperature", 
                         "Macroalgae-temperature", "Seagrass-temperature","Turf algae-temperature", "Zooplankton-temperature")
GP.2_post <- melt(GP.2_post)
colnames(GP.2_post) <- c("Parameter", "Value")
GP.2_post <-subset(GP.2_post, Parameter !="Intercept")
GP.2_post <- droplevels(GP.2_post)
GP.2_post$Group <- c(rep("Diet", 80000), rep("Temperature",80000))

GP.2_postdiet <- subset(GP.2_post, Group == "Diet") #separate into diet and temperature*diet interaction effects because scales very different
GP.2_posttemperature <- subset(GP.2_post, Group == "Temperature")

fig_2b <- ggplot(data = GP.2_postdiet, aes(y = Value, x = Parameter)) +
  stat_pointinterval(orientation = "vertical", .width = 0.95) + 
  geom_hline(yintercept = 0, linetype = "longdash", size =1 ) +
  theme_bw(base_size = 18) + 
  ylab("Diet\neffect size") + xlab("") + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) + scale_y_continuous(label = label_number(accuracy = .01))
fig_2b

GP.2_posttemperature$Parameter <- recode(GP.2_posttemperature$Parameter, "Diatoms-temperature" = "Diatoms", 
                                         "Macroalgae-temperature" = "Macroalgae", "Seagrass-temperature" = "Seagrass",
                                         "Turf algae-temperature" = "Turf algae","Zooplankton-temperature" = "Zooplankton")

fig_2c <- ggplot(data = GP.2_posttemperature, aes(y = Value, x = Parameter)) +
  stat_pointinterval(orientation = "vertical", .width = 0.95) + 
  geom_hline(yintercept = 0, linetype = "longdash", size =1 ) +
  theme_bw(base_size = 18) + 
  ylab("Temperature*diet\neffect size") + xlab("") + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
fig_2c

#Put all plots together
ggarrange(fig_2a, 
          ggarrange(fig_2b, fig_2c, ncol = 2, labels = c("b", "c")),
          nrow = 2, heights = c(450, 300), labels ="a")

#############################
# 4. Absorption efficiency  #
#############################

AE.dat <- read.csv("data_absorption_efficiency.csv")

############################################################################
# In-text values from Absorption efficiency section of results (Main text) #
############################################################################

aggregate(AE.dat$absorption.efficiency, by = list(AE.dat$component), mean)
aggregate(AE.dat$absorption.efficiency, by = list(AE.dat$component), sd)

#Setting invertebrate diets and total absorption efficiency as the intercepts for ease of comparison

AE.dat$diet <- as.factor(AE.dat$diet)
AE.dat <- within(AE.dat, diet <- relevel(diet, ref = "invertebrate")) 
AE.dat$component <- as.factor(AE.dat$component)
AE.dat <- within(AE.dat, component <- relevel(component, ref = "total"))

#Models are saved so that they don't need to be rerun for every session

AE.0 <- readRDS("AE0.rds")
AE.1 <- readRDS("AE1.rds")
AE.2 <- readRDS("AE2.rds")
AE.3 <- readRDS("AE3.rds")

##############
# Table S4.2 #
##############

AE.priors <- c(prior(normal(0, 1), class = "Intercept"), #Corresponds to logit-transformed %50 absorption efficiency
               prior(normal(0, 1), class = "b"), #For fixed effects
               prior(exponential(1), class = "sd"), #For random effects
               prior(exponential(1),class = "phi")) #Precision parameter

##################
### Model AE.0 ###
##################

#Prior predictive check

AE.0_prior <- brm(absorption.efficiency ~ component + method + (1|species) + (1|study.id),
                data = AE.dat, 
                family = Beta(link = "logit"),
                prior = AE.priors,
                cores = 4,
                sample_prior = "only",
                control = list(adapt_delta = 0.999, max_treedepth = 15))

AE.0_prior #See GP.0 for explanation of divergent effects in priors

#Note that outcomes of prior predictive checks when nsamples = 11 may vary widely depending on the set.seed()

pp_check(AE.0_prior, type = "hist", nsamples = 11) #Priors look fine

#Build model 
AE.0 <- brm(absorption.efficiency ~ component + method + (1|species) + (1|study.id), 
            data = AE.dat,
            family = Beta(link = "logit"),
            prior = AE.priors,
            chains = 4, cores = 4, 
            iter = 8000, warmup = 4000, 
            control = list(adapt_delta = 0.99))


        ##############
AE.0   # Table S4.4 #
        ############## 

#Run LOO-PSIS cross-validation

AE.0_PK <- loo(AE.0) 
AE.0_PK

AE.0 <- add_criterion(AE.0, "loo", reloo = TRUE)

                     ##############
AE.0$criteria$loo    # Table S4.3 #
                     ##############

#saveRDS(AE.0, file = "AE0.rds")  

#Posterior predictive checks
pp_check(AE.0, type="scatter_avg",nsamples=100)                #looks good                             
pp_check(AE.0, type = "loo_pit_overlay",nsamples = NULL)       #not perfect
pp_check(AE.0, type = "dens_overlay",nsamples = 99)            #looks good
pp_check(AE.0, type = "hist", nsamples = 11)                   #looks ok


##################
### Model AE.1 ###
##################

#Prior predictive check

AE.1_prior<-brm(absorption.efficiency ~ diet + component + method + (1|species) + (1|study.id),
                data = AE.dat, 
                family = Beta(link = "logit"),
                prior = AE.priors,
                cores = 4,
                sample_prior = "only",
                control = list(adapt_delta = 0.999, max_treedepth = 15))
AE.1_prior #See GP.0 for explanation of divergent effects in priors

#Note that outcomes of prior predictive checks when nsamples = 11 may vary widely depending on the set.seed()

pp_check(AE.1_prior, type = "hist", nsamples = 11) #Priors look fine

#Build model 
AE.1 <- brm(absorption.efficiency ~ diet + component + method + (1|species) + (1|study.id), 
            data = AE.dat,
            family = Beta(link = "logit"),
            prior = AE.priors,
            chains = 4, cores = 4, 
            iter = 8000, warmup = 4000, 
            control = list(adapt_delta = 0.99))

          ##############
AE.1      # Table S4.4 #
          ##############

#Run LOO-PSIS cross-validation

AE.1_PK <- loo(AE.1)
AE.1_PK

AE.1 <- add_criterion(AE.1, "loo", reloo = TRUE)

                     ##############
AE.1$criteria$loo    # Table S4.3 #
                     ##############

#saveRDS(AE.1, file = "AE1.rds")  

#Posterior predictive checks
pp_check(AE.1, type="scatter_avg",nsamples=100)                #looks good                             
pp_check(AE.1, type = "loo_pit_overlay",nsamples = NULL)       #looks good
pp_check(AE.1, type = "dens_overlay",nsamples = 99)            #looks good
pp_check(AE.1, type = "hist", nsamples = 11)                   #looks good


##################
### Model AE.2 ###
##################

#Prior predictive check

AE.2_prior<-brm(absorption.efficiency ~ temperature.c + diet + component + method + (1|species) + (1|study.id),
                data = AE.dat, 
                family = Beta(link = "logit"),
                prior = AE.priors,
                cores = 4,
                sample_prior = "only",
                control = list(adapt_delta = 0.999, max_treedepth = 15))

#Note that outcomes of prior predictive checks when nsamples = 11 may vary widely depending on the set.seed()
AE.2_prior

pp_check(AE.2_prior, type = "hist", nsamples = 11) #Priors look fine


#Build model 
AE.2 <- brm(absorption.efficiency ~ temperature.c + diet + component + method + (1|species) + (1|study.id), 
            data = AE.dat,
            family = Beta(link = "logit"),
            prior = AE.priors,
            chains = 4, cores = 4, 
            iter = 8000, warmup = 4000, 
            control = list(adapt_delta = 0.99))

           ##############
AE.2       # Table S4.4 #
           ##############

#Run LOO-PSIS cross-validation

AE.2_PK <- loo(AE.2) 
AE.2_PK

AE.2 <- add_criterion(AE.2, "loo", reloo = TRUE)

                     ##############
AE.2$criteria$loo    # Table S4.3 #
                     ##############

#saveRDS(AE.2, file = "AE2.rds")  


#Posterior predictive checks
pp_check(AE.2, type="scatter_avg",nsamples=100)                #looks good                             
pp_check(AE.2, type = "loo_pit_overlay",nsamples = NULL)       #looks good
pp_check(AE.2, type = "dens_overlay",nsamples = 99)            #looks good
pp_check(AE.2, type = "hist", nsamples = 11)                   #looks good

###########


AE.3_prior<-brm(absorption.efficiency ~ temperature.c * diet + component + method + (1|species) + (1|study.id),
                data = AE.dat, 
                family = Beta(link = "logit"),
                prior = AE.priors,
                cores = 4,
                sample_prior = "only",
                control = list(adapt_delta = 0.999, max_treedepth = 15))
AE.3_prior

#Note that outcomes of prior predictive checks when nsamples = 11 may vary widely depending on the set.seed()
pp_check(AE.3_prior, type = "hist", nsamples = 11) #Priors look fine


#Build model 
AE.3 <- brm(absorption.efficiency ~ temperature.c * diet + component + method + (1|species) + (1|study.id), 
            data = AE.dat,
            family = Beta(link = "logit"),
            prior = AE.priors,
            chains = 4, cores = 4, 
            iter = 8000, warmup = 4000, 
            control = list(adapt_delta = 0.99))

          ############## 
AE.3      # Table S4.4 #
          ##############

#Run LOO-PSIS cross-validation

AE.3_PK <- loo(AE.3)
AE.3_PK

AE.3 <- add_criterion(AE.3, "loo", reloo = TRUE)


                     ##############
AE.3$criteria$loo    # Table S4.3 #
                     ##############

#saveRDS(AE.3, file = "AE3.rds")  

#Posterior predictive checks
pp_check(AE.3, type="scatter_avg",nsamples=100)                #looks good                             
pp_check(AE.3, type = "loo_pit_overlay",nsamples = NULL)       #looks ok
pp_check(AE.3, type = "dens_overlay",nsamples = 99)            #looks good
pp_check(AE.3, type = "hist", nsamples = 11)                   #looks good

######################################################
# Absorption efficiency model selection (Table S4.3) #
######################################################

#Note: Because models with relatively small sample sizes can produce variable model weights, we reran models and model selection 10 times with a random seed 
#to quantify the variability around model weight estimates.  This was done for gut passage time, absorption efficiency and gut length models.

model_weights(AE.0, AE.1, AE.2, AE.3, weights = "loo") #Compare model performance using LOO-PSIS as criterion

model_weights(AE.0, AE.1, AE.2, AE.3, weights = "stacking") #Stacked model weights

##################
# Build Figure 3 #
##################

#First take posterior draws of all parameter estimates and reformat so ready for ggplot

AE.1_post <- posterior_samples(AE.1)
AE.1_post <- AE.1_post[, 2:14]
colnames(AE.1_post) <- c("Diatoms", "Fish", "Macroalgae", "Seagrass", "Turf algae", "Zooplankton", 
                         "Carbohydrate", "Carbon", "Energy", "Lipid", "Nitrogen", "Organic", "Protein")
AE.1_post <- melt(AE.1_post)
AE.1_post$Group <- c(rep("Diet", 96000), rep("Component", 112000))
colnames(AE.1_post) <- c("Parameter", "Value", "Group")

#Set colours for plot

model.colours_data <- c("herbivore" = "steelblue2", "carnivore" = "red2")

fig_3a <- ggplot(data = AE.dat, aes(x = component, y = absorption.efficiency)) + geom_boxplot(aes(fill = trophic.group)) +
  theme_bw(base_size = 18) + 
  scale_fill_manual(values = model.colours_data) + 
  xlab("") + ylab("Absorption efficiency") +
  theme(legend.position = c(0.9,0.15)) + theme(legend.background = element_blank()) + labs(fill = "") +
  theme(plot.margin = unit(c(0.2, 0.1, 0, 0.1), "cm"))
fig_3a

AE.1_postdiet <- subset(AE.1_post, Group == "Diet")
AE.1_postcomponent <- subset(AE.1_post, Group == "Component")

fig_3b <- ggplot(data = AE.1_postdiet, aes(y = Value, x = Parameter)) +
  stat_pointinterval(orientation = "vertical", .width = 0.95) + 
  geom_hline(yintercept = 0, linetype = "longdash") +
  theme_bw(base_size = 18) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ylab("Diet effect size") + xlab("") +
  theme(plot.margin = unit(c(0.1, 0.1, -0.3, 0.1), "cm")) +
  theme(legend.position = "none") + scale_y_continuous(label = label_number(accuracy = .1))
fig_3b


fig_3c <- ggplot(data = AE.1_postcomponent, aes(y = Value, x = Parameter)) +
  stat_pointinterval(orientation = "vertical", .width = 0.95) + 
  geom_hline(yintercept = 0, linetype = "longdash") +
  theme_bw(base_size = 18) + theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  ylab("Component effect size") + xlab("") +
  theme(plot.margin = unit(c(0.1, 0.1, -0.3, 0.1), "cm")) +
  theme(legend.position = "none")
fig_3c

ggarrange(fig_3a, 
          ggarrange(fig_3b, fig_3c, ncol = 2, labels = c("b", "c")),
                    nrow = 2, labels = "a")
###############
# Figure S4.1 #
###############

xlab <- expression("Temperature " ( degree*C))

AE.dat$diet <- factor(AE.dat$diet, levels=c("diatoms", "macroalgae", "seagrass", "turf algae", "fish", "invertebrate", "zooplankton"))

fig.s4.1 <- ggplot(AE.dat, aes(x = temperature.c, y = absorption.efficiency)) + 
            geom_point(size = 3, alpha = 0.5, aes(colour = trophic.group)) +
            facet_grid(rows = vars(component), cols = vars(diet)) + theme_bw(base_size = 15) + xlab(xlab) + ylab("Absorption efficiency") +
            guides(colour=guide_legend(title = "Trophic group"))
fig.s4.1


##########################
# 5. Gut length analyses #
##########################

GL.dat<-read.csv("data_gut_length.csv")

GL.dat <- subset(GL.dat, diet != "seagrass") #not enough data to include seagrass consumers as a separate category
GL.dat<-subset(GL.dat,diet!="")              #couldn't classify specific diet for two species
GL.dat <- subset(GL.dat, preferred == "Y")   #if a species had multiple observations of gut length, I excluded worse observations (e.g., measurements of the 
                                             #intestine and not the entire gut, no temperature observation in paper) and then randomly selected from the best/
                                             #most complete observations (next line of code)
GL.dat <- distinct(GL.dat, species, .keep_all = TRUE) 

GL.dat$diet <- as.factor(GL.dat$diet) #set fish consuming invertebrates as reference level (i.e., intercept)
GL.dat <- within(GL.dat, diet <- relevel(diet, ref = "invertebrates"))


#################################################################
# In-text values from Gut length section of results (Main text) #
#################################################################

aggregate(GL.dat$relativelength.mean, by = list(GL.dat$trophic.group), mean)
aggregate(GL.dat$relativelength.mean, by = list(GL.dat$trophic.group), sd)

aggregate(GL.dat$relativelength.mean, by = list(GL.dat$diet), mean)
aggregate(GL.dat$relativelength.mean, by = list(GL.dat$diet), sd)

#Models are saved so that they don't need to be rerun for every session.
GL.0 <- readRDS("GL0.rds")
GL.1 <- readRDS("GL1.rds")
GL.2 <- readRDS("GL2.rds")
GL.3 <- readRDS("GL3.rds")


##############
# Table S4.2 #
##############

GL.priors <- c(prior(normal(0, 2), class = "Intercept"), 
               prior(normal(0, 1), class = "b"), #For fixed effects
               prior(exponential(1), class = "sd"), #For random effects
               prior(exponential(1), class = "sigma")) #Variance

##############
# Model GL.0 #
##############

#Prior predictive checks

GL.0_prior <- brm(gutlength.mean.mm ~ ln.bodylength + (1|order) + (1|study.id), 
                      data = GL.dat,
                      family = lognormal(), 
                      prior = GL.priors,
                      cores = 4, 
                      sample_prior = "only", control = list(adapt_delta = 0.9, max_treedepth = 15))
GL.0_prior

pp_check(GL.0_prior, type = "hist", nsamples = 11)
#looks fine. priors allow for a distribution with a long tail, but that doesn't seem to affect the final models

#Build model 

GL.0 <- brm(gutlength.mean.mm ~ ln.bodylength + (1|order) + (1|study.id), 
            data = GL.dat,
            family = lognormal(), 
            prior = GL.priors,
            chains = 4, cores = 4, 
            iter = 8000, warmup = 4000, 
            control = list(adapt_delta = 0.9))

          ##############
GL.0      # Table S5.4 #
          ##############

#Run LOO-PSIS cross-validation
          
GL.0_PK <- loo(GL.0)     
GL.0_PK     

GL.0 <- add_criterion(GL.0, "loo", reloo = TRUE)

                    ##############
GL.0$criteria$loo   # Table S5.3 #
                    ##############

#saveRDS(GL.0, file = "GL0.rds")  

#Posterior predictive checks                                    
pp_check(GL.0, type="scatter_avg", nsamples=100)             #looks good         
pp_check(GL.0, type = "loo_pit_overlay",nsamples = NULL)     #looks good                                
pp_check(GL.0, type = "dens_overlay",nsamples = 99)          #looks good
pp_check(GL.0, type = "hist", nsamples = 11)                 #looks good

##############
# Model GL.1 #
##############

#Prior predictive checks

GL.1_prior <-brm(gutlength.mean.mm ~ diet + ln.bodylength + (1|order) + (1|study.id),
                 data = GL.dat,
                 family = lognormal(),
                 prior = GL.priors,
                 sample_prior = "only",
                 control = list(adapt_delta = 0.999))
GL.1_prior 

pp_check(GL.1_prior, type = "hist", nsamples = 11)
#looks fine. priors allow for a distribution with a long tail, but that doesn't seem to affect the final models

#Build model

GL.1 <- brm(gutlength.mean.mm ~ diet + ln.bodylength + (1|order) + (1|study.id),
            data = GL.dat,
            family = lognormal(),
            prior = GL.priors,
            chains = 4, cores = 4,
            iter = 8000, warmup = 4000, 
            control = list(adapt_delta = 0.9),
            save_all_pars = TRUE)

         ##############
GL.1     # Table S5.4 #
         ##############

#Run LOO-PSIS cross-validation

GL.1_PK <- loo(GL.1)            
GL.1_PK

GL.1 <- add_criterion(GL.1, "loo", reloo = TRUE)

                     ##############
GL.1$criteria$loo    # Table S5.3 #
                     ##############

#saveRDS(GL.1, file = "GL1.rds")  

#Posterior predictive checks                                    
pp_check(GL.1, type="scatter_avg", nsamples=100)             #looks good         
pp_check(GL.1, type = "loo_pit_overlay",nsamples = NULL)     #looks good                                
pp_check(GL.1, type = "dens_overlay",nsamples = 99)          #looks good
pp_check(GL.1, type = "hist", nsamples = 11)                 #looks good

##############
# Model GL.2 #
##############

#Prior predictive checks

GL.2_prior <-brm(gutlength.mean.mm ~ temperature.c + diet + ln.bodylength + (1|order) + (1|study.id),
                 data = GL.dat,
                 family = lognormal(),
                 prior = GL.priors,
                 sample_prior = "only",
                 control = list(adapt_delta = 0.999))
GL.2_prior

pp_check(GL.2_prior, type = "hist", nsamples = 11)
#looks fine. priors allow for a distribution with a long tail, but that doesn't seem to affect the final models


GL.2 <- brm(gutlength.mean.mm ~ temperature.c + diet + ln.bodylength + (1|order) + (1|study.id),
            data = GL.dat, 
            family = lognormal(),
            prior = GL.priors,
            chains = 4, cores = 4, 
            iter = 8000, warmup = 4000,
            control = list(adapt_delta = 0.9),
            save_all_pars = TRUE)

         ##############
GL.2     # Table S5.4 #
         ##############

#Run LOO-PSIS cross-validation

GL.2_PK <- loo(GL.2)    
GL.2_PK

GL.2 <- add_criterion(GL.2, "loo", reloo = TRUE)

                     ##############
GL.2$criteria$loo    # Table S5.3 #
                     ##############

#saveRDS(GL.2, file = "GL2.rds")  

#Posterior predictive checks                                    
pp_check(GL.2, type="scatter_avg", nsamples=100)             #looks good         
pp_check(GL.2, type = "loo_pit_overlay",nsamples = NULL)     #looks good                                
pp_check(GL.2, type = "dens_overlay",nsamples = 99)          #looks good
pp_check(GL.2, type = "hist", nsamples = 11)                 #looks good

##############
# Model GL.3 #
##############

GL.3_prior <-brm(gutlength.mean.mm ~ temperature.c * diet + ln.bodylength + (1|order) + (1|study.id),
                 data = GL.dat,
                 family = lognormal(),
                 prior = GL.priors,
                 sample_prior = "only",
                 control = list(adapt_delta = 0.999))
GL.3_prior

pp_check(GL.3_prior, type = "hist", nsamples = 11)
#looks fine. priors allow for a distribution with a long tail, but that doesn't seem to affect the final models

#Build model

GL.3 <- brm(gutlength.mean.mm ~ temperature.c * diet + ln.bodylength + (1|order) + (1|study.id),
            data = GL.dat, 
            family = lognormal(),
            prior = GL.priors,
            chains = 4, cores = 4, 
            iter = 8000, warmup = 4000,
            control = list(adapt_delta = 0.9),
            save_all_pars = TRUE)

         ##############
GL.3     # Table S5.4 #
         ##############

#Run PSIS-LOO cross-validation

GL.3_PK <- loo(GL.3) 
GL.3_PK

GL.3 <- add_criterion(GL.3, "loo", reloo = TRUE)

                     ##############
GL.3$criteria$loo    # Table S5.3 #
                     ##############


#saveRDS(GL.3, file = "GL3.rds")  

#Posterior predictive checks                                    
pp_check(GL.3, type="scatter_avg", nsamples=100)             #looks good         
pp_check(GL.3, type = "loo_pit_overlay",nsamples = NULL)     #looks good                                
pp_check(GL.3, type = "dens_overlay",nsamples = 99)          #looks good
pp_check(GL.3, type = "hist", nsamples = 11)                 #looks good

###########################################
# Gut length model selection (Table S5.3) #
###########################################

#Note: Because models with relatively small sample sizes can produce variable model weights, we reran models and model selection 10 times with a random seed 
#to quantify the variability around model weight estimates.  This was done for gut passage time, absorption efficiency and gut length models.

model_weights(GL.0, GL.1, GL.2, GL.3, weights = "loo")   #Compare model performance using LOO-PSIS as criterion
model_weights(GL.0, GL.1, GL.2, GL.3, weights = "stacking")  #Compare stacked model weights


#####################
# Assemble Figure 4 #
#####################

#Generate predictions for effect of temperature on gut length
GL.2_temperature <- conditional_effects(GL.2, effects = c("temperature.c"), re_formula = NULL, method = "predict", transform = log)
GL.2_eff <- GL.2_temperature[[1]]


GL_colours <- c("carnivore" = "red2", "herbivore" = "steelblue2", "omnivore" = "purple4")

xlab <- expression("Temperature " (degree*C))


fig.4a <- ggplot(data = GL.dat, aes(x = temperature.c, y = ln.gutlength.mean)) +
  geom_smooth(data = GL.2_eff, aes(ymin = lower__, ymax = upper__, y = estimate__ ), colour = "black", fill = "snow3", stat = "identity") +
  theme_bw(base_size = 18) + scale_colour_manual(values = GL_colours) +
  geom_point(data = GL.dat, aes(x = temperature.c, y = ln.gutlength.mean, colour = trophic.group), size = 5, alpha = 0.7, inherit.aes = FALSE) +
  xlab(xlab) + ylab("\nln(gut length)") + theme(plot.margin = unit(c(0.1,0.1,0.3,0.2), "cm")) +
  scale_y_continuous(label = label_number(accuracy = .1)) +
  theme(legend.position = c(0.07, 0.9)) + theme(legend.background = element_blank()) +labs(colour = "")
fig.4a

GL.dat$diet<- factor(GL.dat$diet, levels = c("invertebrates", "invertebrates and fish", "fish", "zooplankton", "omnivore","diatoms", "herbivore-detritivore", "macroalgae", "turf algae"))

fig.4b <- ggplot(GL.dat, aes(x = diet, y = relativelength.mean)) + geom_boxplot(aes(fill = trophic.group)) + ylab("\nRGL") + 
  theme_bw(base_size = 18) + scale_fill_manual(values = GL_colours) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(plot.margin = unit(c(0.1,0.1,0.1,0.2), "cm")) + labs(colour = "") + xlab("") + theme(legend.position = "none")
fig.4b

#GL_habitat.plot <-ggplot(GL.dat, aes(x = habitat2, y = RGL)) + geom_boxplot(aes()) + ylab("RGL") + 
#  theme_bw(base_size =14) + scale_fill_manual(values = GL_colours) + xlab("") +
#  theme(plot.margin = unit(c(0.1,0.1,0.1,0.1), "cm")) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
#GL_habitat.plot

#Plot estimated effect sizes using posterior draws

GL.2_post <- posterior_samples(GL.2)
GL.2_post <- GL.2_post[,3:10]

colnames(GL.2_post) <- c("diatoms", "fish", "herbivore-detritivore", "invertebrates and fish", "macroalgae", 
                         "omnivore","turf algae", "zooplankton")
GL.2_post <- melt(GL.2_post)

colnames(GL.2_post) <- c("Parameter", "Value")

GL.2_post$Parameter <- ordered(GL.2_post$Parameter, levels = c("invertebrates and fish", "fish", "zooplankton", "omnivore", "diatoms", "herbivore-detritivore",
                                                                       "macroalgae", "turf algae"))

fig.4c <-ggplot(data = GL.2_post, aes(x = Parameter, y = Value)) +
  stat_pointinterval(orientation = "vertical", .width = 0.95) + 
  theme_bw(base_size = 18) + geom_hline(yintercept = 0, linetype = "longdash") +
  xlab("") + ylab("Diet\neffect size") + theme(plot.margin = unit(c(0.1,0.1,0.1,0.2), "cm")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
fig.4c

#Combine all three plots
ggarrange(fig.4a, fig.4b, fig.4c, labels = c("a", "b", "c"),
          nrow = 3)

#####################################################
# Estimated model weights (Tables S3.3, S4.3, S5.3) #
#####################################################

#The models in the preceding code were run 10 times using a random seed for each run.  This was done because the relatively small sample sizes of the models
#can lead to variability in the model weight estimates (https://projecteuclid.org/euclid.ba/1516093227).  Here we calculate the mean and standard error of
#those estimates.

model_weights <- read.csv("model_weights.csv")

aggregate(model_weights$LOO.weight, by = list(model_weights$Model), mean)
aggregate(model_weights$LOO.weight, by = list(model_weights$Model), std.error)


aggregate(model_weights$Stacking.weight, by = list(model_weights$Model), mean)
aggregate(model_weights$Stacking.weight, by = list(model_weights$Model), std.error)

###########################
# 6. Nutrient composition #
###########################

NUT.dat <- read.csv("data_nutrients.csv")

NUT.dat <- subset(NUT.dat, diet.details != "unspecified") #exclude food items that don't have a clear identification
NUT.dat <- subset(NUT.dat, diet.details != "gut contents") #exclude food items that don't have a clear identification


NUT.priors <- c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(gamma(0.01, 0.01), class = phi))


NUT.carbohydrate <- readRDS("nutrient-carbohydrate.rds")
NUT.carbon <- readRDS("nutrient-carbon.rds")
NUT.energy <- readRDS("nutrient-energy.rds")
NUT.lipid <- readRDS("nutrient-lipid.rds")
NUT.nitrogen <- readRDS("nutrient-nitrogen.rds")
NUT.organic <- readRDS("nutrient-organic.rds")
NUT.protein <- readRDS("nutrient-protein.rds")

#Prior predictive checks (will be the same for all models except energy since sample_prior = "only")
#Note that outcomes of prior predictive checks when nsamples = 11 may vary widely depending on the set.seed()

NUT.prior.test <-brm(protein ~ diet.type, 
                     data = NUT.dat, 
                     family=Beta(link = "logit"),
                     chains = 4, cores = 4,
                     prior = NUT.priors,
                     sample_prior = "only",
                     iter = 8000, warmup = 4000,
                     control = list(adapt_delta = 0.999, max_treedepth = 15))

pp_check(NUT.prior.test, type = "hist", nsamples = 11)  #Priors look fine

#Build models

NUT.carbohydrate <- brm(carbohydrate ~ diet.type, 
                        data = NUT.dat, 
                        family=Beta(link = "logit"),
                        chains = 4, cores = 4,
                        prior = NUT.priors,
                        iter = 8000, warmup = 4000)

LOO.carbohydrate <- loo(NUT.carbohydrate)
LOO.carbohydrate

#saveRDS(NUT.carbohydrate, file = "nutrient-carbohydrate.rds")

NUT.carbon <- brm(carbon ~ diet.type, 
                  data = NUT.dat, 
                  family=Beta(link = "logit"),
                  chains = 4, cores = 4,
                  prior = NUT.priors,
                  iter = 8000, warmup = 4000)

LOO.carbon <- loo(NUT.carbon)
LOO.carbon

#saveRDS(NUT.carbon, file = "nutrient-carbon.rds")

NUT.lipid <- brm(lipid ~ diet.type,
                 data = NUT.dat, 
                 family = Beta(link = "logit"),
                 chains = 4, cores = 4,
                 prior = NUT.priors,
                 iter = 8000, warmup = 4000)

LOO.lipid <- loo(NUT.lipid)
LOO.lipid

#saveRDS(NUT.lipid, file = "nutrient-lipid.rds")

NUT.organic <- brm(organic ~ diet.type,
                   data = NUT.dat, 
                   family = Beta(link = "logit"),
                   chains = 4, cores = 4,
                   prior = NUT.priors,
                   iter = 8000, warmup = 4000)

LOO.organic <- loo(NUT.organic)
LOO.organic

#saveRDS(NUT.organic, file = "nutrient-organic.rds")

NUT.protein <- brm(protein ~ diet.type, 
                  data = NUT.dat, 
                  family=Beta(link = "logit"),
                  chains = 4, cores = 4,
                  prior = NUT.priors,
                  iter = 8000, warmup = 4000)

LOO.protein <- loo(NUT.protein)
LOO.protein

#saveRDS(NUT.protein, file = "nutrient-protein.rds")

NUT.nitrogen <- brm(nitrogen ~ diet.type, 
                   data = NUT.dat, 
                   family=Beta(link = "logit"),
                   chains = 4, cores = 4,
                   prior = NUT.priors,
                   iter = 8000, warmup = 4000)

LOO.nitrogen <- loo(NUT.nitrogen)
LOO.nitrogen

#saveRDS(NUT.nitrogen, file = "nutrient-nitrogen.rds")

#Energy priors/model different because energy wasn't measured as a percentage

EN.priors <- c(prior(normal(13.9, 5), class = Intercept),      ##############
               prior(normal(0,5), class = b),                  # Table S6.1 #
               prior(exponential(1), class = sigma))           ##############


EN.prior.test <-brm(energy.kJg ~ diet.type, 
                     data = NUT.dat, 
                     family=gaussian(),
                     chains = 4, cores = 4,
                     prior = EN.priors,
                     sample_prior = "only",
                     iter = 8000, warmup = 4000)

#Note that outcomes of prior predictive checks when nsamples = 11 may vary widely depending on the set.seed()
pp_check(EN.prior.test, type = "hist", nsamples = 11)  #Priors look fine

NUT.energy<-brm(energy.kJg ~ diet.type, 
                data = NUT.dat, 
                family = gaussian(), 
                prior=EN.priors,
                chains = 4, cores = 4, 
                iter = 8000, warmup = 4000)

#saveRDS(NUT.energy, file = "nutrient-energy.rds")

##############
# Table S5.4 #
##############

NUT.carbohydrate
NUT.carbon
NUT.energy
NUT.lipid
NUT.nitrogen
NUT.organic
NUT.protein


NUT.melt <- melt(NUT.dat, id.vars = c("study.id","diet.type","diet.details","reference"))
NUT.melt$value <- as.numeric(NUT.melt$value)

EN.data <- subset(NUT.melt, variable == "energy.kJg")
EN.data$variable <- "energy"

NUT.melt <-subset(NUT.melt, variable != "energy.kJg")
NUT.melt <-subset(NUT.melt, variable != "energy.reported")
NUT.melt <-subset(NUT.melt, variable != "energy.reported.units")
NUT.melt <-subset(NUT.melt, variable != "ash")
NUT.melt <-subset(NUT.melt, variable != "Notes")


##############################################################################
# In-text values from Nutrient concentrations section of results (Main text) #
##############################################################################

aggregate(NUT.melt$value, by = list(NUT.melt$variable, NUT.melt$diet.type), mean, na.rm = TRUE)
aggregate(NUT.melt$value, by = list(NUT.melt$variable, NUT.melt$diet.type), sd, na.rm = TRUE)

aggregate(EN.data$value, by = list(EN.data$diet.type), mean, na.rm = TRUE)
aggregate(EN.data$value, by = list(EN.data$diet.type), sd, na.rm = TRUE)

#Pull posterior draws from each model to plot estimated effect size

protein_post <- posterior_samples(NUT.protein)
protein_post <- protein_post[2]
protein_post$component <- "protein"

lipid_post <- posterior_samples(NUT.lipid)
lipid_post <- lipid_post[2]
lipid_post$component <- "lipid"

carbohydrate_post <- posterior_samples(NUT.carbohydrate)
carbohydrate_post <- carbohydrate_post[2]
carbohydrate_post$component <- "carbohydrate"

nitrogen_post <- posterior_samples(NUT.nitrogen)
nitrogen_post <- nitrogen_post[2]
nitrogen_post$component <- "nitrogen"

carbon_post <- posterior_samples(NUT.carbon)
carbon_post <- carbon_post[2]
carbon_post$component <- "carbon"

organic_post <- posterior_samples(NUT.organic)
organic_post <- organic_post[2]
organic_post$component <- "organic"

energy_post <- posterior_samples(NUT.energy)
energy_post <- energy_post[2]
energy_post$component <- "energy"

NUT.est <- rbind(protein_post, lipid_post, carbohydrate_post, nitrogen_post, carbon_post, organic_post)
EN.est <- energy_post

colnames(NUT.est) <- c("Diet.estimate", "Component")
colnames(EN.est) <- c("Diet.estimate", "Component")

fig5c <-ggplot(data = NUT.est, aes(x = Component, y = Diet.estimate)) +
  stat_pointinterval(orientation = "vertical", .width = 0.95) + 
  theme_bw(base_size = 18) + geom_hline(yintercept = 0, linetype = "longdash") +
  xlab("") + ylab("Diet effect size") + theme(plot.margin = unit(c(0.1,0.1,0.1,0.2), "cm")) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))+ scale_y_continuous(label = label_number(accuracy = .1)) 
fig5c

fig5d <-ggplot(data = EN.est, aes(x = Component, y = Diet.estimate)) +
  stat_pointinterval(orientation = "vertical", .width = 0.95) + 
  theme_bw(base_size = 18) + geom_hline(yintercept = 0, linetype = "longdash") +
  xlab("") + ylab("") + theme(plot.margin = unit(c(0.1,0.1,0.1,0.2), "cm")) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
fig5d


NUT_fill <- c("animal" = "red2", "plant" = "steelblue2")

NUT.melt$variable <- factor(NUT.melt$variable, levels=c("carbohydrate", "carbon", "lipid", "nitrogen", "organic", "protein"))

fig5a <- ggplot(data = NUT.melt, aes(x = variable, y = value)) + geom_boxplot(aes(fill = diet.type)) +
  theme_bw(base_size =18) + scale_fill_manual(values = NUT_fill) + xlab("") + ylab("Proportion of total dry mass") +
  theme(legend.position = c(0.12, 0.95)) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  theme(legend.background = element_blank()) + labs(fill = "")

fig5a

fig5b <-ggplot(data = EN.data, aes(x = variable, y = value)) + geom_boxplot(aes(fill = diet.type)) +
  theme_bw(base_size = 18) + scale_fill_manual(values = NUT_fill) + xlab("") + ylab("kJ/g dry mass") + 
  theme(legend.position = "NONE") + theme(axis.text.x = element_text(angle = 20, hjust = 1)) 
fig5b

ggarrange(fig5a, fig5b, fig5c, fig5d, labels = c("a", "b", "c", "d"),
          nrow = 2, ncol = 2, widths = c(2, 1))

