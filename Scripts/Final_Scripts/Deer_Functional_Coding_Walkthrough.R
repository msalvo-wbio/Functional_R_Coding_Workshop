
########## An Excericse in Functional Coding #############

# I believe functional coding has three main components:

# 1) Identify the task you are trying to accomplish
# 2) Identify the subtasks which need to be accomplished to complete the main task
# 3) Find the code needed to complete each subtask
# 4) Combine everything into a function


# Our example problem here will be that we have GPS collar data for white-tailed deer in Missouri.
# We just joined this project which was starting up, and the year before you joined and the 
# project started, Missouri Department of Conservation put a gps collar (same model) on a deer
# to begin collecting sample data.

# The broad research question which is being asked here is how are deer selecting for habitat in different
# parts of the state based on differing landscape configurations.

# In the overarching context of this, you can imagine a few tasks which may need to be accomplished.
# We will need to:

# 1) build home ranges for each deer
# 2) randomly generate available sites within the home range
# 3) Extract covariates at each point
# 4) Build multiple a priori models
# 5) summarise models

# You're job at the moment is to build code based on this one deer which can be used at the end of the
# project to do all of these things for all deer.

# load functions
source("Scripts/Final_Scripts/99_Deer_Functions.R")

# load libraries
ipak(c("terra", "sf", "sp", "adehabitatHR", "jagsUI", "reshape2"))

#load data
deer <- st_read("data/DeerN15002_Points/DeerN15002_Points.shp")
deer <- deer %>% dplyr::select(-tm_dffr, -anml_p_, -t, -id_yr)

head(deer)

nlcd <- rast('data/NLCD/Annual_NLCD_LndCov_2015_CU_C1V0_j1zCLfSQlTb8jJpQfScA.TIFF')
roads <- rast('data/distance2road.TIFF')

######### Task 1: Build Home Ranges #############

UDs.deer <- estimate.UD(deer, 'id', 95)

######## Task 2: Generate Random Points ###########

deer.use.avail <- generate.available.points(UDs.deer, 'id', deer)

####### Task 3: Covariate extraction ############

#nlcd
deer.final <- extract.custom.nlcd(deer.use.avail, nlcd)

#distance to road
deer.final$road.dist <- extract(scale(roads), deer.final)[,2]

head(deer.final)

###### Task 4: Build Models ################
deer.final.df <- as.data.frame(deer.final)

#look at data briefly
ggplot(deer.final.df, aes(road.dist)) + geom_histogram() + facet_wrap(~used)
ggplot(deer.final.df, aes(nlcd_factor)) + geom_bar() + facet_wrap(~used)

deer.final.df %>% group_by(nlcd_factor, used) %>% count()

# Jags model
sink(file = "data/deer_rsf.jags")
cat("
    model {
    for(j in 1:7){
      beta[j] ~ dnorm(0, 0.1)
    }
    
    for(i in 1:n){
      logit(use[i]) <- beta[1] + beta[2] * crop[i] + beta[3] * developed[i] +
                       beta[4] * grassland[i] + beta[5] * water[i] + beta[6] * wetland[i] +
                       beta[7] * road.dist[i]
                      
      y[i] ~ dbern(use[i])
    }
    
    }
", fill = TRUE)
sink()

# set up data
jags.data <- list(n = nrow(deer.final.df),
                  y = deer.final.df$used,
                  crop = deer.final.df$Crop,
                  developed = deer.final.df$Developed,
                  grassland = deer.final.df$Grassland,
                  water = deer.final.df$Water,
                  wetland = deer.final.df$Wetlands,
                  road.dist = deer.final.df$road.dist)

inits <- function(){list()}  

parameters <- c('beta')

n.chains <- 3
n.thin <- 10
n.iter <- 20000
n.burnin <- 10000

#run model
model <- jags(jags.data, 
          inits, 
          parameters, 
          "data/deer_rsf.jags", 
          parallel = T, 
          n.chains = n.chains,
          n.thin = n.thin, 
          n.iter = n.iter, 
          n.burnin = n.burnin)

model
summary(model)

##### Task 5: Model Predictions and Plots

#roads

m.iters <- length(model$sims.list$beta)
m.intercepts <- model$sims.list$beta[,1]
m.slopes <- list(model$sims.list$beta[,2],
                 model$sims.list$beta[,3],
                 model$sims.list$beta[,4],
                 model$sims.list$beta[,5],
                 model$sims.list$beta[,6],
                 model$sims.list$beta[,7])
m.covs <- list(0,0,0,0,0,
               c(min(deer.final.df$road.dist), max(deer.final.df$road.dist)))

m.road.predict <-
jags.predict(intercepts = m.intercepts,
             slope.coefficients = m.slopes,
             exp.cov.values = m.covs,
             link.function = "log")

plot.jags.predictions(plotting.data = m.road.predict$plotting.data,
                      expected.covariate.values = m.road.predict$expected.cov.values[[6]],
                      quantiles = m.road.predict$expected.quantiles,
                      x.lab = "Distance from Road",
                      y.lab = "Probability of Selection")


#landcover

m.covs <- list("category", 
               "category", 
               "category", 
               "category", 
               "category", 
               mean(deer.final.df$road.dist))

category.names <- c("Forest", "Crop", "Developed", "Grassland", "Water", "Wetland")

m.land.predict <-
  jags.predict(intercepts = m.intercepts,
               slope.coefficients = m.slopes,
               exp.cov.values = m.covs,
               link.function = "log",
               covariate.type = "categorical")

plot.df <- 
  data.frame(category = category.names,
             cov = 1:nrow(m.land.predict$expected.quantiles),
             low.CI = m.land.predict$expected.quantiles[,1],
             median = m.land.predict$expected.quantiles[,2],
             high.CI = m.land.predict$expected.quantiles[,3])

ggplot() + 
  geom_violin(data = m.land.predict$plotting.data, aes(as.factor(cov), response), bw = 0.01, fill = "gray90") + 
  coord_cartesian(ylim = c(0,1)) +
  geom_linerange(data = plot.df, aes(x = as.factor(cov), ymin = low.CI, ymax = high.CI), size = 0.7) +
  geom_point(data = plot.df, aes(as.factor(cov), median), color = "red", size = 2.5) +
  theme_bw() +
  scale_x_discrete(name = waiver(), labels = category.names) +
  ylab("Probability of Selection") +
  xlab("")




