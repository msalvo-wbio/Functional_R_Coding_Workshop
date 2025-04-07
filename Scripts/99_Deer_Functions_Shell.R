#### Deer Functions

###### Data Prep #########

#function to install and load required packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

#load libraries
#ipak()

#read in file dependencies
nlcd.legend <- read.csv('data/NLCD/NLCD_landcover_legend_2018_12_17_j1zCLfSQlTb8jJpQfScA.csv')

###### Data Setup #########

#estimate a home range for each indivdual animal
estimate.UD <- function(sf.data,   # an sf data frame
                        id.column, #character  with animal id column name
                        percentUD = 95) # percent kde home range extracted numeric from 0-100
{
  
  #concert sf to sp
  sp.data <- as(sf.data, "Spatial")
  
  # remove everything other than the id column
  sp.data <- sp.data[,c(id.column)]
  
  #build home range for deer
  myUD <- kernelUD(sp.data, grid = 30, extent = 0.5, same4all = FALSE)
  
  #volume of each pixel
  my.volumeUD <- getvolumeUD(myUD)
  
  #extract polygon boundaries from uD
  UD.polys <- getverticeshr(my.volumeUD, percent = percentUD)
  
  #return to sf
  sf.ud.polys <- st_as_sf(UD.polys)
  #set crs to source crs
  sf.ud.polys <- st_transform(sf.ud.polys, st_crs(sf.data))
  
  return(sf.ud.polys)
  
}


generate.available.points <- function(UD.polygons,
                                      id.column,
                                      use.data,
                                      use.avail.ratio = 10){
  
  
  # indicate used points were used
  use.data$used <- 1
  
  #For each animal's home range
  for(ibg in 1:nrow(UD.polygons)){
    
    #select the polygon
    polygon.UD <- UD.polygons[ibg,]
    #find the name of the animal ID
    animal.id <- polygon.UD[[id.column]]
    #get the number of used points for the animal
    n.points <- nrow(use.data[use.data[[id.column]] == animal.id,])
    #generate random points
    random.points <- st_sample(polygon.UD,
                               size = n.points * use.avail.ratio,
                               type = 'random')
    
    #convert points to a dataframe
    rp.sf <- data.frame(geometry = random.points)
    
    #get the column names of the data which aren't geometry
    data.col.names <- colnames(use.data)[colnames(use.data) != 'geometry']
    
    #set the column values of the random points to the animal data
    rp.sf[data.col.names] <- as.data.frame(use.data[use.data[[id.column]] == animal.id,][1,]) %>% dplyr::select(-geometry)
    
    #set used to 0
    rp.sf$used <- 0
    
    if(ibg == 1){
      output <- rbind(rp.sf, use.data)
    } else{
      output <- rbind(output, rp.sf)
    }
  }
  
  output <- st_as_sf(output)
  
  if(nrow(output) != nrow(use.data) * (use.avail.ratio + 1)){
    stop('Did not generate the right number of available points')
  } else {
    return(output)
  }
  
  
}


#extract nlcd for a seires of points and transform as needed
extract.custom.nlcd <- function(points, #the sf dataframe containing the points
                                raster, #a terra SpatRaster object of nlcd
                                raster.legend = nlcd.legend) #the legend of the nlcd raster, loaded in the functions script
{
  
  #re-project the raster to the point epsg
  raster.proj <- project(raster, paste0('EPSG:', st_crs(points)$epsg))
  
  # extract the points
  points$nlcd <- extract(raster.proj, points)[,2]
  points$nlcd <- round(points$nlcd) #round just in case
  
  #convert nlcd values to characters
  points <- 
    points %>% 
    left_join(raster.legend, by = join_by(nlcd == Value)) %>% #combine the nlcd values with the legend
    rename(nlcd_factor = Legend) #rename legend to nlcd_factor
  
  #combine groups
  points$nlcd_factor[points$nlcd_factor %in% c("Developed, Open Space", 
                                               "Developed, Low Intensity",
                                               "Developed, Medium Intensity",
                                               "Developed, High Intensity")] <- "Developed"
  
  points$nlcd_factor[points$nlcd_factor %in% c("Deciduous Forest",
                                               "Evergreen Forest",
                                               "Mixed Forest")] <- "Forest"
  
  points$nlcd_factor[points$nlcd_factor %in% c("Shrub/Scrub",
                                               "Herbaceuous")] <- "Grassland"
  
  points$nlcd_factor[points$nlcd_factor %in% c("Hay/Pasture",
                                               "Cultivated Crops")] <- "Crop"
  
  points$nlcd_factor[points$nlcd_factor %in% c("Woody Wetlands",
                                               "Emergent Herbaceuous Wetlands")] <- "Wetlands"
  
  points$nlcd_factor[points$nlcd_factor %in% c("Open Water")] <- "Water"
  
  points$nlcd_factor[points$nlcd_factor %in% c("Unclassified", "Barren Land", "Perennial Snow/Ice")] <- "Other"
  
  #pivot wider so each land cover has a 0-1 column
  points <-
    points %>% 
    mutate(value = as.numeric(1),
           nlcd_factor2 = nlcd_factor) %>%
    pivot_wider(names_from = nlcd_factor2, 
                values_from = value,
                values_fill = list(value = 0),
                values_fn = list(value = max))
  
  return(points)
}


###### Model Predictions ########

# predict from jags output
jags.predict <- function(intercepts,
                         slope.coefficients,
                         exp.cov.values,
                         credible.intervals = c(0.025, 0.5, 0.975),
                         res = 100,
                         link.function,
                         lognormal.sd,
                         covariate.type = "continuous"
){
  #extract the number of iterations from your simulated intercepts
  iters <- length(intercepts)
  
  ### Setup Data:
  
  #contin
  if(covariate.type == "continuous"){
    #create sequences of your covariate values
    for(ixi in 1:length(exp.cov.values)){ #for every covariate value
      if("numeric" %in% class(exp.cov.values[[ixi]])){ # if it is numeric:
        if(length(exp.cov.values[[ixi]]) == 1){ #if expected value is length 1, repeat it up the times of res
          exp.cov.values[[ixi]] <- rep(exp.cov.values[[ixi]], times = res)
        } else if(length(exp.cov.values[[ixi]]) == 2){ #if expected values are length 2 (a range), create a sequence of times of res
          exp.cov.values[[ixi]] <- seq(exp.cov.values[[ixi]][1], 
                                       exp.cov.values[[ixi]][2],
                                       length.out = res)
        } else { #if a covariate value is not length 1 or 2, give an error
          stop("Expected covariate values must be of length 1 or 2")
        }
        
      }
    }
  } else if(covariate.type == "categorical"){
    
    #check if the name of each list item is "category"
    is.category <- sapply(exp.cov.values, function(x) x == "category")
    
    #for every covariate value
    for(ixi in 1:length(exp.cov.values)){ 
      if(is.category[ixi] == T){ #if it is categorical
        exp.cov.values[[ixi]] <- rep(0, times = length(exp.cov.values[is.category]) + 1) #create a sequence of 0's which is 1+length of the number of categories
      } else if(is.category[ixi] == F){ #else if it's continuous
        if(length(exp.cov.values[[ixi]]) == 1){
          exp.cov.values[[ixi]] <- rep(exp.cov.values[[ixi]], 
                                       times = length(exp.cov.values[is.category]) + 1) #repeat the means the same number of times as the number of categories
        } else { #no varying values within a categorical prediction (yet)
          stop("Current version does not account for interactions, set continuous covariates to their mean value")
        }
      }
    }
    
    #Fill in 1's in every 0 sequence
    for(i in seq_along(exp.cov.values)){
      if(i == 1){pos <- 2} #set position index to 2 for the first value
      
      if (all(exp.cov.values[[i]] == 0)) {  # Check if the vector is all zeros
        if(i != 1){pos <- pos + 1} #add a 1 to the position index
        exp.cov.values[[i]][pos] <- 1 #set the position value to 1
      }
    }
    
    #set res
    res <- length(exp.cov.values[is.category]) + 1
    
  } else {stop("covariate.type must be one of 'continuous' or 'categorical'.")}
  
  
  #create a matrix for:
  pred.mat <- matrix(NA, iters, res) #predicted values
  pred.cred <- matrix(NA, res, length(credible.intervals)) #predicted median and credible intervals
  
  # Get your estimates
  for(jlv in 1:res){ # for all resolutions:
    
    parameter.values <- list() #initialize list to store parameter values in
    for(ixi in 1:length(exp.cov.values)){ #for all covariates in the model:
      
      if("numeric" %in% class(exp.cov.values[[ixi]])){ #different indexing for vectors and matrices
        parameter.values[[ixi]] <- slope.coefficients[[ixi]] * exp.cov.values[[ixi]][jlv]#get the product of the slope and the covariate value
      } else if("matrix" %in% class(exp.cov.values[[ixi]])){
        parameter.values[[ixi]] <- slope.coefficients[[ixi]] * exp.cov.values[[ixi]][,jlv]#get the product of the slope and the covariate value
      } else {
        stop("expected values must be a numeric vector or a expected.response.matrix matrix derived from this function")
      }
      
      
      if(ixi == 1){ #for the first value, add the intercepts and the first covariate term
        pred.mat[,jlv] <- intercepts + parameter.values[[ixi]]
      } else { #for all other values, add them to the current matrix column to get the complete expected value based on a normal distribution
        pred.mat[,jlv] <- pred.mat[,jlv] + parameter.values[[ixi]]
      } 
      
    }
    
    #If the model is not following a normal distribution:
    if(link.function == "log"){ #for log links
      pred.mat[,jlv] <- exp(pred.mat[,jlv])
    } else if(link.function == "lognormal"){ #for log-normal links
      pred.mat[,jlv] <- rlnorm(iters, pred.mat[,jlv], lognormal.sd)
    }
    
    #get predicted credible intervals
    pred.cred[jlv,] <- quantile(pred.mat[,jlv], credible.intervals)
  }
  
  #melted observations for graphing
  melted.obs <- melt(pred.mat); names(melted.obs) <- c('iteration','cov','response')
  
  return(list(expected.cov.values = exp.cov.values,
              expected.response.matrix = pred.mat,
              expected.quantiles = pred.cred,
              plotting.data = melted.obs))
  
}


# plot jags predictions
plot.jags.predictions <- 
  function(plotting.data, 
           expected.covariate.values,
           quantiles,
           x.lab,
           y.lab){
    
    smoothScatter(plotting.data$response ~ expected.covariate.values[plotting.data$cov],
                  las = 1, nrpoints = 0,
                  ylab = y.lab,
                  xlab = x.lab)
    lines(quantiles[,2] ~ expected.covariate.values, lty = 1, lwd = 3, col = 'white')
    lines(quantiles[,1] ~ expected.covariate.values, lty = 2, lwd = 3, col = 'white')
    lines(quantiles[,3] ~ expected.covariate.values, lty = 2, lwd = 3, col = 'white')
    
  }