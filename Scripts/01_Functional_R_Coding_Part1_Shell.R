
#load packages

  #function to install and load required packages: Our first example of custom functions!
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

ipak(c('tidyverse', 'datasets'))


##### Function Basics #########

print("Hellow World")

?print

View(print)

View(lm)

data(iris)

model <- lm(Sepal.Length ~ Species, data = iris)

summary(model)


##### Creating Custom Functions ########

sd(iris$Sepal.Length)/sqrt(length(iris$Sepal.Length))

2^2

squared <- function(x){
  x.sq <- x^2
  return(x.sq)
}

squared <- function(x){
  return(x^2)
}


squared(2)


SE <- function(values){
  fun.sd <- sd(values)
  fun.n <- length(values)
  se <- fun.sd/sqrt(fun.n)
  return(se)
}

SE(iris$Sepal.Length)
SE(iris$Petal.Length)


mean_se_summary <- function(df) # a dataframe
{ 
  
  #initialize vectors
  means <- c()
  SEs <- c()
  my.names <- c()
  
  #iterate through every column
  for(i in 1:ncol(df)){
    
    #filter to only numeric columns
    if(class(df[,i]) == "numeric"){
      
      #mean
      means <- c(means, mean(df[,i]))
      #SE
      SEs <- c(SEs, SE(df[,i]))
      #name
      my.names <- c(my.names, colnames(df)[i])
      
    } else { #if the column isn't numeric
      next() #skip it
    }
    
  }
  
  #create output data frame
  output <- data.frame(value = my.names,
                       mean = means,
                       se = SEs)
  
  return(output)
}


sum.iris <- mean_se_summary(iris)


ggplot(sum.iris, aes(value, mean)) + geom_point()




