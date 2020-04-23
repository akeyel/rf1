# Code for building package (Not Run)
#devtools::document()
#devtools::check()

#' @importFrom caret RMSE
#' @importFrom psych pairs.panels
#' @importFrom randomForest randomForest
#' @importFrom quantregForest quantregForest
#' @importFrom grDevices dev.off tiff
#' @importFrom graphics barplot par plot segments text
#' @importFrom stats aggregate as.formula cor median na.exclude optimize predict qchisq runif uniroot rnorm
#' @importFrom utils read.csv write.table
#' @importFrom usethis use_data
#' @importFrom dplyr '%>%'
NULL

#' Run the rf1 model
#'
#' The RF1 model first fits a random forest, then excludes all variables with
#' importance scores below the mean importance (NOTE: This will be problematic
#' if all input variables are relevant) It then further pares down the variable
#' list using variance partitioning. It retains only variables that contribute
#' uniquely to explaining variation in the unmeasured year (via
#' leave-one-year-out cross-validation). The RF1 model uses the randomForest
#' package by Liaw & Wiener 2002, that implements the Random Forest method
#' developed by Breiman 2001. MLE calcualtions from MLE_IR.R were written by
#' Williams and Moffit 2005 (reformatted by A. Keyel). Note that an arbitrary
#' starting seed is set, to ensure that results are repeatable.
#'
#' Citations
#'    Breiman, L. 2001. Random forests. Machine Learning 45: 5- 32
#'    Keyel, A.C. et al. 2019 PLOS ONE 14(6): e0217854. https://doi.org/10.1371/journal.pone.0217854
#'    Liaw & Wiener 2002. Classification and Regression by randomForest. R News 2: 18-22
#'    Williams, C and C. Moffitt 2005. Estimation of pathogen prevalence in pooled samples
#'    using maximum likelihood methods and open source software. Journal of Aquatic Animal Health 17: 386-391
#'
#' @param forecast.targets A vector containing options of what to forecast.
#' 'annual.human.cases' generates human case predictions, while 'seasonal.mosquito.MLE'
#' provides options for mosquito predictions
#' @param human.data Data on human cases of the disease. Must be formatted with
#'   two columns: location and date. The location column contains the spatial
#'   unit (typically county), while the date corresponds to the date of the
#'   onset of symptoms for the human case.
#' @param mosq.data Data on mosquito pools tested for the disease. Must be
#'   formatted with 4 columns: location (the spatial unit, e.g. county),
#'   col_date: the date the mosquitoes were tested, wnv_result: whether or not
#'   the pool was positive, pool_size: the number of mosquitoes tested in the
#'   pool. A fifth column species is optional but is not used by the code
#' @param weather.data Data on weather variables to be included in the analysis.
#'   See the read.weather.data function for details about data format. The
#'   read.weather.data function from ArboMAP is a useful way to process one or
#'   more data files downloaded via Google Earth Engine.
#' @param weekinquestion The focal week for the forecast. For the Random Forest
#'   model, this will be the last day used for making the forecast
#' @param rf1.inputs Inputs specific to the RF1 model, see
#'   \code{\link{rf1.inputs}}. If this model is not included, this should be set
#'   to 'none' or omitted from the function call #**# LINK TO AN OBJECT WITH
#'   MORE DETAILS
#' @param results.path The base path in which to place the modeling results.
#'   Some models will create sub-folders for model specific results
#' @param id.string An id to use for labeling the aggregations across all locations
#' @param break.type The temporal frequency to use for the data. The default is
#'   'seasonal' which breaks the environmental data into January, February,
#'   March; April, May, June; July, August, September; October, November,
#'   December. Other options may be supported in the future.
#' @param response.type Whether data should be treated as continuous (mosquito
#'   rates, number of cases) or binary (0 or 1).
#' @param quantile.model Whether (1) or not (0) to use a quantile random forest
#' for the final model output. All other calculations and model fitting use the
#' standard randomForest package.
#' @param n.draws The number of random realizations to draw for the RF1.distributions object
#' @param bins Bin break points for the CDC forecast challenge
#' @param use.testing.objects An indicator. If TRUE, the analysis will not run,
#' but will load previously saved outputs in order to expedite testing the formatting
#' of the code outputs.
#'
#' @return Four outputs are generated: The Results dataframe, the Distributions dataframe, the Bins dataframe, and the model object results
#'
#' @export rf1
rf1 = function(forecast.targets, human.data, mosq.data, weather.data,
               weekinquestion, rf1.inputs, results.path, id.string, break.type = "seasonal", response.type = "continuous",
               quantile.model = 1, n.draws = 1000, bins = c(0,seq(1,51,1),101,151,201,1000),
               use.testing.objects = FALSE){
  
  out = FormatDataForRF1(human.data, mosq.data, weekinquestion, weather.data, rf1.inputs, results.path, break.type)
  my.data = out[[1]]
  independent.vars = out[[2]]
  
  # Split the data set into the forecast year and the historical data
  forecast.year = as.numeric(substr(as.character(weekinquestion), 1, 4))
  
  forecast.data = my.data[my.data$year == forecast.year, ]
  historical.data = my.data[my.data$year < forecast.year, ] #**# This will prevent later years from informing hindcasts of earlier years

  message(sprintf("nrow(historical.data) = %s", nrow(historical.data)))
  #message(paste(colnames(my.data), collapse = ', '))
  message("Independent Vars:")
  message(paste(independent.vars, collapse = ', '))
  
  #message(sprintf("Forecast data rows: %s", nrow(forecast.data)))
  #message(forecast.year)
  if (nrow(forecast.data) == 0){ stop(sprintf("Forecast subset has no data for forecast year %s. Please ensure that all temporally-merged data sets reach the final year.", forecast.year))}
  #stop("Need to fix stuff, maybe")
  
    
  #**# How do I go from a MIR to positive location weeks? Isn't the MIR more useful?
  # First question: Could estimate number of trap nights, number of mosquitoes sampled, and then use the infection rate to get an estimate of the number of positive location-weeks.
  # For second question - does seem like a lot of work to get at a number that tells us what? In ArboMAP it is used to estimate human cases, but we do that directly.
  #**# So for now, maybe the mosquito results aren't relevant? Or phrased better - are not relevant to the specific comparison question being asked.
  annual.positive.location.weeks = NA

  # Create output objects
  RF1.results = data.frame(forecast.target = NA, location = NA, value = NA)
  RF1.distributions = data.frame(forecast.target = NA, location = NA)
  for (i in 1:n.draws){
    new.col = sprintf("DRAW%s", i)
    RF1.distributions[[new.col]] = NA}
  RF1.bins = data.frame(location = NA, target = 'Total WNV neuroinvasive disease cases',
                        type = 'bin', unit = 'cases', bin_start_inclusive = NA, bin_end_notinclusive = NA )
    
  #**# WHAT AM I DOING WITH INFECTION RATE? PREVIOUSLY I ESTIMATED IT OVER THE ENTIRE YEAR, BUT HERE, THE DATA MIGHT BE MISSING OR FOR PART OF THE YEAR.
  #**# I could use the forecast mosquito infection rate? Would that be more useful than just leaving it out?
  # NOTE: the mosquito-only analysis is not directly connected to the prediction of human cases #**# Is there a way to do that, if human cases are not directly estimated?
  if ('seasonal.mosquito.MLE' %in% forecast.targets){
    # First forecast mosquito infection rates
    mosq.model = rf1.inputs[[6]]
    
    # If the mosquito model needs to be estimated, estimate it
    if (is.na(mosq.model)){
      dep.var = "IR"
      # drop IR variable from independent.vars
      m.independent.vars = independent.vars[independent.vars != "IR"]
      mosquito.results.path = sprintf("%s/mosquitoes", results.path)
      m.label = "mosquitoes" #**# THINK ABOUT THIS AND WHETHER THIS SHOULD BE AN ENTRY in rf1.inputs
      
      # Load previously run data for some testing situations
      if (use.testing.objects == 1){
        mosquito.results = rf1::mosquito.results
      }else{
        mosquito.results = do.rf(historical.data, dep.var, m.independent.vars, mosquito.results.path,
                                 response.type = response.type, label = m.label, quantile.model = quantile.model) #do.spatial = 0, create.test.set = 0, create.ci.null = 0
      }
      mosq.model = mosquito.results[[1]]
      kept.vars = mosquito.results[[6]]
    }
    
    # Update RF1.results
    if (quantile.model == 1){
      data.subset = forecast.data[ ,kept.vars] #**# HOPE THIS WORKS
      predictions = predict(mosq.model, data.subset, what = c(0.5)) # Make a single prediction at the 50% quantile value
      dist.predictions = predict(mosq.model, data.subset, what = runif(n.draws)) # select random quantile values. #**# Should talk with someone if this is a legitimate way to sample a distribution
    }else{
      predictions = predict(mosq.model, newdata = forecast.data)
      dist.predictions = rep(predictions, n.draws)
      dist.predictions = matrix(dist.predictions, ncol = n.draws) #**# Confirm the order is correct on this. Should be - should fill by column first.
    }
    
    locations = forecast.data$location
    targets = rep('seasonal.mosquito.MLE', length(predictions))
    mosq.records = data.frame(forecast.target = targets, location = locations, value = predictions)
    RF1.results = rbind(RF1.results, mosq.records)

    # Add statewide result
    seasonal.mosquito.MLE = mean(predictions, na.rm = TRUE) #**# NOTE: perhaps locations should be weighted in this estimate
    statewide.record = c('seasonal.mosquito.MLE', sprintf("%s-STATEWIDE", id.string), seasonal.mosquito.MLE)
    RF1.results = rbind(RF1.results, statewide.record)
    
    # Update RF1.distributions & add statewide result
    mosq.distributions = cbind(targets, locations, dist.predictions)
    colnames(mosq.distributions) = c('forecast.target', 'location', sprintf("DRAW%s", seq(1,n.draws)))
    RF1.distributions = rbind(RF1.distributions, mosq.distributions)
    statewide.prediction = apply(dist.predictions, c(2), mean, na.rm = TRUE) # Should take average MLE across all counties. #**# Has some obvious statistical issues, but not sure how else to do it.
    statewide.distribution = c('seasonal.mosquito.MLE', sprintf("%s-STATEWIDE", id.string), statewide.prediction) 
    names(statewide.distribution) =  c('forecast.target', 'location', sprintf("DRAW%s", seq(1,n.draws)))
    RF1.distributions = rbind(RF1.distributions, statewide.distribution)
        
    # Update RF1.bins
    #**# NOT YET SCRIPTED
  }
  
  # Run the human analysis (or not)
  if ('annual.human.cases' %in% forecast.targets){
    # Second, forecast human cases
    human.model = rf1.inputs[[7]]
    
    # If the human model needs to be estimated, estimate it
    if (is.na(human.model)){
      dep.var = "Cases"
      human.results.path = sprintf("%s/humans", results.path)
      h.label = "humans"
      
      # Load previously run data for some testing situations
      if (use.testing.objects == 1){
        human.results = rf1::human.results
      }else{
        human.results = do.rf(historical.data, dep.var, independent.vars, human.results.path,
                              label = h.label, response.type = response.type,
                              quantile.model = quantile.model)
      }
      human.model = human.results[[1]]
      kept.vars = human.results[[6]]
    }
    
    
    if (quantile.model == 1){
      data.subset = forecast.data[ , kept.vars] #**# Will this line work?
      predictions = predict(human.model, data.subset, what = 0.5)
      dist.predictions = predict(human.model, data.subset, what = runif(n.draws))
    }else{
      predictions = predict(human.model, newdata = forecast.data)
      dist.predictions = rep(predictions, n.draws)
      dist.predictions = matrix(dist.predictions, ncol = n.draws) #**# Confirm correct matrix formation
    }

    locations = forecast.data$location
    targets = rep('annual.human.cases', length(predictions))
    human.records = data.frame(forecast.target = targets, location = locations, value = predictions)
    RF1.results = rbind(RF1.results, human.records)
    
    statewide.cases = sum(predictions)
    statewide.record = c('annual.human.cases', sprintf("%s-STATEWIDE", id.string), statewide.cases)
    RF1.results = rbind(RF1.results, statewide.record)
    
    # Update forecast distributions
    human.distributions = cbind(targets, locations, dist.predictions)
    colnames(human.distributions) = c('forecast.target', 'location', sprintf("DRAW%s", seq(1,n.draws)))
    RF1.distributions = rbind(RF1.distributions, human.distributions)
    statewide.prediction = apply(dist.predictions, c(2), sum, na.rm = TRUE) # Should take average MLE across all counties. #**# Has some obvious statistical issues, but not sure how else to do it.
    statewide.distribution = c('annual.human.cases', sprintf("%s-STATEWIDE", id.string), statewide.prediction) 
    names(statewide.distribution) =  c('forecast.target', 'location', sprintf("DRAW%s", seq(1,n.draws)))
    RF1.distributions = rbind(RF1.distributions, statewide.distribution)

    # Update forecast bins
    #**# NOT SCRIPTED
    
  }
  
  # Remove the leading NA row
  RF1.results = RF1.results[2:nrow(RF1.results), ]
  RF1.distributions = RF1.distributions[2:nrow(RF1.distributions), ]
  #RF1.bins = RF1.bins[2:nrow(RF1.bins), ]

  # Check that fields come out as appropriate format
  RF1.results$value = as.numeric(as.character(RF1.results$value))
  for (i in 1:n.draws){
    RF1.distributions[ , (i + 2)] = as.numeric(as.character(RF1.distributions[ , (i + 2)]))
  }

  # Set up other outputs
  RF1.other = list()
  if ('annual.human.cases' %in% forecast.targets){
    RF1.other$human.results =  human.results
  }
  if ('seasonal.mosquito.MLE' %in% forecast.targets){
    RF1.other$mosquito.results = mosquito.results
  }
    
  return(list(RF1.results, RF1.distributions, RF1.bins, RF1.other = RF1.other))
}

#' rf1.inputs
#'
#' The model-specific inputs required to run the RF1 model
#'
#' @param files.to.add A list of file names of other data sources to be
#' included in the Random Forest model. If no additional data is to be added,
#' this should be an empty list (i.e. list()) or NA.
#' @param merge.type.vec A vector identifying how the files should be joined to
#' the mosquito and human data. Options are spatial_temporal where merges will
#' be performed on county and year (e.g., climate data); state_temporal where
#' merges will be performed on state and year (e.g. BBS data at state level),
#' and spatial, where merges will be performed on county only (e.g., landcover
#' and census data). This must have the same length as files.to.add, and if no
#' additional data is to be added, this should be an empty vector (i.e. c()) or NA.
#' @param analysis.counties A list of counties included in the analysis (this
#' is for ensuring that a county and year that does not have a human case is
#' treated as a 0)
#' @param analysis.years A list of years included in the analysis (this is for
#' ensuring that a county and year that does not have a human case is treated
#' as a 0 rather than as missing) #**# This is problematic if a specific county
#' or year is missing.
#' @param user.drop.vars A list of independent variables that should be
#' excluded from the analysis. If no variables should be excluded, this should
#' be an empty vector c() or NA.
#' @param mosq.model The Random Forest model to use for forecasting mosquito
#' infection rates. If this is to be fitted from the empirical data, this
#' should be set to NA
#' @param human.model The Random Forest model to use for forecasting human
#' cases. If this is to be fitted from the empirical data, this should be set
#' to NA.
#' @param no.data.exceptions #**# COMING SOON A list of county-years within the
#' range of analysis.years and analysis.counties that are missing data rather
#' than 0's
#'
#' @name rf1.inputs
NULL
#**# NOTE: DUPLICATED IN DOCUMENTATION IN DFMIP. NEED TO FIGURE OUT HOW TO COORDINATE


#' Random Forest 1 Model Outputs
#'
#' @description other.outputs$rf1 gives a list of 2 objects. The first is the
#' mosquito results, the second is the human results. Each of these objects are
#' lists of objects. The first $MODEL is the Random Forest Model object. This
#' can be used to generate novel predictions. The second $TEMPORAL.ACCURACY is
#' a summary of the cross-validation accuracy, with $OVERALL giving the overall
#' accuracy. Other accuracy results include $RMSE, $SPEARMAN, $N, $ERROR.DF
#' The third is $ SPATIAL.ACCURACY which performs the same cross-validation,
#' except on spatial units instead of on years. This has the same structure as
#' $TEMPORAL.ACCURACY. The fourth is $TEMPORAL.NULL, which contains two null
#' model results for comparison ($NULL.MEAN and $NULL.CI) each with the same
#' structure as $TEMPORAL.ACCURACY. The fifth is $SPATIAL.NULL, which contains
#' null model results for spatial structure, following the same structure as
#' $TEMPORAL.NULL. The sixth object is $RETAINED.VARS, which includes a list of
#' independent variables that were retained in the final model. (NOTE: These
#' can also be extracted from the $MODEL object)
#'
#' @name rf1.outputs
NULL
#**# Link to documentation from main RF1 description
#**# NOTE: DUPLICATED IN DOCUMENTATION IN DFMIP. NEED TO FIGURE OUT HOW TO
#**# USE A LINK TO JUST HAVE ONE COPY OF DOCUMENTATION
#**# SWITCH DESCRIPTION TO USE A TABLE/LIST FORMAT

#' A random forest based disease forecasting approach
#'
#'@description The do.rf function runs the random forest portion of the code
#'
#'@details Majority of code was writen in wnv_hlpr.R; January - June 2018.
#'Transferred to this file & modified beginning in March 2019. Code adapted to
#'use common inputs with the ArboMAP code of Davis & Wimberly
#' 
#' @param trap.data A data frame containing human cases, mosquito infection rate (if applicable),
#' and independent variables for prediction of human cases and/or mosquito infection rates
#' @param dep.var The field name containing the dependent variable to analyze
#' @param independent.vars The field names of the independent variables to include in the analysis
#'   Variables not listed here but present in the data frame will be excluded from the analysis.
#' @param results.path The base path in which to place the modeling results.
#'   Some models will create sub-folders for model specific results
#' @param do.spatial Whether or not to do spatial crossvalidation (can be time-consuming)
#' @param create.ci.null Whether or not to generate null model results
#' @param label A label for the analysis run
#' @param response.type Whether data should be treated as continuous (mosquito
#'   rates, number of cases) or binary (0 or 1).
#' @param exploratory Whether to identify the best model (exploratory = TRUE),
#'   or whether to run the random forest with all input independent variables (exploratory = FALSE)
#' @param input.seed A seed for the random number generator to ensure that the results are repeatable
#' @param temporal.field The field containing the temporal units
#' @param spatial.field The field containing the spatial units
#' @param quantile.model Whether (1) or not (0) to use a quantile random forest
#'   for the final model output. All other calculations and model fitting use the
#'   standard randomForest package.
#' @param display.messages Whether or not update messages should be output
#' 
#' @return A list containing:\tabular{ll}{
#' MODEL \tab the random forest prediction model\cr
#' TEMPORAL.ACCURACY \tab  an accuracy assessment using leave one year out cross-validation\cr
#' SPATIAL.ACCURACY \tab an accuracy assessment based on leave one location out cross-validation\cr
#' TEMPORAL.NULL \tab Null model results based on estimates across time\cr
#' SPATIAL.NULL \tab Null model results based on estimates across space\cr
#' RETAINED.VARS \tab Variables retained in the final prediction model\cr} 
#' 
#' @export do.rf
#'
do.rf = function(trap.data, dep.var, independent.vars, results.path, do.spatial = 0,
                 create.ci.null = 0, label = "", response.type = "continuous", exploratory = TRUE,
                 input.seed = 20180830, temporal.field = "year", spatial.field = "location",
                 quantile.model = 1, display.messages = 1){
  
  #require(psych)
  
  response.types = c("continuous", "binary")
  if (!response.type %in% response.types){ stop(sprintf("Response type must be %s, and is case sensitive", paste(response.types, collapse = ", "))) }
  
  # Reset the seed to ensure that the random elements in the code are repeataable
  set.seed(input.seed)
  
  message(nrow(trap.data))
  message(dep.var)
  message(names(trap.data))
  
  # Subset data to just non-NA values
  trap.data = trap.data[!is.na(trap.data[[dep.var]]), ]
  
  # Add MERGE_ID #**# FIX THIS TO BE MORE GENERAL / MORE ROBUST
  trap.data$MERGE_ID = trap.data$location_year
  
  #**# Also formally output these somewhere  Report after test data set is removed
  # Report sample size
  N = nrow(trap.data)
  NO.IR.N = nrow(trap.data[is.na(trap.data$IR), ]) # Get number of NA values for Infection Rate
  IR.not.zero = nrow(trap.data[trap.data$IR > 0, ])
  WNV.N = IR.not.zero + NO.IR.N
  if (display.messages == 1){
    message(sprintf("Total N: %s\nWNV+ N: %s\nNo IR N: %s", N, WNV.N, NO.IR.N))
    message("These numbers are not valid for the human case analyses")
  }
  
  
  # Set up the dependent variable
  if (response.type == "continuous"){ trap.data[[dep.var]] = as.numeric(as.character(trap.data[[dep.var]])) }
  if (response.type == "binary"){ trap.data[[dep.var]] = as.factor(trap.data[[dep.var]]) }
  
  message(nrow(trap.data))
  message(paste(independent.vars, collapse = ', '))
  
  ### EXAMINE PREDICTORS FOR COLLINEARITY
  # Run for Everything #**# Put someplace that makes more sense
  correlation.path = sprintf("%s/ModelResults/correlations", results.path)
  dir.create(correlation.path, showWarnings = FALSE, recursive = TRUE)
  tiff(filename = sprintf("%s/pairspanel_%s.tif", correlation.path, label), height = 3000, width = 3000)
  #message(paste(independent.vars, collapse = ', '))
  psych::pairs.panels(trap.data[ , independent.vars])
  dev.off()
  
  if (exploratory == TRUE){
    ### RUN MODEL FOR ALL PREDICTOR VARIABLES
    f = as.formula(paste(dep.var, ' ~ ', paste(independent.vars, collapse = "+")))
    
    best.m = get.best.m(f, independent.vars, trap.data, response.type)
    #best.m = 3 #**#
    
    # Regenerate best model using more trees
    if (display.messages == 1){  message("Creating the model with the best m using 5000 trees")  }
    rf.model = randomForest(f, data = trap.data, na.action = na.exclude, stringsAsFactors = TRUE, importance = TRUE, mtry = best.m, ntree = 5000)
    
    important.stuff = sort(rf.model$importance[,1], decreasing = TRUE)
    
    # Show plot of mean importance, with a line indicating the mean importance cut-off
    importance.path = sprintf("%s/ModelResults/importance", results.path)
    dir.create(importance.path, showWarnings = FALSE, recursive = TRUE)
    tiff(filename = sprintf("%s/importance_%s.tif", importance.path, label))
    plot(seq(1,length(important.stuff)), important.stuff)
    segments(0,mean(important.stuff), length(independent.vars), mean(important.stuff))
    dev.off()
    
    # Keep any variables above the mean importance
    best.vars = names(important.stuff)[important.stuff > mean(important.stuff)]
    
    if(display.messages == 1){  message("Creating a refined model that only includes variables of high importance")  }
    f2 = as.formula(paste(dep.var, ' ~ ', paste(best.vars, collapse = "+")))
    best.m = get.best.m(f2, best.vars, trap.data, response.type)# Recreate best.m variable for the subset. This should remove one of the warnings R was giving me
    #best.m = 3 #**#
    rf.model2 = randomForest(f2, data = trap.data, na.action = na.exclude, stringsAsFactors = TRUE, importance = TRUE, mtry = best.m, ntree = 5000)
    
    ## Refine variable set based on variance partitioning results
    # Calculate a preliminary temporal accuracy to provide a baseline
    temporal.accuracy.prelim = systematic.validation(trap.data, dep.var, f2, best.m, "year", response.type, display.messages)
    temporal.accuracy.R2 = temporal.accuracy.prelim$OVERALL[["R2.sasha"]]
    correlation.threshold = "" # 0.5 #**# not using this aspect initially
    drop.thresholds = c(0, 0.001, 0.002, 0.005, 0.01) #**# Try a vector of thresholds. Use the threshold that corresponds to the fewest variables with no more than a 5% decrease in prediction accuracy
    
    # Run variance partitioning
    kept.vars = do.variance.partitioning(trap.data, dep.var, best.vars, best.m, temporal.accuracy.R2, drop.thresholds,
                                         correlation.threshold, response.type, do.spatial, results.path,
                                         spatial.field, temporal.field, label, temporal.field, display.messages)
    
    kept.vars = as.character(kept.vars) # Ensure they are not treated as factors
    # Re-run the model with just those that are kept after the variance partitioning process
    f3 = as.formula(paste(dep.var, ' ~ ', paste(kept.vars, collapse = "+")))
    best.m = get.best.m(f3, kept.vars, trap.data, response.type)# Recreate best.m variable for the subset. This should remove one of the warnings R was giving me
    #best.m = 3 #**#
    if (quantile.model == 1){
      y = trap.data[[dep.var]]
      #message(y[1:5])
      #message(paste(kept.vars, collapse = ', '))
      x.df = trap.data[ , kept.vars]
      #message(head(x.df))
      rf.model3 = quantregForest(x.df, y, na.action = na.exclude, stringsAsFactors = TRUE, importance = TRUE, mtry = best.m, ntree = 5000)
    }else{
      rf.model3 = randomForest(f3, data = trap.data, na.action = na.exclude, stringsAsFactors = TRUE, importance = TRUE, mtry = best.m, ntree = 5000)
    }
    
  }
  
  # Just run the model for the input independent variables
  if (exploratory == FALSE){
    best.m = get.best.m(f3, independent.vars, trap.data, response.type)# Recreate best.m variable for the subset. This should remove one of the warnings R was giving me
    f3 = as.formula(paste(dep.var, ' ~ ', paste(independent.vars, collapse = "+")))
    kept.vars = as.character(independent.vars) # For results output
    
    if (quantile.model == 1){
      y = trap.data[[dep.var]]
      x.df = trap.data[ , kept.vars]
      rf.model3 = quantregForest(x.df, y, na.action = na.exclude, stringsAsFactors = TRUE, importance = TRUE, mtry = best.m, ntree = 5000)
    }else{
      rf.model3 = randomForest(f3, data = trap.data, na.action = na.exclude, stringsAsFactors = TRUE, importance = TRUE, mtry = best.m, ntree = 5000)
    }
  }
  
  ## Temporal validation step
  if (display.messages == 1){  message("Calculating temporal accuracy statistics")  }
  temporal.accuracy = systematic.validation(trap.data, dep.var, f3, best.m, temporal.field, response.type, display.messages) # formerly "TEMPORAL" instead of being temporal.field
  
  if (do.spatial == 1){
    ## Spatial validation step
    if (display.messages == 1){ message("Calculating spatial accuracy statistics") }
    spatial.accuracy = systematic.validation(trap.data, dep.var, f3, best.m, spatial.field, response.type, display.messages) #**# formerly SPATIAL instead of location
  }else{
    if (display.messages == 1){message("Skipping spatial accuracy assessment (may take > 17 hours for some analyses)")}
    spatial.accuracy = NA
  }
  
  if (display.messages == 1){message("Creating null models")}
  temporal.null = null.validation(trap.data, dep.var, temporal.field, create.ci.null, display.messages)
  spatial.null = null.validation(trap.data, dep.var, spatial.field, create.ci.null, display.messages)
  
  if (quantile.model == 0){
    # Write overall model observations and predictions to file #**# But these will be predictions from the same data used to generate the model. Yes. But the accuracy was assessed separately, and those will be the accuracy metrics reported.
    if (display.messages == 1){message("Writing predictions")}
    write.predictions(trap.data, dep.var, rf.model3, results.path, spatial.field, temporal.field, label)
  }
  
  #**# MAKE OUTPUTS #**# CLEAN THIS UP
  #make.maps(model.results, spatial.resolution, temporal.resolution, results.path)
  #plot.residuals(model.results, spatial.resolution, temporal.resolution)
  if (do.spatial == 1){
    if(display.messages == 1){ message("Creating barplots")}
    spatial.temporal.barplots(temporal.accuracy, spatial.accuracy, spatial.field, temporal.field, results.path, label)
  }
  
  return(list(MODEL = rf.model3, TEMPORAL.ACCURACY = temporal.accuracy,
              SPATIAL.ACCURACY = spatial.accuracy, TEMPORAL.NULL = temporal.null,
              SPATIAL.NULL = spatial.null, RETAINED.VARS = kept.vars))
}


#' Convert data in the format for ArboMAP to a format that can be used by the
#' RF1 model
#'
#' @param human.data human data from the main input, but read into memory
#' @param mosq.data mosquito data from main input, but read into memory
#' @param weekinquestion #**# ADD
#' @param weather.data non-temporal data from main input
#' @param rf1.inputs see \code{\link{rf1.inputs}}
#' @param results.path #**# ADD
#' @param break.type #**# ADD
#'
#' @noRd
#'
FormatDataForRF1 = function(human.data, mosq.data, weekinquestion, weather.data, rf1.inputs, results.path, break.type){
  
  # If it is NA or any other single entry, do not process the mosquito data
  if (length(mosq.data) != 1){
    # Need to convert mosquito data into mosquito infection rates (uses MLE code, which is not mine. Need to wait for response from Moffit. Perhaps a phone call? If no email response.)
    temporal.resolution = "annual"

    #Switch to a location_year naming convention & add year field from col_date field
    mosq.data$year = mapply(substr, mosq.data$col_date, nchar(as.character(mosq.data$col_date)) - 3, nchar(as.character(mosq.data$col_date)))
    mosq.data$year = as.numeric(mosq.data$year)
    mosq.data = district.to.location(mosq.data, "mosq.data")
    
    md.data = calculate.MLE.v2(mosq.data, temporal.resolution)
  }
  
  # Need to convert human data to cases per county per year, and add 0's for county years without cases
  all.counties = rf1.inputs[[3]] # Get list of counties
  all.years = rf1.inputs[[4]] # Get list of all years of analysis
  # Ensure human data is in data frame format, and not just a file name
  if (typeof(human.data) == "character"){  human.data = read.csv(human.data)  }
  
  hd.data = convert.human.data(human.data, all.counties, all.years)
  #message(max(hd.data$year, na.rm = TRUE))
  #message(weekinquestion)
  #message(break.type)
  breaks = assign.breaks(weekinquestion, break.type)
  #message(paste(breaks, collapse = ", "))
  
  # Need to aggregate the climate data by season & merge with static.data
  # Only process if weather.data is not NA
  if (length(weather.data) != 1){
    weather.data = district.to.location(weather.data)
    env.data = convert.env.data(weather.data, all.counties, breaks) #**# Could save this for faster re-use. How should I do that?
    #**# Watch for problem where location_year contains lower-case entries, while district has been converted to upper case for merging purposes.
  }else{
    #**# Watch input field names. These will need to be standardized, and the standardized names listed in the documentation
    # Otherwise, initialize an empty env.data with only the location_year, county, and year fields
    in.counties = rep(all.counties, length(all.years))
    in.years = sort(rep(all.years, length(all.counties)))
    env.data = data.frame(location = in.counties, year = in.years, location_year = sprintf("%s_%s", in.counties, in.years))
  }
  
  #message(max(env.data$year))
  
  env.data = add.rf1.inputs(env.data, rf1.inputs, breaks)
  #message(max(env.data$year))
  
  if (length(mosq.data) == 1){  my.data = hd.data
  }else{
    # Merge on county.year. Keep only records that have mosquito data
    my.data = merge(md.data, hd.data, by = "location_year", all.x = TRUE) #, all = TRUE
  }
  
  # Merge to environmental data.
  my.data = merge(my.data, env.data, by = "location_year") # , all = TRUE
  my.data = cleanup.garbage(my.data)
  
  #message(max(my.data$year))
  
  # Pull independent variable names from the climate data. Add those from the mosquito data
  independent.vars = colnames(my.data)
  
  #**# maybe independent.vars should be an input? That might be easier, and give the user more control. But how will the user know what fields there are?
  
  # Remove non-analysis variables #**# This needs an upgrade to allow a user input and user control!
  user.drop.vars = rf1.inputs[[5]]
  drop.vars = c(user.drop.vars, "location_year", "GROUP", "CI.lower", "CI.upper", "COUNTY", "Cases", "YEAR", "year", "county", "district", "breaks", "location") #**# I have a names problem! Should really clean up the name usage!
  
  for (drop.var in drop.vars){
    if (drop.var %in% independent.vars){
      independent.vars = independent.vars[independent.vars != drop.var]
    }
  }
  
  # Ensure variables come out as numeric, unless they are specified as factors
  for (var in (independent.vars)){
    # If it is not a factor, ensure it is 
    if (!typeof(my.data[[var]]) == 'factor'){
      my.data[[var]] = as.numeric(as.character(my.data[[var]]))
    }
  }
  
  #message(max(my.data$year))
  
  out = list(my.data, independent.vars)
  return(out)
}

#' Convert human data
#'
#' Simple little function to aggregate human case data by year and county
#'
#' @noRd
#'
convert.human.data = function(hd, all.counties, all.years){
  #**# The steps creating the year were already done in forecast_NYCT.R and could be required as part of the input.
  hd$year = mapply(substr, hd$date, nchar(as.character(hd$date)) - 3, nchar(as.character(hd$date)))
  hd$year = as.numeric(hd$year)
  # Check if data sets are input with district and district_year instead of location and location_year
  hd = district.to.location(hd)
  hd$location = as.character(hd$location)
  #hd$location_year = sprintf("%s_%s", hd$location, hd$year)
  hd$count = 1 # One case per entry
  hd.data = aggregate(hd$count, by = list(hd$location_year), "sum")
  colnames(hd.data) = c("location_year", "Cases")
  hd.data$location = sapply(hd.data$location_year, splitter, "_", 1, 1)
  hd.data$year =  sapply(hd.data$location_year, splitter, "_", 2, 0)
  
  # Make sure there is a record for every year and county included in the data set
  for (county in all.counties){
    for (year in all.years){
      location_year = sprintf("%s_%s", county, year)
      if (!location_year %in% hd.data$location_year){
        # Add 0 cases for missing location_years
        #test.type(c(county_year, 0, county, year), 'L1480')
        hd.data = rbind(hd.data, c(location_year, 0, county, year))
      }
    }
  }
  
  # For some reason, this field is being converted to text following the rbind. Odd
  hd.data$Cases = as.numeric(hd.data$Cases)
  
  return(hd.data)
}

#' Assign Breaks
#'
#' @noRd
#'
assign.breaks = function(weekinquestion, break.type){
  # Get the date of the forecast. Only include breaks prior to forecast.doy
  forecast.year=  as.numeric(substr(as.character(weekinquestion), 1, 4))
  forecast.month = as.numeric(substr(as.character(weekinquestion), 6,7))
  forecast.day = as.numeric(substr(as.character(weekinquestion), 9,10))
  forecast.doy = get.DOY(forecast.year, forecast.month, forecast.day)
  
  # Assign all breaks
  if (break.type == "seasonal"){
    # Breaks on 3/31, 6/30, 9/30, 12/31
    breaks = c(get.DOY(forecast.year,3,31), get.DOY(forecast.year, 6, 30), get.DOY(forecast.year,9,30), get.DOY(forecast.year, 12, 31))
  }
  
  if (break.type != "seasonal"){ stop(sprintf("%s needs to be scripted out", break.type)) }
  
  # Drop any breaks that go past the forecast.doy - these cannot be used to make a prediction as they have not been observed yet!
  breaks = breaks[breaks < forecast.doy] # Assume we have data up to the day before the forecast.doy #**# May not be a good assumption
  
  return(breaks)
}

#' Convert environmental data to the seasonal format used by the RF1 model
#'
#' @noRd
#'
convert.env.data = function(weather.data, all.counties, season.breaks){
  
  # These are column names, but the dplyr syntax makes them flag up as undefined global variables
  # So I define them here to make R happy, even though these are not actually used anywhere in the following
  # code, as the entries below refer to column names in a data object
  pr = rmean = tmaxc = tmeanc = tminc = vpd = wbreaks = NA
  
  #Save processing time by restricting to just locations of interest
  weather.data = weather.data[weather.data$location %in% all.counties, ]
  
  # Drop weather data beyond the last break
  last.break = season.breaks[length(season.breaks)]
  weather.data = weather.data[weather.data$doy <= last.break, ]
  
  # Add grouping factor based on day of year
  weather.data$wbreaks = sapply(weather.data$doy, assign.groups, season.breaks)
  weather.data$location_year = sprintf("%s_%s", weather.data$location, weather.data$year)
  
  # This looks like a job for dplyr; https://datacarpentry.org/dc_zurich/R-ecology/04-dplyr
  #**# This is hard-coded to specific variables. Can dplyr take a more general input?
  env.data.pre = weather.data %>%
    dplyr::group_by(location_year, wbreaks) %>%
    dplyr::summarize(TMINC = mean(tminc), TMEANC = mean(tmeanc), TMAXC = mean(tmaxc), PR = mean(pr),
                     RMEAN = mean(rmean), VPD = mean(vpd))
  
  
  #**# May be a better way to do this - this is not optimized for speed
  # Reformat env.data to have a separate column per break
  env.data = data.frame(location_year = unique(env.data.pre$location_year))
  
  # Add variables to the env.data data frame
  vars = c("TMINC", "TMEANC", "TMAXC", "PR", "RMEAN", "VPD") #**# HARD CODED, BECAUSE dplyr step is hard coded
  for (var in vars){
    for (b in unique(env.data.pre$wbreaks)){
      this.var = sprintf("%s_%s", var, b)
      env.data[[this.var]] = rep(NA, nrow(env.data))
      
      # Subset data frame to just this break
      env.data.subset = env.data.pre[env.data.pre$wbreaks == b, ]
      # Loop through subset and extract break values for this variable
      for (i in 1:nrow(env.data.subset)){
        location_year = env.data.subset$location_year[i]
        this.value = env.data.subset[[var]][env.data.subset$location_year == location_year]
        env.data[[this.var]][env.data$location_year == location_year] = this.value
      }
    }
  }
  
  env.data$location_year = as.character(env.data$location_year) # Somehow this was turning into a factor (!)
  env.data$location = sapply(env.data$location_year, splitter, "_", 1, 1)
  env.data$year = sapply(env.data$location_year, splitter, "_", 2, 0)
  
  env.data$location = toupper(env.data$location) #**# Is this going to break things badly?
  
  return(env.data)
}


#' Assign break grouping for aggregating weather data
#'
#' @noRd
#'
assign.groups = function(doy, breaks){
  
  out.group = 1
  
  # Increment the out.group until the appropriate break period is found
  # Starts at 1
  for (b in breaks){
    if (doy > b){
      out.group = out.group + 1
    }
  }
  
  if (out.group > length(breaks)){
    stop("An input day of year (doy) exceeded the last break")
  }
  
  # Simple manual version
  # Use of get.DOY function adjusts based on whether the year is a leap year
  #if (doy <= breaks[1]){  out.group = 1  }
  #if (doy > breaks[1] & doy <= breaks[2]){ out.group = 2  }
  #if (doy > breaks[2] & doy <= breaks[3]){ out.group = 3  }
  #if (doy > breaks[3]){ out.group = 4 }
  
  return(out.group)
}

#' Add data layers specific for the Random Forest model
#'
#' Modified from the add.other function in wnv_hlpr.R
#'
#' @param env.data The environmental data object
#' @param files.to.add A list containing the full paths for .csv files or the data objects to merge into the env.data object
#' @param merge.type.vec A vector containing the merge type for each of the files to be added
#'
#' @noRd
#'
add.rf1.inputs = function(env.data, rf1.inputs, breaks){
  
  # unpack rf1.inputs
  files.to.add = rf1.inputs[[1]]
  merge.type.vec = rf1.inputs[[2]]

  # env.data will be coming in with location_year, location, and year fields. Goal is to use consistent terminology
  
  # If other inputs are present, merge them in
  if (length(files.to.add) > 0){
    
      # Check that files.to.add has same length as merge.type
    if (length(files.to.add) != length(merge.type.vec)){  stop("files.to.add must have the same number of elements as merge.type.vec")  }
    
    for (i in 1:length(files.to.add)){
      
      # Add covariate information
      this.file = files.to.add[[i]]
      
      run.file = 1
      # Check if it is NA # More complicated because this.file could be a data set with length > 1, and is.na() just checks first element and makes R cranky
      if (length(this.file) == 1){
        if (is.na(this.file)){ run.file = 0 }
      }
      
      if(run.file == 1){ 
        
        this.merge = merge.type.vec[i]
        
        env.data = add.data(env.data, this.file, this.merge, breaks)
        env.data = cleanup.garbage(env.data) #**# This step may need modification in the context of the changes made to dfmip
      }
    }
    #}
  }
  
  # If rf1.inputs are NULL, do nothing, and just return the env.data object
  
  return(env.data)
}

#' Cleanup Garbage
#'
#' Taken from wnv_hlpr.R, rs.data is env.data object, but with the name used in wnv_hlpr.R
#'
#' Quick function to restore original field names after the merge
#' 
#' @noRd
#' 
cleanup.garbage = function(rs.data){
  # Clean up merge garbage to restore original variable names
  if (length(rs.data$year) == 0){
    rs.data$year = rs.data$year.x
    rs.data$year.x = NULL
    rs.data$year.y = NULL
  }
  
  #**# Using year as the temporal field for now
  #if (length(rs.data$TEMPORAL) == 0){
  #  rs.data$TEMPORAL = rs.data$TEMPORAL.x
  #  rs.data$TEMPORAL.x = NULL
  #  rs.data$TEMPORAL.y = NULL
  #}
  
  if (length(rs.data$location) == 0){
    rs.data$location = rs.data$location.x
    rs.data$location.x = NULL
    rs.data$location.y = NULL
  }
  
  return(rs.data)
}

#' Function to add data to the env.data object
#'
#' Modified from wnv_hlpr.
#' 
#' @param rs.data env.data, the object to have data appended to
#' @param this.file the file containing the data to append
#' @param this.merge the structure of the data for the merge - i.e. "spatial_temporal" for spatial and temporal,
#' "spatial" for spatial only or "state_temporal" for statewide data that varies by year (i.e. the BBS data)
#' temporal only would be theoretically possible, but has not been scripted.
#' For now, it assumes the spatial is county, and temporal is year. Need to generalize this to spatial and temporal.
#'
#' @noRd
#'
add.data = function(rs.data, this.file, this.merge, breaks){

  n.breaks = length(breaks)
  
  # If a file, load the file. If an R data object, merge the data object
  if (typeof(this.file) == 'character'){
    these.data = read.csv(this.file)
  }else{ these.data = this.file }
  
  # Merge by location and year
  if (this.merge == "spatial_temporal"){
    
    if (length(these.data$location_year) == 0){
      #**# ideally the stop will be in the external loop, and will give information about which data set is missing the required field
      #**# or we can add an rf1.inputs check tool that just checks this stuff right off the bat in the model
      stop("location_year field is missing from one or more of the spatail_temporal merge input data sets")
    }
    
    # # Format for my climate and anomaly data was variable_break. As this is being done for a specific forecast day, we need to drop any that exceed the forecast day
    # # Drop any fields that have a numeric ending beyond the last break
    # these.names = colnames(these.data)
    # keep.vec = c()
    # #**# This is a tricky split, as the variable names might have variable numbers of "_", but if we ever move past seasons to months, this will require a double-digit extraction
    # 
    # if (n.breaks < 10){
    #   # For now, do the simple approach and assume a single-digit extraction
    #   for (name in these.names){
    #     keep = FALSE
    #     # Just get last digit
    #     season = as.numeric(as.character(substr(name,length(name) - 1, length(name))))
    #     
    #     # Only keep if this variable corresponds to a season (i.e. ends in a number)
    #     if (season %in% seq(1,9)){
    #       if (season < n.breaks){
    #         keep = TRUE
    #       }
    #     }
    #     keep.vec = c(keep.vec, keep)
    #   }
    # }
    # 
    # if (max(breaks) > 9){ stop("Need to upgrade the add.data function to handle more than 9 breaks")  }
    # 
    # keep.vars = these.names[keep.vec] #NOTE: 0 and 1 do not seem to work here.
    
    #these.data$county_year = these.data$COUNTY_YEAR
    rs.data = merge(rs.data, these.data, by = "location_year")
  }
  
  # Merge by unit and year. Needs updating, but currently only BBS would be applicable here, and the county-scale merge could be done separately
  if (this.merge == "state_temporal"){
    stop("Need to add state as an input to the rf1 function")
    # Merge by State and year
    if (length(rs.data$year) > 0){ rs.data$STATE_YEAR = sprintf("%s_%s", rs.data$STATE, rs.data$year) }
    rs.data = merge(rs.data, these.data, by = "STATE_YEAR") #, all.x = TRUE)
  }
  
  # Merge by location only
  if (this.merge == "spatial"){
    
    # check that location field is present
    if (length(these.data$location) == 0){
      #**# ideally the stop will be in the external loop, and will give information about which data set is missing the required field
      #**# or we can add an rf1.inputs check tool that just checks this stuff right off the bat in the model
      stop("location_year field is missing from one or more of the spatial_temporal merge input data sets")
    }
    
    
    # Merge by spatial unit - no temporal resolution
    #these.data$location = toupper(these.data$SPATIAL) # Ensure data is upper-case for proper join
    #these.data$SPATIAL = NULL # Remove this field, otherwise it will later get converted to a SPATIAL.x if more than one gets merged. Can always re-assign from location if it is needed later
    
    rs.data = merge(rs.data, these.data, by = "location")
  }
  
  #message(sprintf("Joined. Max year is %s. Max data year is %s", max(rs.data$year), max(these.data$year)))
  
  return(rs.data)
}




#### HELPER FUNCTIONS FOR do.rf FUNCTION CALL ####
#### MODIFIED (OR NOT) FROM wnv_hlpr.R ####


#' Function to identify the m value that produces the best results
#'
#' @noRd
#'
get.best.m = function(f, independent.vars, trap.data, response.type){
  best.r2 = 0
  best.m = 0
  # Try different values of "m" # With 37 iterations of 5000 trees, this is going slow.
  # Demo with 500 trees - probably too few.
  # Getting m of 3-4. Seems about right
  for (m in 1:length(independent.vars)){
    #message(sprintf("Trying m = %s", m))
    rf.model = randomForest(f, data = trap.data, na.action = na.exclude, stringsAsFactors = TRUE, importance = TRUE, mtry = m, ntree = 1000)
    
    if (response.type == "continuous"){
      if(median(rf.model$rsq) > best.r2){
        best.m = m
        best.r2 = median(rf.model$rsq)
      }
    }
    if (response.type == "binary"){
      #message("R2 used for choice of m is actually the out-of-box error rate, not the R2")
      if(median(rf.model$err.rate[ ,1]) < best.r2){
        best.m = m
        best.r2 = median(rf.model$err.rate[ ,1])
      }
      
    }
    
  }
  return(best.m)
}

#' Systematic Validation
#'
#' @noRd
#'
systematic.validation = function(trap.data, dep.var, f2, best.m, drop.field, response.type = "continuous", display.messages = 1){
  #require(caret)
  
  # Get list of unique values from the field to have one unit sequentially dropped
  drop.vals = unique(trap.data[[drop.field]])
  #message(paste(names(trap.data), collapse = ", "))
  #message(drop.field)
  #message(drop.vals)
  
  # Create a data frame to hold the results
  errors = data.frame(observed = NA, predicted = NA, unit = NA, MERGE_ID = NA)
  rmse.errors = c()
  spearmans = c()
  samples = c()
  to.exclude = c()
  
  for (drop.val in drop.vals){
    # Subset the data to exclude the field in question
    sub.data = trap.data[trap.data[[drop.field]] != drop.val, ]
    this.sample = nrow(sub.data)
    samples = c(samples, this.sample)
    
    # Check to make sure there is variation to explain with a random forest
    if (length(unique(sub.data[[dep.var]])) > 0){
      # Fit the random forest
      this.result = randomForest(f2, data = sub.data, na.action = na.exclude, stringsAsFactors = TRUE, importance = FALSE, mtry = best.m, ntree = 5000)
      
      # Get validation subset
      validation.data = trap.data[trap.data[[drop.field]] == drop.val, ]
      outcomes = predict(this.result, newdata = validation.data)
      
      # Get RMSE for this unit
      this.rmse = caret::RMSE(validation.data[[dep.var]], outcomes, na.rm = TRUE)
      rmse.errors = c(rmse.errors, this.rmse)
      
      # Get scaled RMSE for this unit #**# Don't need this for every unit - just do once overall
      #scaled.rmse = this.rmse / (mean(validation.data[[dep.var]]))
      
      if (response.type == "continuous"){
        # Get Spearman correlations for this unit
        this.spearman = suppressWarnings(cor(validation.data[[dep.var]], outcomes, use = "complete.obs"))
        spearmans = c(spearmans, this.spearman)
      }
      if (response.type == "binary"){
        spearmans = NA
      }
      
      # Add results to data frame
      stuff = cbind(validation.data[[dep.var]], outcomes, rep(as.character(drop.val), length(outcomes)), validation.data$MERGE_ID)
      colnames(stuff) = c("observed", "predicted", "unit", "MERGE_ID")
      #test.type(stuff, 'L1991')
      errors = rbind(errors, stuff)
    }else{
      if(display.messages == 1){message(sprintf("%s had no variation in the dependent variable %s", drop.val, dep.var))}
      to.exclude = c(to.exclude, drop.val)
    }
    
  }
  
  # Exclude any variables where the Random Forest could not be tested due to lack of variation
  for (ex in to.exclude){
    if (ex %in% drop.vals){
      drop.vals = drop.vals[drop.vals != ex]
    }
  }
  
  # Drop NAs in first row added during data frame initialization
  errors = errors[2:nrow(errors), ]
  errors$observed = as.numeric(as.character(errors$observed))
  errors$predicted = as.numeric(as.character(errors$predicted))
  overall.rmse = caret::RMSE(errors$observed, errors$predicted, na.rm = TRUE)
  #0.0039 So, our average error is 3.9 cases per 1000.
  # Our average infection rate is 3.0 cases per 1000. So... on average we are just very wrong.
  min.error = min(abs(errors$observed - errors$predicted))
  max.error = max(abs(errors$observed - errors$predicted))
  #0.0231 (!) That is a HUGE mistake!
  
  custom.r2 = estimate.custom.r2(errors$observed, errors$predicted)
  #0.26 (!) # Better than I expected!
  
  # get median rmse of those calculated for each unit
  median.rmse = median(rmse.errors)
  # The median is 0.0029. So better than the mean overall RMSE.
  # Are the years with the most WNV the ones with the worst predictions?
  
  # Divide the root mean squared error by the mean value in the observed dataset. How big is the typical error relative to the mean?
  # One citation said a good rule of thumb is aim for an RMSE that is 10% of the target value or less
  scaled.rmse = overall.rmse / mean(errors$observed)
  
  #errors.by.year = cbind(as.numeric(as.character(drop.vals)), rmse.errors)
  #errors.by.year = errors.by.year[order(errors.by.year[ ,2], decreasing = TRUE), ]
  
  if (response.type == "continuous"){
    # Look at rank correlation ? See how well it does at ranking the magnitude of WNV?
    spearman = cor(errors[,c(1,2)], method = "spearman")[1,2]
    # 0.58 spearman correlation. So... there is a correlation between what the model predicts is a bad year and how bad the year is.
    # But it's not a great correlation.
    # And this is at the county level. Ouch.
    pearson = cor(errors[,c(1,2)], method = "pearson")[1,2]
  }
  
  if (response.type == "binary"){
    spearman = NA
    pearson = NA
  }
  
  overall.results = c(overall.rmse, median.rmse, scaled.rmse, custom.r2, spearman, pearson, min.error, max.error)
  #message(overall.rmse)
  #message(median.rmse)
  #message(scaled.rmse)
  #message(custom.r2)
  #message("break")
  #message(spearman)
  #message(pearson)
  #message(min.error)
  #message(max.error)
  names(overall.results) = c("Overall.RMSE", "Median.RMSE", "Scaled.RMSE", "R2.sasha", "Spearman", "Pearson", "Min.Error", "Max.Error")
  #**# 2012 was a bad year - totally underpredicted the WNV infection rates
  #**# So... this is limited as an early warning system (!)
  
  # Add unit informatoin to the rmse.errors & spearman correlations
  names(rmse.errors) = as.character(drop.vals)
  if (response.type == "continuous"){ names(spearmans) = as.character(drop.vals) }
  names(samples) = as.character(drop.vals)
  
  #hist(temporal.accuracy$SPEARMAN, breaks = seq(-1,1,0.1))
  
  return(list(OVERALL = overall.results, RMSE = rmse.errors, SPEARMAN = spearmans, N = samples, ERROR.DF = errors))
}


#' Similar to systematic, but for the production of the null comparisons
#' 
#' @noRd
#' 
null.validation = function(trap.data, dep.var, drop.field, create.ci.null, display.messages){
  #require(caret)
  
  # Add a message to warn users that some functionality has been turned off - to ensure it is intentional
  if (create.ci.null != 1 & display.messages == 1){ message("Not creating a Null model based on the confidence intervals of the original data")}
  
  # Get list of unique values from the field to have one unit sequentially dropped
  drop.vals = unique(trap.data[[drop.field]])
  
  # Create a data frame to hold the results
  errors = data.frame(observed = NA, mean.predicted = NA, ci.predicted = NA, unit = NA, MERGE_ID = NA)
  mean.rmse.errors = c()
  mean.spearmans = c()
  
  if (create.ci.null == 1){
    ci.rmse.errors = c()
    ci.spearmans = c()
  }
  samples = c()
  
  for (drop.val in drop.vals){
    # Subset the data to exclude the field in question
    sub.data = trap.data[trap.data[[drop.field]] != drop.val, ]
    this.sample = nrow(sub.data)
    samples = c(samples, this.sample)
    
    # Fit the random forest
    mean.result = mean(sub.data[[dep.var]])
    
    # Get validation subset
    validation.data = trap.data[trap.data[[drop.field]] == drop.val, ]
    mean.outcomes = rep(mean.result, nrow(validation.data))
    
    # Get RMSE for this unit
    this.mean.rmse = caret::RMSE(validation.data[[dep.var]], mean.outcomes)
    mean.rmse.errors = c(mean.rmse.errors, this.mean.rmse)
    
    if (create.ci.null == 1){
      # Repeat for CI Null
      ci.outcomes = mapply(runif, 1, validation.data$CI.lower, validation.data$CI.upper)
      this.ci.rmse = caret::RMSE(validation.data[[dep.var]], ci.outcomes)
      ci.rmse.errors = c(ci.rmse.errors, this.ci.rmse)
      
      # Get Spearman correlations for this unit
      # correlations for the mean values are meaningless! they're all the same predictions!
      #this.mean.spearman = cor(validation.data[[dep.var]], mean.outcomes)
      #mean.spearmans = c(mean.spearmans, this.mean.spearman)
      
      this.ci.spearman = cor(validation.data[[dep.var]], ci.outcomes)
      ci.spearmans = c(ci.spearmans, this.ci.spearman)
    }else{ ci.outcomes = rep(NA, length(mean.outcomes)) }
    
    # Add results to data frame
    stuff = cbind(validation.data[[dep.var]], mean.outcomes, ci.outcomes, rep(as.character(drop.val), length(mean.outcomes)), validation.data$MERGE_ID)
    colnames(stuff) = c("observed", "mean.predicted", "ci.predicted", "unit", "MERGE_ID")
    #test.type(stuff, 'L2128')
    errors = rbind(errors, stuff)
  }
  
  # Drop NAs in first row added during data frame initialization
  errors = errors[2:nrow(errors), ]
  errors$observed = as.numeric(as.character(errors$observed))
  errors$mean.predicted = as.numeric(as.character(errors$mean.predicted))
  
  mean.overall.rmse = RMSE(errors$observed, errors$mean.predicted)
  
  # Scale by mean to put into the context of the data
  mean.scaled.rmse = mean.overall.rmse / mean(errors$observed)
  
  mean.min.error = min(abs(errors$observed - errors$mean.predicted))
  mean.max.error = max(abs(errors$observed - errors$mean.predicted))
  mean.custom.r2 = estimate.custom.r2(errors$observed, errors$mean.predicted)
  
  # Look at median of the individual RMSE's
  mean.median.rmse = median(mean.rmse.errors)
  
  mean.overall.results = c(mean.overall.rmse, mean.median.rmse, mean.scaled.rmse, mean.custom.r2, NA, NA, mean.min.error, mean.max.error)
  names(mean.overall.results) = c("Overall.RMSE", "Median.RMSE", "Scaled.RMSE", "R2.sasha", "Spearman", "Pearson", "Min.Error", "Max.Error")
  
  # Add unit informatoin to the rmse.errors & spearman correlations
  names(mean.rmse.errors) = as.character(drop.vals)
  #names(mean.spearmans) = as.character(drop.vals)
  names(samples) = as.character(drop.vals)
  
  mean.results = list(OVERALL = mean.overall.results, RMSE = mean.rmse.errors, SPEARMAN = NA, N = samples, ERROR.DF = errors)
  #hist(temporal.accuracy$SPEARMAN, breaks = seq(-1,1,0.1))
  
  # Repeat process for CI results, if applicable
  ci.results = NA
  if (create.ci.null == 1){
    errors$ci.predicted = as.numeric(as.character(errors$ci.predicted))
    ci.overall.rmse = caret::RMSE(errors$observed, errors$ci.predicted)
    ci.scaled.rmse = ci.overall.rmse / mean(errors$observed)
    
    ci.min.error = min(abs(errors$observed - errors$ci.predicted))
    ci.max.error = max(abs(errors$observed - errors$ci.predicted))
    ci.custom.r2 = estimate.custom.r2(errors$observed, errors$ci.predicted)
    
    ci.median.rmse = median(ci.rmse.errors)
    # Look at rank correlation ? See how well it does at ranking the magnitude of WNV?
    #mean.spearman = cor(errors[,c(1,2)], method = "spearman")[1,2]
    #mean.pearson = cor(errors[,c(1,2)], method = "pearson")[1,2]
    ci.spearman = cor(errors[,c(1,3)], method = 'spearman')[1,2]
    ci.pearson = cor(errors[,c(1,3)], method = 'pearson')[1,2]
    
    ci.overall.results = c(ci.overall.rmse, ci.median.rmse, ci.scaled.rmse, ci.custom.r2, ci.spearman, ci.pearson, ci.min.error, ci.max.error)
    names(ci.overall.results) = c("Overall.RMSE", "Median.RMSE", "Scaled.RMSE", "R2.sasha", "Spearman", "Pearson", "Min.Error", "Max.Error")
    
    names(ci.rmse.errors) = as.character(drop.vals)
    names(ci.spearmans) = as.character(drop.vals)
    
    ci.results = list(OVERALL = ci.overall.results, RMSE = ci.rmse.errors, SPEARMAN = ci.spearmans, N = samples, ERROR.DF = errors)
    
  }
  
  return(list(NULL.MEAN = mean.results, NULL.CI = ci.results))
}

#' Calculate R2 estimate of variation explained
#'
#' This R2 is not bounded as zero, as the model is not derived from the original data
#' The model is derived from a novel data set, consquently, it could fail to explain even
#' as much variation as the data mean. In these cases, negative R2's will be obtained
#'
#' Consequently, a negative R2 might be better than no knowledge at all, but it's not apt to be very good (depends on how variable the data are from the mean!)
#'  Basic assumption is that there is some variability around the mean, and that typically, a model should be able to explain that variation better than a static mean.
#'
#' @noRd
#'
estimate.custom.r2 = function(observed, predicted){
  
  deviation = observed - predicted
  
  squared.deviation = deviation ^ 2
  
  # Calculate prediction R2 (this is not PRESS, which is based on leaving one out. This is an entirely new validation data set, and so the R2 could be negative)
  sum.squared.deviation = sum(squared.deviation)
  
  data.mean = mean(observed, na.rm = TRUE)
  deviation.from.mean = data.mean - observed
  deviation.from.mean.squared = deviation.from.mean ^ 2
  sum.dfm.squared = sum(deviation.from.mean.squared)
  
  # Scale the error remaining after the model by the error obtained from just using the mean
  # Subtract this from one.
  R2.pred = 1 - (sum.squared.deviation / sum.dfm.squared)
  
  return(R2.pred)
}

#' Determine the contribution of each variable to the total variation explained
#' 
#' @noRd
#' 
do.variance.partitioning = function(trap.data, dep.var, best.vars, best.m, temporal.accuracy.R2,
                                    drop.thresholds, correlation.threshold, response.type,
                                    do.spatial, results.path,
                                    spatial.resolution, temporal.resolution, this.label,
                                    temporal.field, display.messages){
  
  if (display.messages == 1){
    message("Running variance partitioning")
    message(sprintf("Baseline R2 is %.3f", temporal.accuracy.R2))
  }
  
  # Check if a spatial partitioning is desired - if so, throw a warning
  if (do.spatial == 1){ warning("Variance partitioning was performed for temporal crossvalidation, NOT spatial crossvalidation. The code could easily be adapted to support the spatial partitioning, but currently it does not.") }
  
  if (correlation.threshold != ""){
    # Identify correlations among variables
    cor.matrix = cor(trap.data[ ,best.vars])
    corr.variation = c()
    n.corr = c()
  }
  
  unique.variation = c()
  
  # Run variance partitioning leaving each variable out in sequence
  for (this.var in best.vars){
    if (display.messages == 1){message(sprintf("Processing %s", this.var))}
    remaining.vars = best.vars[best.vars != this.var] # Drop the variable in question
    out.R2 = run.partitioning(trap.data, dep.var, remaining.vars, best.m, temporal.accuracy.R2,
                              response.type, temporal.field, display.messages)
    unique.variation = c(unique.variation, out.R2)
    
    if (correlation.threshold != ""){
      # Run variance partitioning leaving out correlated suites of variables (greater than a threshold)
      cor.row = cor.matrix[this.var, ]
      cor.index = sapply(cor.row, simple.threshold, correlation.threshold)
      #cor.index = matrix(cor.index, ncol = 1)
      remaining.vars2 = best.vars[!cor.index]
      out.cor.R2 = run.partitioning(trap.data, dep.var, remaining.vars2, best.m, temporal.accuracy.R2,
                                    response.type, temporal.field, display.messages)
      corr.variation = c(corr.variation, out.cor.R2)
      n.corr = c(n.corr, sum(cor.index)) # This works because TRUE evaluates to 1 and FALSE evaluates to 0
    }
  }
  
  out.df = data.frame(VARIABLES = best.vars, UNIQUE.R2 = unique.variation)
  
  if (correlation.threshold != ""){
    out.df$CORR.R2 = corr.variation
    out.df$N.CORRELATED = n.corr
  }
  
  best.threshold = NA
  best.R2.diff = 0
  kept.vars = best.vars
  
  for (threshold in drop.thresholds){
    # Try dropping all variables that contribute less than 0.1% uniquely to the overall result #**# 0.1% is also arbitrary
    remaining.vars3 = out.df[out.df$UNIQUE.R2 > threshold, 1] # 1 indicates the first column
    
    if (length(remaining.vars3) == 0){
      if(display.messages == 1){message(sprintf("No variables remaining at threshold %s", threshold))}
    }else{
      this.R2.diff = run.partitioning(trap.data, dep.var, remaining.vars3, best.m, temporal.accuracy.R2,
                                      response.type, temporal.field, display.messages)
      if(display.messages == 1){message(sprintf("R2 difference for threshold %s is %.3f", threshold, this.R2.diff))}
      
      # Right now, this JUST accepts models that constitute an improvement.
      # We likely want to add a penalty term that favors simpler models with similar explanatory power.
      if (this.R2.diff < best.R2.diff){
        best.R2.diff = this.R2.diff
        best.threshold = threshold
        kept.vars = remaining.vars3
      }
    }
  }
  
  # Add a threshold column to the data set for reference as to what was used to identify the final model variables
  out.df$best.threshold = best.threshold
  
  # Write variance partitioning results to an output file
  out.path = sprintf("%s/ModelResults/VariancePartitions", results.path)
  dir.create(out.path, showWarnings = FALSE, recursive = TRUE) # showWarnings = FALSE suppresses a warning if the path already exists
  out.file = sprintf("%s/varpar_%s_%s_%s.csv", out.path, spatial.resolution, temporal.resolution, this.label)
  write.table(out.df, out.file, sep = ',', row.names = FALSE, col.names = TRUE)
  
  return(kept.vars)
}

#' Calculate the variance partitioning run
#' 
#' @noRd
#' 
run.partitioning = function(trap.data, dep.var, remaining.vars, best.m, temporal.accuracy.R2, response.type,
                            temporal.field = "TEMPORAL", display.messages = 1){
  
  my.f = as.formula(paste(dep.var, ' ~ ', paste(remaining.vars, collapse = "+")))
  
  # Redefine best.m if it is no longer applicable when the variable set is reduced
  if (best.m > length(remaining.vars)){ best.m = length(remaining.vars) }
  
  temporal.accuracy = systematic.validation(trap.data, dep.var, my.f, best.m, temporal.field, response.type, display.messages)
  
  # Calculate difference between the main run with all variables and this run
  R2.difference = temporal.accuracy.R2 - temporal.accuracy$OVERALL[["R2.sasha"]]
  
  return(R2.difference)
}

#' A simple function for using a threshold
#' 
#' @noRd
#' 
simple.threshold = function(x, threshold){
  y = NA
  if (x >= threshold){ y = TRUE }
  if (x < threshold){ y = FALSE }
  return(y)
}

#' Write model predictions
#'
#' @noRd
#'
write.predictions = function(trap.data, dep.var, rf.model, results.path, spatial.resolution, temporal.resolution, label){
  # Code used to patch the human analyses prior to addition of dep.var
  # DiseasePath = "C:/Users/ak697777/University at Albany - SUNY/Elison Timm, Oliver - CCEID/RESULTS/WNV_statistical/ModelResults"
  # spatial.resolution = "county"; temporal.resolution = "annual"; analysis.label = label = "human" #"human_subset"
  # trap.data = read.csv(sprintf("%s/%s_%s_%s.csv", DiseasePath, spatial.resolution, temporal.resolution, analysis.label))
  # load(sprintf("C:/Users/ak697777/University at Albany - SUNY/Elison Timm, Oliver - CCEID/RESULTS/WNV_statistical/ModelRuns/Results_%s_%s_%s.RData", spatial.resolution, temporal.resolution, analysis.label))
  # rf.model = model.results$MODEL
  trap.data$OBSERVED = trap.data[[dep.var]]
  trap.data$PREDICTED = predict(rf.model) #**# More commands needed?
  #output = trap.data[ , c("OBSERVED", "PREDICTED", "SPATIAL", "TEMPORAL")] #**# Do we need more than this? Yes for the visualization
  
  out.file = sprintf("%s/ModelResults/%s_%s_%s.csv", results.path, spatial.resolution, temporal.resolution, label)
  
  write.table(trap.data, file = out.file, sep = ',', row.names = FALSE)
  
}

#' Create spatial/temporal barplots
#' 
#' @noRd
#' 
spatial.temporal.barplots = function(temporal.accuracy, spatial.accuracy, spatial.resolution, temporal.resolution, results.path, label){
  
  outfile = sprintf("%s/RMSE_Spearman_%s_%s%s.tif", results.path, spatial.resolution, temporal.resolution, label)
  tiff(filename = outfile, compression = c("lzw"), height = 2400, width = 2400, res = 300)
  
  par(mfrow = c(2,2))
  par(mar = c(8,6,2,0))
  rmse.y.lim = c(0, 10)
  sp.y.lim = c(-1,1)
  
  x.axis = "s"
  x.offset = 0.8
  x.adj1 = x.adj2 = 0
  if (spatial.resolution == "point"){
    x.axis = "n"
    x.offset = 0.8
    x.adj1 = 20
    x.adj2 = 4.2
  }
  
  
  barplot(sort(temporal.accuracy$RMSE * 1000), ylab = "RMSE\n(# infected per 1000 mosquitoes)" , las = 2, ylim = rmse.y.lim, xpd = NA) #main = "RMSE"
  text(par("usr")[1] + x.offset, 9.8, "a.", cex = 1, xpd = NA)
  barplot(sort(spatial.accuracy$RMSE * 1000), ylab = 'RMSE\n(# infected per 1000 mosquitoes)', las = 2, ylim = rmse.y.lim, xaxt = x.axis, xpd = NA) #main = "RMSE"
  text(par("usr")[1] + x.offset + x.adj1, 9.8, "b.", cex = 1, xpd = NA)
  barplot(sort(temporal.accuracy$SPEARMAN), ylab = "Spearman coefficient", las = 2, ylim = sp.y.lim)
  text(par("usr")[1] + x.offset, 0.95, "c.", cex = 1, xpd = NA)
  barplot(sort(spatial.accuracy$SPEARMAN), ylab = "Spearman coefficient", las = 2, ylim = sp.y.lim, xaxt = x.axis)
  text(par("usr")[1] + x.offset + x.adj2, 0.95, "d.", cex = 1, xpd = NA)
  
  dev.off()
}

#' Calculate MLE for a specified subset of the data using R
#'
#' Modified from wnv_hlpr.R
#'
#' @param md Data on mosquito pools tested for the disease. Must be formatted with 4 columns: location (the spatial unit, e.g. county), col_date: the date the mosquitoes were tested, wnv_result: whether or not the pool was positive, pool_size: the number of mosquitoes tested in the pool. A fifth column species is optional but is not used by the code
#' @param temporal.resolution Must be set to 'annual' for now. Ideally support for finer-scale resolution will be added.
#'
#' @return Mosquito data compiled to give estimated mosquito infection rates by spatial unit and temporal resolution
#'
#' @export calculate.MLE.v2
calculate.MLE.v2 = function(md, temporal.resolution = "annual"){
  
  # Remove any NA rows from the calculations, report their presence to the user
  virus.nas = length(is.na(md$wnv_result))
  pool.nas = length(is.na(md$pool_size))
  md = md[!is.na(md$wnv_result), ]
  md = md[!is.na(md$pool_size), ]
  
  if (virus.nas != 0 | pool.nas != 0){
    warning(sprintf("NAs present in the input mosquito data: %s for wnv_result and %s for pool_size. These NA records have been removed from further analysis", virus.nas, pool.nas))
  }
  
  # Check that expected column names are present, if not, give an informative error
  expected.names = c('location', "col_date", "wnv_result", "pool_size") #'species' is not actually required at this point
  is.error = 0
  missing.vec = c()
  for (e.name in expected.names){
    if (!e.name %in% names(md)){
      is.error = 1
      missing.vec = c(missing.vec, e.name)
    }
  }
  if (is.error == 1){
    m = "One or more required fields is missing. Inputs are case sensitive. Missing fields are:"
    
    stop(sprintf("%s %s. Input fields were: %s", m, paste(missing.vec, collapse = ", "), paste(names(md), collapse = ', ')))
  }
  
  # Call dprev function from MLE_IR.R
  if (!exists("dprev")){ stop(" Load the dprev function from Williams & MOffit 2005 into computer memory. It is required for the calculation of MLE's") }
  # Make sure mosq.data is a data frame and not a file string
  if (typeof(md) == "character"){  md = read.csv(md)  }
  
  if (length(unique(md$species)) > 1){ warning(sprintf("Pooling across %s species", length(unique(md$species))))  }
  if (temporal.resolution != "annual"){ stop("Need to code other temporal resolution options")  }
  
  #**# THIS IS SPECIFIC TO THE FORMAT MIKE & JUSTIN USED FOR ARBOMAP - REVISIT TO SEE IF A MORE GENERAL APPROACH IS DESIRABLE
  # Get year from date
  md$year = mapply(substr, md$col_date, nchar(as.character(md$col_date)) - 3, nchar(as.character(md$col_date)))
  md$year = as.numeric(md$year)
  
  # Add VIRUS and ABUND fields #**# CONDSIDER WHETHER THESE ARE THE FIELD NAMES YOU REALLY WANT TO USE
  
  if (length(md$wnv_result) == 0){
    stop("A field with WNV result must be present, and must be labeled 'wnv_result'")
  }
  
  md$VIRUS = md$wnv_result
  md$wnv_result = NULL # Prevent an overabundance of redundant fields
  
  # Check that VIRUS field has only positives, negatives, and NA values
  md.positives = md$VIRUS[md$VIRUS == 1] 
  md.negatives = md$VIRUS[md$VIRUS == 0]
  #md.nas = md$VIRUS[is.na(md$VIRUS)]
  md.check = length(md.positives) + length(md.negatives) # + length(md.nas)
  #message(md.check)
  #message(length(md.positives))
  #message(length(md.negatives))
  #message(length(md.nas))
  #message(nrow(md))
  if (md.check != nrow(md)){
    m1 = "The VIRUS field (wnv_result) can only contain 0, 1, or NA values."
    m2 = sprintf("Entered values were: %s", paste(unique(md$VIRUS), collapse = ', '))
    m3 = sprintf("There were %s positive pools, %s negative pools, and %s pools total")
    m4 = sprintf("Thus %s pools are missing, as there should be %s pools", (nrow(md) - md.check), nrow(md))
    stop(sprintf("%s %s %s %s", m1, m2, m3, m4))
  }
  
  if (length(md$pool_size) == 0){
    stop("A field with mosquito pool sizes must be present, and must be labeled 'pool_size'")
  }
  
  md$ABUND = md$pool_size
  md$pool_size = NULL # Prevent redundnat fields
  
  # Drop any missing abundance field #**# Perhaps want to do this PRIOR to this function to ensure quality data inputs?
  md = md[!is.na(md$ABUND), ]
  
  # Aggregate by spatial field and year
  md$location_year = sprintf("%s_%s", md$location, md$year)
  
  # Get group ID's
  group.ids = unique(md$location_year)
  
  # Create a new data frame with group ID information
  md.data = data.frame(GROUP = group.ids, CI.lower = NA, CI.upper = NA,
                       IR = NA, COUNTY = NA)
  
  # Loop through each group and calculate MLE
  for (i in 1:length(group.ids)){
    group.id = group.ids[i]
    # Subset to just this group.id
    group.data = md[md$location_year == group.id, ]
    county = as.character(group.data$location[1]) # Should all be the same, just use first value
    this.year = group.data$year[1] # Should all be same, just use first value #**# Modify if temporal resolution is not annual
    
    # Get number of positive and negative pools for this group
    positives.index = group.data$VIRUS == 1
    negatives.index = group.data$VIRUS == 0
    #Explicitly finding negatives to exclude any NA values (i.e. from NY)
    
    # Get positive and negative pools including pool sizes
    positives = group.data$ABUND[positives.index]
    negatives = group.data$ABUND[negatives.index]
    abundance = sum(positives) + sum(negatives) # Get total number of mosquitoes tested. In CT, this is true abundance, in NY, pool size is capped, so this will underestimate abundance
    density = abundance / (length(positives) + length(negatives)) # Get number of mosquitoes divided by number of pools tested
    
    # If no pools, assign a value of -999 to distinguish from undefined IR's
    if (length(positives) == 0 & length(negatives) == 0){
      warning(sprintf("No pools for group %s. Not sure how that happened. Values of -999 assigned", group.id))
      IR = -999
      CI.upper = -999
      CI.lower = -999
      abundance = 0
      density = 0
    }else{
      # Proceed with checking for missing positives or negatives alone
      # If no positive pools, assign an IR of 0 and CI's of 0 to 0
      if (length(positives) == 0){
        IR = 0
        CI.lower = 0
        CI.upper = 0
      }
      
      # If no negative pools, assign an IR of NA, and CI's of NA
      if (length(negatives) == 0){
        IR = NA
        CI.lower = NA
        CI.upper = NA
      }
      
      if (length(positives) != 0 & length(negatives) != 0){
        # Calculate the infection rate (IR)
        output = dprev(positives, negatives, disp = 'n')
        CI.lower = output[1]
        CI.upper = output[3]
        IR = output[2]
      }
    } # END OF ELSE STATEMENT
    
    # update the data frame
    md.data$CI.lower[i] = CI.lower
    md.data$CI.upper[i] = CI.upper
    md.data$IR[i] = IR
    md.data$abundance[i] = abundance
    md.data$density[i] = density
    
    # Update the COUNTY entry
    md.data$COUNTY[i] = county
    
    # Ensure there is a clearly delineated year field. Adjust as appropriate for other temporal resolutions
    if (temporal.resolution == "annual") { md.data$year[i] = this.year }
  } # END OF LOOP OVER GROUPS
  
  md.data$location_year = sprintf("%s_%s", md.data$COUNTY, md.data$year)
  
  return(md.data)
}

#**# The below functions are redefined in dfmip. Need to think about the flow of code and whether an extra package is needed.

#' Get Day of Year, given year, month, and day
#'
#' COPIED FROM wnv_hlpr.R
#'
#' @param year Input year
#' @param month Input month
#' @param day Input day
#'
#' @return Return ordinal day of the year
#'
get.DOY = function(year, month, day){
  
  # Determine if it is a leap year
  days = get.days(year)
  
  # Sum up days from months
  Jan = 31
  Feb = 28
  Mar = 31
  Apr = 30
  May = 31
  Jun = 30
  Jul = 31
  Aug = 31
  Sep = 30
  Oct = 31
  Nov = 30
  Dec = 31
  
  if (days == 366){  Feb = 29  }
  months = c(Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec)
  
  month.days = 0
  if (month != 1){
    month.days = sum(months[1:(month - 1)])
  }
  
  # Add days from days
  doy = month.days + day
  
  return(doy)
}

#' Function to get number of days in a year (copied from wnv_hlpr.R)
#'
#' Does not work for years 1900 and before or after 2100 (i.e. it does not handle the Century cases)
#'
#' @param year The year to examine
#'
#' @return The number of days in the year (365 for non-leap years, 366 for a leap year)
#'
get.days = function(year){
  
  # Add check - not designed for 1900 leap year, and won't work past 2400
  if (year <= 1900 | year >= 2100){ stop("Not designed years <= 1900 or >= 2100")}
  
  # Days are 365 unless it is a leap year
  days = 365
  if (year %% 4 == 0){
    days = 366
  }
  return(days)
}


#' Split strings using strsplit in a way that can easily be used within sapply
#'
#' @param string The string to be split
#' @param delimiter What to use to split the string
#' @param position Which part of the string to retain. Only one piece can be retained
#' To retain multiple pieces, do a separate function call for each piece
#' @param as.string 0 Result will be numeric, 1 Result will be string
#'
splitter = function(string, delimiter, position, as.string = 0){
  parts = strsplit(string, delimiter)[[1]]
  out = parts[position]
  
  if (as.string == 0){
    out = as.numeric(as.character(out)) # Make sure it is in number format and drop any leading 0's
  }
  if (as.string == 1){  out = as.character(out)  }
  
  return(out)
}


#' District to location
#' 
#' Could be generalized even
#' Copy in dfmip as well
#' 
#' @noRd
district.to.location = function(in.data, data.label, old.name = 'district', new.name = 'location'){
  
  fields = colnames(in.data)
  new.name.regex = sprintf("\\b%s\\b", new.name)
  old.name.regex = sprintf("\\b%s\\b", old.name)
  new.pos = grep(new.name.regex, fields)
  old.pos = grep(old.name.regex, fields)
  
  if (length(new.pos) == 0){
    # If the old name is present and the new name is absent, rename the field
    if (length(old.pos) == 1){
      colnames(in.data)[old.pos] = new.name
      warning(sprintf("%s field missing. Substituting values from %s field", new.name, old.name))
    }
    if (length(old.pos) == 0){
      stop(sprintf("Required 'location' field is missing. Field names are %s", paste(fields, collapse = ', ')))
    }
    if (length(old.pos) > 1){
      stop(sprintf("More than one field returned for %s. Probably an error with the regular expressions. Ideally just include a %s field", old.name, new.name))
    }
  }
  
  new.name_year = sprintf("%s_year", new.name)
  new.name_year.regex = sprintf("\\b%s\\b", new.name_year)
  new.name_year.pos = grep(new.name_year.regex, fields)
  
  if (length(new.name_year.pos) == 0){
    # Search for variables associated with the old year
    old.name_year = sprintf("%s_year", old.name)
    old.name_year.regex = sprintf("\\b%s\\b", old.name_year)
    old.name_year.pos = grep(old.name_year.regex, fields)

    # This check must precede the next, otherwise the new.name_year field will exist!
    if (length(old.name_year.pos) == 0){
      # If the field is not there, generate the location_year field from the location and year fields
      in.data[[new.name_year]] = sprintf("%s_%s", in.data[[new.name]], in.data[["year"]])
      
      # Check that this was successful
      if (length(in.data[[new.name_year]]) == 0){
        stop(sprintf("Something went wrong with generation of the %s_year field for data set %s", new.name, data.label))
      }
    }
    
    if (length(old.name_year.pos) == 1){
      warning(sprintf("%s field missing. Subsituting values from %s field", new.name_year, old.name_year))
      colnames(in.data)[old.name_year.pos] = new.name_year
    }
    
    if (length(old.name_year.pos) > 1){
      stop(sprintf("More than one field returned for %s. Probably an error with the regular expressions. Ideally just include a %s field", old.name_year, new.name_year))
    }
  }
 
  return(in.data) 
}
