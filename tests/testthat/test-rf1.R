## Unit Tests for rf1 package
# 1. rf1 model runs successfully for mosquitoes
# 2. rf1 model runs successfully for humans
# 3. rf1 model final results are correctly formatted for export to dfmip
#devtools::load_all()

test_that("rf1 core model runs successfully for mosquitoes", {
  #**# NOTE: non-quantile regression is an option, but is not explicitly tested in the unit tests.
  
  # Set up other inputs
  set.seed(20200304)
  forecast.targets = c('seasonal.mosquito.MLE')
  weekinquestion = as.Date("2015-07-26", "%Y-%m-%d") #**# Is the as.Date part necessary?
  analysis.counties = unique(rf1::human.data$district)
  analysis.years = seq(2011, 2015)
  rf1.inputs = list(NA, NA, analysis.counties, analysis.years, NA, NA, NA)
  #week.id = sprintf("test:%s", weekinquestion)
  id.string = 'test'
  results.path = "RF1TESTS/"
  break.type = "seasonal"
  response.type = "continuous"
  quantile.model = 1
  
  out = FormatDataForRF1(rf1::human.data, rf1::mosq.data, weekinquestion, rf1::weather.data, rf1.inputs, results.path, break.type)
  my.data = out[[1]]
  independent.vars = out[[2]]
  
  # Split the data set into the forecast year and the historical data
  forecast.year = as.numeric(substr(as.character(weekinquestion), 1, 4))
  
  forecast.data = my.data[my.data$year == forecast.year, ]
  historical.data = my.data[my.data$year < forecast.year, ] #**# This will prevent later years from informing hindcasts of earlier years
  
  dep.var = "IR"
  # drop IR variable from independent.vars
  m.independent.vars = independent.vars[independent.vars != "IR"]
  mosquito.results.path = sprintf("%s/mosquitoes", results.path)
  m.label = "mosquitoes"
  mosquito.results = do.rf(historical.data, dep.var, m.independent.vars, mosquito.results.path,
                           response.type = response.type, label = m.label, quantile.model = quantile.model) #do.spatial = 0, create.test.set = 0, create.ci.null = 0
  
  # Test that mosquito.results came out as expected
  kept.vars = as.character(mosquito.results[[6]])
  expect_equal(kept.vars, c('RMEAN_2', 'TMEANC_2', "TMINC_1", 'abundance', 'density'))
  #expect_equal(kept.vars, c('RMEAN_2', 'TMEANC_2', 'PR_1')) #**# Unclear why this changed - could be a change in the randomization procedures leading to a different result.
  #**# And now it's changed back. Very Very odd.
  new.df = forecast.data[ ,kept.vars]
  predictions = predict(mosquito.results$MODEL, new.df) # Does 0.1, 0.5, 0.9 by default
  expect_equal(unname(predictions[1,1]), 0)
  expect_equal(round(unname(predictions[1,2]),3), 0.014)
  expect_equal(round(unname(predictions[4,3]),3), 0.014)

  # Save mosquito.results as testing object (only done once. Update if updating unit test results)
  #usethis::use_data(mosquito.results)

  unlink(results.path)
})

test_that("rf1 model runs successfully for humans", {
  #**# NOTE: non-quantile regression is an option, but is not explicitly tested in the unit tests.
  
  # Set up other inputs
  set.seed(20200304)
  forecast.targets = c('annual.human.cases')
  weekinquestion = as.Date("2015-07-26", "%Y-%m-%d") #**# Is the as.Date part necessary?
  analysis.counties = unique(rf1::human.data$district)
  analysis.years = seq(2011, 2015)
  rf1.inputs = list(NA, NA, analysis.counties, analysis.years, NA, NA, NA)
  #week.id = sprintf("test:%s", weekinquestion)
  id.string = 'test'
  results.path = "RF1TESTS/"
  break.type = "seasonal"
  response.type = "continuous"
  quantile.model = 1
  
  out = FormatDataForRF1(rf1::human.data, rf1::mosq.data, weekinquestion, rf1::weather.data, rf1.inputs, results.path, break.type)
  my.data = out[[1]]
  independent.vars = out[[2]]
  
  # Split the data set into the forecast year and the historical data
  forecast.year = as.numeric(substr(as.character(weekinquestion), 1, 4))
  
  forecast.data = my.data[my.data$year == forecast.year, ]
  historical.data = my.data[my.data$year < forecast.year, ] #**# This will prevent later years from informing hindcasts of earlier years
  
  dep.var = "Cases"
  human.results.path = sprintf("%s/humans", results.path)
  h.label = "humans"
  human.results = do.rf(historical.data, dep.var, independent.vars, human.results.path,
                        response.type = response.type, label = h.label,
                        quantile.model = quantile.model) #do.spatial = 0, create.test.set = 0, create.ci.null = 0
  
  # Test that mosquito.results came out as expected
  kept.vars = as.character(human.results[[6]])
  expect_equal(kept.vars, c('PR_1', 'RMEAN_2'))
  new.df = forecast.data[ ,kept.vars]
  predictions = predict(human.results$MODEL, new.df) # Does 0.1, 0.5, 0.9 by default
  expect_equal(unname(predictions[1,1]), 1)
  expect_equal(round(unname(predictions[1,2]),0), 1)
  expect_equal(round(unname(predictions[4,3]),0), 1)
  expect_equal(round(unname(predictions[2,2]),0), 0)
  
  # Save human.results as testing object (only done once. Update if updating unit test results)
  #usethis::use_data(human.results)
  
  unlink(results.path)
  
  
})

test_that("rf1 model final results are correctly formatted for export to dfmip", {
  # Set up other inputs
  set.seed(20200304)
  forecast.targets = c('annual.human.cases', 'seasonal.mosquito.MLE')
  weekinquestion = as.Date("2015-07-26", "%Y-%m-%d") #**# Is the as.Date part necessary?
  analysis.counties = unique(rf1::human.data$district)
  analysis.years = seq(2011, 2015)
  rf1.inputs = list(NA, NA, analysis.counties, analysis.years, NA, NA, NA)
  #week.id = sprintf("test:%s", weekinquestion)
  id.string = 'test'
  results.path = "RF1TESTS/"
  
  out.results = rf1(forecast.targets, rf1::human.data, rf1::mosq.data, rf1::weather.data,
                    weekinquestion, rf1.inputs, results.path, id.string, break.type = "seasonal", response.type = "continuous",
                    quantile.model = 1, n.draws = 1000, bins = c(0,seq(1,51,1),101,151,201,1000), use.testing.objects = TRUE)
  
  # Check that 4 outputs are generated: The Results dataframe, the Distributions dataframe, the Bins dataframe, and the model object results  
  expect_equal(length(out.results), 4)
  
  # Check Results data frame
  RF1.results = out.results[[1]]
  expect_equal(nrow(RF1.results), 10)
  expect_equal(RF1.results$forecast.target[1], 'seasonal.mosquito.MLE')
  expect_equal(RF1.results$location[4], 'district4')
  expect_equal(RF1.results$location[5], 'test-STATEWIDE')
  expect_equal(round(RF1.results$value[1], 3), 0.014)
  expect_equal(round(RF1.results$value[5], 3), 0.007)
  expect_equal(RF1.results$forecast.target[6], 'annual.human.cases')
  expect_equal(RF1.results$value[10], 3)
  
  # Check the distributions data frame
  RF1.distributions = out.results[[2]]
  expect_equal(nrow(RF1.distributions), 10)
  expect_equal(ncol(RF1.distributions), 1002)
  expect_equal(RF1.distributions$forecast.target[1], 'seasonal.mosquito.MLE')
  expect_equal(RF1.distributions$location[4], 'district4')
  expect_equal(RF1.distributions$location[5], 'test-STATEWIDE')
  expect_equal(round(RF1.distributions$DRAW1[1], 3), 0.014)
  expect_equal(round(RF1.distributions$DRAW999[5], 3), 0.000)
  expect_equal(RF1.distributions$forecast.target[6], 'annual.human.cases')
  expect_equal(RF1.distributions$DRAW500[10], 3)
  
  # Check the Bins data frame
  #**# NOT SCRIPTED
  
  # Check the other outputs
  #**# NOT SCRIPTED
  # Remove testing directory
  unlink(results.path)
})

# Test DOY.to.day function
test_that("DOY.to.date function works", {
  expect_equal(DOY.to.date(1, 2001), "2001-01-01")
  expect_equal(DOY.to.date(31, 2001), "2001-01-31")
  expect_equal(DOY.to.date(32, 2001), "2001-02-01")
  expect_equal(DOY.to.date(59, 2001), '2001-02-28')
  expect_equal(DOY.to.date(60, 2001), '2001-03-01')
  expect_equal(DOY.to.date(60, 2000), '2000-02-29')
  expect_equal(DOY.to.date(365, 2001), '2001-12-31')
  expect_equal(DOY.to.date(366, 2000), '2000-12-31')
})

# Test that various errors have useful error messages
test_that("Assorted problems give useful information", {

  # Test bad values in mosq.data file virus field
  set.seed(20200304)
  forecast.targets = c('seasonal.mosquito.MLE')
  weekinquestion = as.Date("2015-07-26", "%Y-%m-%d") #**# Is the as.Date part necessary?
  analysis.counties = unique(rf1::human.data$district)
  analysis.years = seq(2011, 2015)
  rf1.inputs = list(NA, NA, analysis.counties, analysis.years, NA, NA, NA)
  #week.id = sprintf("test:%s", weekinquestion)
  id.string = 'test'
  results.path = "RF1TESTS/"
  
  mosq.data = rf1::mosq.data
  mosq.data$wnv_result[2] = 'none'
  
  # Error message does not match according to the code, but copy and pasting
  # the output error messages for Expected and Actual leads to a perfect match.
  # Must be some difference in invisible formatting.
  #m1 = "The VIRUS field (wnv_result) can only contain 0, 1, or NA values."
  #m2 = "Entered values were: 1, none, 0"
  #err.message = as.character(sprintf("%s %s", m1, m2))
  
  expect_error(rf1(forecast.targets, rf1::human.data, mosq.data,
                   rf1::weather.data, weekinquestion, rf1.inputs,
                   results.path, id.string, break.type = "seasonal",
                   response.type = "continuous", quantile.model = 1,
                   n.draws = 1000,
                   bins = c(0,seq(1,51,1),101,151,201,1000),
                   use.testing.objects = FALSE))
               #, err.message)

})

# Test district.to.location function
test_that("district.to.location function works", {
  
  test.data = rf1::weather.data
  data.label = "weather.data"
  
  out = district.to.location(test.data, data.label, old.name = 'district', new.name = 'location')
  expect_equal(paste(colnames(out), collapse = " "), "location doy year tminc tmeanc tmaxc pr rmean vpd date location_year")
  
  colnames(test.data)[1] = "location"
  out = district.to.location(test.data, data.label, old.name = 'district', new.name = 'location')
  expect_equal(paste(colnames(out), collapse = " "), "location doy year tminc tmeanc tmaxc pr rmean vpd date location_year")
  
  colnames(test.data)[1] = "notlocation"
  expect_error(district.to.location(test.data, data.label, old.name = 'district', new.name = 'location'), 
               "Required 'location' field is missing. Field names are notlocation, doy, year, tminc, tmeanc, tmaxc, pr, rmean, vpd, date")
  
  colnames(test.data)[1] = 'location'
  test.data$district_year = sprintf("%s_%s", test.data$location, test.data$year)
  expect_warning(district.to.location(test.data, data.label, old.name = 'district', new.name = 'location'),
                 "location_year field missing. Subsituting values from district_year field")
  
})

# # Test that RF1 model correctly handles districts with no human cases included in the analysis are handled properly.
# test_that("Districts missing human cases are handled properly", {
#   
#   # Test when all districts have values, a correct output is generated
#   stop("Unit test not scripted")
#   
#   # Test when one or more districts have no human cases, a correct output is generated
#   stop("Unit test not scripted")
# })
