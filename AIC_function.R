library(MuMIn) # to run the dredge function
library(tidyverse) # for wrangling
library(MASS) # negative binomial tests

aic_selection <- function(y_var, data) {
  # List predictors for each subset
  predictors_overall <- c("log_TotalSub", "log_tdr", "log_tdh")
  predictors_forest <- c("tree100m", "log_TotalSub", "log_tdh", "log_patch", "log_tdr")
  predictors_center <- c("log_TotalSub", "log_tdr", "log_tdh", "spec_rad", "log_patch", "road_length_m")

  # create subsets of data
  data_overall <- data
  data_forest <- data %>% 
    filter(Location == "Urban Forest")
  data_center <- data %>% 
    filter(Location == "Urban Center")
  
  # Performing glm and dredge
  run_analysis <- function(predictors, y_var, data_subset) {
    
    formula <- as.formula(paste(y_var, "~", paste(predictors, collapse = " + ")))
    
    # fit the global model
    global <- glm.nb(formula, data = data_subset)
    
    # Run dredge for model selection
    options(na.action = "na.fail") # Required for dredge to run
    model_dredge <- dredge(global, beta = "none", evaluate = TRUE, rank = AICc)
    options(na.action = "na.omit") # Set back to default
    
    # Models within 2 delta AIC
    model_selection <- model.sel(model_dredge)
    model_selection <- model_selection %>% 
      filter(delta < 2)
    
    # Get the top model
    top_model <- get.models(model_dredge, subset = 1)[[1]]
    
    # Extract the summary of the top model
    model_summary <- summary(top_model)
    
    # Extract coefficient table (including standard errors)
    info <- coef(model_summary)
    
    # Extract exponentiated coefficients (and standard errors) for all predictors
    coeffs <- list()
    for (var in predictors) {
      if (var %in% rownames(info)) {
        coeffs[[var]] <- c(exp(info[var, 1]), exp(info[var, 2]))  # Exponentiated coefficients
      } else {
        coeffs[[var]] <- NA  # Handle case where the predictor is not in the top model
      }
    }
    
    # Calculate R squared value
    r_squared <- with(model_summary, 1 - deviance / null.deviance)
    
    # Check VIF for multicollinearity
    # Check VIF for multicollinearity, only if there are more than 1 predictor
    vif_values <- if (length(coef(top_model)) > 2) {
      car::vif(top_model)
    } else {
      NA  # VIF cannot be computed for models with fewer than 2 predictors
    }    
    # return the results as a list
    return(list(
      model_summary = model_summary,
      r_squared = r_squared,
      model_selection = model_selection,
      vif_values = vif_values,
      coefficients = coeffs
    ))
  }
  
  # Run the analysis for each subset
  overall_results <- run_analysis(predictors_overall, y_var, data_overall)
  forest_results <- run_analysis(predictors_forest, y_var, data_forest)
  center_results <- run_analysis(predictors_center, y_var, data_center)
  
  # Return results for all three analyses
  return(list(
    overall = overall_results,
    forest = forest_results,
    center = center_results
  ))
  
}

