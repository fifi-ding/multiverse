library(dplyr)
library(survival)
library(ROCR)

# =============================================================================
# NC Multiverse Analysis with Reduced Universe Sampling
# =============================================================================
# This script performs multiverse analysis on NC prisoner data with the following features:
# - Randomly samples 500 universes from each part (Part 1: 1-3000, Part 2: 3001-6000, Part 3: 6001-9600)
# - Total of 1,500 universes processed instead of 9,600 (83% reduction for faster computation)
# - Calculates recidivism probability, accuracy, and AUC for each universe
# - Uses reproducible random sampling with different seeds for each part
# - Configurable universe count via max_universes_per_part parameter (default: 500)
# - Progress tracking shows sampling information for transparency
# =============================================================================

# ======== NC Combinations =======
generate_nc_garden <- function(data, scaling_parameters, convert_followup_parameters, age_category_parameters, split_parameters, imbalancing_parameters, model_parameters, predictor_parameters, define_recid_parameters, profile_age = 25, profile_gender = "male", profile_race = "caucasian", max_universes_per_part = 500) {
  
  # each row in the resulting data frame will be a unique "universe."
  all_universes <- expand.grid(
    scaling = scaling_parameters,
    convert_followup = convert_followup_parameters,
    split = split_parameters,
    age_category = age_category_parameters,
    imbalancing = imbalancing_parameters,
    model = model_parameters,
    predictor = predictor_parameters,
    define_recid = define_recid_parameters,
    stringsAsFactors = FALSE
  )
  
  # Profile base - use parameters passed from Python
  cat("DEBUG Part 1: Profile parameters received - age:", profile_age, "gender:", profile_gender, "race:", profile_race, "\n")
  high_risk_base <- list(
    AGE = profile_age * 12, # Convert years to months
    TSERVD = 100,
    PRIORS = 15,
    WHITE = ifelse(tolower(profile_race) == "caucasian", 1, 0),
    FELON = 1,
    ALCHY = 1,
    JUNKY = 0,
    PROPTY = 1,
    MALE = ifelse(tolower(profile_gender) == "male", 1, 0),
    RULE = 5,
    MARRIED = 1,
    SCHOOL = 12,
    WORKREL = 0,
    PERSON = 0,
    SUPER = 1
  )
  cat("DEBUG Part 1: High risk base created - AGE:", high_risk_base$AGE, "MALE:", high_risk_base$MALE, "WHITE:", high_risk_base$WHITE, "\n")
  
  # initialize a new column to store the results.
  all_universes$recidivism_prob <- NA
  all_universes$accuracy <- NA
  all_universes$auc <- NA
  
  # Randomly sample universes from Part 1 (indices 1-3000)
  set.seed(42)  # For reproducibility
  part1_indices <- sample(1:3000, size = min(max_universes_per_part, 3000), replace = FALSE)
  universes_part1 <- all_universes[part1_indices, ]
  num_universes_part1 <- nrow(universes_part1)
  
  # Write initial progress message
  progress_file <- "progress_nc_garden.txt"
  total_possible <- nrow(all_universes)
  total_sampled <- max_universes_per_part * 3
  write(paste("Starting analysis: Sampling", total_sampled, "universes from", total_possible, "total possible combinations"), file = progress_file, append = FALSE)
  cat("Starting analysis: Sampling", total_sampled, "universes from", total_possible, "total possible combinations\n")
  
  # Loop through each universe and perform the analysis
  for (i in 1:num_universes_part1) {
    
    # Extract the current universe's procedural choices
    universe <- universes_part1[i, ]
    
    # Check if we've reached the end of valid universes
    if (is.na(universe$scaling)) {
      cat("\n--- Stopping Loop: Reached end of generated universes. ---\n")
      break
    }
    
    # Print a header for the current universe to track progress
    cat("\n--- Running Universe", i, "of", num_universes_part1, "---\n")
    cat("  Scaling PV:", universe$scaling, "\n")
    cat("  Convert Follow up PV:", universe$convert_followup, "\n")
    cat("  Split PV:", universe$split, "\n")
    cat("  Age Category PV:", universe$age_category, "\n")
    cat("  Imbalancing PV:", universe$imbalancing, "\n")
    cat("  Modeling PV:", universe$model, "\n")
    cat("  Predictor PV:", universe$predictor, "\n")
    cat("  Define Recidivism PV:", universe$define_recid, "\n")
    
    # Progress update every 25 universes - write to file for Python to read
    if (i %% 25 == 0 || i == num_universes_part1) {
      progress_file <- "progress_nc_garden.txt"
      total_universes <- max_universes_per_part * 3  # Part1 + Part2 + Part3
      write(paste("Universe", i, "of", total_universes, "completed (sampled from", nrow(all_universes), "total universes)"), file = progress_file, append = FALSE)
      cat("PROGRESS_UPDATE: Universe", i, "of", total_universes, "completed (sampled from", nrow(all_universes), "total universes)\n")
    }
    
    current_data <- data
    high_risk_df <- data.frame(high_risk_base)
    
    # ============= SCALING ==============
    if (universe$scaling == "schmidt") {
      cat("     Applying schmidt's scaling procedure...\n")
      high_risk_df <- high_risk_df |> mutate(
        TSERVD_100 = TSERVD / 100,
        PRIORS_10 = PRIORS / 10,
        SCHOOL_10 = SCHOOL / 10,
        RULE_100 = RULE / 100
      )
      
      current_data <- current_data |> mutate(
        TSERVD_100 = TSERVD / 100,
        PRIORS_10 = PRIORS / 10,
        SCHOOL_10 = SCHOOL / 10,
        RULE_100 = RULE / 100,
      )
    } else if (universe$scaling == "no_scaling") {
      cat("     No scaling applied...\n")
      # Need to create the predictor columns even if not scaled
      high_risk_df <- high_risk_df |> mutate(
        TSERVD_100 = TSERVD,
        PRIORS_10 = PRIORS,
        SCHOOL_10 = SCHOOL,
        RULE_100 = RULE
      )
      current_data <- current_data |> mutate(
        TSERVD_100 = TSERVD,
        PRIORS_10 = PRIORS,
        SCHOOL_10 = SCHOOL,
        RULE_100 = RULE
      )
    }
    
    # ============= CONVERT FOLLOW UP ==============
    if(universe$convert_followup == "strict_def"){
      cat("     Apply Convert Follow Up...\n")
      current_data <- current_data |> mutate(
        TIME = ifelse(RECID == 0 & TIME == 0, FOLLOW, TIME),
        TIME = ifelse(TIME == 0, 0.5, TIME),
        RECID = as.numeric(RECID)
      )
    } else if (universe$convert_followup == "general_def"){
      cat("     Apply Convert Follow Up...\n")
      current_data <- current_data |> mutate(
        TIME = ifelse(RECID == 0, FOLLOW, TIME),
        TIME = ifelse(TIME == 0, 0.5, TIME),
        RECID = as.numeric(RECID)
      )
    }
    
    # ============= DEFINE RECID LOGISTIC ==============
    cat("     Apply LOGISTIC RECID DEFINITION...\n")
    current_data <- current_data |> mutate(
      RECID_1 = ifelse(RECID == 1 & (TIME/12) <  1, 1, 0),
      RECID_2 = ifelse(RECID == 1 & (TIME/12) <  2, 1, 0),
      RECID_3 = ifelse(RECID == 1 & (TIME/12) <  3, 1, 0),
      RECID_4 = ifelse(RECID == 1 & (TIME/12) <  4, 1, 0),
      RECID_5 = ifelse(RECID == 1 & (TIME/12) <  5, 1, 0)
    )
    
    # ============= AGE CATEGORY ==============
    if (universe$age_category == "year"){
      cat("     Apply AGE_YEAR...\n")
      high_risk_df <- high_risk_df |> mutate(AGE_YEAR = round(AGE / 12))
      current_data <- current_data |> mutate(AGE_YEAR = round(AGE / 12))
      
    } else if (universe$age_category == "schmidt"){
      cat("     Apply AGE_1000...\n")
      high_risk_df <- high_risk_df |> mutate(AGE_1000 = AGE / 1000)
      current_data <- current_data |> mutate(AGE_1000 = AGE / 1000)
      
    } else if (universe$age_category == "compas"){
      cat("     Apply AGE_COMPAS...\n")
      high_risk_df <- high_risk_df |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_COMPAS = case_when(
            AGE_YEAR < 25 ~ "less than 25",
            AGE_YEAR >= 25 & AGE_YEAR <= 45 ~ "25~45",
            AGE_YEAR > 45 ~ "over 45"
          )
        )
      
      current_data <- current_data |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_COMPAS = case_when(
          AGE_YEAR < 25 ~ "less than 25",
          AGE_YEAR >= 25 & AGE_YEAR <= 45 ~ "25~45",
          AGE_YEAR > 45 ~ "over 45"
          )
        )
      
    } else if (universe$age_category == "nij"){
      cat("     Apply AGE_NIJ...\n")
      high_risk_df <- high_risk_df |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_NIJ = case_when(
            AGE_YEAR >= 18 & AGE_YEAR <= 22 ~ "18~22",
            AGE_YEAR >= 23 & AGE_YEAR <= 27 ~ "23~27",
            AGE_YEAR >= 28 & AGE_YEAR <= 32 ~ "28~32",
            AGE_YEAR >= 33 & AGE_YEAR <= 37 ~ "33~37",
            AGE_YEAR >= 38 & AGE_YEAR <= 42 ~ "38~42",
            AGE_YEAR >= 43 & AGE_YEAR <= 47 ~ "43~47",
            AGE_YEAR > 48 ~ "over 48"
          )
        )
      current_data <- current_data |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_NIJ = case_when(
            AGE_YEAR >= 18 & AGE_YEAR <= 22 ~ "18~22",
            AGE_YEAR >= 23 & AGE_YEAR <= 27 ~ "23~27",
            AGE_YEAR >= 28 & AGE_YEAR <= 32 ~ "28~32",
            AGE_YEAR >= 33 & AGE_YEAR <= 37 ~ "33~37",
            AGE_YEAR >= 38 & AGE_YEAR <= 42 ~ "38~42",
            AGE_YEAR >= 43 & AGE_YEAR <= 47 ~ "43~47",
            AGE_YEAR > 48 ~ "over 48"
          )
        )
    }
    
    # ============= PREDICTORS ============
    if (universe$scaling == "schmidt") {
      if (universe$predictor == "full") {
        predictor_string <- "TSERVD_100 + PRIORS_10 + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE + RULE_100 + MARRIED + SCHOOL_10 + WORKREL + PERSON + SUPER"
      } else if (universe$predictor == "schmidt") {
        predictor_string <- "TSERVD_100 + PRIORS_10 + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE"
      } else if (universe$predictor == "fair") {
        predictor_string <- "TSERVD_100 + PRIORS_10 + FELON + ALCHY + JUNKY + PROPTY"
      }
    } else { # no_scaling
      if (universe$predictor == "full") {
        predictor_string <- "TSERVD + PRIORS + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE + RULE + MARRIED + SCHOOL + WORKREL + PERSON + SUPER"
      } else if (universe$predictor == "schmidt") {
        predictor_string <- "TSERVD + PRIORS + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE"
      } else if (universe$predictor == "fair") {
        predictor_string <- "TSERVD + PRIORS + FELON + ALCHY + JUNKY + PROPTY"
      }
    }
    
    # Add the appropriate age variable dynamically to the string
    if (universe$age_category == "year") {
      predictor_string <- paste(predictor_string, "AGE_YEAR", sep = " + ")
    } else if (universe$age_category == "schmidt") {
      predictor_string <- paste(predictor_string, "AGE_1000", sep = " + ")
    } else if (universe$age_category == "compas") {
      predictor_string <- paste(predictor_string, "AGE_COMPAS", sep = " + ")
    } else if (universe$age_category == "nij") {
      predictor_string <- paste(predictor_string, "AGE_NIJ", sep = " + ")
    }
    cat("     Final predictor string: ", predictor_string, "\n")
    
    # Final Profile Characteristics
    model_vars <- unlist(strsplit(predictor_string, " \\+ "))
    cat("     MODEL VARS: ", model_vars, "\n")
    high_risk <- high_risk_df[, model_vars, drop = FALSE]
        
    # ============= SPLIT ==============
    if(universe$split == "1:2"){
      analysis_1978 <- current_data |> filter(FILE == 1)
      validation_1978 <- current_data |> filter(FILE == 2)
    } else if (universe$split == "6:4"){
      set.seed(123)
      n <- nrow(current_data)
      train_indices <- sample(seq_len(n), size = 0.6 * n)
      analysis_1978 <- current_data[train_indices, ]
      validation_1978 <- current_data[-train_indices, ]
    } else if (universe$split == "7:3"){
      set.seed(123)
      n <- nrow(current_data)
      train_indices <- sample(seq_len(n), size = 0.7 * n)
      analysis_1978 <- current_data[train_indices, ]
      validation_1978 <- current_data[-train_indices, ]
    } else if (universe$split == "8:2"){
      set.seed(123)
      n <- nrow(current_data)
      train_indices <- sample(seq_len(n), size = 0.8 * n)
      analysis_1978 <- current_data[train_indices, ]
      validation_1978 <- current_data[-train_indices, ]
    } 
    
    # ============= IMBALANCING ==============
    if (universe$imbalancing == "over") {
      cat("     Applying Oversampling method...\n")
      male_sample <- analysis_1978 |> filter(MALE == 1)
      female_sample <- analysis_1978 |> filter(MALE == 0)
      set.seed(123)
      females_oversampled <- female_sample |> sample_n(nrow(male_sample), replace = TRUE)
      analysis_1978 <- bind_rows(male_sample, females_oversampled)
      cat("     Balance check after oversampling:\n")
      #print(analysis_1978 |> count(MALE))
      
    } else if (universe$imbalancing == "under"){
      cat("     Applying Undersampling method...\n")
      male_sample <- analysis_1978 |> filter(MALE == 1)
      female_sample <- analysis_1978 |> filter(MALE == 0)
      # Downsample males to match number of females
      set.seed(123)
      males_undersampled <- male_sample |> sample_n(nrow(female_sample), replace = FALSE)
      # Combine balanced dataset
      analysis_1978 <- bind_rows(males_undersampled, female_sample)
      cat("     Balance check after undersampling:\n")
      #print(current_data |> count(MALE))
      
    } else if (universe$imbalancing == "male only"){
      cat("     Applying Male Only...\n")
      analysis_1978 <- analysis_1978 |> filter(MALE == 1)
      model_vars <- model_vars[model_vars != "MALE"]
    } else if (universe$imbalancing == "female only"){
      cat("     Applying Female Only...\n")
      analysis_1978 <- analysis_1978 |> filter(MALE == 0)
      model_vars <- model_vars[model_vars != "MALE"]
    } else if (universe$imbalancing == "weight"){
      cat("     Applying Weighting...\n")
      analysis_1978$weights <- ifelse(analysis_1978$MALE == 0, 4348 / 270, 1)
    }
    
    # =============  MODEL ==============
    # After filtering, check for and remove constant variables
    if (length(model_vars) > 0) {
      constant_vars <- sapply(analysis_1978[, model_vars, drop = FALSE], function(x) length(unique(x[!is.na(x)])) <= 1)
      model_vars <- model_vars[!constant_vars]
    }
    
    # Rebuild the predictor string from the corrected model_vars
    predictor_string <- paste(model_vars, collapse = " + ")
    
    # Ensure we have at least one predictor
    if (nchar(trimws(predictor_string)) == 0 || length(model_vars) == 0) {
      predictor_string <- "1" # Intercept-only model
    }
    
    #if (universe$imbalancing == "female only" || universe$imbalancing == "male only") {
      # Remove the 'MALE' predictor from the string
      #predictor_string <- gsub("\\+ MALE", "", predictor_string)
    #}
    
    if (universe$model == "survival") {
      # Add the survival response variable to the formula
      cox_formula <- as.formula(paste("Surv(TIME, RECID)", predictor_string, sep = " ~ "))
      
      if(universe$imbalancing == "weight"){
        cat("     Applying Cox Weight Model...\n")
        cox_model <- coxph(cox_formula, data = analysis_1978, weights = analysis_1978$weights)
        surv_pred <- survfit(cox_model, newdata = high_risk)
        
      } else {
        
        # Run the survival model
        cat("     Applying Cox Model...\n")
        #print(cox_formula)
    cox_model <- coxph(cox_formula, data = analysis_1978)
        surv_pred <- survfit(cox_model, newdata = high_risk)
      }
      
      # ========== DEFINE RECID SURVIVAL=======
      time_point <- case_when(
        universe$define_recid == "1yr" ~ 12,
        universe$define_recid == "2yr" ~ 24,
        universe$define_recid == "3yr" ~ 36,
        universe$define_recid == "4yr" ~ 48,
        universe$define_recid == "5yr" ~ 60
      )
      
      #print(universe$define_recid)
      #print(time_point)
      
      # Find the index of the time point closest to the desired time
      time_index <- which.min(abs(surv_pred$time - time_point))
      
      # Extract the survival probability at that index
      S_t <- surv_pred$surv[time_index]
      
      # Calculate the probability of recidivism
      result_prob <- 1 - S_t
      #print(result_prob)
    
    # =============  MODEL ==============
    } else if (universe$model == "logistic") {
      cat("     Applying Logistic Model...\n")
      # Get the name of the response variable dynamically
      recidivism_response_var <- case_when(
        universe$define_recid == "1yr" ~ "RECID_1",
        universe$define_recid == "2yr" ~ "RECID_2",
        universe$define_recid == "3yr" ~ "RECID_3",
        universe$define_recid == "4yr" ~ "RECID_4",
        universe$define_recid == "5yr" ~ "RECID_5",
      )
      
      # Ensure the response variable is a factor with 0 and 1 levels
      analysis_1978[[recidivism_response_var]] <- as.factor(as.integer(analysis_1978[[recidivism_response_var]]))

      # Construct the formula using the dynamic response variable
      logistic_formula <- as.formula(paste(recidivism_response_var, "~", predictor_string))

      # Add weights if specified in the imbalancing method
      if (universe$imbalancing == "weight") {
        logistic_model <- glm(logistic_formula, data = analysis_1978, family = binomial, weights = analysis_1978$weights)
      } else {
        logistic_model <- glm(logistic_formula, data = analysis_1978, family = binomial)
      }
      
      # Predict the recidivism probability for the high-risk profile
      result_prob <- predict(logistic_model, newdata = high_risk, type = "response")
    }
    
    # ================= Model Evaluation on Test Set =================
    accuracy_result <- NA
    auc_result <- NA
    
    if (!is.na(result_prob)) {
      cat("     Evaluating model on test set...\n")
      
      # Apply the same preprocessing to test set
      test_processed <- validation_1978
      
      # Apply the same factor level synchronization
      if (length(model_vars) > 0) {
        test_processed <- droplevels(test_processed)
        
        # Update gender_factor in test set if it was removed due to "only" filtering
        if (universe$imbalancing %in% c("male only", "female only")) {
          test_processed$MALE <- droplevels(as.factor(test_processed$MALE), exclude = if (universe$imbalancing == "male only") {0} else {1})
        }
      }
      
      # Debug: Check test set size after factor alignment
      cat(paste0("     DEBUG: Test set size after factor alignment: ", nrow(test_processed), "\n"))
      
      if (nrow(test_processed) == 0) {
        cat("     WARNING: Test set is empty after factor alignment. Skipping evaluation.\n")
        accuracy_result <- NA
        auc_result <- NA
      } else {
        if (universe$model == "survival") {
          # For survival models, we need to predict survival probability
          if (exists("cox_model") && !is.null(cox_model)) {
            tryCatch({
              # Predict survival probability at the time point
              surv_pred_test <- survfit(cox_model, newdata = test_processed)
              
              # Calculate recidivism probability for each test case
              time_point <- case_when(
                universe$define_recid == "1yr" ~ 12,
                universe$define_recid == "2yr" ~ 24,
                universe$define_recid == "3yr" ~ 36,
                universe$define_recid == "4yr" ~ 48,
                universe$define_recid == "5yr" ~ 60
              )
              
              # Get survival probability at the time point
              surv_summary_test <- summary(surv_pred_test, times = time_point)
              
              if (is.null(surv_summary_test$surv) || length(surv_summary_test$surv) == 0) {
                # Fallback: get the last survival probability
                last_surv <- summary(surv_pred_test)$surv
                if (length(last_surv) > 0) {
                  test_predictions <- rep(1 - last_surv[length(last_surv)], nrow(test_processed))
                } else {
                  test_predictions <- rep(NA, nrow(test_processed))
                }
              } else {
                # Ensure we have the right number of predictions
                n_test <- nrow(test_processed)
                n_pred <- length(surv_summary_test$surv)
                
                if (n_pred == n_test) {
                  test_predictions <- 1 - surv_summary_test$surv
                } else if (n_pred == 1) {
                  # Single prediction for all test cases
                  test_predictions <- rep(1 - surv_summary_test$surv, n_test)
                } else {
                  # Mismatch - use the first prediction for all
                  test_predictions <- rep(1 - surv_summary_test$surv[1], n_test)
                }
              }
              
              # Convert to binary predictions (threshold = 0.5)
              binary_predictions <- ifelse(test_predictions > 0.5, 1, 0)
              
              # Get true labels based on recidivism definition
              recidivism_response_var <- case_when(
                universe$define_recid == "1yr" ~ "RECID_1",
                universe$define_recid == "2yr" ~ "RECID_2",
                universe$define_recid == "3yr" ~ "RECID_3",
                universe$define_recid == "4yr" ~ "RECID_4",
                universe$define_recid == "5yr" ~ "RECID_5"
              )
              true_labels <- test_processed[[recidivism_response_var]]
              
              # Ensure all vectors have the same length
              min_length <- min(length(binary_predictions), length(true_labels))
              binary_predictions <- binary_predictions[1:min_length]
              true_labels <- true_labels[1:min_length]
              test_predictions <- test_predictions[1:min_length]
              
              # Calculate accuracy
              accuracy_result <- mean(binary_predictions == true_labels, na.rm = TRUE)
              
              # Calculate AUC
              auc_result <- tryCatch({
                roc_obj <- ROCR::prediction(test_predictions, true_labels)
                ROCR::performance(roc_obj, "auc")@y.values[[1]]
              }, error = function(e) NA)
              
            }, error = function(e) {
              cat("     Error in survival model evaluation:", conditionMessage(e), "\n")
            })
          }
          
        } else if (universe$model == "logistic") {
          # For logistic models
          if (exists("logistic_model") && !is.null(logistic_model)) {
            tryCatch({
              test_predictions <- predict(logistic_model, newdata = test_processed, type = "response")
              binary_predictions <- ifelse(test_predictions > 0.5, 1, 0)
              
              # Get true labels based on recidivism definition
              recidivism_response_var <- case_when(
                universe$define_recid == "1yr" ~ "RECID_1",
                universe$define_recid == "2yr" ~ "RECID_2",
                universe$define_recid == "3yr" ~ "RECID_3",
                universe$define_recid == "4yr" ~ "RECID_4",
                universe$define_recid == "5yr" ~ "RECID_5"
              )
              true_labels <- test_processed[[recidivism_response_var]]
              
              # Calculate accuracy
              accuracy_result <- mean(binary_predictions == true_labels, na.rm = TRUE)
              
              # Calculate AUC
              auc_result <- tryCatch({
                roc_obj <- ROCR::prediction(test_predictions, true_labels)
                ROCR::performance(roc_obj, "auc")@y.values[[1]]
              }, error = function(e) NA)
              
            }, error = function(e) {
              cat("     Error in logistic model evaluation:", conditionMessage(e), "\n")
            })
          }
        }
      }
      
      cat("     Test Accuracy:", round(accuracy_result, 4), "| AUC:", round(auc_result, 4), "\n")
    }
    
    cat(result_prob)
    universes_part1$recidivism_prob[i] <- result_prob
    universes_part1$accuracy[i] <- accuracy_result
    universes_part1$auc[i] <- auc_result
  }
  
  # Write completion message
  write(paste("Part 1 completed: Processed", num_universes_part1, "of", total_sampled, "total sampled universes"), file = progress_file, append = TRUE)
  cat("Part 1 completed: Processed", num_universes_part1, "of", total_sampled, "total sampled universes\n")
  
  # Save the results of Part 1 to a file
  save(universes_part1, file = "results_part1.RData")
  
  return(universes_part1)
}

# Legacy function for backward compatibility
run_multiverse_analysis <- function(data, preprocessing_methods, split_methods, age_categories, imbalancing_methods, predictor_methods, define_recid_methods, profile_age, profile_gender, profile_race) {
  # Convert old parameter names to new ones
  scaling_parameters = c("schmidt", "no_scaling")
  convert_followup_parameters = c("strict_def", "general_def")
  age_category_parameters = c("year", "schmidt", "compas", "nij")
  split_parameters = split_methods
  imbalancing_parameters = c("under", "over", "male only", "female only", "weight")
  model_parameters = c("survival", "logistic")
  predictor_parameters = c("full", "schmidt", "fair")
  define_recid_parameters = define_recid_methods
  
  return(generate_nc_garden(data, scaling_parameters, convert_followup_parameters, age_category_parameters, split_parameters, imbalancing_parameters, model_parameters, predictor_parameters, define_recid_parameters))
}

# ======== NC Combinations Part 2 =======
generate_nc_garden_part2 <- function(data, scaling_parameters, convert_followup_parameters, age_category_parameters, split_parameters, imbalancing_parameters, model_parameters, predictor_parameters, define_recid_parameters, profile_age = 25, profile_gender = "male", profile_race = "caucasian", max_universes_per_part = 500) {
  
  # each row in the resulting data frame will be a unique "universe."
  all_universes <- expand.grid(
    scaling = scaling_parameters,
    convert_followup = convert_followup_parameters,
    split = split_parameters,
    age_category = age_category_parameters,
    imbalancing = imbalancing_parameters,
    model = model_parameters,
    predictor = predictor_parameters,
    define_recid = define_recid_parameters,
    stringsAsFactors = FALSE
  )
  
  # Profile base - use parameters passed from Python
  cat("DEBUG Part 2: Profile parameters received - age:", profile_age, "gender:", profile_gender, "race:", profile_race, "\n")
  high_risk_base <- list(
    AGE = profile_age * 12, # Convert years to months
    TSERVD = 100,
    PRIORS = 15,
    WHITE = ifelse(tolower(profile_race) == "caucasian", 1, 0),
    FELON = 1,
    ALCHY = 1,
    JUNKY = 0,
    PROPTY = 1,
    MALE = ifelse(tolower(profile_gender) == "male", 1, 0),
    RULE = 5,
    MARRIED = 1,
    SCHOOL = 12,
    WORKREL = 0,
    PERSON = 0,
    SUPER = 1
  )
  cat("DEBUG Part 2: High risk base created - AGE:", high_risk_base$AGE, "MALE:", high_risk_base$MALE, "WHITE:", high_risk_base$WHITE, "\n")
  
  # initialize a new column to store the results.
  all_universes$recidivism_prob <- NA
  all_universes$accuracy <- NA
  all_universes$auc <- NA
  
  # Randomly sample universes from Part 2 (indices 3001-6000)
  set.seed(43)  # Different seed for Part 2
  part2_indices <- sample(3001:6000, size = min(max_universes_per_part, 3000), replace = FALSE)
  universes_part2 <- all_universes[part2_indices, ]
  num_universes_part2 <- nrow(universes_part2)
  
  # Loop through each universe and perform the analysis
  for (i in 1:num_universes_part2) {
    
    # Extract the current universe's procedural choices
    universe <- universes_part2[i, ]
    
    # Print a header for the current universe to track progress
    cat("\n--- Running Universe", i, "of", num_universes_part2, "---\n")
    cat("  Scaling PV:", universe$scaling, "\n")
    cat("  Convert Follow up PV:", universe$convert_followup, "\n")
    cat("  Split PV:", universe$split, "\n")
    cat("  Age Category PV:", universe$age_category, "\n")
    cat("  Imbalancing PV:", universe$imbalancing, "\n")
    cat("  Modeling PV:", universe$model, "\n")
    cat("  Predictor PV:", universe$predictor, "\n")
    cat("  Define Recidivism PV:", universe$define_recid, "\n")
    
    # Progress update every 25 universes - write to file for Python to read
    if (i %% 25 == 0 || i == num_universes_part2) {
      progress_file <- "progress_nc_garden.txt"
      total_universes <- max_universes_per_part * 3  # Part1 + Part2 + Part3
      write(paste("Universe", i + max_universes_per_part, "of", total_universes, "completed (sampled from", nrow(all_universes), "total universes)"), file = progress_file, append = FALSE)
      cat("PROGRESS_UPDATE: Universe", i + max_universes_per_part, "of", total_universes, "completed (sampled from", nrow(all_universes), "total universes)\n")
    }
    
    current_data <- data
    high_risk_df <- data.frame(high_risk_base)
    
    # ============= SCALING ==============
    if (universe$scaling == "schmidt") {
      cat("     Applying schmidt's scaling procedure...\n")
      high_risk_df <- high_risk_df |> mutate(
        TSERVD_100 = TSERVD / 100,
        PRIORS_10 = PRIORS / 10,
        SCHOOL_10 = SCHOOL / 10,
        RULE_100 = RULE / 100
      )
      
      current_data <- current_data |> mutate(
        TSERVD_100 = TSERVD / 100,
        PRIORS_10 = PRIORS / 10,
        SCHOOL_10 = SCHOOL / 10,
        RULE_100 = RULE / 100,
      )
    } else if (universe$scaling == "no_scaling") {
      cat("     No scaling applied...\n")
      # Need to create the predictor columns even if not scaled
      high_risk_df <- high_risk_df |> mutate(
        TSERVD_100 = TSERVD,
        PRIORS_10 = PRIORS,
        SCHOOL_10 = SCHOOL,
        RULE_100 = RULE
      )
      current_data <- current_data |> mutate(
        TSERVD_100 = TSERVD,
        PRIORS_10 = PRIORS,
        SCHOOL_10 = SCHOOL,
        RULE_100 = RULE
      )
    }
    
    # ============= CONVERT FOLLOW UP ==============
    if(universe$convert_followup == "strict_def"){
      cat("     Apply Convert Follow Up...\n")
      current_data <- current_data |> mutate(
        TIME = ifelse(RECID == 0 & TIME == 0, FOLLOW, TIME),
        TIME = ifelse(TIME == 0, 0.5, TIME),
        RECID = as.numeric(RECID)
      )
    } else if (universe$convert_followup == "general_def"){
      cat("     Apply Convert Follow Up...\n")
      current_data <- current_data |> mutate(
        TIME = ifelse(RECID == 0, FOLLOW, TIME),
        TIME = ifelse(TIME == 0, 0.5, TIME),
        RECID = as.numeric(RECID)
      )
    }
    
    # ============= DEFINE RECID LOGISTIC ==============
    cat("     Apply LOGISTIC RECID DEFINITION...\n")
    current_data <- current_data |> mutate(
      RECID_1 = ifelse(RECID == 1 & (TIME/12) <  1, 1, 0),
      RECID_2 = ifelse(RECID == 1 & (TIME/12) <  2, 1, 0),
      RECID_3 = ifelse(RECID == 1 & (TIME/12) <  3, 1, 0),
      RECID_4 = ifelse(RECID == 1 & (TIME/12) <  4, 1, 0),
      RECID_5 = ifelse(RECID == 1 & (TIME/12) <  5, 1, 0)
    )
    
    # ============= AGE CATEGORY ==============
    if (universe$age_category == "year"){
      cat("     Apply AGE_YEAR...\n")
      high_risk_df <- high_risk_df |> mutate(AGE_YEAR = round(AGE / 12))
      current_data <- current_data |> mutate(AGE_YEAR = round(AGE / 12))
      
    } else if (universe$age_category == "schmidt"){
      cat("     Apply AGE_1000...\n")
      high_risk_df <- high_risk_df |> mutate(AGE_1000 = AGE / 1000)
      current_data <- current_data |> mutate(AGE_1000 = AGE / 1000)
      
    } else if (universe$age_category == "compas"){
      cat("     Apply AGE_COMPAS...\n")
      high_risk_df <- high_risk_df |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_COMPAS = case_when(
            AGE_YEAR < 25 ~ "less than 25",
            AGE_YEAR >= 25 & AGE_YEAR <= 45 ~ "25~45",
            AGE_YEAR > 45 ~ "over 45"
          )
        )
      
      current_data <- current_data |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_COMPAS = case_when(
          AGE_YEAR < 25 ~ "less than 25",
          AGE_YEAR >= 25 & AGE_YEAR <= 45 ~ "25~45",
          AGE_YEAR > 45 ~ "over 45"
          )
        )
      
    } else if (universe$age_category == "nij"){
      cat("     Apply AGE_NIJ...\n")
      high_risk_df <- high_risk_df |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_NIJ = case_when(
            AGE_YEAR >= 18 & AGE_YEAR <= 22 ~ "18~22",
            AGE_YEAR >= 23 & AGE_YEAR <= 27 ~ "23~27",
            AGE_YEAR >= 28 & AGE_YEAR <= 32 ~ "28~32",
            AGE_YEAR >= 33 & AGE_YEAR <= 37 ~ "33~37",
            AGE_YEAR >= 38 & AGE_YEAR <= 42 ~ "38~42",
            AGE_YEAR >= 43 & AGE_YEAR <= 47 ~ "43~47",
            AGE_YEAR > 48 ~ "over 48"
          )
        )
      current_data <- current_data |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_NIJ = case_when(
            AGE_YEAR >= 18 & AGE_YEAR <= 22 ~ "18~22",
            AGE_YEAR >= 23 & AGE_YEAR <= 27 ~ "23~27",
            AGE_YEAR >= 28 & AGE_YEAR <= 32 ~ "28~32",
            AGE_YEAR >= 33 & AGE_YEAR <= 37 ~ "33~37",
            AGE_YEAR >= 38 & AGE_YEAR <= 42 ~ "38~42",
            AGE_YEAR >= 43 & AGE_YEAR <= 47 ~ "43~47",
            AGE_YEAR > 48 ~ "over 48"
          )
        )
    }
    
    # ============= PREDICTORS ============
    if (universe$scaling == "schmidt") {
      if (universe$predictor == "full") {
        predictor_string <- "TSERVD_100 + PRIORS_10 + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE + RULE_100 + MARRIED + SCHOOL_10 + WORKREL + PERSON + SUPER"
      } else if (universe$predictor == "schmidt") {
        predictor_string <- "TSERVD_100 + PRIORS_10 + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE"
      } else if (universe$predictor == "fair") {
        predictor_string <- "TSERVD_100 + PRIORS_10 + FELON + ALCHY + JUNKY + PROPTY"
      }
    } else { # no_scaling
      if (universe$predictor == "full") {
        predictor_string <- "TSERVD + PRIORS + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE + RULE + MARRIED + SCHOOL + WORKREL + PERSON + SUPER"
      } else if (universe$predictor == "schmidt") {
        predictor_string <- "TSERVD + PRIORS + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE"
      } else if (universe$predictor == "fair") {
        predictor_string <- "TSERVD + PRIORS + FELON + ALCHY + JUNKY + PROPTY"
      }
    }
    
    # Add the appropriate age variable dynamically to the string
    if (universe$age_category == "year") {
      predictor_string <- paste(predictor_string, "AGE_YEAR", sep = " + ")
    } else if (universe$age_category == "schmidt") {
      predictor_string <- paste(predictor_string, "AGE_1000", sep = " + ")
    } else if (universe$age_category == "compas") {
      predictor_string <- paste(predictor_string, "AGE_COMPAS", sep = " + ")
    } else if (universe$age_category == "nij") {
      predictor_string <- paste(predictor_string, "AGE_NIJ", sep = " + ")
    }
    cat("     Final predictor string: ", predictor_string, "\n")
    
    # Final Profile Characteristics
    model_vars <- unlist(strsplit(predictor_string, " \\+ "))
    cat("     MODEL VARS: ", model_vars, "\n")
    high_risk <- high_risk_df[, model_vars, drop = FALSE]
        
    # ============= SPLIT ==============
    if(universe$split == "1:2"){
      analysis_1978 <- current_data |> filter(FILE == 1)
      validation_1978 <- current_data |> filter(FILE == 2)
    } else if (universe$split == "6:4"){
      set.seed(123)
      n <- nrow(current_data)
      train_indices <- sample(seq_len(n), size = 0.6 * n)
      analysis_1978 <- current_data[train_indices, ]
      validation_1978 <- current_data[-train_indices, ]
    } else if (universe$split == "7:3"){
      set.seed(123)
      n <- nrow(current_data)
      train_indices <- sample(seq_len(n), size = 0.7 * n)
      analysis_1978 <- current_data[train_indices, ]
      validation_1978 <- current_data[-train_indices, ]
    } else if (universe$split == "8:2"){
      set.seed(123)
      n <- nrow(current_data)
      train_indices <- sample(seq_len(n), size = 0.8 * n)
      analysis_1978 <- current_data[train_indices, ]
      validation_1978 <- current_data[-train_indices, ]
    } 
    
    # ============= IMBALANCING ==============
    if (universe$imbalancing == "over") {
      cat("     Applying Oversampling method...\n")
      male_sample <- analysis_1978 |> filter(MALE == 1)
      female_sample <- analysis_1978 |> filter(MALE == 0)
      set.seed(123)
      females_oversampled <- female_sample |> sample_n(nrow(male_sample), replace = TRUE)
      analysis_1978 <- bind_rows(male_sample, females_oversampled)
      
    } else if (universe$imbalancing == "under"){
      cat("     Applying Undersampling method...\n")
      male_sample <- analysis_1978 |> filter(MALE == 1)
      female_sample <- analysis_1978 |> filter(MALE == 0)
      set.seed(123)
      males_undersampled <- male_sample |> sample_n(nrow(female_sample), replace = FALSE)
      analysis_1978 <- bind_rows(males_undersampled, female_sample)
      
    } else if (universe$imbalancing == "male only"){
      cat("     Applying Male Only...\n")
      analysis_1978 <- analysis_1978 |> filter(MALE == 1)
      model_vars <- model_vars[model_vars != "MALE"]
    } else if (universe$imbalancing == "female only"){
      cat("     Applying Female Only...\n")
      analysis_1978 <- analysis_1978 |> filter(MALE == 0)
      model_vars <- model_vars[model_vars != "MALE"]
    } else if (universe$imbalancing == "weight"){
      cat("     Applying Weighting...\n")
      analysis_1978$weights <- ifelse(analysis_1978$MALE == 0, 4348 / 270, 1)
    }
    
    # After filtering, check for and remove constant variables
    if (length(model_vars) > 0) {
      constant_vars <- sapply(analysis_1978[, model_vars, drop = FALSE], function(x) length(unique(x[!is.na(x)])) <= 1)
      model_vars <- model_vars[!constant_vars]
    }
    
    # Rebuild the predictor string from the corrected model_vars
    predictor_string <- paste(model_vars, collapse = " + ")
    
    # Ensure we have at least one predictor
    if (nchar(trimws(predictor_string)) == 0 || length(model_vars) == 0) {
      predictor_string <- "1" # Intercept-only model
    }
    
    # =============  MODEL ==============
    if (universe$model == "survival") {
      # Add the survival response variable to the formula
      cox_formula <- as.formula(paste("Surv(TIME, RECID)", predictor_string, sep = " ~ "))
      
      if(universe$imbalancing == "weight"){
        cat("     Applying Cox Weight Model...\n")
        cox_model <- coxph(cox_formula, data = analysis_1978, weights = analysis_1978$weights)
        surv_pred <- survfit(cox_model, newdata = high_risk)
        
      } else {
        
        # Run the survival model
        cat("     Applying Cox Model...\n")
        cox_model <- coxph(cox_formula, data = analysis_1978)
        surv_pred <- survfit(cox_model, newdata = high_risk)
      }
      
      # ========== DEFINE RECID SURVIVAL=======
      time_point <- case_when(
        universe$define_recid == "1yr" ~ 12,
        universe$define_recid == "2yr" ~ 24,
        universe$define_recid == "3yr" ~ 36,
        universe$define_recid == "4yr" ~ 48,
        universe$define_recid == "5yr" ~ 60
      )
      
      # Find the index of the time point closest to the desired time
      time_index <- which.min(abs(surv_pred$time - time_point))
      
      # Extract the survival probability at that index
      S_t <- surv_pred$surv[time_index]
      
      # Calculate the probability of recidivism
      result_prob <- 1 - S_t
    
    # =============  MODEL ==============
    } else if (universe$model == "logistic") {
      cat("     Applying Logistic Model...\n")
      # Get the name of the response variable dynamically
      recidivism_response_var <- case_when(
        universe$define_recid == "1yr" ~ "RECID_1",
        universe$define_recid == "2yr" ~ "RECID_2",
        universe$define_recid == "3yr" ~ "RECID_3",
        universe$define_recid == "4yr" ~ "RECID_4",
        universe$define_recid == "5yr" ~ "RECID_5",
      )
      
      # Ensure the response variable is a factor with 0 and 1 levels
      analysis_1978[[recidivism_response_var]] <- as.factor(as.integer(analysis_1978[[recidivism_response_var]]))

      # Construct the formula using the dynamic response variable
      logistic_formula <- as.formula(paste(recidivism_response_var, "~", predictor_string))

      # Add weights if specified in the imbalancing method
      if (universe$imbalancing == "weight") {
        logistic_model <- glm(logistic_formula, data = analysis_1978, family = binomial, weights = analysis_1978$weights)
      } else {
        logistic_model <- glm(logistic_formula, data = analysis_1978, family = binomial)
      }
      
      # Predict the recidivism probability for the high-risk profile
      result_prob <- predict(logistic_model, newdata = high_risk, type = "response")
    }
    
    # ================= Model Evaluation on Test Set =================
    accuracy_result <- NA
    auc_result <- NA
    
    if (!is.na(result_prob)) {
      cat("     Evaluating model on test set...\n")
      
      # Apply the same preprocessing to test set
      test_processed <- validation_1978
      
      # Apply the same factor level synchronization
      if (length(model_vars) > 0) {
        test_processed <- droplevels(test_processed)
        
        # Update gender_factor in test set if it was removed due to "only" filtering
        if (universe$imbalancing %in% c("male only", "female only")) {
          test_processed$MALE <- droplevels(as.factor(test_processed$MALE), exclude = if (universe$imbalancing == "male only") {0} else {1})
        }
      }
      
      # Debug: Check test set size after factor alignment
      cat(paste0("     DEBUG: Test set size after factor alignment: ", nrow(test_processed), "\n"))
      
      if (nrow(test_processed) == 0) {
        cat("     WARNING: Test set is empty after factor alignment. Skipping evaluation.\n")
        accuracy_result <- NA
        auc_result <- NA
      } else {
        if (universe$model == "survival") {
          # For survival models, we need to predict survival probability
          if (exists("cox_model") && !is.null(cox_model)) {
            tryCatch({
              # Predict survival probability at the time point
              surv_pred_test <- survfit(cox_model, newdata = test_processed)
              
              # Calculate recidivism probability for each test case
              time_point <- case_when(
                universe$define_recid == "1yr" ~ 12,
                universe$define_recid == "2yr" ~ 24,
                universe$define_recid == "3yr" ~ 36,
                universe$define_recid == "4yr" ~ 48,
                universe$define_recid == "5yr" ~ 60
              )
              
              # Get survival probability at the time point
              surv_summary_test <- summary(surv_pred_test, times = time_point)
              
              if (is.null(surv_summary_test$surv) || length(surv_summary_test$surv) == 0) {
                # Fallback: get the last survival probability
                last_surv <- summary(surv_pred_test)$surv
                if (length(last_surv) > 0) {
                  test_predictions <- rep(1 - last_surv[length(last_surv)], nrow(test_processed))
                } else {
                  test_predictions <- rep(NA, nrow(test_processed))
                }
              } else {
                # Ensure we have the right number of predictions
                n_test <- nrow(test_processed)
                n_pred <- length(surv_summary_test$surv)
                
                if (n_pred == n_test) {
                  test_predictions <- 1 - surv_summary_test$surv
                } else if (n_pred == 1) {
                  # Single prediction for all test cases
                  test_predictions <- rep(1 - surv_summary_test$surv, n_test)
                } else {
                  # Mismatch - use the first prediction for all
                  test_predictions <- rep(1 - surv_summary_test$surv[1], n_test)
                }
              }
              
              # Convert to binary predictions (threshold = 0.5)
              binary_predictions <- ifelse(test_predictions > 0.5, 1, 0)
              
              # Get true labels based on recidivism definition
              recidivism_response_var <- case_when(
                universe$define_recid == "1yr" ~ "RECID_1",
                universe$define_recid == "2yr" ~ "RECID_2",
                universe$define_recid == "3yr" ~ "RECID_3",
                universe$define_recid == "4yr" ~ "RECID_4",
                universe$define_recid == "5yr" ~ "RECID_5"
              )
              true_labels <- test_processed[[recidivism_response_var]]
              
              # Ensure all vectors have the same length
              min_length <- min(length(binary_predictions), length(true_labels))
              binary_predictions <- binary_predictions[1:min_length]
              true_labels <- true_labels[1:min_length]
              test_predictions <- test_predictions[1:min_length]
              
              # Calculate accuracy
              accuracy_result <- mean(binary_predictions == true_labels, na.rm = TRUE)
              
              # Calculate AUC
              auc_result <- tryCatch({
                roc_obj <- ROCR::prediction(test_predictions, true_labels)
                ROCR::performance(roc_obj, "auc")@y.values[[1]]
              }, error = function(e) NA)
              
            }, error = function(e) {
              cat("     Error in survival model evaluation:", conditionMessage(e), "\n")
            })
          }
          
        } else if (universe$model == "logistic") {
          # For logistic models
          if (exists("logistic_model") && !is.null(logistic_model)) {
            tryCatch({
              test_predictions <- predict(logistic_model, newdata = test_processed, type = "response")
              binary_predictions <- ifelse(test_predictions > 0.5, 1, 0)
              
              # Get true labels based on recidivism definition
              recidivism_response_var <- case_when(
                universe$define_recid == "1yr" ~ "RECID_1",
                universe$define_recid == "2yr" ~ "RECID_2",
                universe$define_recid == "3yr" ~ "RECID_3",
                universe$define_recid == "4yr" ~ "RECID_4",
                universe$define_recid == "5yr" ~ "RECID_5"
              )
              true_labels <- test_processed[[recidivism_response_var]]
              
              # Calculate accuracy
              accuracy_result <- mean(binary_predictions == true_labels, na.rm = TRUE)
              
              # Calculate AUC
              auc_result <- tryCatch({
                roc_obj <- ROCR::prediction(test_predictions, true_labels)
                ROCR::performance(roc_obj, "auc")@y.values[[1]]
              }, error = function(e) NA)
              
            }, error = function(e) {
              cat("     Error in logistic model evaluation:", conditionMessage(e), "\n")
            })
          }
        }
      }
      
      cat("     Test Accuracy:", round(accuracy_result, 4), "| AUC:", round(auc_result, 4), "\n")
    }
    
    cat(result_prob)
    universes_part2$recidivism_prob[i] <- result_prob
    universes_part2$accuracy[i] <- accuracy_result
    universes_part2$auc[i] <- auc_result
  }
  
  # Write completion message
  write(paste("Part 2 completed: Processed", num_universes_part2, "of", max_universes_per_part * 3, "total sampled universes"), file = "progress_nc_garden.txt", append = TRUE)
  cat("Part 2 completed: Processed", num_universes_part2, "of", max_universes_per_part * 3, "total sampled universes\n")
  
  # Save the results of Part 2 to a file
  save(universes_part2, file = "results_part2.RData")
  
  return(universes_part2)
}

# ======== NC Combinations Part 3 =======
generate_nc_garden_part3 <- function(data, scaling_parameters, convert_followup_parameters, age_category_parameters, split_parameters, imbalancing_parameters, model_parameters, predictor_parameters, define_recid_parameters, profile_age = 25, profile_gender = "male", profile_race = "caucasian", max_universes_per_part = 500) {
  
  # each row in the resulting data frame will be a unique "universe."
  all_universes <- expand.grid(
    scaling = scaling_parameters,
    convert_followup = convert_followup_parameters,
    split = split_parameters,
    age_category = age_category_parameters,
    imbalancing = imbalancing_parameters,
    model = model_parameters,
    predictor = predictor_parameters,
    define_recid = define_recid_parameters,
    stringsAsFactors = FALSE
  )
  
  # Profile base - use parameters passed from Python
  cat("DEBUG Part 3: Profile parameters received - age:", profile_age, "gender:", profile_gender, "race:", profile_race, "\n")
  high_risk_base <- list(
    AGE = profile_age * 12, # Convert years to months
    TSERVD = 100,
    PRIORS = 15,
    WHITE = ifelse(tolower(profile_race) == "caucasian", 1, 0),
    FELON = 1,
    ALCHY = 1,
    JUNKY = 0,
    PROPTY = 1,
    MALE = ifelse(tolower(profile_gender) == "male", 1, 0),
    RULE = 5,
    MARRIED = 1,
    SCHOOL = 12,
    WORKREL = 0,
    PERSON = 0,
    SUPER = 1
  )
  cat("DEBUG Part 3: High risk base created - AGE:", high_risk_base$AGE, "MALE:", high_risk_base$MALE, "WHITE:", high_risk_base$WHITE, "\n")
  
  # initialize a new column to store the results.
  all_universes$recidivism_prob <- NA
  all_universes$accuracy <- NA
  all_universes$auc <- NA
  
  # Randomly sample universes from Part 3 (indices 6001-9600)
  set.seed(44)  # Different seed for Part 3
  part3_indices <- sample(6001:9600, size = min(max_universes_per_part, 3600), replace = FALSE)
  universes_part3 <- all_universes[part3_indices, ]
  num_universes_part3 <- nrow(universes_part3)
  
  # Loop through each universe and perform the analysis
  for (i in 1:num_universes_part3) {
    
    # Extract the current universe's procedural choices
    universe <- universes_part3[i, ]
    
    # Print a header for the current universe to track progress
    cat("\n--- Running Universe", i, "of", num_universes_part3, "---\n")
    cat("  Scaling PV:", universe$scaling, "\n")
    cat("  Convert Follow up PV:", universe$convert_followup, "\n")
    cat("  Split PV:", universe$split, "\n")
    cat("  Age Category PV:", universe$age_category, "\n")
    cat("  Imbalancing PV:", universe$imbalancing, "\n")
    cat("  Modeling PV:", universe$model, "\n")
    cat("  Predictor PV:", universe$predictor, "\n")
    cat("  Define Recidivism PV:", universe$define_recid, "\n")
    
    # Progress update every 25 universes - write to file for Python to read
    if (i %% 25 == 0 || i == num_universes_part3) {
      progress_file <- "progress_nc_garden.txt"
      total_universes <- max_universes_per_part * 3  # Part1 + Part2 + Part3
      write(paste("Universe", i + (max_universes_per_part * 2), "of", total_universes, "completed (sampled from", nrow(all_universes), "total universes)"), file = progress_file, append = FALSE)
      cat("PROGRESS_UPDATE: Universe", i + (max_universes_per_part * 2), "of", total_universes, "completed (sampled from", nrow(all_universes), "total universes)\n")
    }
    
    current_data <- data
    high_risk_df <- data.frame(high_risk_base)
    
    # ============= SCALING ==============
    if (universe$scaling == "schmidt") {
      cat("     Applying schmidt's scaling procedure...\n")
      high_risk_df <- high_risk_df |> mutate(
        TSERVD_100 = TSERVD / 100,
        PRIORS_10 = PRIORS / 10,
        SCHOOL_10 = SCHOOL / 10,
        RULE_100 = RULE / 100
      )
      
      current_data <- current_data |> mutate(
        TSERVD_100 = TSERVD / 100,
        PRIORS_10 = PRIORS / 10,
        SCHOOL_10 = SCHOOL / 10,
        RULE_100 = RULE / 100,
      )
    } else if (universe$scaling == "no_scaling") {
      cat("     No scaling applied...\n")
      # Need to create the predictor columns even if not scaled
      high_risk_df <- high_risk_df |> mutate(
        TSERVD_100 = TSERVD,
        PRIORS_10 = PRIORS,
        SCHOOL_10 = SCHOOL,
        RULE_100 = RULE
      )
      current_data <- current_data |> mutate(
        TSERVD_100 = TSERVD,
        PRIORS_10 = PRIORS,
        SCHOOL_10 = SCHOOL,
        RULE_100 = RULE
      )
    }
    
    # ============= CONVERT FOLLOW UP ==============
    if(universe$convert_followup == "strict_def"){
      cat("     Apply Convert Follow Up...\n")
      current_data <- current_data |> mutate(
        TIME = ifelse(RECID == 0 & TIME == 0, FOLLOW, TIME),
        TIME = ifelse(TIME == 0, 0.5, TIME),
        RECID = as.numeric(RECID)
      )
    } else if (universe$convert_followup == "general_def"){
      cat("     Apply Convert Follow Up...\n")
      current_data <- current_data |> mutate(
        TIME = ifelse(RECID == 0, FOLLOW, TIME),
        TIME = ifelse(TIME == 0, 0.5, TIME),
        RECID = as.numeric(RECID)
      )
    }
    
    # ============= DEFINE RECID LOGISTIC ==============
    cat("     Apply LOGISTIC RECID DEFINITION...\n")
    current_data <- current_data |> mutate(
      RECID_1 = ifelse(RECID == 1 & (TIME/12) <  1, 1, 0),
      RECID_2 = ifelse(RECID == 1 & (TIME/12) <  2, 1, 0),
      RECID_3 = ifelse(RECID == 1 & (TIME/12) <  3, 1, 0),
      RECID_4 = ifelse(RECID == 1 & (TIME/12) <  4, 1, 0),
      RECID_5 = ifelse(RECID == 1 & (TIME/12) <  5, 1, 0)
    )
    
    # ============= AGE CATEGORY ==============
    if (universe$age_category == "year"){
      cat("     Apply AGE_YEAR...\n")
      high_risk_df <- high_risk_df |> mutate(AGE_YEAR = round(AGE / 12))
      current_data <- current_data |> mutate(AGE_YEAR = round(AGE / 12))
      
    } else if (universe$age_category == "schmidt"){
      cat("     Apply AGE_1000...\n")
      high_risk_df <- high_risk_df |> mutate(AGE_1000 = AGE / 1000)
      current_data <- current_data |> mutate(AGE_1000 = AGE / 1000)
      
    } else if (universe$age_category == "compas"){
      cat("     Apply AGE_COMPAS...\n")
      high_risk_df <- high_risk_df |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_COMPAS = case_when(
            AGE_YEAR < 25 ~ "less than 25",
            AGE_YEAR >= 25 & AGE_YEAR <= 45 ~ "25~45",
            AGE_YEAR > 45 ~ "over 45"
          )
        )
      
      current_data <- current_data |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_COMPAS = case_when(
          AGE_YEAR < 25 ~ "less than 25",
          AGE_YEAR >= 25 & AGE_YEAR <= 45 ~ "25~45",
          AGE_YEAR > 45 ~ "over 45"
          )
        )
      
    } else if (universe$age_category == "nij"){
      cat("     Apply AGE_NIJ...\n")
      high_risk_df <- high_risk_df |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_NIJ = case_when(
            AGE_YEAR >= 18 & AGE_YEAR <= 22 ~ "18~22",
            AGE_YEAR >= 23 & AGE_YEAR <= 27 ~ "23~27",
            AGE_YEAR >= 28 & AGE_YEAR <= 32 ~ "28~32",
            AGE_YEAR >= 33 & AGE_YEAR <= 37 ~ "33~37",
            AGE_YEAR >= 38 & AGE_YEAR <= 42 ~ "38~42",
            AGE_YEAR >= 43 & AGE_YEAR <= 47 ~ "43~47",
            AGE_YEAR > 48 ~ "over 48"
          )
        )
      current_data <- current_data |>
        mutate(AGE_YEAR = round(AGE / 12)) |>
        mutate(
          AGE_NIJ = case_when(
            AGE_YEAR >= 18 & AGE_YEAR <= 22 ~ "18~22",
            AGE_YEAR >= 23 & AGE_YEAR <= 27 ~ "23~27",
            AGE_YEAR >= 28 & AGE_YEAR <= 32 ~ "28~32",
            AGE_YEAR >= 33 & AGE_YEAR <= 37 ~ "33~37",
            AGE_YEAR >= 38 & AGE_YEAR <= 42 ~ "38~42",
            AGE_YEAR >= 43 & AGE_YEAR <= 47 ~ "43~47",
            AGE_YEAR > 48 ~ "over 48"
          )
        )
    }
    
    # ============= PREDICTORS ============
    if (universe$scaling == "schmidt") {
      if (universe$predictor == "full") {
        predictor_string <- "TSERVD_100 + PRIORS_10 + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE + RULE_100 + MARRIED + SCHOOL_10 + WORKREL + PERSON + SUPER"
      } else if (universe$predictor == "schmidt") {
        predictor_string <- "TSERVD_100 + PRIORS_10 + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE"
      } else if (universe$predictor == "fair") {
        predictor_string <- "TSERVD_100 + PRIORS_10 + FELON + ALCHY + JUNKY + PROPTY"
      }
    } else { # no_scaling
      if (universe$predictor == "full") {
        predictor_string <- "TSERVD + PRIORS + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE + RULE + MARRIED + SCHOOL + WORKREL + PERSON + SUPER"
      } else if (universe$predictor == "schmidt") {
        predictor_string <- "TSERVD + PRIORS + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE"
      } else if (universe$predictor == "fair") {
        predictor_string <- "TSERVD + PRIORS + FELON + ALCHY + JUNKY + PROPTY"
      }
    }
    
    # Add the appropriate age variable dynamically to the string
    if (universe$age_category == "year") {
      predictor_string <- paste(predictor_string, "AGE_YEAR", sep = " + ")
    } else if (universe$age_category == "schmidt") {
      predictor_string <- paste(predictor_string, "AGE_1000", sep = " + ")
    } else if (universe$age_category == "compas") {
      predictor_string <- paste(predictor_string, "AGE_COMPAS", sep = " + ")
    } else if (universe$age_category == "nij") {
      predictor_string <- paste(predictor_string, "AGE_NIJ", sep = " + ")
    }
    cat("     Final predictor string: ", predictor_string, "\n")
    
    # Final Profile Characteristics
    model_vars <- unlist(strsplit(predictor_string, " \\+ "))
    cat("     MODEL VARS: ", model_vars, "\n")
    high_risk <- high_risk_df[, model_vars, drop = FALSE]
        
    # ============= SPLIT ==============
    if(universe$split == "1:2"){
      analysis_1978 <- current_data |> filter(FILE == 1)
      validation_1978 <- current_data |> filter(FILE == 2)
    } else if (universe$split == "6:4"){
      set.seed(123)
      n <- nrow(current_data)
      train_indices <- sample(seq_len(n), size = 0.6 * n)
      analysis_1978 <- current_data[train_indices, ]
      validation_1978 <- current_data[-train_indices, ]
    } else if (universe$split == "7:3"){
      set.seed(123)
      n <- nrow(current_data)
      train_indices <- sample(seq_len(n), size = 0.7 * n)
      analysis_1978 <- current_data[train_indices, ]
      validation_1978 <- current_data[-train_indices, ]
    } else if (universe$split == "8:2"){
      set.seed(123)
      n <- nrow(current_data)
      train_indices <- sample(seq_len(n), size = 0.8 * n)
      analysis_1978 <- current_data[train_indices, ]
      validation_1978 <- current_data[-train_indices, ]
    } 
    
    # ============= IMBALANCING ==============
    if (universe$imbalancing == "over") {
      cat("     Applying Oversampling method...\n")
      male_sample <- analysis_1978 |> filter(MALE == 1)
      female_sample <- analysis_1978 |> filter(MALE == 0)
      set.seed(123)
      females_oversampled <- female_sample |> sample_n(nrow(male_sample), replace = TRUE)
      analysis_1978 <- bind_rows(male_sample, females_oversampled)
      
    } else if (universe$imbalancing == "under"){
      cat("     Applying Undersampling method...\n")
      male_sample <- analysis_1978 |> filter(MALE == 1)
      female_sample <- analysis_1978 |> filter(MALE == 0)
      set.seed(123)
      males_undersampled <- male_sample |> sample_n(nrow(female_sample), replace = FALSE)
      analysis_1978 <- bind_rows(males_undersampled, female_sample)
      
    } else if (universe$imbalancing == "male only"){
      cat("     Applying Male Only...\n")
      analysis_1978 <- analysis_1978 |> filter(MALE == 1)
      model_vars <- model_vars[model_vars != "MALE"]
    } else if (universe$imbalancing == "female only"){
      cat("     Applying Female Only...\n")
      analysis_1978 <- analysis_1978 |> filter(MALE == 0)
      model_vars <- model_vars[model_vars != "MALE"]
    } else if (universe$imbalancing == "weight"){
      cat("     Applying Weighting...\n")
      analysis_1978$weights <- ifelse(analysis_1978$MALE == 0, 4348 / 270, 1)
    }
    
    # After filtering, check for and remove constant variables
    if (length(model_vars) > 0) {
      constant_vars <- sapply(analysis_1978[, model_vars, drop = FALSE], function(x) length(unique(x[!is.na(x)])) <= 1)
      model_vars <- model_vars[!constant_vars]
    }
    
    # Rebuild the predictor string from the corrected model_vars
    predictor_string <- paste(model_vars, collapse = " + ")
    
    # Ensure we have at least one predictor
    if (nchar(trimws(predictor_string)) == 0 || length(model_vars) == 0) {
      predictor_string <- "1" # Intercept-only model
    }
    
    # =============  MODEL ==============
    if (universe$model == "survival") {
      # Add the survival response variable to the formula
      cox_formula <- as.formula(paste("Surv(TIME, RECID)", predictor_string, sep = " ~ "))
      
      if(universe$imbalancing == "weight"){
        cat("     Applying Cox Weight Model...\n")
        cox_model <- coxph(cox_formula, data = analysis_1978, weights = analysis_1978$weights)
        surv_pred <- survfit(cox_model, newdata = high_risk)
        
      } else {
        
        # Run the survival model
        cat("     Applying Cox Model...\n")
        cox_model <- coxph(cox_formula, data = analysis_1978)
        surv_pred <- survfit(cox_model, newdata = high_risk)
      }
      
      # ========== DEFINE RECID SURVIVAL=======
      time_point <- case_when(
        universe$define_recid == "1yr" ~ 12,
        universe$define_recid == "2yr" ~ 24,
        universe$define_recid == "3yr" ~ 36,
        universe$define_recid == "4yr" ~ 48,
        universe$define_recid == "5yr" ~ 60
      )
      
      # Find the index of the time point closest to the desired time
      time_index <- which.min(abs(surv_pred$time - time_point))
      
      # Extract the survival probability at that index
      S_t <- surv_pred$surv[time_index]
      
      # Calculate the probability of recidivism
      result_prob <- 1 - S_t
    
    # =============  MODEL ==============
    } else if (universe$model == "logistic") {
      cat("     Applying Logistic Model...\n")
      # Get the name of the response variable dynamically
      recidivism_response_var <- case_when(
        universe$define_recid == "1yr" ~ "RECID_1",
        universe$define_recid == "2yr" ~ "RECID_2",
        universe$define_recid == "3yr" ~ "RECID_3",
        universe$define_recid == "4yr" ~ "RECID_4",
        universe$define_recid == "5yr" ~ "RECID_5",
      )
      
      # Ensure the response variable is a factor with 0 and 1 levels
      analysis_1978[[recidivism_response_var]] <- as.factor(as.integer(analysis_1978[[recidivism_response_var]]))

      # Construct the formula using the dynamic response variable
      logistic_formula <- as.formula(paste(recidivism_response_var, "~", predictor_string))

      # Add weights if specified in the imbalancing method
      if (universe$imbalancing == "weight") {
        logistic_model <- glm(logistic_formula, data = analysis_1978, family = binomial, weights = analysis_1978$weights)
      } else {
        logistic_model <- glm(logistic_formula, data = analysis_1978, family = binomial)
      }
      
      # Predict the recidivism probability for the high-risk profile
      result_prob <- predict(logistic_model, newdata = high_risk, type = "response")
    }
    
    # ================= Model Evaluation on Test Set =================
    accuracy_result <- NA
    auc_result <- NA
    
    if (!is.na(result_prob)) {
      cat("     Evaluating model on test set...\n")
      
      # Apply the same preprocessing to test set
      test_processed <- validation_1978
      
      # Apply the same factor level synchronization
      if (length(model_vars) > 0) {
        test_processed <- droplevels(test_processed)
        
        # Update gender_factor in test set if it was removed due to "only" filtering
        if (universe$imbalancing %in% c("male only", "female only")) {
          test_processed$MALE <- droplevels(as.factor(test_processed$MALE), exclude = if (universe$imbalancing == "male only") {0} else {1})
        }
      }
      
      # Debug: Check test set size after factor alignment
      cat(paste0("     DEBUG: Test set size after factor alignment: ", nrow(test_processed), "\n"))
      
      if (nrow(test_processed) == 0) {
        cat("     WARNING: Test set is empty after factor alignment. Skipping evaluation.\n")
        accuracy_result <- NA
        auc_result <- NA
      } else {
        if (universe$model == "survival") {
          # For survival models, we need to predict survival probability
          if (exists("cox_model") && !is.null(cox_model)) {
            tryCatch({
              # Predict survival probability at the time point
              surv_pred_test <- survfit(cox_model, newdata = test_processed)
              
              # Calculate recidivism probability for each test case
              time_point <- case_when(
                universe$define_recid == "1yr" ~ 12,
                universe$define_recid == "2yr" ~ 24,
                universe$define_recid == "3yr" ~ 36,
                universe$define_recid == "4yr" ~ 48,
                universe$define_recid == "5yr" ~ 60
              )
              
              # Get survival probability at the time point
              surv_summary_test <- summary(surv_pred_test, times = time_point)
              
              if (is.null(surv_summary_test$surv) || length(surv_summary_test$surv) == 0) {
                # Fallback: get the last survival probability
                last_surv <- summary(surv_pred_test)$surv
                if (length(last_surv) > 0) {
                  test_predictions <- rep(1 - last_surv[length(last_surv)], nrow(test_processed))
                } else {
                  test_predictions <- rep(NA, nrow(test_processed))
                }
              } else {
                # Ensure we have the right number of predictions
                n_test <- nrow(test_processed)
                n_pred <- length(surv_summary_test$surv)
                
                if (n_pred == n_test) {
                  test_predictions <- 1 - surv_summary_test$surv
                } else if (n_pred == 1) {
                  # Single prediction for all test cases
                  test_predictions <- rep(1 - surv_summary_test$surv, n_test)
                } else {
                  # Mismatch - use the first prediction for all
                  test_predictions <- rep(1 - surv_summary_test$surv[1], n_test)
                }
              }
              
              # Convert to binary predictions (threshold = 0.5)
              binary_predictions <- ifelse(test_predictions > 0.5, 1, 0)
              
              # Get true labels based on recidivism definition
              recidivism_response_var <- case_when(
                universe$define_recid == "1yr" ~ "RECID_1",
                universe$define_recid == "2yr" ~ "RECID_2",
                universe$define_recid == "3yr" ~ "RECID_3",
                universe$define_recid == "4yr" ~ "RECID_4",
                universe$define_recid == "5yr" ~ "RECID_5"
              )
              true_labels <- test_processed[[recidivism_response_var]]
              
              # Ensure all vectors have the same length
              min_length <- min(length(binary_predictions), length(true_labels))
              binary_predictions <- binary_predictions[1:min_length]
              true_labels <- true_labels[1:min_length]
              test_predictions <- test_predictions[1:min_length]
              
              # Calculate accuracy
              accuracy_result <- mean(binary_predictions == true_labels, na.rm = TRUE)
              
              # Calculate AUC
              auc_result <- tryCatch({
                roc_obj <- ROCR::prediction(test_predictions, true_labels)
                ROCR::performance(roc_obj, "auc")@y.values[[1]]
              }, error = function(e) NA)
              
            }, error = function(e) {
              cat("     Error in survival model evaluation:", conditionMessage(e), "\n")
            })
          }
          
        } else if (universe$model == "logistic") {
          # For logistic models
          if (exists("logistic_model") && !is.null(logistic_model)) {
            tryCatch({
              test_predictions <- predict(logistic_model, newdata = test_processed, type = "response")
              binary_predictions <- ifelse(test_predictions > 0.5, 1, 0)
              
              # Get true labels based on recidivism definition
              recidivism_response_var <- case_when(
                universe$define_recid == "1yr" ~ "RECID_1",
                universe$define_recid == "2yr" ~ "RECID_2",
                universe$define_recid == "3yr" ~ "RECID_3",
                universe$define_recid == "4yr" ~ "RECID_4",
                universe$define_recid == "5yr" ~ "RECID_5"
              )
              true_labels <- test_processed[[recidivism_response_var]]
              
              # Calculate accuracy
              accuracy_result <- mean(binary_predictions == true_labels, na.rm = TRUE)
              
              # Calculate AUC
              auc_result <- tryCatch({
                roc_obj <- ROCR::prediction(test_predictions, true_labels)
                ROCR::performance(roc_obj, "auc")@y.values[[1]]
              }, error = function(e) NA)
              
            }, error = function(e) {
              cat("     Error in logistic model evaluation:", conditionMessage(e), "\n")
            })
          }
        }
      }
      
      cat("     Test Accuracy:", round(accuracy_result, 4), "| AUC:", round(auc_result, 4), "\n")
    }
    
    cat(result_prob)
    universes_part3$recidivism_prob[i] <- result_prob
    universes_part3$accuracy[i] <- accuracy_result
    universes_part3$auc[i] <- auc_result
  }
  
  # Write completion message
  write(paste("Part 3 completed: Processed", num_universes_part3, "of", max_universes_per_part * 3, "total sampled universes"), file = "progress_nc_garden.txt", append = TRUE)
  cat("Part 3 completed: Processed", num_universes_part3, "of", max_universes_per_part * 3, "total sampled universes\n")
  
  # Save the results of Part 3 to a file
  save(universes_part3, file = "results_part3.RData")
  
  return(universes_part3)
}

# Define parameter vectors for the new NC garden approach
scaling_parameters = c("schmidt", "no_scaling")
convert_followup_parameters = c("strict_def", "general_def")
age_category_parameters = c("year", "schmidt", "compas", "nij")
split_parameters = c("1:2", "6:4", "7:3", "8:2")
imbalancing_parameters = c("under", "over", "male only", "female only", "weight")
model_parameters = c("survival", "logistic")
predictor_parameters = c("full", "schmidt", "fair")
define_recid_parameters = c("1yr", "2yr", "3yr", "4yr", "5yr")