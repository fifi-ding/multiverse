library(dplyr)
library(survival)
library(glmnet)
library(xgboost) 
library(tibble)
library(readr)
library(ROCR)
library(caret)

# Load the NIJ GA data
nij_full <- read.csv('nij_full.csv')

# Remove the problematic columns
nij_full <- nij_full |> select(-c(X_v2, X_v3, X_v4))

# Print column names to verify
cat("Available columns:\n")
print(colnames(nij_full))
cat("\nData dimensions:", nrow(nij_full), "rows,", ncol(nij_full), "columns\n")

generate_ga_garden <- function(data, collect_data_parameters, split_parameters, age_category_parameters, race_category_parameters, gender_imbalancing_parameters, model_parameters, predictor_parameters, define_recidivism_parameters, profile_age = 25, profile_gender = "male", profile_race = "white") {
    
    # --- 1. GENERATE UNIVERSE GRID ---
    all_universes <- expand.grid(
        collect_data = collect_data_parameters,
        split = split_parameters,
        age_category = age_category_parameters,
        race_category = race_category_parameters,
        gender_imbalancing = gender_imbalancing_parameters,
        model = model_parameters,
        predictor = predictor_parameters,
        define_recidivism = define_recidivism_parameters,
        stringsAsFactors = TRUE
    )
    
    # --- 2. SETUP HIGH-RISK PROFILE FOR PREDICTION (Test Data) ---
    high_risk_df <- data.frame(
        Gender = ifelse(tolower(profile_gender) == "male", "M", "F"),
        Race = ifelse(tolower(profile_race) == "white", "WHITE", "BLACK"),
        Age_at_Release = "23-27",  # Default age category
        age_numeric = 25,  # For continuous age
        Supervision_Risk_Score_First = 7,  # High risk score
        Education_Level = "Less than HS diploma",
        Dependents = "3 or more",
        Prison_Offense = "Violent/Non-Sex",
        Prison_Years = "More than 3 years",
        Prior_Arrest_Episodes_Felony = factor(6, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10 or more")),
        Prior_Arrest_Episodes_Misd = factor("6 or more", levels = c("0", "1", "2", "3", "4", "5", "6 or more")),
        Prior_Arrest_Episodes_Violent = factor("3 or more", levels = c("0", "1", "2", "3 or more")),
        Prior_Arrest_Episodes_Property = factor("5 or more", levels = c("0", "1", "2", "3", "4", "5 or more")),
        Prior_Arrest_Episodes_Drug = factor("5 or more", levels = c("0", "1", "2", "3", "4", "5 or more")),
        Gang_Affiliated = "Yes",
        Supervision_Level_First = "High",
        race_binary = ifelse(tolower(profile_race) == "white", "White", "Non-White"),
        stringsAsFactors = FALSE
    )
    
    all_universes$recidivism_prob <- NA
    all_universes$accuracy <- NA
    all_universes$auc <- NA
    
    num_universes <- nrow(all_universes)
    
    # Define universes_part1 here, inside the function
    universes_part1 <- all_universes[seq_len(min(100, nrow(all_universes))), ]  # Test with first 100 universes
    num_universes_part1 <- nrow(universes_part1)

    # --- 3. LOOP THROUGH EACH UNIVERSE ---
    for (i in seq_len(num_universes_part1)) {
        
        universe <- universes_part1[i, ]
        
        if (is.na(universe$split)) {
             cat("\n--- Stopping Loop: Reached end of generated universes. ---\n")
             break
        }
        
        cat("\n--- Running Universe", i, "of", num_universes_part1, "---\n")
        cat("  Split:", universe$split, "\n")
        cat("  Age Category:", universe$age_category, "\n")
        cat("  Race Category:", universe$race_category, "\n")
        cat("  Gender Imbalancing:", universe$gender_imbalancing, "\n")
        cat("  Model:", universe$model, "\n")
        cat("  Predictor:", universe$predictor, "\n")
        cat("  Define Recidivism:", universe$define_recidivism, "\n")
        
        current_data <- data
        
        # Convert character variables to factors and numeric as needed
        current_data$Gender <- as.factor(current_data$Gender)
        current_data$Race <- as.factor(current_data$Race)
        current_data$Age_at_Release <- as.factor(current_data$Age_at_Release)
        current_data$Education_Level <- as.factor(current_data$Education_Level)
        current_data$Dependents <- as.factor(current_data$Dependents)
        current_data$Prison_Offense <- as.factor(current_data$Prison_Offense)
        current_data$Prison_Years <- as.factor(current_data$Prison_Years)
        current_data$Gang_Affiliated <- as.factor(current_data$Gang_Affiliated)
        current_data$Supervision_Level_First <- as.factor(current_data$Supervision_Level_First)
        
        # Convert Prior_Arrest_Episodes variables to factors
        current_data$Prior_Arrest_Episodes_Felony <- as.factor(current_data$Prior_Arrest_Episodes_Felony)
        current_data$Prior_Arrest_Episodes_Misd <- as.factor(current_data$Prior_Arrest_Episodes_Misd)
        current_data$Prior_Arrest_Episodes_Violent <- as.factor(current_data$Prior_Arrest_Episodes_Violent)
        current_data$Prior_Arrest_Episodes_Property <- as.factor(current_data$Prior_Arrest_Episodes_Property)
        current_data$Prior_Arrest_Episodes_Drug <- as.factor(current_data$Prior_Arrest_Episodes_Drug)
        
        # Convert recidivism variables to numeric (Yes=1, No=0)
        current_data$Recidivism_Within_3years <- ifelse(current_data$Recidivism_Within_3years == "Yes", 1, 0)
        current_data$Recidivism_Arrest_Year1 <- ifelse(current_data$Recidivism_Arrest_Year1 == "Yes", 1, 0)
        current_data$Recidivism_Arrest_Year2 <- ifelse(current_data$Recidivism_Arrest_Year2 == "Yes", 1, 0)
        current_data$Recidivism_Arrest_Year3 <- ifelse(current_data$Recidivism_Arrest_Year3 == "Yes", 1, 0)
        
        # ================= Age Category Processing =================
        if (universe$age_category == "continuous"){
            cat("    Using continuous age...\n")
            # Extract numeric age from age categories
            current_data <- current_data %>%
                mutate(age_numeric = case_when(
                    Age_at_Release == "18-22" ~ 20,
                    Age_at_Release == "23-27" ~ 25,
                    Age_at_Release == "28-32" ~ 30,
                    Age_at_Release == "33-37" ~ 35,
                    Age_at_Release == "38-42" ~ 40,
                    Age_at_Release == "43-47" ~ 45,
                    Age_at_Release == "48 or older" ~ 50
                ))
        } else if (universe$age_category == "categorical") {
            cat("    Using categorical age...\n")
            # Keep Age_at_Release as is
        }
        
        # ================= Race Category Processing =================
        if (universe$race_category == "binary"){
            cat("    Using binary race (White/Non-White)...\n")
            current_data <- current_data %>%
                mutate(race_binary = ifelse(Race == "WHITE", "White", "Non-White"))
            current_data$race_binary <- as.factor(current_data$race_binary)
        } else if (universe$race_category == "original") {
            cat("    Using original race categories...\n")
            # Keep Race as is
        }
        
        # ================= Predictor Selection =================
        if (universe$predictor == "all_ga"){
            cat("    Using all GA predictors...\n")
            base_predictor_string <- "Gender + Race + Age_at_Release + Supervision_Risk_Score_First + Education_Level + Dependents + Prison_Offense + Prison_Years + Prior_Arrest_Episodes_Felony + Prior_Arrest_Episodes_Misd + Prior_Arrest_Episodes_Violent + Prior_Arrest_Episodes_Property + Prior_Arrest_Episodes_Drug + Gang_Affiliated + Supervision_Level_First"
        } else if (universe$predictor == "fair") { 
            cat("    Using fair predictors only...\n")
            base_predictor_string <- "Age_at_Release + Education_Level + Dependents + Prison_Offense + Prison_Years + Prior_Arrest_Episodes_Felony + Prior_Arrest_Episodes_Misd + Prior_Arrest_Episodes_Violent + Prior_Arrest_Episodes_Property + Prior_Arrest_Episodes_Drug + Gang_Affiliated + Supervision_Level_First"
        } else if (universe$predictor == "lasso") {
            cat("    Using all predictors with lasso selection...\n")
            base_predictor_string <- "Gender + Race + Age_at_Release + Supervision_Risk_Score_First + Education_Level + Dependents + Prison_Offense + Prison_Years + Prior_Arrest_Episodes_Felony + Prior_Arrest_Episodes_Misd + Prior_Arrest_Episodes_Violent + Prior_Arrest_Episodes_Property + Prior_Arrest_Episodes_Drug + Gang_Affiliated + Supervision_Level_First"
        } else {
            base_predictor_string <- "Supervision_Risk_Score_First"
        }
        
        # Add age variable based on age category
        if (universe$age_category == "continuous" && universe$predictor != "fair") {
            base_predictor_string <- paste(base_predictor_string, "age_numeric", sep = " + ")
        } else if (universe$age_category == "categorical" && universe$predictor != "fair") {
            base_predictor_string <- paste(base_predictor_string, "Age_at_Release", sep = " + ")
        }
        
        # Add race variable based on race category
        if (universe$race_category == "binary" && universe$predictor != "fair") {
            base_predictor_string <- paste(base_predictor_string, "race_binary", sep = " + ")
        } else if (universe$race_category == "original" && universe$predictor != "fair") {
            base_predictor_string <- paste(base_predictor_string, "Race", sep = " + ")
        }
        
        cat("    Predictor string: ", base_predictor_string, "\n")
        
        lasso_base_vars <- unlist(strsplit(base_predictor_string, " \\+ "))
        
        # ================= Stratified Split & Imbalancing =================
        set.seed(123)
        
        # Create stratification variable based on key categorical variables
        stratification_vars <- c()
        for (var in lasso_base_vars) {
            if (is.factor(current_data[[var]]) && length(unique(current_data[[var]])) > 1) {
                stratification_vars <- c(stratification_vars, var)
            }
        }
        
        # Also include critical variables that might cause factor level issues
        critical_vars <- c("Gender", "Race", "Age_at_Release", "Education_Level", "Prison_Offense")
        for (var in critical_vars) {
            if (var %in% names(current_data)) {
                if (!is.factor(current_data[[var]])) {
                    current_data[[var]] <- as.factor(current_data[[var]])
                }
                if (length(unique(current_data[[var]])) > 1) {
                    stratification_vars <- c(stratification_vars, var)
                }
            }
        }
        
        # Remove duplicates and ensure we have valid stratification variables
        stratification_vars <- unique(stratification_vars)
        stratification_vars <- stratification_vars[stratification_vars %in% names(current_data)]
        
        cat(paste0("        Available stratification variables: ", paste(stratification_vars, collapse = ", "), "\n"))
        
        # Create stratification groups
        if (length(stratification_vars) > 0) {
            # Use up to 3 key variables for stratification
            key_vars <- stratification_vars[seq_len(min(3, length(stratification_vars)))]
            current_data$strat_group <- interaction(current_data[key_vars], drop = TRUE)
            
            cat(paste0("        Using stratified sampling with variables: ", paste(key_vars, collapse = ", "), "\n"))
            
            # Use stratified sampling
            split_ratio <- switch(as.character(universe$split), "6:4" = 0.6, "7:3" = 0.7, "8:2" = 0.8)
            train_indices <- createDataPartition(current_data$strat_group, p = split_ratio, list = FALSE)
            train <- current_data[train_indices, ]
            test <- current_data[-train_indices, ]
            
            # Remove stratification variable
            train$strat_group <- NULL
            test$strat_group <- NULL
        } else {
            # Fallback to simple random sampling
            n <- nrow(current_data)
            split_ratio <- switch(as.character(universe$split), "6:4" = 0.6, "7:3" = 0.7, "8:2" = 0.8) 
            train_indices <- sample(seq_len(n), size = split_ratio * n)
            train <- current_data[train_indices, ]
            test <- current_data[-train_indices, ]
        }
        
        # Ensure factor levels are consistent between train and test
        for (var in lasso_base_vars) {
            if (is.factor(train[[var]]) && is.factor(test[[var]])) {
                all_levels <- unique(c(levels(train[[var]]), levels(test[[var]])))
                train[[var]] <- factor(train[[var]], levels = all_levels)
                test[[var]] <- factor(test[[var]], levels = all_levels)
            }
        }

        # ================= Gender Imbalancing =================
        if (universe$gender_imbalancing == "over") {
            male_sample <- train %>% filter(Gender == "M")
            female_sample <- train %>% filter(Gender == "F")
            set.seed(123)
            females_oversampled <- female_sample %>% sample_n(nrow(male_sample), replace = TRUE)
            train <- bind_rows(male_sample, females_oversampled)
        } else if (universe$gender_imbalancing == "under"){
            male_sample <- train %>% filter(Gender == "M")
            female_sample <- train %>% filter(Gender == "F")
            set.seed(123)
            males_undersampled <- male_sample %>% sample_n(nrow(female_sample), replace = FALSE)
            train <- bind_rows(males_undersampled, female_sample)
        } else if (universe$gender_imbalancing == "male only"){
            train <- train %>% filter(Gender == "M")
            lasso_base_vars <- lasso_base_vars[lasso_base_vars != "Gender"]
        } else if (universe$gender_imbalancing == "female only"){
            train <- train %>% filter(Gender == "F")
            lasso_base_vars <- lasso_base_vars[lasso_base_vars != "Gender"]
        } else if (universe$gender_imbalancing == "weight"){
            train$weights <- ifelse(train$Gender == "F", sum(train$Gender == "M") / sum(train$Gender == "F"), 1)
        }
        
        # ================= Data Cleaning =================
        initial_train_rows <- nrow(train)
        response_var_name <- if (universe$model == "survival") "Recidivism_Within_3years" else "Recidivism_Within_3years"
        
        # Clean the Response Variable
        train <- train[!is.na(train[[response_var_name]]), ]
        
        if (nrow(train) == 0) {
             cat("        FATAL ERROR: Training set reduced to 0 rows after removing NA in response variable. Skipping.\n")
             universes_part1$recidivism_prob[i] <- NA
             next
        }
        
        # Clean ALL Predictors
        required_vars_for_model <- lasso_base_vars[lasso_base_vars %in% names(train)]
        
        if (length(required_vars_for_model) > 0) {
             train_before_predictor_omit <- nrow(train)
             train <- train %>%
                 na.omit(subset = required_vars_for_model)
             
             if (train_before_predictor_omit != nrow(train)) {
                 cat(paste0("        NOTE: Removed ", train_before_predictor_omit - nrow(train), " rows due to NA in required predictors.\n"))
             }
        }
        
        # Final Row Count Check
        if (initial_train_rows != nrow(train)) {
             cat(paste0("        NOTE: Total ", initial_train_rows - nrow(train), " rows removed due to NA.\n"))
        }
        
        cat(paste("        DEBUG: Rows in training set after NA filter:", nrow(train), "\n"))
        cat(paste("        DEBUG: NA count in", response_var_name, "column:", sum(is.na(train[[response_var_name]])), "\n"))
        
        if (sum(is.na(train[[response_var_name]])) > 0) {
            cat("        FATAL DEBUG CHECK: The response variable is STILL dirty after all cleaning. Skipping.\n")
            universes_part1$recidivism_prob[i] <- NA
            next
        }

        # ================= Lasso Predictor Selection =================
        if (universe$predictor == "lasso" && length(lasso_base_vars) > 0) {
            cat("    Applying Lasso Variable Selection...\n")
            
            # Use droplevels to remove unused factor levels
            train <- droplevels(train)
            
            # Re-check for constant variables after droplevels
            constant_vars_logical <- sapply(train[, lasso_base_vars], function(x) {
                 length(unique(x[!is.na(x)])) <= 1
            })
            lasso_base_vars <- lasso_base_vars[!constant_vars_logical]
            
            current_predictor_string <- paste(lasso_base_vars, collapse = " + ")
            if (nchar(trimws(current_predictor_string)) == 0) {
                current_predictor_string <- "1" # Intercept-only
            }
            
            if (current_predictor_string != "1") {
                lasso_formula_full <- as.formula(paste("Recidivism_Within_3years ~", current_predictor_string))
                
                lasso_prep_result <- tryCatch({
                    model_frame_lasso <- model.frame(lasso_formula_full, data = train, na.action = na.omit)
                    x_train <- model.matrix(lasso_formula_full, data = model_frame_lasso)[, -1, drop = FALSE]
                    y_train <- model.response(model_frame_lasso)
                    
                    if(any(is.na(x_train))) stop("NA detected in the predictor matrix.")
                    if(any(is.na(y_train))) stop("NA detected in the response vector.") 
                    
                    list(x = x_train, y = y_train)
                }, error = function(e) {
                    cat("        ERROR during Lasso data preparation:", conditionMessage(e), "\n")
                    return(NULL)
                })

                if (!is.null(lasso_prep_result)) {
                    x_train <- lasso_prep_result$x
                    y_train <- lasso_prep_result$y
                    
                    if (ncol(x_train) >= 2) {
                        set.seed(42)
                        cv_lasso_result <- tryCatch({
                            cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, type.measure = "auc", nfolds = 10)
                        }, error = function(e) {
                            cat("        Error in cv.glmnet: ", conditionMessage(e), "\n")
                            return(NULL)
                        })

                        if (!is.null(cv_lasso_result)) {
                            lasso_coef <- coef(cv_lasso_result, s = "lambda.min")
                            coef_vector <- as.vector(lasso_coef)
                            
                            if (is.null(coef_vector) || all(is.na(coef_vector))) {
                                retained_dummies <- c()
                            } else {
                                retained_dummies <- rownames(lasso_coef)[abs(coef_vector) > 1e-6]
                                retained_dummies <- retained_dummies[retained_dummies != "(Intercept)"]
                            }
                            
                            original_predictors_retained <- c()
                            
                            if (length(retained_dummies) > 0) {
                                for (orig_var in lasso_base_vars) { 
                                    if (any(grepl(paste0("^", orig_var, "|", orig_var), retained_dummies, fixed = FALSE))) {
                                        original_predictors_retained <- c(original_predictors_retained, orig_var)
                                    }
                                }
                            }
                            
                            if (length(original_predictors_retained) > 0) {
                                current_predictor_string <- paste(original_predictors_retained, collapse = " + ")
                            } else {
                                current_predictor_string <- "1" 
                            }
                            cat("        Lasso selected predictors: ", current_predictor_string, "\n")
                        }
                    }
                }
            }
        } else {
            current_predictor_string <- base_predictor_string
        }
        
        cat("    Final model predictor string: ", current_predictor_string, "\n")

        # ================= Final Factor Level Check =================
        # Remove any factors that have only one level after all processing
        final_predictor_vars <- unlist(strsplit(current_predictor_string, " \\+ "))
        single_level_vars <- c()
        
        for (var in final_predictor_vars) {
            if (var %in% names(train) && is.factor(train[[var]])) {
                if (length(levels(train[[var]])) <= 1) {
                    single_level_vars <- c(single_level_vars, var)
                    cat("        Removing single-level factor:", var, "\n")
                }
            }
        }
        
        # Update predictor string to remove single-level factors
        if (length(single_level_vars) > 0) {
            final_predictor_vars <- final_predictor_vars[!final_predictor_vars %in% single_level_vars]
            if (length(final_predictor_vars) > 0) {
                current_predictor_string <- paste(final_predictor_vars, collapse = " + ")
            } else {
                current_predictor_string <- "1"  # Intercept-only model
            }
            cat("        Updated predictor string after removing single-level factors:", current_predictor_string, "\n")
        }
        
        # Check if we have any predictors left after removing single-level factors
        if (current_predictor_string == "1") {
            cat("        WARNING: Only intercept model available. Skipping this universe.\n")
            universes_part1$recidivism_prob[i] <- NA
            universes_part1$accuracy[i] <- NA
            universes_part1$auc[i] <- NA
            next
        }

        # ================= Final Model Validation =================
        # Check for single-level factors right before model fitting
        final_predictor_vars <- unlist(strsplit(current_predictor_string, " \\+ "))
        single_level_vars <- c()
        
        for (var in final_predictor_vars) {
            if (var %in% names(train) && is.factor(train[[var]])) {
                if (length(levels(train[[var]])) <= 1) {
                    single_level_vars <- c(single_level_vars, var)
                    cat("        Removing single-level factor before model fitting:", var, "\n")
                }
            }
        }
        
        # Update predictor string to remove single-level factors
        if (length(single_level_vars) > 0) {
            final_predictor_vars <- final_predictor_vars[!final_predictor_vars %in% single_level_vars]
            if (length(final_predictor_vars) > 0) {
                current_predictor_string <- paste(final_predictor_vars, collapse = " + ")
            } else {
                current_predictor_string <- "1"  # Intercept-only model
            }
            cat("        Updated predictor string after removing single-level factors:", current_predictor_string, "\n")
        }
        
        # Check if we have any predictors left after removing single-level factors
        if (current_predictor_string == "1") {
            cat("        WARNING: Only intercept model available. Skipping this universe.\n")
            universes_part1$recidivism_prob[i] <- NA
            universes_part1$accuracy[i] <- NA
            universes_part1$auc[i] <- NA
            next
        }

        # ================= Model Fitting =================
        result_prob <- NA
        
        # --- Survival Model ---
        if (universe$model == "survival") {
             cat("    Applying Cox Model...\n")
             
             cox_formula <- as.formula(paste('Recidivism_Within_3years ~', current_predictor_string))
             cox_model <- NULL
             
             cox_model <- tryCatch({
                 if(universe$gender_imbalancing == "weight"){ 
                     coxph(cox_formula, data = train, weights = train$weights, control=coxph.control(iter.max=20))
                 } else {
                     coxph(cox_formula, data = train, control=coxph.control(iter.max=20))
                 }
             }, warning = function(w) {
                     cat("        Warning detected:", conditionMessage(w), "Attempting robust fit...\n")
                     if(universe$gender_imbalancing == "weight"){ 
                         coxph(cox_formula, data = train, weights = train$weights, control=coxph.control(iter.max=5))
                     } else {
                         coxph(cox_formula, data = train, control=coxph.control(iter.max=5))
                     }
             }, error = function(e) {
                     if (grepl("exp overflow", conditionMessage(e))) {
                         cat("        ERROR: Complete separation (exp overflow) detected. Skipping prediction.\n")
                         return(NULL)
                     } else {
                         cat("        ERROR in coxph:", conditionMessage(e), "Skipping prediction.\n")
                         return(NULL)
                     }
             })

             
             if (!is.null(cox_model)) {
                 # For survival models, we predict the probability of recidivism within 3 years
                 # Since we're using Recidivism_Within_3years as the outcome, we can use logistic regression approach
                 result_prob <- predict(cox_model, newdata = high_risk_df, type = "response")
             } else {
                 result_prob <- NA 
             }
             
        # --- Logistic Model ---
        } else if (universe$model == "logistic") {
             cat("    Applying Logistic Model...\n")
             
             logistic_formula <- as.formula(paste("Recidivism_Within_3years ~", current_predictor_string))
             
            logistic_model <- tryCatch({
                 if (universe$gender_imbalancing == "weight") { 
                     glm(logistic_formula, data = train, family = binomial, weights = train$weights)
                 } else {
                     glm(logistic_formula, data = train, family = binomial)
                 }
            }, error = function(e) {
                 if (grepl("fitted probabilities numerically 0 or 1 occurred", conditionMessage(e))) {
                     cat("        ERROR: Complete separation in GLM detected. Skipping prediction.\n")
                     return(NULL)
                 } else if (grepl("contrasts can be applied only to factors with 2 or more levels", conditionMessage(e))) {
                     cat("        ERROR: Factor became constant during model fitting. Skipping prediction.\n")
                     return(NULL)
                 } else {
                     stop(e) 
                 }
            })
             
             if (!is.null(logistic_model)) {
                 result_prob <- predict(logistic_model, newdata = high_risk_df, type = "response")[1]
             } else {
                  result_prob <- NA
             }

        # --- XGBoost Model ---
        } else if (universe$model == "xgboost") {
            cat("    Applying XGBoost Model...\n")

            if (current_predictor_string == "1") {
                cat("        Warning: Intercept-only model not suitable for XGBoost. Skipping.\n")
                result_prob <- NA
            } else {
                
                # Factor Level Synchronization and Constant Factor Removal
                predictor_vars <- unlist(strsplit(current_predictor_string, " \\+ "))
                train_for_mm <- droplevels(train)
                # Ensure response variable is included in train_for_mm
                train_for_mm$Recidivism_Within_3years <- train$Recidivism_Within_3years
                high_risk_for_mm <- high_risk_df
                
                safe_predictor_vars <- c()
                for (var in predictor_vars) {
                    if (is.factor(train_for_mm[[var]]) || is.character(train_for_mm[[var]])) {
                        
                        if(length(levels(train_for_mm[[var]])) <= 1 || length(unique(train_for_mm[[var]][!is.na(train_for_mm[[var]])])) <= 1) {
                            cat(paste0("        Removed constant factor '", var, "' from XGBoost formula.\n"))
                            next 
                        }
                        
                        current_levels <- levels(train_for_mm[[var]])
                        if (is.factor(high_risk_for_mm[[var]])) {
                            high_risk_for_mm[[var]] <- factor(high_risk_for_mm[[var]], levels = current_levels)
                        } else {
                            high_risk_for_mm[[var]] <- factor(high_risk_for_mm[[var]], levels = current_levels)
                        }
                    }
                    safe_predictor_vars <- c(safe_predictor_vars, var) 
                }
                
                current_xgb_predictor_string <- paste(safe_predictor_vars, collapse = " + ")
                if (nchar(trimws(current_xgb_predictor_string)) == 0) {
                    cat("        Error: No predictors remaining after constant check. Skipping.\n")
                    result_prob <- NA
                } else {
                    # Data Preparation
                    xgb_formula <- as.formula(paste("Recidivism_Within_3years ~", current_xgb_predictor_string))
                    
                    tryCatch({
                        # Create Model Frame/Matrix with proper NA handling
                        model_frame_train <- model.frame(xgb_formula, data = train_for_mm, na.action = na.omit) 
                        x_train_xgb <- model.matrix(xgb_formula, data = model_frame_train)
                        x_train_xgb <- x_train_xgb[, colnames(x_train_xgb) != "(Intercept)", drop = FALSE]

                        # Get the response variable from the same model frame to ensure perfect alignment
                        y_train_xgb <- model.response(model_frame_train)
                        
                        # Debug information
                        cat(paste0("        DEBUG XGBoost: Model matrix rows: ", nrow(x_train_xgb), ", Response length: ", length(y_train_xgb), "\n"))
                        cat(paste0("        DEBUG XGBoost: NA count in response: ", sum(is.na(y_train_xgb)), "\n"))
                        
                        # Additional safety checks
                        if (length(y_train_xgb) != nrow(x_train_xgb)) {
                            stop("Mismatch between model matrix rows and response variable length")
                        }
                        if (any(is.na(y_train_xgb))) {
                             stop("Response variable contains NA/NaN after model matrix creation, preventing DMatrix creation.")
                        }

                        train_weights <- if (universe$gender_imbalancing == "weight") {
                            train_for_mm$weights[as.numeric(row.names(model_frame_train))] 
                        } else { NULL }

                        # Create DMatrix
                        dtrain <- xgb.DMatrix(
                            data = x_train_xgb, 
                            label = y_train_xgb, 
                            weight = train_weights
                        )
                        
                        # Train the model
                        params <- list(
                            booster = "gbtree", 
                            objective = "binary:logistic", 
                            eval_metric = "auc",
                            eta = 0.1, max_depth = 4, subsample = 0.8, colsample_bytree = 0.8
                        )
                        
                        xgboost_model <- tryCatch({
                            xgboost::xgb.train(params = params, data = dtrain, nrounds = 50, verbose = 0)
                        }, error = function(e) {
                            cat("        ERROR in xgboost::xgb.train:", conditionMessage(e), "\n")
                            return(NULL)
                        })

                        # Prediction
                        if (!is.null(xgboost_model)) {
                            model_frame_test <- model.frame(xgb_formula, data = high_risk_for_mm, na.action = na.pass)
                            x_test_xgb <- model.matrix(xgb_formula, data = model_frame_test)
                            x_test_xgb <- x_test_xgb[, colnames(x_test_xgb) != "(Intercept)", drop = FALSE]

                            missing_in_test <- setdiff(colnames(x_train_xgb), colnames(x_test_xgb))
                            
                            if (length(missing_in_test) > 0) {
                                x_test_xgb <- cbind(x_test_xgb, matrix(0, nrow = nrow(x_test_xgb), ncol = length(missing_in_test)))
                                colnames(x_test_xgb)[(ncol(x_test_xgb) - length(missing_in_test) + 1):ncol(x_test_xgb)] <- missing_in_test
                            }
                            
                            x_test_xgb <- x_test_xgb[, colnames(x_train_xgb), drop = FALSE]

                            dtest <- xgb.DMatrix(data = x_test_xgb)
                            result_prob <- predict(xgboost_model, dtest)[1]
                        } else {
                            result_prob <- NA
                        }
                    }, error = function(e) {
                        if (grepl("contrasts can be applied only to factors with 2 or more levels", conditionMessage(e))) {
                            cat("        FATAL ERROR: Factor became constant during model matrix creation. Skipping.\n")
                            result_prob <- NA
                        } else if (grepl("Label contains NaN|Response variable", conditionMessage(e))) {
                            cat("        FATAL ERROR: Label contains NaN/NA. Skipping.\n")
                            result_prob <- NA
                        } else {
                            cat("        FATAL ERROR in data prep:", conditionMessage(e), "\n")
                            result_prob <- NA
                        }
                    })
                }
            }
        }
        
        # ================= Model Evaluation on Test Set =================
        accuracy_result <- NA
        auc_result <- NA
        
        if (!is.na(result_prob)) {
            cat("    Evaluating model on test set...\n")
            
            # Apply the same preprocessing to test set
            test_processed <- test
            
            # Apply the same factor level synchronization
            if (length(lasso_base_vars) > 0) {
                test_processed <- droplevels(test_processed)
                
                # Update Gender in test set if it was removed due to "only" filtering
                if (universe$gender_imbalancing %in% c("male only", "female only")) {
                    test_processed$Gender <- droplevels(test_processed$Gender, exclude = if (universe$gender_imbalancing == "male only") {"F"} else {"M"})
                }
            }
            
            # Get the final predictor variables used in the model
            final_predictor_vars <- unlist(strsplit(current_predictor_string, " \\+ "))
            
            # Debug: Check test set size after factor alignment
            cat(paste0("        DEBUG: Test set size after factor alignment: ", nrow(test_processed), "\n"))
            
            if (nrow(test_processed) == 0) {
                cat("        WARNING: Test set is empty after factor alignment. Skipping evaluation.\n")
                accuracy_result <- NA
                auc_result <- NA
            } else {
                if (universe$model == "logistic") {
                    # For logistic models
                    if (!is.null(logistic_model)) {
                        tryCatch({
                            test_predictions <- predict(logistic_model, newdata = test_processed, type = "response")
                            binary_predictions <- ifelse(test_predictions > 0.5, 1, 0)
                            true_labels <- test_processed$Recidivism_Within_3years
                            
                            # Calculate accuracy
                            accuracy_result <- mean(binary_predictions == true_labels, na.rm = TRUE)
                            
                            # Calculate AUC
                            if (requireNamespace("pROC", quietly = TRUE)) {
                                auc_result <- pROC::auc(true_labels, test_predictions, quiet = TRUE)
                            } else {
                                auc_result <- tryCatch({
                                    roc_obj <- ROCR::prediction(test_predictions, true_labels)
                                    ROCR::performance(roc_obj, "auc")@y.values[[1]]
                                }, error = function(e) NA)
                            }
                            
                        }, error = function(e) {
                            cat("        Error in logistic model evaluation:", conditionMessage(e), "\n")
                        })
                    }
                    
                } else if (universe$model == "xgboost") {
                    # For XGBoost models - only evaluate if we have the training matrix
                    if (!is.null(xgboost_model) && exists("x_train_xgb")) {
                        tryCatch({
                        # Prepare test data the same way as training data
                        test_for_mm <- droplevels(test_processed)
                        # Ensure response variable is included in test_for_mm
                        test_for_mm$Recidivism_Within_3years <- test_processed$Recidivism_Within_3years
                        high_risk_for_mm <- high_risk_df
                            
                            safe_predictor_vars <- c()
                            for (var in final_predictor_vars) {
                                if (is.factor(test_for_mm[[var]]) || is.character(test_for_mm[[var]])) {
                                    if(length(levels(test_for_mm[[var]])) <= 1 || length(unique(test_for_mm[[var]][!is.na(test_for_mm[[var]])])) <= 1) {
                                        next 
                                    }
                                    
                                    current_levels <- levels(test_for_mm[[var]])
                                    if (is.factor(high_risk_for_mm[[var]])) {
                                        high_risk_for_mm[[var]] <- factor(high_risk_for_mm[[var]], levels = current_levels)
                                    } else {
                                        high_risk_for_mm[[var]] <- factor(high_risk_for_mm[[var]], levels = current_levels)
                                    }
                                }
                                safe_predictor_vars <- c(safe_predictor_vars, var) 
                            }
                            
                            current_xgb_predictor_string <- paste(safe_predictor_vars, collapse = " + ")
                            
                        if (nchar(trimws(current_xgb_predictor_string)) > 0) {
                            xgb_formula <- as.formula(paste("Recidivism_Within_3years ~", current_xgb_predictor_string))
                            
                            model_frame_test <- model.frame(xgb_formula, data = test_for_mm, na.action = na.omit)
                            x_test_xgb <- model.matrix(xgb_formula, data = model_frame_test)
                            x_test_xgb <- x_test_xgb[, colnames(x_test_xgb) != "(Intercept)", drop = FALSE]
                                
                                # Align with training features
                                missing_in_test <- setdiff(colnames(x_train_xgb), colnames(x_test_xgb))
                                if (length(missing_in_test) > 0) {
                                    x_test_xgb <- cbind(x_test_xgb, matrix(0, nrow = nrow(x_test_xgb), ncol = length(missing_in_test)))
                                    colnames(x_test_xgb)[(ncol(x_test_xgb) - length(missing_in_test) + 1):ncol(x_test_xgb)] <- missing_in_test
                                }
                                x_test_xgb <- x_test_xgb[, colnames(x_train_xgb), drop = FALSE]
                                
                                dtest <- xgb.DMatrix(data = x_test_xgb)
                                test_predictions <- predict(xgboost_model, dtest)
                                
                            # Get true labels aligned with test predictions
                            true_labels <- model.response(model_frame_test)
                                
                                binary_predictions <- ifelse(test_predictions > 0.5, 1, 0)
                                
                                # Calculate accuracy
                                accuracy_result <- mean(binary_predictions == true_labels, na.rm = TRUE)
                                
                                # Calculate AUC
                                if (requireNamespace("pROC", quietly = TRUE)) {
                                    auc_result <- pROC::auc(true_labels, test_predictions, quiet = TRUE)
                                } else {
                                    auc_result <- tryCatch({
                                        roc_obj <- ROCR::prediction(test_predictions, true_labels)
                                        ROCR::performance(roc_obj, "auc")@y.values[[1]]
                                    }, error = function(e) NA)
                                }
                            }
                            
                        }, error = function(e) {
                            cat("        Error in XGBoost model evaluation:", conditionMessage(e), "\n")
                        })
                    }
                }
            } # Close the else block for test_processed check
            
            cat("        Test Accuracy:", round(accuracy_result, 4), "| AUC:", round(auc_result, 4), "\n")
        }
        
        cat(" Recidivism Probability:", result_prob, "\n")
        
        universes_part1$recidivism_prob[i] <- result_prob
        universes_part1$accuracy[i] <- accuracy_result
        universes_part1$auc[i] <- auc_result
    }
    
    save(universes_part1, file = "ga_results_part1.RData")
    return(universes_part1)
}

# Define parameter vectors for the GA garden approach
collect_data_parameters = c("ga")
split_parameters = c("6:4", "7:3", "8:2")
age_category_parameters = c("continuous", "categorical")
race_category_parameters = c("binary", "original")
gender_imbalancing_parameters = c("under", "over", "male only", "female only", "weight")
model_parameters = c("logistic", "xgboost")  # Note: survival model needs special handling for this dataset
predictor_parameters = c("all_ga", "fair", "lasso")
define_recidivism_parameters = c("3yr")

# Example usage
# Call the function (Ensure 'nij_full' is present)
multiverse_ga <- generate_ga_garden(
     nij_full,
     collect_data_parameters,
     split_parameters, 
     age_category_parameters, 
     race_category_parameters, 
     gender_imbalancing_parameters, 
     model_parameters, 
     predictor_parameters,
     define_recidivism_parameters,
     profile_age = 25,
     profile_gender = "male", 
     profile_race = "white"
)

# Print summary of results
cat("Analysis completed!\n")
cat("Number of universes processed:", nrow(multiverse_ga), "\n")
cat("Average accuracy:", round(mean(multiverse_ga$accuracy, na.rm = TRUE), 4), "\n")
cat("Average AUC:", round(mean(multiverse_ga$auc, na.rm = TRUE), 4), "\n")
cat("Results saved to: ga_results_part1.RData\n")