library(dplyr)
library(survival)
library(glmnet)
library(xgboost) 
library(tibble)
library(readr)
library(ROCR)
library(caret)

# Load the COMPAS data
compas_2yr_recid <- read_csv("compas.csv")

# Rename duplicate columns to more meaningful names
colnames(compas_2yr_recid)[colnames(compas_2yr_recid) == "priors_count...15"] <- "priors_count"
colnames(compas_2yr_recid)[colnames(compas_2yr_recid) == "priors_count...49"] <- "priors_count_2"
colnames(compas_2yr_recid)[colnames(compas_2yr_recid) == "decile_score...12"] <- "decile_score"
colnames(compas_2yr_recid)[colnames(compas_2yr_recid) == "decile_score...40"] <- "decile_score_2"

# Print column names to verify
cat("Available columns:\n")
print(colnames(compas_2yr_recid))
cat("\nData dimensions:", nrow(compas_2yr_recid), "rows,", ncol(compas_2yr_recid), "columns\n") 

generate_compas_garden <- function(data, collect_data_parameters, time_partition_parameters, age_category_parameters, race_category_parameters, split_parameters, gender_imbalancing_parameters, model_parameters, predictor_parameters) {
    
    # --- 1. GENERATE UNIVERSE GRID ---
    all_universes <- expand.grid(
        collect_data = collect_data_parameters,
        time_partition = time_partition_parameters,
        age_category = age_category_parameters,
        race_category = race_category_parameters,
        split = split_parameters,
        gender_imbalancing = gender_imbalancing_parameters,
        model = model_parameters,
        predictor = predictor_parameters,
        stringsAsFactors = TRUE
    )
    
    # --- 2. SETUP HIGH-RISK PROFILE FOR PREDICTION (Test Data) ---
    high_risk_df <- data.frame(
        age = 25, priors_count = 10, sex = "Male", race = "African-American",
        c_charge_degree = "F", score_text = "High",
        juv_misd_count = 1, 
        juv_other_count = 1,
        r_charge_degree = "(F3)", 
        gender_factor = factor("Male", levels = c("Female", "Male")),
        crime_factor = factor("F", levels = c("F", "M")),  
        score_factor = factor("HighScore", levels = c("LowScore", "HighScore")),
        compas_age = factor("Less than 25", levels = c("Less than 25", "25 - 45", "Greater than 45")), 
        nij_age = "23~27", year_age = 25, 
        nc_race = "BLACK", nij_race = "Black",
        stringsAsFactors = FALSE
    )
    
    all_universes$recidivism_prob <- NA
    all_universes$accuracy <- NA
    all_universes$auc <- NA
    
    num_universes <- nrow(all_universes)
    
    # Define universes_part1 here, inside the function
    universes_part1 <- all_universes[1600:1620, ]  # Test with first 50 universes
    num_universes_part1 <- nrow(universes_part1)

    # --- 3. LOOP THROUGH EACH UNIVERSE ---
    for (i in 1:num_universes_part1) {
        
        universe <- universes_part1[i, ]
        
        if (is.na(universe$time_partition)) {
             cat("\n--- Stopping Loop: Reached end of generated universes. ---\n")
             break
        }
        
        cat("\n--- Running Universe", i, "of", num_universes_part1, "---\n")
        cat("  Time Partition:", universe$time_partition, "\n")
        cat("  Age Category:", universe$age_category, "\n")
        cat("  Race Category:", universe$race_category, "\n")
        cat("  Split:", universe$split, "\n")
        cat("  Gender Imbalancing:", universe$gender_imbalancing, "\n")
        cat("  Model:", universe$model, "\n")
        cat("  Predictor:", universe$predictor, "\n")
        
        current_data <- data
        current_data$start <- as.numeric(current_data$start)
        current_data$end <- as.numeric(current_data$end)
        current_data$event <- as.numeric(current_data$event)
        current_data$two_year_recid <- as.numeric(current_data$two_year_recid)
        
        # ================= Time Partition, Preprocessing, Age/Race Category =================
        if (universe$time_partition == "barenstein"){
            cat("    Applying Barenstein Time Partition...\n")
            current_data$compas_date_numeric <- as.numeric(as.Date(current_data$compas_screening_date))
            current_data$post_april_2014 <- ifelse(current_data$compas_date_numeric > 16161, 1, 0)
            current_data <- current_data[which(current_data$post_april_2014==0),]
        } else if (universe$time_partition == "compas") {
            cat("    Applying Propublica Time Partition...\n")
        }
        
        cat("    Applying Propublica Preprocessing...\n")
        current_data <-  dplyr::select(current_data, age, c_charge_degree, r_charge_degree, race, age_cat, 
                                     score_text, sex, priors_count, days_b_screening_arrest, juv_misd_count, juv_other_count,
                                     decile_score, is_recid, two_year_recid, c_jail_in, c_jail_out, start, end, event) %>% 
             filter(days_b_screening_arrest <= 30) %>% filter(days_b_screening_arrest >= -30) %>%
             filter(is_recid != -1) %>% filter(c_charge_degree != "O") %>% filter(score_text != 'N/A') %>%
             filter(end > start) %>% 
             mutate(
                 length_of_stay = as.numeric(as.Date(c_jail_out)-as.Date(c_jail_in)),
                 crime_factor = factor(c_charge_degree),
                 age_factor = as.factor(age_cat)
             ) %>%
             within(age_factor <- relevel(age_factor, ref = 1)) %>%
             mutate(race_factor = factor(race, labels = c("African-American", "Asian", "Caucasian", "Hispanic", "Native American", "Other"))) %>%
             within(race_factor <- relevel(race_factor, ref = 3)) %>% 
             mutate(gender_factor = factor(sex, labels= c("Female","Male"))) %>%
             within(gender_factor <- relevel(gender_factor, ref = 2)) %>% 
             mutate(score_factor = factor(score_text != "Low", labels = c("LowScore","HighScore")))

        if (universe$age_category == "compas"){
            cat("    Applying COMPAS Age Category...\n")
            current_data <- current_data %>% mutate( compas_age = age_factor)
        } else if (universe$age_category == "nij") {
            cat("    Applying NIJ Age Category...\n")
            current_data <- current_data %>%
                 mutate( nij_age = case_when(
                     age >= 18 & age <= 22 ~ "18~22", age >= 23 & age <= 27 ~ "23~27",
                     age >= 28 & age <= 32 ~ "28~32", age >= 33 & age <= 37 ~ "33~37",
                     age >= 38 & age <= 42 ~ "38~42", age >= 43 & age <= 47 ~ "43~47",
                     age > 48 ~ "over 48",
                     .default = NA_character_
                 ))
             current_data$nij_age <- factor(current_data$nij_age) # Ensure NIJ Age is a factor
        }else if (universe$age_category == "year") {
            current_data <- current_data %>% mutate( year_age = age)
        }

        if (universe$race_category == "nc"){ 
            cat("    Applying NC Race Category...\n")
            current_data <- current_data %>%
                 mutate( nc_race = case_when(
                     race == "African-American" ~ "BLACK",
                     .default = "WHITE"
                 ))
            current_data$nc_race <- factor(current_data$nc_race) # Ensure NC Race is a factor
        } else if (universe$race_category == "nij") { 
            cat("    Applying NIJ Race Category...\n")
            current_data <- current_data %>%
                 mutate( nij_race = case_when(
                     race == "Caucasian" ~ "White",
                     race == "African-American" ~ "Black",
                     .default = NA_character_
                 ))
             current_data$nij_race <- factor(current_data$nij_race) # Ensure NIJ Race is a factor
        }
        
        # ================= Predictor Selection - Initial String =================
        if (universe$predictor == "compas" || universe$predictor == "lasso"){
          cat("    All predictors then apply lasso")
             base_predictor_string <- "sex + juv_misd_count + juv_other_count + c_charge_degree + r_charge_degree + priors_count + crime_factor + score_factor"
        } else if (universe$predictor == "fair") { 
             cat("    Only Fair Predictors")
             base_predictor_string <- "juv_misd_count + juv_other_count + c_charge_degree + r_charge_degree + priors_count + crime_factor + score_factor"
        } else {
             base_predictor_string <- "priors_count"
        }
        
        predictor_string <- base_predictor_string
        
        # Applying the conditional age/race additions
        cat("    Before age and race predictor string: ", predictor_string, "\n")
        if (universe$age_category == "year" && universe$predictor != "fair") { predictor_string <- paste(predictor_string, "year_age", sep = " + ")} 
        else if (universe$age_category == "compas" && universe$predictor != "fair") { predictor_string <- paste(predictor_string, "compas_age", sep = " + ")} 
        else if (universe$age_category == "nij" && universe$predictor != "fair") { predictor_string <- paste(predictor_string, "nij_age", sep = " + ")}
        
        if (universe$race_category == "nc" && universe$predictor != "fair") { predictor_string <- paste(predictor_string, "nc_race", sep = " + ")} 
        else if (universe$race_category == "nij" && universe$predictor != "fair") { predictor_string <- paste(predictor_string, "nij_race", sep = " + ")}
        
        cat("    After age and race predictor string: ", predictor_string, "\n")
        
        lasso_base_vars <- unlist(strsplit(predictor_string, " \\+ "))
        
        # ================= Stratified Split & Imbalancing =================
        set.seed(123)
        
        # Create stratification variable based on key categorical variables
        # This ensures both train and test sets have similar distributions
        stratification_vars <- c()
        for (var in lasso_base_vars) {
            if (is.factor(current_data[[var]]) && length(unique(current_data[[var]])) > 1) {
                stratification_vars <- c(stratification_vars, var)
            }
        }
        
        # Also include ALL categorical variables that might cause factor level issues
        # These are critical for preventing "new levels" errors
        critical_vars <- c("r_charge_degree", "c_charge_degree", "sex", "race", "age_cat")
        for (var in critical_vars) {
            if (var %in% names(current_data)) {
                # Convert to factor if not already, then check if it has multiple levels
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
            # Prioritize critical variables that cause factor level issues
            critical_vars <- c("r_charge_degree", "c_charge_degree", "sex", "race", "age_cat")
            priority_vars <- stratification_vars[stratification_vars %in% critical_vars]
            other_vars <- stratification_vars[!stratification_vars %in% critical_vars]
            
            # Use up to 3 key variables, prioritizing critical ones
            key_vars <- c(priority_vars, other_vars)[1:min(3, length(stratification_vars))]
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
            # Fallback to simple random sampling if no categorical variables
            n <- nrow(current_data)
            split_ratio <- switch(as.character(universe$split), "6:4" = 0.6, "7:3" = 0.7, "8:2" = 0.8) 
            train_indices <- sample(seq_len(n), size = split_ratio * n)
            train <- current_data[train_indices, ]
            test <- current_data[-train_indices, ]
        }
        
        # Ensure factor levels are consistent between train and test
        for (var in lasso_base_vars) {
            if (is.factor(train[[var]]) && is.factor(test[[var]])) {
                # Get all unique levels from both train and test
                all_levels <- unique(c(levels(train[[var]]), levels(test[[var]])))
                
                # Set both train and test to have the same levels
                train[[var]] <- factor(train[[var]], levels = all_levels)
                test[[var]] <- factor(test[[var]], levels = all_levels)
            }
        }

        if (universe$gender_imbalancing == "over") {
            male_sample <- train %>% filter(gender_factor == "Male")
            female_sample <- train %>% filter(gender_factor == "Female")
            set.seed(123)
            females_oversampled <- female_sample %>% sample_n(nrow(male_sample), replace = TRUE)
            train <- bind_rows(male_sample, females_oversampled)
        } else if (universe$gender_imbalancing == "under"){
            male_sample <- train %>% filter(gender_factor == "Male")
            female_sample <- train %>% filter(gender_factor == "Female")
            set.seed(123)
            males_undersampled <- male_sample %>% sample_n(nrow(female_sample), replace = FALSE)
            train <- bind_rows(males_undersampled, female_sample)
        } else if (universe$gender_imbalancing == "male only"){
            train <- train %>% filter(gender_factor == "Male")
            lasso_base_vars <- lasso_base_vars[lasso_base_vars != "gender_factor"]
            lasso_base_vars <- lasso_base_vars[lasso_base_vars != "sex"]
        } else if (universe$gender_imbalancing == "female only"){
            train <- train %>% filter(gender_factor == "Female")
            lasso_base_vars <- lasso_base_vars[lasso_base_vars != "gender_factor"]
            lasso_base_vars <- lasso_base_vars[lasso_base_vars != "sex"]
        } else if (universe$gender_imbalancing == "weight"){
            train$weights <- ifelse(train$gender_factor == "Female", sum(train$gender_factor == "Male") / sum(train$gender_factor == "Female"), 1)
        }
        
        # ----------------------------------------------------------------------------------
        # **CRITICAL FIX BLOCK: Robust NA handling for Label/Response Variable**
        # ----------------------------------------------------------------------------------
        
        initial_train_rows <- nrow(train)
        response_var_name <- if (universe$model == "survival") "event" else "two_year_recid"
        
        # 1. **FIRST PRIORITY: Clean the Response Variable (Label)**
        # Use base R indexing for maximum robustness against NA in the response variable.
        train <- train[!is.na(train[[response_var_name]]), ]
        
        if (nrow(train) == 0) {
             cat("        FATAL ERROR: Training set reduced to 0 rows after removing NA in response variable. Skipping.\n")
             universes_part1$recidivism_prob[i] <- NA
             next
        }
        
        # 2. **SECOND PRIORITY: Clean ALL Predictors (Current and Final)**
        # List of all variables that are currently in the base model string
        required_vars_for_model <- lasso_base_vars[lasso_base_vars %in% names(train)]
        
        if (length(required_vars_for_model) > 0) {
             train_before_predictor_omit <- nrow(train)
             # Crucially, remove NAs from ALL predictors that WILL be used for model.matrix
             train <- train %>%
                 na.omit(subset = required_vars_for_model)
             
             if (train_before_predictor_omit != nrow(train)) {
                 cat(paste0("        NOTE: Removed ", train_before_predictor_omit - nrow(train), " rows due to NA in required predictors.\n"))
             }
        }
        
        # Final Row Count Check/Report
        if (initial_train_rows != nrow(train)) {
             cat(paste0("        NOTE: Total ", initial_train_rows - nrow(train), " rows removed due to NA.\n"))
        }
        
        # --- DEBUG CHECK: Response Variable Status ---
        cat(paste("        DEBUG: Rows in training set after NA filter:", nrow(train), "\n"))
        cat(paste("        DEBUG: NA count in", response_var_name, "column:", sum(is.na(train[[response_var_name]])), "\n"))
        cat(paste("        DEBUG: Unique values in", response_var_name, "column:", paste(sort(unique(train[[response_var_name]])), collapse = ", "), "\n"))
        if (sum(is.na(train[[response_var_name]])) > 0) {
            cat("        FATAL DEBUG CHECK: The response variable is STILL dirty after all cleaning. Skipping.\n")
            universes_part1$recidivism_prob[i] <- NA
            next
        }
        # ----------------------------------------------------------------------------------

        # 3. **Lasso/Model Predictor Cleanup (Factor Levels)**
        if (length(lasso_base_vars) > 0) {
            # 1. Use droplevels to remove unused factor levels
            train <- droplevels(train)
            
            # 2. Re-check for constant variables after droplevels/subsetting
            constant_vars_logical <- sapply(train[, lasso_base_vars], function(x) {
                 length(unique(x[!is.na(x)])) <= 1
            })
            lasso_base_vars <- lasso_base_vars[!constant_vars_logical]
            
            # 3. Update gender_factor in high_risk_df if it was removed due to "only" filtering
            if (universe$gender_imbalancing %in% c("male only", "female only")) {
                high_risk_df$gender_factor <- droplevels(high_risk_df$gender_factor, exclude = if (universe$gender_imbalancing == "male only") {"Female"} else {"Male"})
            }
        }
        
        current_predictor_string <- paste(lasso_base_vars, collapse = " + ")
        if (nchar(trimws(current_predictor_string)) == 0) {
            current_predictor_string <- "1" # Intercept-only
        }
        cat("    Pre-Lasso model predictor string: ", current_predictor_string, "\n")
        
        # ----------------------------------------------------------------------------------
        # End of CRITICAL FIX BLOCK
        # ----------------------------------------------------------------------------------

        # ================= Lasso Predictor Selection (FIXED BLOCK) =================
        if (universe$predictor == "lasso" && current_predictor_string != "1") {
            cat("    Applying Lasso Variable Selection...\n")
            
            # The full formula MUST now include the response variable (two_year_recid)
            lasso_formula_full <- as.formula(paste("two_year_recid ~", current_predictor_string))
            
            lasso_prep_result <- tryCatch({
                
                # 1. Create a model frame and use na.action = na.omit for final, perfect alignment
                model_frame_lasso <- model.frame(lasso_formula_full, data = train, na.action = na.omit)
                
                # 2. Extract the Predictor Matrix (X)
                # The -1 removes the intercept term, as glmnet handles it implicitly
                x_train <- model.matrix(lasso_formula_full, data = model_frame_lasso)[, -1, drop = FALSE]
                
                # 3. Extract the Response Vector (Y) directly from the model frame (guaranteed alignment)
                y_train <- model.response(model_frame_lasso)
                
                if(any(is.na(x_train))) stop("NA detected in the predictor matrix.")
                if(any(is.na(y_train))) stop("NA detected in the response vector.") 
                
                list(x = x_train, y = y_train)
            }, error = function(e) {
                # This error block should now only catch non-NA related issues
                cat("        ERROR during Lasso data preparation (model.matrix):", conditionMessage(e), "\n")
                return(NULL)
            })

            if (is.null(lasso_prep_result)) {
                cat("        Lasso data preparation failed, retaining original predictors.\n")
            } else {
                x_train <- lasso_prep_result$x
                y_train <- lasso_prep_result$y
                
                if (ncol(x_train) < 2) {
                    cat("        Warning: Too few predictors remaining for Lasso selection, skipping.\n")
                } else {
                    set.seed(42)
                    cv_lasso_result <- tryCatch({
                        # cv.glmnet is now safe because y_train is guaranteed NA-free and aligned
                        cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, type.measure = "auc", nfolds = 10)
                    }, error = function(e) {
                        cat("        Error in cv.glmnet: ", conditionMessage(e), "\n")
                        return(NULL)
                    })

                    if (is.null(cv_lasso_result)) {
                        cat("        Lasso failed, retaining all original predictors for robustness.\n")
                    } else {
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
                        cat("        Lasso selected original predictors: ", current_predictor_string, "\n")
                    }
                }
            }
        }
        
        cat("    Final model predictor string: ", current_predictor_string, "\n")

        # ================= Model Fitting =================
        result_prob <- NA
        
        # --- Survival Model ---
        if (universe$model == "survival") {
             cat("    Applying Cox Model...\n")
            
            cox_formula <- as.formula(paste('Surv(start, end, event, type="counting")', current_predictor_string, sep = " ~ "))
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
                surv_pred <- survfit(cox_model, newdata = high_risk_df)
                surv_summary <- summary(surv_pred, times=730)
                
                if (is.null(surv_summary$surv) || length(surv_summary$surv) == 0) {
                     last_surv <- summary(surv_pred)$surv
                     if (length(last_surv) > 0) {
                          result_prob <- 1 - last_surv[length(last_surv)]
                     }
                } else {
                     result_prob <- 1 - surv_summary$surv
                }
            } else {
                result_prob <- NA 
            }
            
        # --- Logistic Model ---
        } else if (universe$model == "logistic") {
             cat("    Applying Logistic Model...\n")
            
            recidivism_response_var <- "two_year_recid" 
            logistic_formula <- as.formula(paste(recidivism_response_var, "~", current_predictor_string))
            
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
        } else if (universe$model == "xg_boost") {
            cat("    Applying XGBoost Model...\n")

            if (current_predictor_string == "1") {
                cat("        Warning: Intercept-only model not suitable for XGBoost. Skipping.\n")
                result_prob <- NA
            } else {
                
                # --- Factor Level Synchronization and Constant Factor Removal ---
                predictor_vars <- unlist(strsplit(current_predictor_string, " \\+ "))
                train_for_mm <- droplevels(train)
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
                    next
                }
                
                # --- Data Preparation ---
                xgb_formula <- as.formula(paste("~", current_xgb_predictor_string))
                
                tryCatch({
                    # 1. Create Model Frame/Matrix with proper NA handling
                    model_frame_train <- model.frame(xgb_formula, data = train_for_mm, na.action = na.omit) 
                    x_train_xgb <- model.matrix(xgb_formula, data = model_frame_train)
                    x_train_xgb <- x_train_xgb[, colnames(x_train_xgb) != "(Intercept)", drop = FALSE]

                    # 2. Get the response variable from the same model frame to ensure perfect alignment
                    y_train_xgb <- train_for_mm$two_year_recid[as.numeric(row.names(model_frame_train))]
                    
                    # Debug information
                    cat(paste0("        DEBUG XGBoost: Model matrix rows: ", nrow(x_train_xgb), ", Response length: ", length(y_train_xgb), "\n"))
                    cat(paste0("        DEBUG XGBoost: NA count in response: ", sum(is.na(y_train_xgb)), "\n"))
                    
                    # Additional safety checks
                    if (length(y_train_xgb) != nrow(x_train_xgb)) {
                        stop("Mismatch between model matrix rows and response variable length")
                    }
                    if (any(is.na(y_train_xgb))) {
                         stop("Response variable 'two_year_recid' contains NA/NaN after model matrix creation, preventing DMatrix creation.")
                    }

                    train_weights <- if (universe$gender_imbalancing == "weight") {
                        train_for_mm$weights[as.numeric(row.names(model_frame_train))] 
                    } else { NULL }


                    # 3. Create DMatrix
                    dtrain <- xgb.DMatrix(
                        data = x_train_xgb, 
                        label = y_train_xgb, 
                        weight = train_weights
                    )
                    
                    # 4. Train the model
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

                    # 5. Prediction
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
                        cat("        FATAL ERROR: Factor became constant during model matrix creation (secondary check failed). Skipping.\n")
                        result_prob <- NA
                    } else if (grepl("Label contains NaN|Response variable 'two_year_recid'", conditionMessage(e))) {
                        cat("        FATAL ERROR: Label contains NaN/NA. Skipping. (Filtered earlier).\n")
                        result_prob <- NA
                    } else {
                        cat("        FATAL ERROR in data prep:", conditionMessage(e), "\n")
                        result_prob <- NA
                    }
                })
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
                
                # Update gender_factor in test set if it was removed due to "only" filtering
                if (universe$gender_imbalancing %in% c("male only", "female only")) {
                    test_processed$gender_factor <- droplevels(test_processed$gender_factor, exclude = if (universe$gender_imbalancing == "male only") {"Female"} else {"Male"})
                }
                
                # Factor levels should already be aligned from earlier processing
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
                if (universe$model == "survival") {
                # For survival models, we need to predict survival probability
                if (!is.null(cox_model)) {
                    tryCatch({
                        # Predict survival probability at 730 days (2 years)
                        surv_pred_test <- survfit(cox_model, newdata = test_processed)
                        surv_summary_test <- summary(surv_pred_test, times = 730)
                        
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
                        true_labels <- test_processed$event
                        
                        # Ensure all vectors have the same length
                        min_length <- min(length(binary_predictions), length(true_labels))
                        binary_predictions <- binary_predictions[1:min_length]
                        true_labels <- true_labels[1:min_length]
                        test_predictions <- test_predictions[1:min_length]
                        
                        # Calculate accuracy
                        accuracy_result <- mean(binary_predictions == true_labels, na.rm = TRUE)
                        
                        # Calculate AUC
                        if (requireNamespace("pROC", quietly = TRUE)) {
                            auc_result <- pROC::auc(true_labels, test_predictions, quiet = TRUE)
                        } else {
                            # Simple AUC calculation
                            auc_result <- tryCatch({
                                roc_obj <- ROCR::prediction(test_predictions, true_labels)
                                ROCR::performance(roc_obj, "auc")@y.values[[1]]
                            }, error = function(e) NA)
                        }
                        
                    }, error = function(e) {
                        cat("        Error in survival model evaluation:", conditionMessage(e), "\n")
                    })
                }
                
            } else if (universe$model == "logistic") {
                # For logistic models
                if (!is.null(logistic_model)) {
                    tryCatch({
                        test_predictions <- predict(logistic_model, newdata = test_processed, type = "response")
                        binary_predictions <- ifelse(test_predictions > 0.5, 1, 0)
                        true_labels <- test_processed$two_year_recid
                        
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
                
            } else if (universe$model == "xg_boost") {
                # For XGBoost models - only evaluate if we have the training matrix
                if (!is.null(xgboost_model) && exists("x_train_xgb")) {
                    tryCatch({
                        # Prepare test data the same way as training data
                        test_for_mm <- droplevels(test_processed)
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
                            xgb_formula <- as.formula(paste("~", current_xgb_predictor_string))
                            
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
                            test_indices <- as.numeric(row.names(model_frame_test))
                            true_labels <- test_for_mm$two_year_recid[test_indices]
                            
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
    
    save(universes_part1, file = "compas_results_part1.RData")
    return(universes_part1)
}
# --------------------------------------------------------------------------------------------------

# Example usage (assuming 'compas_2yr_recid' data object exists)
# Define your vectors of procedural choices
collect_data_parameters = c("compas")
time_partition_parameters = c("compas", "barenstein")
age_category_parameters = c("compas", "nij", "year")
race_category_parameters = c("nc", "nij")
split_parameters = c("6:4", "7:3", "8:2")
gender_imbalancing_parameters = c("under", "over", "male only", "female only", "weight")
model_parameters = c("survival", "logistic", "xg_boost") # <--- ADDED XGBOOST
predictor_parameters = c("compas","fair", "lasso")
define_recid_parameters = c("2yr")

# Call the function (Ensure 'compas_2yr_recid' is present)
multiverse_compas <- generate_compas_garden(
     compas_2yr_recid,
     collect_data_parameters,
     time_partition_parameters, 
     age_category_parameters, 
     race_category_parameters, 
     split_parameters, 
     gender_imbalancing_parameters, 
     model_parameters, 
     predictor_parameters
)

# Print summary of results
cat("Analysis completed!\n")
cat("Number of universes processed:", nrow(multiverse_compas), "\n")
cat("Average accuracy:", round(mean(multiverse_compas$accuracy, na.rm = TRUE), 4), "\n")
cat("Average AUC:", round(mean(multiverse_compas$auc, na.rm = TRUE), 4), "\n")
cat("Results saved to: compas_results_part1.RData\n")