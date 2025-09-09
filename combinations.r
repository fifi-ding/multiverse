#high risk
run_multiverse_analysis <- function(data, preprocessing_methods, split_methods, age_categories, imbalancing_methods, predictor_methods, define_recid_methods, profile_age, profile_gender, profile_race) {
  
    calculate_recid_surv <- function(surv_object, time_point) {
        # Find the index of the time point closest to the desired time
        time_index <- which.min(abs(surv_object$time - time_point))
        
        # Extract the survival probability at that index
        S_t <- surv_object$surv[time_index]
        
        # Calculate the probability of recidivism
        P_recid <- 1 - S_t
        
        # return survival probability and recidivism probability
        return(list(survival_prob = S_t, recidivism_prob = P_recid))
    }

    preprocess_A <- function(data){
        # recidivate during follow up period
        data <- data |>
        mutate(
            TIME = ifelse(RECID == 0 & TIME == 0, FOLLOW, TIME),
            TIME = ifelse(TIME == 0, 0.5, TIME),
            RECID = as.numeric(RECID)
        )

        # rescale
        data <- data |>
            mutate(
            TSERVD_100 = TSERVD / 100,
            AGE_1000 = AGE / 1000, 
            PRIORS_10 = PRIORS / 10,
            SCHOOL_10 = SCHOOL / 10,
            RULE_100 = RULE / 100,
            # convert month to year
            AGE_YEAR = round(AGE / 12)
        )
    }

    preprocess_B <- function(data){
        # recidivate during follow up period
        data <- data |>
        mutate(
            TIME = ifelse(RECID == 0, FOLLOW, TIME),
            TIME = ifelse(TIME == 0, 0.5, TIME),
            RECID = as.numeric(RECID)
        )

        # rescale
        data <- data |>
            mutate(
            TSERVD_100 = TSERVD / 100,
            AGE_1000 = AGE / 1000, 
            PRIORS_10 = PRIORS / 10,
            SCHOOL_10 = SCHOOL / 10,
            RULE_100 = RULE / 100,
            # convert month to year
            AGE_YEAR = round(AGE / 12)
        )
    }

    # no rescaling
    preprocess_C <- function(data){
        # recidivate during follow up period
        data <- data |>
        mutate(
            TIME = ifelse(RECID == 0, FOLLOW, TIME),
            TIME = ifelse(TIME == 0, 0.5, TIME),
            RECID = as.numeric(RECID)
        )
        
        # rescale
        data <- data |>
            mutate(
            TSERVD_100 = TSERVD / 100,
            AGE_1000 = AGE / 1000, 
            PRIORS_10 = PRIORS / 10,
            SCHOOL_10 = SCHOOL / 10,
            RULE_100 = RULE / 100,
            # convert month to year
            AGE_YEAR = round(AGE / 12)
        )
    }

  # Each row in the resulting data frame will be a unique "universe."
  all_universes <- expand.grid(
    preprocessing = preprocessing_methods,
    split = split_methods,
    age_category = age_categories,
    imbalancing_method = imbalancing_methods,
    predictor_method = predictor_methods,
    define_recid_method = define_recid_methods,
    stringsAsFactors = TRUE
  )
  
  # initialize a new column to store the results.
  all_universes$recidivism_prob <- NA
  
  num_universes <- nrow(all_universes)
  
  # Loop through each universe and perform the analysis
  for (i in 1:num_universes) {
    
    # Extract the current universe's procedural choices
    universe <- all_universes[i, ]
    
    # Print a header for the current universe to track progress
    cat("\n--- Running Universe", i, "of", num_universes, "---\n")
    cat("  Preprocessing Method:", universe$preprocessing, "\n")
    cat("  Split Method:", universe$split, "\n")
    cat("  Age Category:", universe$age_category, "\n")
    cat("  Imbalancing Method:", universe$imbalancing_method, "\n")
    cat("  Predictor Method:", universe$predictor_method, "\n")
    cat("  Define Recidivism:", universe$define_recid_method, "\n")
    
    # Progress update every 50 universes - write to file for Python to read
    if (i %% 50 == 0 || i == num_universes) {
      progress_file <- paste0("progress_", profile_age, "_", profile_gender, "_", profile_race, ".txt")
      write(paste("Universe", i, "of", num_universes, "completed"), file = progress_file, append = TRUE)
      cat("PROGRESS_UPDATE: Universe", i, "of", num_universes, "completed\n")
    }
    
    # ------------------------------------------------------------------------
    
   current_data <- data
    
    # Preprocessing method
    if (universe$preprocessing == "Method A") {
      cat("    Applying Preprocess A...\n")
      current_data <- preprocess_A(current_data)
    } else if (universe$preprocessing == "Method B") {
      cat("    Applying Preprocess B...\n")
      current_data <- preprocess_B(current_data)
    }else{
      cat("    Applying Preprocess C...\n")
      current_data <- preprocess_C(current_data)
    }
    
    # Age Category
    if (universe$age_category == "raw_age_year"){
      current_data <- current_data |>
        mutate(
          AGE_YEAR = round(AGE / 12)
      )
    } else if (universe$age_category == "age_cat_compas"){
      current_data <- current_data |>
        mutate(
          AGE_CAT_COMPAS = case_when(
          AGE_YEAR < 25 ~ "less than 25",
          AGE_YEAR >= 25 & AGE_YEAR <= 45 ~ "25~45",
          AGE_YEAR > 45 ~ "over 45"
          )
        )
    } else if (universe$age_category == "age_cat_nij"){
      current_data <- current_data |>
        mutate(
          AGE_CAT_NIJ = case_when(
            # age under 17 isn't considered, there are missing cells
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
    
    # Split Method
    if(universe$split == "1:2"){
      # analysis file
      analysis_1978 <- current_data |> filter(FILE == 1)
      # validation file
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
    
    # Imbalancing method
    if (universe$imbalancing_method == "Oversampling") {
      cat("    Applying Oversampling method...\n")
      male_sample <- analysis_1978 |> filter(MALE == 1)
      female_sample <- analysis_1978 |> filter(MALE == 0)
      set.seed(123)
      females_oversampled <- female_sample |> sample_n(nrow(male_sample), replace = TRUE)
      analysis_1978 <- bind_rows(male_sample, females_oversampled)
      # cat("    Balance check after oversampling:\n")
      # print(analysis_1978 |> count(MALE))
    } else if (universe$imbalancing_method == "Undersampling"){
      cat("    Applying Undersampling method...\n")
      male_sample <- analysis_1978 |> filter(MALE == 1)
      female_sample <- analysis_1978 |> filter(MALE == 0)
      # Downsample males to match number of females
      set.seed(123)
      males_undersampled <- male_sample |> sample_n(nrow(female_sample), replace = FALSE)
      # Combine balanced dataset
      analysis_1978 <- bind_rows(males_undersampled, female_sample)
      # cat("    Balance check after undersampling:\n")
      # print(current_data |> count(MALE))
    } else if (universe$imbalancing_method == "Male Only"){
      cat("    Applying Male Only...\n")
      analysis_1978 <- analysis_1978 |> filter(MALE == 1)
    } else if (universe$imbalancing_method == "Female Only"){
      cat("    Applying Female Only...\n")
      analysis_1978 <- analysis_1978 |> filter(MALE == 0)
    } else {
      cat("    Applying Weighting...\n")
      analysis_1978$weights <- ifelse(analysis_1978$MALE == 0, 4348 / 270, 1)
    }
    
    # Predictor
    if (universe$predictor_method == "full") {
      cox_formula <- Surv(TIME, RECID) ~ TSERVD_100 + AGE_1000 + PRIORS_10 + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE + RULE_100 + MARRIED + SCHOOL_10 + WORKREL + PERSON + SUPER
    } else if (universe$predictor_method == "schmidtl"){
      cox_formula <- Surv(TIME, RECID) ~ TSERVD_100 + AGE_1000 + PRIORS_10 + WHITE + FELON + ALCHY + JUNKY + PROPTY + MALE
    } else if (universe$predictor_method == "protected") {
      cox_formula <- Surv(TIME, RECID) ~ TSERVD_100 + PRIORS_10 + FELON + ALCHY + JUNKY + PROPTY
    }
    
    # Create profile data frame for prediction
    profile_data <- data.frame(
      TSERVD_100 = 0,  # Default values - you may need to adjust these
      AGE_1000 = profile_age / 1000,
      PRIORS_10 = 0,
      WHITE = ifelse(profile_race == "caucasian", 1, 0),
      FELON = 0,
      ALCHY = 0,
      JUNKY = 0,
      PROPTY = 0,
      MALE = ifelse(profile_gender == "male", 1, 0),
      RULE_100 = 0,
      MARRIED = 0,
      SCHOOL_10 = 0,
      WORKREL = 0,
      PERSON = 0,
      SUPER = 0
    )
    
    # Model
    cox_model <- coxph(cox_formula, data = analysis_1978)
    surv_pred <- survfit(cox_model, newdata = profile_data)
    
    # Define Recidivism
    if (universe$define_recid_method == "1yr") {
      # within 1 year
      result <- calculate_recid_surv(surv_pred, 12)
    } else if (universe$define_recid_method == "2yr"){
      # within 2 years
      result <- calculate_recid_surv(surv_pred, 24)
    } else if (universe$define_recid_method == "3yr") {
      # within 3 years
      result <- calculate_recid_surv(surv_pred, 36)
    } else if (universe$define_recid_method == "4yr") {
      # within 4 years
      result <- calculate_recid_surv(surv_pred, 48)
    } else if (universe$define_recid_method == "5yr") {
      # within 5 years
      result <- calculate_recid_surv(surv_pred, 60)
    }
    
    all_universes$recidivism_prob[i] <- result$recidivism_prob
    cat("  Predicted Recidivism Probability:", round(all_universes$recidivism_prob[i], 3), "\n")
  }
  
  return(all_universes)
}

# Define your vectors of procedural choices
preprocessing_options <- c("Method_A", "Method_B", "Method_C")
split_options <- c("1:2", "6:4", "7:3", "8:2")
age_cat_options <- c("raw_age_year","age_cat_compas", "age_cat_nij")
imbalancing_options <- c("Undersampling", "Oversampling", "Male Only", "Female Only", "Weighting")
predictor_options <- c("full", "schmidt", "protected")
define_recid_options <- c("1yr", "2yr", "3yr", "4yr", "5yr")

# Function will be called from Python with profile parameters