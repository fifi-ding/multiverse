import dash
from dash import dcc, html, Input, Output, State
import plotly.graph_objects as go
import dash_bootstrap_components as dbc
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Initialize empty dataframes that will be populated by the R analysis
df_nc = pd.DataFrame()  # Will hold focal profile results
df_low_risk = pd.DataFrame()  # Will hold counterfactual profile results

# Import R packages for decision tree analysis
base_pkg = importr('base')
stats = importr('stats')
rpart_pkg = importr('rpart')
rpart_plot_pkg = importr('rpart.plot')

# R code for decision tree analysis
r_tree_code = """
# Helper function to safely reset row names
reset_rownames <- function(data) {
    if (is.data.frame(data)) {
        row.names(data) <- NULL
    }
    return(data)
}

# Function to create decision tree and get variable importance
create_decision_tree <- function(data) {
    # Reset row names to avoid conversion issues
    data <- reset_rownames(data)
    
    # Create the decision tree
    tree <- rpart(recidivism_prob ~ scaling + convert_followup + split + age_category + 
                  imbalancing + model + predictor + define_recid,
                  data = data,
                  method = "anova",
                  control = rpart.control(minsplit = 2, cp = 0))
    
    # Get variable importance
    var_importance <- tree$variable.importance
    
    # Handle case where variable.importance is NULL (no splits in tree)
    if (is.null(var_importance)) {
        # Create a default variable importance with all variables set to 0
        all_vars <- c("scaling", "convert_followup", "split", "age_category", 
                     "imbalancing", "model", "predictor", "define_recid")
        var_importance <- rep(0, length(all_vars))
        names(var_importance) <- all_vars
    }
    
    # Ensure variable importance has proper names
    if (is.null(names(var_importance)) || length(names(var_importance)) == 0) {
        # Get variable names from the tree frame
        tree_frame <- tree$frame
        var_names <- as.character(tree_frame$var)
        used_vars <- unique(var_names[var_names != "<leaf>"])
        
        # Assign names to variable importance
        if (length(used_vars) >= length(var_importance)) {
            names(var_importance) <- used_vars[1:length(var_importance)]
        } else {
            names(var_importance) <- paste0("Var_", 1:length(var_importance))
        }
    }
    
    # Return both tree and variable importance
    return(list(tree = tree, var_importance = var_importance))
}

# Function to get variable importance for a dataset
get_variable_importance <- function(data) {
    tryCatch({
        result <- create_decision_tree(data)
        var_importance <- result$var_importance
        
        # Debug: print what we get from variable.importance
        cat("Variable importance result:\n")
        print(var_importance)
        cat("Names:", names(var_importance), "\n")
        
        # Sort by importance values (descending) to ensure proper ordering
        var_importance_sorted <- sort(var_importance, decreasing = TRUE)
        
        cat("Sorted variable importance:\n")
        print(var_importance_sorted)
        
        # Return sorted result
        return(var_importance_sorted)
    }, error = function(e) {
        cat("Error in get_variable_importance:", e$message, "\n")
        # Return default variable importance with all variables set to 0
        all_vars <- c("scaling", "convert_followup", "split", "age_category", 
                     "imbalancing", "model", "predictor", "define_recid")
        default_importance <- rep(0, length(all_vars))
        names(default_importance) <- all_vars
        return(default_importance)
    })
}

#' Get and rank the split rules from an rpart tree by their position in the tree.
#'
#' @param tree An rpart object.
#' @return A data.frame with details for each split, ranked by top-down position.
get_split_rules_by_position <- function(tree) {
  # Check if the input is an rpart object
  if (!inherits(tree, "rpart")) {
    stop("The provided object is not an rpart tree.")
  }

  # The frame and splits contain all the necessary info
  frame <- tree$frame
  splits <- tree$splits

  # 1. Filter for internal nodes only (where splits occur)
  internal_nodes_frame <- frame[frame$var != "<leaf>", ]
  if (nrow(internal_nodes_frame) == 0 || is.null(splits)) {
    return(data.frame()) # Return empty frame if no splits
  }
  
  node_ids <- as.numeric(rownames(internal_nodes_frame))

  # 2. Get the paths to the direct children of these nodes
  left_children <- node_ids * 2
  right_children <- node_ids * 2 + 1
  
  # Safely get child paths
  tryCatch({
    child_paths <- path.rpart(tree, nodes = c(left_children, right_children), print.it = FALSE)
    
    left_rules <- sapply(child_paths[1:length(node_ids)], tail, 1)
    right_rules <- sapply(child_paths[(length(node_ids) + 1):length(child_paths)], tail, 1)
  }, error = function(e) {
    # If path.rpart fails, create simple rules
    left_rules <- rep("", length(node_ids))
    right_rules <- rep("", length(node_ids))
    child_paths <- list()
  })

  # 3. Combine everything into a data frame
  split_rules_df <- data.frame(
    node = node_ids,
    variable = internal_nodes_frame$var,
    n_obs_split = internal_nodes_frame$n,
    improvement = splits[1:nrow(internal_nodes_frame), "improve"],
    left_rule = left_rules,
    right_rule = right_rules
  )
  
  # --- THIS IS THE CORRECTED LINE ---
  # We now order by the node number to match the tree's visual structure.
  split_rules_df <- split_rules_df[order(split_rules_df$node), ]
  split_rules_df$rank <- 1:nrow(split_rules_df)
  
  # Reorder columns for clarity
  split_rules_df <- split_rules_df[, c("rank", "node", "variable", "improvement", "n_obs_split", "left_rule", "right_rule")]

  return(split_rules_df)
}

# Function to get regression tree nodes and split rules
get_tree_nodes <- function(data) {
    tryCatch({
        result <- create_decision_tree(data)
        tree <- result$tree
        
        # Use the new function to get split rules by position
        ranked_splits_by_position <- get_split_rules_by_position(tree)
        
        # Convert to the expected format for backward compatibility
        node_rules <- list()
        if (nrow(ranked_splits_by_position) > 0) {
        for (i in 1:nrow(ranked_splits_by_position)) {
            node_rules[[length(node_rules) + 1]] <- list(
                rank = ranked_splits_by_position$rank[i],
                node = ranked_splits_by_position$node[i],
                variable = as.character(ranked_splits_by_position$variable[i]),
                rule = as.character(ranked_splits_by_position$left_rule[i])
            )
        }
    }
    
    return(node_rules)
    }, error = function(e) {
        cat("Error in get_tree_nodes:", e$message, "\n")
        # Return empty list if analysis fails
        return(list())
    })
}
"""

# Execute the R code
robjects.r(r_tree_code)

# Load the combinations.r file
with open('combinations.r', 'r') as f:
    combinations_r_code = f.read()

# Load required R packages and data
r_setup_code = """
# Load required packages
library(survival)
library(dplyr)

# Load the analysis data
# You'll need to adjust this path to your actual data file
analysis_1978 <- read.csv('nc_prisoner_1978.csv')

# Fix row names to avoid conversion issues with rpy2
# Reset row names to simple integer sequence
rownames(analysis_1978) <- NULL

# Define helper functions that might be missing
calculate_recid_surv <- function(surv_pred, time_months) {
  # Extract survival probability at the specified time
  time_years <- time_months / 12
  surv_prob <- summary(surv_pred, times = time_years)$surv[1]
  recid_prob <- 1 - surv_prob
  return(list(recidivism_prob = recid_prob))
}

# Preprocessing functions are now defined in combinations.r
"""

# Execute the setup code
robjects.r(r_setup_code)

# Execute the combinations R code
robjects.r(combinations_r_code)

# Function to check if two profiles have identical characteristics
def profiles_are_identical(age1, gender1, race1, age2, gender2, race2):
    """Check if two profiles have identical demographic characteristics"""
    return age1 == age2 and gender1 == gender2 and race1 == race2

# Function to run multiverse analysis with profile parameters
def run_multiverse_analysis_with_profile(profile_age, profile_gender, profile_race, is_shared_analysis=False):
    """Run multiverse analysis for a specific demographic profile"""
    global progress_updates, progress_counts
    
    # Determine profile key based on parameters
    if profile_age == 25 and profile_gender == "male" and profile_race == "caucasian":
        profile_key = "focal"
        profile_name = "Focal"
    else:
        profile_key = "counterfactual"
        profile_name = "Counterfactual"
    
    print(f"=== run_multiverse_analysis_with_profile called ===")
    print(f"Parameters: age={profile_age}, gender={profile_gender}, race={profile_race}")
    print(f"Detected profile: {profile_key} ({profile_name})")
    print(f"Is shared analysis: {is_shared_analysis}")
    print(f"Current progress counts: focal={progress_counts['focal']}, counterfactual={progress_counts['counterfactual']}")
    
    try:
        # Only reset progress if not already at completion (meaning it's being reused)
        if progress_counts[profile_key] < total_universes:
            print(f"Resetting progress for {profile_key}")
            progress_updates[profile_key] = []
            progress_counts[profile_key] = 0
        else:
            print(f"NOT resetting progress for {profile_key} - already at {progress_counts[profile_key]}")
            # If it's a reused profile, just update the message
            if len(progress_updates[profile_key]) == 0:
                progress_updates[profile_key].append(f"{profile_name} profile results reused from previous analysis")
        
        # If this is a shared analysis, also reset the other profile's progress
        if is_shared_analysis:
            other_profile_key = "counterfactual" if profile_key == "focal" else "focal"
            progress_updates[other_profile_key] = []
            progress_counts[other_profile_key] = 0
        
        # Update progress
        progress_updates[profile_key].append(f"Starting {profile_name} analysis...")
        progress_updates[profile_key].append(f"Generating 9600 universes in 3 parts...")
        
        # If shared analysis, update both profiles
        if is_shared_analysis:
            other_profile_name = "Counterfactual" if profile_name == "Focal" else "Focal"
            progress_updates[other_profile_key].append(f"Starting {other_profile_name} analysis...")
            progress_updates[other_profile_key].append(f"Generating 9600 universes in 3 parts...")
        
        with localconverter(robjects.default_converter + pandas2ri.converter):
            # Get the R functions - use new NC garden approach for all three parts
            r_run_multiverse_part1 = robjects.globalenv['generate_nc_garden']
            r_run_multiverse_part2 = robjects.globalenv['generate_nc_garden_part2']
            r_run_multiverse_part3 = robjects.globalenv['generate_nc_garden_part3']
            
            # Define the procedural choice options for new NC garden approach
            scaling_options = robjects.StrVector(["schmidt_scale", "no_scaling"])
            convert_followup_options = robjects.StrVector(["strict_def", "general_def"])
            split_options = robjects.StrVector(["1:2", "6:4", "7:3", "8:2"])
            age_cat_options = robjects.StrVector(["year", "schmidt_age", "compas", "nij"])
            imbalancing_options = robjects.StrVector(["under", "over", "male only", "female only", "weight"])
            model_options = robjects.StrVector(["survival", "logistic"])
            predictor_options = robjects.StrVector(["full", "schmidt", "fair"])
            define_recid_options = robjects.StrVector(["1yr", "2yr", "3yr", "4yr", "5yr"])
            
            # Profile parameters will be passed directly to R functions
            # They expect: profile_age (numeric), profile_gender (string), profile_race (string)
            
            # Load the analysis data with proper row name handling
            try:
                # First, reset row names in R to avoid conversion issues
                robjects.r('row.names(analysis_1978) <- NULL')
                analysis_data = robjects.globalenv['analysis_1978']
            except Exception as e:
                print(f"Warning: Issue loading analysis data: {e}")
                # Try to reload the data from file
                robjects.r("analysis_1978 <- read.csv('nc_prisoner_1978.csv')")
                robjects.r('row.names(analysis_1978) <- NULL')
                analysis_data = robjects.globalenv['analysis_1978']
            
            print(f"Running {profile_name} analysis: Age={profile_age}, Gender={profile_gender}, Race={profile_race}")
            
            # Start progress monitoring in a separate thread
            def monitor_progress():
                import os
                import time
                progress_file = "progress_nc_garden.txt"
                
                # Clear any existing progress file
                if os.path.exists(progress_file):
                    os.remove(progress_file)
                
                # Monitor the progress file
                while True:
                    if os.path.exists(progress_file):
                        try:
                            with open(progress_file, 'r') as f:
                                lines = f.readlines()
                                if lines:
                                    # Get the last line to find current progress
                                    last_line = lines[-1].strip()
                                    if "Universe" in last_line and "of" in last_line:
                                        # Extract universe number
                                        parts = last_line.split()
                                        universe_num = int(parts[1])
                                        total_universes_from_file = int(parts[3])
                                        
                                        # Only update progress if not already at completion (for reused profiles)
                                        if progress_counts[profile_key] < total_universes or universe_num > progress_counts[profile_key]:
                                            progress_counts[profile_key] = universe_num
                                        
                                        # If shared analysis, update both profiles' progress
                                        if is_shared_analysis:
                                            other_profile_key = "counterfactual" if profile_key == "focal" else "focal"
                                            if progress_counts[other_profile_key] < total_universes or universe_num > progress_counts[other_profile_key]:
                                                progress_counts[other_profile_key] = universe_num
                                        
                                        # Update progress text every 100 universes
                                        if universe_num % 100 == 0:
                                            progress_updates[profile_key].append(f"Processing Universe {universe_num}/{total_universes}...")
                                            if is_shared_analysis:
                                                other_profile_name = "Counterfactual" if profile_name == "Focal" else "Focal"
                                                progress_updates[other_profile_key].append(f"Processing Universe {universe_num}/{total_universes}...")
                                        
                                        # Check if analysis is complete
                                        if universe_num >= total_universes:
                                            print(f"Progress monitoring: Detected completion at {universe_num}/{total_universes}")
                                            break
                        except Exception as e:
                            # Log errors but continue monitoring
                            print(f"Progress monitoring error: {e}")
                            pass
                    time.sleep(0.5)  # Check every 0.5 seconds
                
                print("Progress monitoring thread exiting")
            
            # Start progress monitoring
            import threading
            progress_thread = threading.Thread(target=monitor_progress)
            progress_thread.daemon = True
            progress_thread.start()
            
            # Run Part 1: universes 1-3000
            progress_updates[profile_key].append(f"Running Part 1 (universes 1-3000)...")
            if is_shared_analysis:
                progress_updates[other_profile_key].append(f"Running Part 1 (universes 1-3000)...")
            
            print(f"Calling R function Part 1 with profile: age={profile_age}, gender={profile_gender}, race={profile_race}")
            try:
                results_part1 = r_run_multiverse_part1(analysis_data, scaling_options, convert_followup_options,
                                                       age_cat_options, split_options, imbalancing_options, 
                                                       model_options, predictor_options, define_recid_options,
                                                       profile_age, profile_gender, profile_race)
                print("Part 1 R function completed, converting to pandas...")
                results_df_part1 = robjects.conversion.rpy2py(results_part1)
                print(f"Part 1 completed: {len(results_df_part1)} universes")
            except Exception as e:
                print(f"ERROR in Part 1: {e}")
                import traceback
                traceback.print_exc()
                raise
            progress_counts[profile_key] = 3000
            if is_shared_analysis:
                progress_counts[other_profile_key] = 3000
            
            # Run Part 2: universes 3001-6000
            progress_updates[profile_key].append(f"Running Part 2 (universes 3001-6000)...")
            if is_shared_analysis:
                progress_updates[other_profile_key].append(f"Running Part 2 (universes 3001-6000)...")
            
            print(f"Calling R function Part 2 with profile: age={profile_age}, gender={profile_gender}, race={profile_race}")
            try:
                results_part2 = r_run_multiverse_part2(analysis_data, scaling_options, convert_followup_options,
                                                       age_cat_options, split_options, imbalancing_options, 
                                                       model_options, predictor_options, define_recid_options,
                                                       profile_age, profile_gender, profile_race)
                print("Part 2 R function completed, converting to pandas...")
                results_df_part2 = robjects.conversion.rpy2py(results_part2)
                print(f"Part 2 completed: {len(results_df_part2)} universes")
            except Exception as e:
                print(f"ERROR in Part 2: {e}")
                import traceback
                traceback.print_exc()
                raise
            progress_counts[profile_key] = 6000
            if is_shared_analysis:
                progress_counts[other_profile_key] = 6000
            
            # Run Part 3: universes 6001-9600
            progress_updates[profile_key].append(f"Running Part 3 (universes 6001-9600)...")
            if is_shared_analysis:
                progress_updates[other_profile_key].append(f"Running Part 3 (universes 6001-9600)...")
            
            print(f"Calling R function Part 3 with profile: age={profile_age}, gender={profile_gender}, race={profile_race}")
            try:
                results_part3 = r_run_multiverse_part3(analysis_data, scaling_options, convert_followup_options,
                                                       age_cat_options, split_options, imbalancing_options, 
                                                       model_options, predictor_options, define_recid_options,
                                                       profile_age, profile_gender, profile_race)
                print("Part 3 R function completed, converting to pandas...")
                results_df_part3 = robjects.conversion.rpy2py(results_part3)
                print(f"Part 3 completed: {len(results_df_part3)} universes")
            except Exception as e:
                print(f"ERROR in Part 3: {e}")
                import traceback
                traceback.print_exc()
                raise
            
            # Combine all three parts into one DataFrame
            print("Combining all three parts into one DataFrame...")
            results_df = pd.concat([results_df_part1, results_df_part2, results_df_part3], ignore_index=True)
            print(f"Combined DataFrame: {len(results_df)} total universes")
            print(f"DataFrame columns: {list(results_df.columns)}")
            print(f"Scaling column unique values: {results_df['scaling'].unique() if 'scaling' in results_df.columns else 'SCALING COLUMN NOT FOUND'}")
            print(f"First few rows:\n{results_df.head()}")
            
            # Convert factor columns to strings if they came as numeric codes from R
            categorical_columns = ['scaling', 'convert_followup', 'split', 'age_category', 
                                  'imbalancing', 'model', 'predictor', 'define_recid']
            for col in categorical_columns:
                if col in results_df.columns:
                    # If column is numeric (factor codes), we have a problem
                    if results_df[col].dtype in ['int64', 'int32', 'float64', 'float32']:
                        print(f"WARNING: {col} column is numeric type {results_df[col].dtype}, R factors may not have converted properly")
                    # Convert to string type to ensure consistency
                    results_df[col] = results_df[col].astype(str)
            
            print(f"After type conversion - Scaling unique values: {results_df['scaling'].unique()}")
            
            # Finalize progress
            progress_counts[profile_key] = 9600
            progress_updates[profile_key].append(f"Completed: Generated all {len(results_df)} universes")
            progress_updates[profile_key].append(f"{profile_name} analysis completed!")
            
            # If shared analysis, finalize both profiles
            if is_shared_analysis:
                other_profile_key = "counterfactual" if profile_key == "focal" else "focal"
                other_profile_name = "Counterfactual" if profile_name == "Focal" else "Focal"
                progress_counts[other_profile_key] = 9600
                progress_updates[other_profile_key].append(f"Completed: Generated all {len(results_df)} universes")
                progress_updates[other_profile_key].append(f"{other_profile_name} analysis completed!")
            
            # Clean up progress file
            import os
            progress_file = "progress_nc_garden.txt"
            if os.path.exists(progress_file):
                os.remove(progress_file)
        
        print(f"{profile_name} analysis completed. Generated {len(results_df)} universes.")
        return results_df
        
    except Exception as e:
        # Update progress with error
        progress_updates[profile_key].append(f"Error in {profile_name} analysis")
        print(f"Error running {profile_name} multiverse analysis: {e}")
        import traceback
        traceback.print_exc()
        return None

# Removed parallel function - using sequential execution with proper progress tracking

# Data loaded and ready for dashboard

# Constants for dropdown options - using the new NC garden parameter values
SCALING_OPTIONS = [{'label': method, 'value': method} for method in ["schmidt_scale", "no_scaling"]]
CONVERT_FOLLOWUP_OPTIONS = [{'label': method, 'value': method} for method in ["strict_def", "general_def"]]
SPLIT_OPTIONS = [{'label': split, 'value': split} for split in ["1:2", "6:4", "7:3", "8:2"]]
AGE_CATEGORY_OPTIONS = [{'label': age_cat, 'value': age_cat} for age_cat in ["year", "schmidt_age", "compas", "nij"]]
IMBALANCING_OPTIONS = [{'label': method, 'value': method} for method in ["under", "over", "male only", "female only", "weight"]]
MODEL_OPTIONS = [{'label': method, 'value': method} for method in ["survival", "logistic"]]
PREDICTOR_OPTIONS = [{'label': method, 'value': method} for method in ["full", "schmidt", "fair"]]
RECID_METHOD_OPTIONS = [{'label': method, 'value': method} for method in ["1yr", "2yr", "3yr", "4yr", "5yr"]]

# Initialize the Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True)
app.title = "Multiverse Dashboard"


# Add state to track which profile's grid is currently displayed
current_grid_profile = 1  # 1 for Profile 1, 2 for Profile 2

def create_layout():
    """Create the main layout of the dashboard"""
    return html.Div([
        # Hidden input to trigger updates
        dcc.Input(id='dummy-input', value='', style={'display': 'none'}),
        
        # Hidden divs to store selected universes
        html.Div(id='selected-universes-1', style={'display': 'none'}),
        html.Div(id='selected-universes-2', style={'display': 'none'}),
        
        # Hidden divs to store highlighted universe for each profile
        html.Div(id='highlighted-universe-1', style={'display': 'none'}),
        html.Div(id='highlighted-universe-2', style={'display': 'none'}),
        
        # Header
        html.H1("Multiverse Dashboard", style={'textAlign': 'center', 'marginBottom': '20px'}),
        
        # Main horizontal flex container
        html.Div([
            # Left Sidebar - Counterfactual Profiles
            html.Div([
                # Step 1 Instruction Block (includes Profiles title and profile cards)
                html.Div([
                    html.Div("STEP 1", style={
                        'fontSize': '12px', 
                        'fontWeight': 'bold', 
                        'color': 'white', 
                        'backgroundColor': '#6c757d', 
                        'padding': '4px 8px', 
                        'borderRadius': '4px',
                        'textAlign': 'center',
                        'marginBottom': '12px'
                    }),
                    html.P("Customize your counterfactual profile by adjusting characteristics.", 
                           style={'fontSize': '11px', 'color': '#333', 'margin': '0 0 15px 0', 'lineHeight': '1.4', 'fontStyle': 'italic'}),
                
                # Profile containers side by side
                html.Div([
                    # Profile 1 - Focal
                    html.Div([
                        html.Div("Focal (Locked)", style={'textAlign': 'center', 'marginBottom': '8px', 'fontSize': '12px', 'color': 'white', 'fontWeight': 'bold', 'backgroundColor': '#d63384', 'padding': '4px 8px', 'borderRadius': '12px'}),
                        
                        html.Div([
                            html.Label("Age:", style={'fontSize': '14px', 'marginBottom': '4px'}),
                            dcc.Slider(
                                id='age-slider-1',
                                min=0,
                                max=100,
                                step=1,
                                value=25,
                                marks={0: '0', 25: '25', 50: '50', 75: '75', 100: '100'},
                                tooltip={'placement': 'bottom', 'always_visible': True},
                                disabled=True  # Lock the focal profile age slider
                            ),
                        ], style={'marginBottom': '12px'}),
                        
                        html.Div([
                            html.Label("Gender:", style={'fontSize': '14px', 'marginBottom': '4px'}),
                            dcc.Dropdown(
                                id='gender-dropdown-1',
                                options=[
                                    {'label': 'Female', 'value': 'female'},
                                    {'label': 'Male', 'value': 'male'}
                                ],
                                value='male',
                                clearable=False,
                                disabled=True  # Lock the focal profile gender dropdown
                            ),
                        ], style={'marginBottom': '12px'}),
                        
                        html.Div([
                            html.Label("Race:", style={'fontSize': '14px', 'marginBottom': '4px'}),
                            dcc.Dropdown(
                                id='race-dropdown-1',
                                options=[
                                    {'label': 'Asian', 'value': 'asian'},
                                    {'label': 'Caucasian', 'value': 'caucasian'},
                                    {'label': 'African American', 'value': 'african_american'},
                                    {'label': 'Native American', 'value': 'native_american'},
                                    {'label': 'Hispanic', 'value': 'hispanic'},
                                    {'label': 'Other', 'value': 'other'}
                                ],
                                value='caucasian',
                                clearable=False,
                                disabled=True  # Lock the focal profile race dropdown
                            ),
                        ], style={'marginBottom': '12px'}),
                        
                        # Add a note that the focal profile is locked
                        html.P("This profile is locked for comparison", 
                               style={'fontSize': '10px', 'color': '#666', 'textAlign': 'center', 'marginTop': '8px', 'fontStyle': 'italic'})
                        
                    ], style={'flex': '1', 'padding': '12px', 'border': '2px solid #d63384', 'borderRadius': '5px', 'marginRight': '8px', 'backgroundColor': '#fff8fa'}),
                    
                    # Profile 2 - Counterfactual
                    html.Div([
                        html.Div("Counterfactual", style={'textAlign': 'center', 'marginBottom': '8px', 'fontSize': '12px', 'color': 'white', 'fontWeight': 'bold', 'backgroundColor': '#0d6efd', 'padding': '4px 8px', 'borderRadius': '12px'}),
                        
                        html.Div([
                            html.Label("Age:", style={'fontSize': '14px', 'marginBottom': '4px'}),
                            dcc.Slider(
                                id='age-slider-2',
                                min=0,
                                max=100,
                                step=1,
                                value=25,
                                marks={0: '0', 25: '25', 50: '50', 75: '75', 100: '100'},
                                tooltip={'placement': 'bottom', 'always_visible': True}
                            ),
                        ], style={'marginBottom': '12px'}),
                        
                        html.Div([
                            html.Label("Gender:", style={'fontSize': '14px', 'marginBottom': '4px'}),
                            dcc.Dropdown(
                                id='gender-dropdown-2',
                                options=[
                                    {'label': 'Female', 'value': 'female'},
                                    {'label': 'Male', 'value': 'male'}
                                ],
                                value='male',
                                clearable=False
                            ),
                        ], style={'marginBottom': '12px'}),
                        
                        html.Div([
                            html.Label("Race:", style={'fontSize': '14px', 'marginBottom': '4px'}),
                            dcc.Dropdown(
                                id='race-dropdown-2',
                                options=[
                                    {'label': 'Asian', 'value': 'asian'},
                                    {'label': 'Caucasian', 'value': 'caucasian'},
                                    {'label': 'African American', 'value': 'african_american'},
                                    {'label': 'Native American', 'value': 'native_american'},
                                    {'label': 'Hispanic', 'value': 'hispanic'},
                                    {'label': 'Other', 'value': 'other'}
                                ],
                                value='caucasian',
                                clearable=False
                            ),
                        ], style={'marginBottom': '12px'}),
                        
                    ], style={'flex': '1', 'padding': '12px', 'border': '2px solid #0d6efd', 'borderRadius': '5px', 'marginLeft': '8px', 'backgroundColor': '#f8fbff'}),
                    
                                ], style={
                    'display': 'flex',
                    'flexDirection': 'row',
                    'gap': '15px',
                    'marginBottom': '15px'
                }),
                
                    # Run Multiverse Analysis Button
                    html.Button(
                        "Run Multiverse Analysis",
                        id='submit-profiles-button',
                        style={
                            'fontSize': '16px',
                            'padding': '12px 24px',
                            'backgroundColor': '#6c757d',
                            'color': 'white',
                            'border': 'none',
                            'borderRadius': '8px',
                            'cursor': 'pointer',
                            'fontWeight': 'bold',
                            'width': '100%',
                            'marginBottom': '10px'
                        }
                    ),
                    html.P("Click to analyze the selected profiles and generate multiverse results", 
                           style={'fontSize': '11px', 'color': '#666', 'textAlign': 'center', 'marginBottom': '0px', 'fontStyle': 'italic'})
                
                ], style={
                    'backgroundColor': '#f8f9fa', 
                    'border': '2px solid #6c757d', 
                    'borderRadius': '8px', 
                    'padding': '15px', 
                    'marginBottom': '15px'
                }),
                
                # Step 2 Instruction Block (includes progress and overview)
                html.Div([
                    html.Div("STEP 2", style={
                        'fontSize': '12px', 
                        'fontWeight': 'bold', 
                        'color': 'white', 
                        'backgroundColor': '#6c757d', 
                        'padding': '4px 8px', 
                        'borderRadius': '4px',
                        'textAlign': 'center',
                        'marginBottom': '12px'
                    }),
                    html.P("Monitor the analysis progress and view method combinations overview. The analysis will generate all possible universe combinations for both profiles.", 
                           style={'fontSize': '11px', 'color': '#333', 'margin': '0 0 15px 0', 'lineHeight': '1.4', 'fontStyle': 'italic'}),
                    
                    # Analysis Progress Section
                    html.H5("Analysis Progress", style={'textAlign': 'center', 'marginBottom': '10px', 'fontSize': '14px', 'color': '#333'}),
                    
                    # Progress bars for each profile
                    html.Div([
                        html.Div([
                            html.P("Focal Profile", style={'fontSize': '11px', 'color': '#d63384', 'marginBottom': '5px', 'fontWeight': 'bold'}),
                            dcc.Interval(
                                id='progress-interval',
                                interval=2000,  # Update every 2 seconds
                                n_intervals=0
                            ),
                            html.Div(
                                id='focal-progress-bar',
                                style={
                                    'width': '0%',
                                    'height': '8px',
                                    'backgroundColor': '#d63384',
                                    'borderRadius': '4px',
                                    'transition': 'width 0.3s ease'
                                }
                            ),
                            html.P(id='focal-progress-text', children="Universe 0/9600", 
                                   style={'fontSize': '10px', 'color': '#666', 'marginTop': '2px'})
                        ], style={'marginBottom': '10px'}),
                        
                        html.Div([
                            html.P("Counterfactual Profile", style={'fontSize': '11px', 'color': '#0d6efd', 'marginBottom': '5px', 'fontWeight': 'bold'}),
                            html.Div(
                                id='counterfactual-progress-bar',
                                style={
                                    'width': '0%',
                                    'height': '8px',
                                    'backgroundColor': '#0d6efd',
                                    'borderRadius': '4px',
                                    'transition': 'width 0.3s ease'
                                }
                            ),
                            html.P(id='counterfactual-progress-text', children="Universe 0/9600", 
                                   style={'fontSize': '10px', 'color': '#666', 'marginTop': '2px'})
                        ], style={'marginBottom': '10px'}),
                        
                        # Add a note about optimizations
                        html.Div([
                            html.P("Note: Identical profiles run once together. Focal profile results are reused when possible.", 
                                   style={'fontSize': '9px', 'color': '#888', 'textAlign': 'center', 'fontStyle': 'italic', 'marginTop': '5px'})
                        ], style={'marginBottom': '10px'})
                    ], style={'marginBottom': '15px'}),
                    
                    # Method Combinations Overview Section
                    html.H5("Method Combinations Overview", style={'textAlign': 'center', 'marginBottom': '10px', 'fontSize': '14px', 'color': '#333'}),
                    html.Div(
                        id='dataset-overview-content',
                        style={
                            'padding': '8px',
                            'border': '1px solid #ddd',
                            'borderRadius': '5px',
                            'backgroundColor': '#f9f9f9',
                            'fontSize': '11px',
                            'color': '#666'
                        }
                    )
                    
                ], style={
                    'backgroundColor': '#f8f9fa', 
                    'border': '2px solid #6c757d', 
                    'borderRadius': '8px', 
                    'padding': '15px', 
                    'marginBottom': '15px'
                }),
                
            ], style={
                'flex': '0 0 400px',
                'padding': '15px',
                'borderRight': '1px solid #ccc',
                'marginRight': '20px',
                'backgroundColor': '#fefefe'
            }),
            
            # Central Charts & Right Cards Container
            html.Div([
                # Step 3 Instruction Block with Clear All button as tab
                html.Div([
                    # Header row with STEP 3 box
                    html.Div([
                        html.Div([
                            html.Div("STEP 3", style={
                                'fontSize': '12px', 
                                'fontWeight': 'bold', 
                                'color': 'white', 
                                'backgroundColor': '#6c757d', 
                                'padding': '4px 8px', 
                                'borderRadius': '4px',
                                'textAlign': 'center',
                                'marginBottom': '4px'
                            }),
                            html.P("Click/drag to select regions • Multiple selections accumulate • Grid shows method combinations", 
                                   style={'fontSize': '10px', 'color': '#666', 'textAlign': 'center', 'margin': '0 0 8px 0', 'fontStyle': 'italic'})
                        ])
                    ], style={'marginBottom': '10px'}),
                    
                    # Combined Specification Grid
                    html.Div([
                        
                        # Selection status indicator (moved above grid)
                        html.Div([
                            html.Div(
                                id='selection-status-indicator',
                                style={
                                    'fontSize': '11px',
                                    'textAlign': 'center',
                                    'marginBottom': '10px'
                                }
                            )
                        ]),
                        
                        dcc.Graph(
                            id='combined-spec-grid',
                            config={
                                'displayModeBar': True,
                                'modeBarButtonsToRemove': ['pan2d'],
                                'displaylogo': False
                            }
                        )
                        
                    ], style={'flex': '1', 'minWidth': '0', 'marginBottom': '20px'}),
                    
                    
                    # Clear buttons above curves panel - full width
                    html.Div([
                        html.Div([
                            html.Button(
                                "Clear Focal Selections",
                                id='clear-focal-btn',
                                n_clicks=0,
                                style={
                                    'fontSize': '11px',
                                    'padding': '8px 0px',
                                    'backgroundColor': '#d63384',
                                    'color': 'white',
                                    'border': 'none',
                                    'borderRadius': '4px 0px 0px 4px',
                                    'cursor': 'pointer',
                                    'width': '100%',
                                    'height': '100%'
                                }
                            )
                        ], style={'flex': '1', 'marginRight': '2px'}),
                        html.Div([
                            html.Button(
                                "Clear All Selections",
                                id='clear-selections-button',
                                style={
                                    'fontSize': '11px',
                                    'padding': '8px 0px',
                                    'backgroundColor': '#dc3545',
                                    'color': 'white',
                                    'border': 'none',
                                    'borderRadius': '0px',
                                    'cursor': 'pointer',
                                    'width': '100%',
                                    'height': '100%'
                                }
                            )
                        ], style={'flex': '0.5', 'marginRight': '2px'}),
                        html.Div([
                            html.Button(
                                "Clear CF Selections",
                                id='clear-cf-btn',
                                n_clicks=0,
                                style={
                                    'fontSize': '11px',
                                    'padding': '8px 0px',
                                    'backgroundColor': '#0d6efd',
                                    'color': 'white',
                                    'border': 'none',
                                    'borderRadius': '0px 4px 4px 0px',
                                    'cursor': 'pointer',
                                    'width': '100%',
                                    'height': '100%'
                                }
                            )
                        ], style={'flex': '1'})
                    ], style={'display': 'flex', 'marginBottom': '15px', 'height': '40px'}),
                    
                    # Two Column Layout with Stacked Components
                    html.Div([
                        # Left Column: Focal Profile
                        html.Div([
                            # Focal Specification Curve
                            dcc.Graph(
                                id='spec-curve-1',
                                config={
                                    'displayModeBar': True,
                                    'modeBarButtonsToRemove': ['pan2d'],
                                    'displaylogo': False
                                }
                            ),
                            
                            # Focal Decision Importance Card (Step 4 Box)
                            html.Div([
                                html.Div(
                                    id='variable-importance-content-1',
                                    style={
                                        'padding': '15px',
                                        'border': '1px solid #ddd',
                                        'borderRadius': '5px',
                                        'minHeight': '0px',
                                        'maxHeight': '200px',
                                        'backgroundColor': '#f9f9f9',
                                        'overflowY': 'auto',
                                        'marginBottom': '15px'
                                    }
                                ),
                            ], style={
                                'backgroundColor': '#f8f9fa', 
                                'border': '2px solid #6c757d', 
                                'borderRadius': '8px', 
                                'padding': '0px', 
                                'marginBottom': '0px'
                            })
                        ], style={'flex': '1', 'minWidth': '0', 'marginRight': '15px'}),
                        
                        # Right Column: Counterfactual Profile
                        html.Div([
                            # Counterfactual Specification Curve
                            dcc.Graph(
                                id='spec-curve-2',
                                config={
                                    'displayModeBar': True,
                                    'modeBarButtonsToRemove': ['pan2d'],
                                    'displaylogo': False
                                }
                            ),
                            
                            # Counterfactual Decision Importance Card (Step 4 Box)
                            html.Div([
                                html.Div(
                                    id='variable-importance-content-2',
                                    style={
                                        'padding': '15px',
                                        'border': '1px solid #ddd',
                                        'borderRadius': '5px',
                                        'minHeight': '0px',
                                        'maxHeight': '200px',
                                        'backgroundColor': '#f9f9f9',
                                        'overflowY': 'auto',
                                        'marginBottom': '15px'
                                    }
                                ),
                            ], style={
                                'backgroundColor': '#f8f9fa', 
                                'border': '2px solid #6c757d', 
                                'borderRadius': '8px', 
                                'padding': '0px', 
                                'marginBottom': '0px'
                            })
                        ], style={'flex': '1', 'minWidth': '0'})
                        
                    ], style={'display': 'flex', 'flexDirection': 'row', 'alignItems': 'flex-start', 'gap': '15px', 'marginBottom': '20px'}),
                    
                ], style={
                    'backgroundColor': '#f8f9fa', 
                    'border': '2px solid #6c757d', 
                    'borderRadius': '8px', 
                    'padding': '15px', 
                    'marginBottom': '20px'
                }),
                
            ], style={'flexGrow': 1, 'display': 'flex', 'flexDirection': 'column', 'gap': '20px', 'maxWidth': '1400px', 'minWidth': '0'}),
            
        ], style={'display': 'flex', 'flexDirection': 'row', 'alignItems': 'flex-start', 'gap': '20px'}),
        
    ], style={'fontFamily': 'Arial, sans-serif', 'padding': '20px', 'maxWidth': '1800px', 'margin': 'auto'})

def create_specification_curve(df, profile, profile_num, previously_selected=None, highlighted_universe=None):
    """Create specification curve showing recidivism probabilities across different specifications for a single profile"""
    # Check if dataframe is empty (no analysis run yet)
    if df.empty:
        # Create empty figure with message
        fig = go.Figure()
        fig.add_annotation(
            text="No analysis data available.<br>Click 'Run Multiverse Analysis' to generate results.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, xanchor='center', yanchor='middle',
            showarrow=False, font=dict(size=16, color='gray')
        )
        fig.update_layout(
            title=f'{"Focal Profile" if profile_num == 1 else "Counterfactual Profile"}',
            xaxis_title='Analysis Pipeline Index',
            yaxis_title='Predicted Risk Probability',
            height=400,
            showlegend=False,
            yaxis=dict(range=[0, 1])
        )
        return fig
    
    # Get data for the profile
    df_plot = df.copy()
    
    # Sort by probability for better visualization
    df_plot = df_plot.sort_values('recidivism_prob')
    
    # Create the specification curve
    fig = go.Figure()
    
    # Add the profile line
    profile_color = '#d63384' if profile_num == 1 else '#0d6efd'  # Pink for Profile 1, Blue for Profile 2
    profile_name = 'Focal Profile' if profile_num == 1 else 'Counterfactual Profile'
    
    # Create sequential x-axis positions (0, 1, 2, 3...) for proper sorting visualization
    x_positions = list(range(len(df_plot)))
    
    # Create scatter plot with all points - use sequential positions for x-axis, actual universe indices for data
    fig.add_trace(go.Scatter(
        x=x_positions,  # Use sequential positions (0, 1, 2, 3...) for proper sorting
        y=df_plot['recidivism_prob'],
        mode='lines+markers',
        name=profile_name,
        line=dict(color=profile_color, width=2),
        marker=dict(size=6),
        customdata=df_plot.index,  # Store actual universe indices for selection and hover
        hovertemplate='<b>Universe %{customdata}</b><br>Probability: %{y:.3f}<extra></extra>',
        showlegend=False,  # Hide from legend
        selected=dict(marker=dict(size=10, color='yellow', opacity=0.8)),  # Highlight selected points
        unselected=dict(marker=dict(size=6, opacity=0.6))  # Dim unselected points
    ))
    
    # Add highlighting for the selected universe
    if highlighted_universe is not None and highlighted_universe < len(df_plot):
        # Add a vertical line to highlight the selected column
        fig.add_shape(
            type="line",
            x0=highlighted_universe,
            x1=highlighted_universe,
            y0=0,
            y1=1,
            line=dict(color="orange", width=3),
            layer="below"  # Place behind the scatter plot
        )
        
        # Add a background rectangle to highlight the selected column
        fig.add_shape(
            type="rect",
            x0=highlighted_universe - 0.5,
            x1=highlighted_universe + 0.5,
            y0=0,
            y1=1,
            fillcolor="rgba(255, 255, 0, 0.2)",  # Semi-transparent yellow
            line=dict(color="orange", width=1),
            layer="below"  # Place behind the scatter plot
        )
    
    
    # Add highlighting for previously selected regions
    if previously_selected and len(previously_selected) > 0:
        # Identify regions within previously selected indices
        regions = identify_regions(previously_selected, df_plot)
        
        # Add a background trace for each previously selected region
        for i, region in enumerate(regions):
            if region:  # Make sure region is not empty
                # Get the x and y coordinates for this region
                region_x = [x_positions[idx] for idx in region if idx < len(x_positions)]
                region_y = [df_plot['recidivism_prob'].iloc[idx] for idx in region if idx < len(df_plot)]
                
                if region_x and region_y:  # Make sure we have valid coordinates
                    # Add a background highlight for this region
                    fig.add_trace(go.Scatter(
                        x=region_x,
                        y=region_y,
                        mode='markers',
                        name=f'Selected Region {i+1}' if len(regions) > 1 else 'Selected Region',
                        marker=dict(
                            size=8,
                            color='rgba(255, 255, 0, 0.3)',  # Semi-transparent yellow
                            line=dict(width=2, color='orange')
                        ),
                        hoverinfo='skip',  # Skip hover for background traces
                        showlegend=False
                    ))
                    
                    # Add region boundary lines
                    if len(region_x) > 1:
                        # Add vertical lines at the start and end of each region
                        start_x = min(region_x)
                        end_x = max(region_x)
                        y_min = min(region_y)
                        y_max = max(region_y)
                        
                        # Start boundary
                        fig.add_shape(
                            type="line",
                            x0=start_x, x1=start_x,
                            y0=y_min, y1=y_max,
                            line=dict(color="white", width=10, dash="dash")
                        )
                        
                        # End boundary
                        fig.add_shape(
                            type="line",
                            x0=end_x, x1=end_x,
                            y0=y_min, y1=y_max,
                            line=dict(color="white", width=10, dash="dash")
                        )
    

    
    # Profile highlighting removed since we now use dynamic data from R analysis
    
    # Add legend for visual indicators
    legend_annotations = []
    if previously_selected and len(previously_selected) > 0:
        regions = identify_regions(previously_selected, df_plot)
        # if len(regions) > 1:
        #     legend_annotations.append(
        #         dict(
        #             x=0.02, y=0.98,
        #             xref='paper', yref='paper',
        #             text="🟡 Previously Selected Regions<br>🟠 Region Boundaries",
        #             showarrow=False,
        #             font=dict(size=10, color='#666'),
        #             bgcolor='rgba(255, 255, 255, 0.8)',
        #             bordercolor='#ccc',
        #             borderwidth=1
        #         )
        #     )
        # else:
        #     legend_annotations.append(
        #         dict(
        #             x=0.02, y=0.98,
        #             xref='paper', yref='paper',
        #             text="🟡 Previously Selected Region<br>🟠 Region Boundaries",
        #             showarrow=False,
        #             font=dict(size=10, color='#666'),
        #             bgcolor='rgba(255, 255, 255, 0.8)',
        #             bordercolor='#ccc',
        #             borderwidth=1
        #         )
        #     )
    
    fig.update_layout(
        title=f'{profile_name}',
        xaxis_title='Analysis Pipeline Index',
        yaxis_title='Predicted Risk Probability',
        height=400,
        showlegend=False,
        hovermode='x unified',
        dragmode='select',  # Enable selection mode by default
        yaxis=dict(
            range=[0, 1],  # Force y-axis range from 0 to 1 for consistent comparison
            tickmode='linear',
            dtick=0.2,  # Show ticks at 0, 0.2, 0.4, 0.6, 0.8, 1.0
            tickformat='.1f'  # Format ticks as 0.0, 0.2, 0.4, etc.
        ),
        annotations=legend_annotations
    )
    
    # Configure x-axis to show sequential positions with universe labels
    if len(df_plot.index) <= 50:  # Show all positions if 50 or fewer
        fig.update_xaxes(
            tickmode='array',
            tickvals=x_positions,
            ticktext=[f"{df_plot.index[i]}" for i in x_positions]
        )
    else:  # Show every nth position for readability
        step = max(1, len(x_positions) // 10)
        selected_positions = x_positions[::step]
        fig.update_xaxes(
            tickmode='array',
            tickvals=selected_positions,
            ticktext=[f"{df_plot.index[i]}" for i in selected_positions]
        )
    

    
    return fig

def create_combined_specification_grid(df1, df2, selected_universes_1=None, selected_universes_2=None, highlighted_universe_1=None, highlighted_universe_2=None):
    """Create combined specification grid showing procedural choices from both profiles"""
    # Check if both dataframes are empty
    if df1.empty and df2.empty:
        # Create empty figure with message
        fig = go.Figure()
        fig.add_annotation(
            text="No analysis data available.<br>Click 'Run Multiverse Analysis' to generate results.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, xanchor='center', yanchor='middle',
            showarrow=False, font=dict(size=16, color='gray')
        )
        fig.update_layout(
            title='Combined Specification Grid',
            xaxis_title='Analysis Pipeline Index',
            yaxis_title='Pipeline Choices',
            height=500,
            showlegend=False
        )
        return fig
    
    # Combine data from both profiles with region tracking
    combined_data = []
    universe_labels = []
    profile_colors = []
    region_info = []  # Track which region each universe belongs to
    
    # Process focal profile data (df1)
    if not df1.empty:
        df1_sorted = df1.sort_values('recidivism_prob')
        if selected_universes_1 is not None and len(selected_universes_1) > 0:
            try:
                selected_universe_ids_1 = [df1_sorted.index[i] for i in selected_universes_1 if isinstance(i, int) and i < len(df1_sorted)]
                if selected_universe_ids_1:
                    df1_display = df1_sorted.loc[selected_universe_ids_1]
                    # Identify regions within focal selections
                    focal_regions = identify_regions(selected_universes_1, df1_sorted)
                else:
                    df1_display = df1_sorted
                    focal_regions = [list(range(len(df1_sorted)))]
            except (IndexError, KeyError, TypeError):
                df1_display = df1_sorted
                focal_regions = [list(range(len(df1_sorted)))]
        else:
            df1_display = df1_sorted
            focal_regions = [list(range(len(df1_sorted)))]
        
        # Add focal profile data with region tracking
        # Collect all selected universes and sort them by recidivism probability to maintain curve ordering
        all_selected_universes = []
        for region in focal_regions:
            all_selected_universes.extend(region)
        
        # Sort all selected universes by their recidivism probability to match curve ordering
        all_selected_universes.sort(key=lambda idx: df1_sorted.iloc[idx]['recidivism_prob'] if idx < len(df1_sorted) else 0)
        
        universe_counter = 0
        for universe_idx in all_selected_universes:
            if universe_idx < len(df1_sorted):
                row = df1_sorted.iloc[universe_idx]
                combined_data.append(row)
                universe_labels.append(f"F{universe_counter}")
                profile_colors.append('#d63384')  # Pink for focal
                # Find which region this universe belongs to
                region_idx = next((i for i, region in enumerate(focal_regions) if universe_idx in region), 0)
                region_info.append(f"Focal Region {region_idx + 1}")
                universe_counter += 1
    
    # Process counterfactual profile data (df2)
    if not df2.empty:
        df2_sorted = df2.sort_values('recidivism_prob')
        if selected_universes_2 is not None and len(selected_universes_2) > 0:
            try:
                selected_universe_ids_2 = [df2_sorted.index[i] for i in selected_universes_2 if isinstance(i, int) and i < len(df2_sorted)]
                if selected_universe_ids_2:
                    df2_display = df2_sorted.loc[selected_universe_ids_2]
                    # Identify regions within counterfactual selections
                    cf_regions = identify_regions(selected_universes_2, df2_sorted)
                else:
                    df2_display = df2_sorted
                    cf_regions = [list(range(len(df2_sorted)))]
            except (IndexError, KeyError, TypeError):
                df2_display = df2_sorted
                cf_regions = [list(range(len(df2_sorted)))]
        else:
            df2_display = df2_sorted
            cf_regions = [list(range(len(df2_sorted)))]
        
        # Add counterfactual profile data with region tracking
        # Collect all selected universes and sort them by recidivism probability to maintain curve ordering
        all_selected_universes = []
        for region in cf_regions:
            all_selected_universes.extend(region)
        
        # Sort all selected universes by their recidivism probability to match curve ordering
        all_selected_universes.sort(key=lambda idx: df2_sorted.iloc[idx]['recidivism_prob'] if idx < len(df2_sorted) else 0)
        
        universe_counter = 0
        for universe_idx in all_selected_universes:
            if universe_idx < len(df2_sorted):
                row = df2_sorted.iloc[universe_idx]
                combined_data.append(row)
                universe_labels.append(f"CF{universe_counter}")
                profile_colors.append('#0d6efd')  # Blue for counterfactual
                # Find which region this universe belongs to
                region_idx = next((i for i, region in enumerate(cf_regions) if universe_idx in region), 0)
                region_info.append(f"CF Region {region_idx + 1}")
                universe_counter += 1
    
    if not combined_data:
        # No data to display
        fig = go.Figure()
        fig.add_annotation(
            text="No universes selected for display.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, xanchor='center', yanchor='middle',
            showarrow=False, font=dict(size=16, color='gray')
        )
        fig.update_layout(
            title='Combined Specification Grid',
            xaxis_title='Analysis Pipeline Index',
            yaxis_title='Pipeline Choices',
            height=500,
            showlegend=False
        )
        return fig
    
    # Create the grid data matrix
    grid_data = []
    y_labels = []
    
    # Define procedural choices and their options for new NC garden approach
    procedural_choices = [
        ('define_recid', ['5yr', '4yr', '3yr', '2yr', '1yr']),
        ('predictor', ['fair', 'schmidt', 'full']),
        ('model', ['logistic', 'survival']),
        ('imbalancing', ['weight', 'female only', 'male only', 'over', 'under']),
        ('age_category', ['nij', 'compas', 'schmidt_age', 'year']),
        ('split', ['8:2', '7:3', '6:4', '1:2']),
        ('convert_followup', ['general_def', 'strict_def']),
        ('scaling', ['schmidt_scale', 'no_scaling'])  # Show both scaling options
    ]
    
    # Create the grid data matrix
    for choice_name, options in procedural_choices:
        column_display_name = choice_name.replace('_', ' ').title()
        
        # Add the options FIRST (they will appear below the header in the visual)
        for option in options:
            y_labels.append(f"    {option}")
            row = []
            
            for data_row in combined_data:
                # Determine if this option is selected for this universe
                if choice_name == 'define_recid':
                    selected = data_row['define_recid'] == option
                elif choice_name == 'predictor':
                    selected = data_row['predictor'] == option
                elif choice_name == 'imbalancing':
                    selected = data_row['imbalancing'] == option
                elif choice_name == 'age_category':
                    selected = data_row['age_category'] == option
                elif choice_name == 'split':
                    selected = data_row['split'] == option
                elif choice_name == 'scaling':
                    selected = data_row['scaling'] == option
                elif choice_name == 'convert_followup':
                    selected = data_row['convert_followup'] == option
                elif choice_name == 'model':
                    selected = data_row['model'] == option
                else:
                    selected = False
                
                # Binary selection: 1 if selected, 0.5 if not
                if selected:
                    row.append(1)
                else:
                    row.append(0.5)
            
            grid_data.append(row)
        
        # Add the category header AFTER the options (will appear above in visual due to y-axis ordering)
        y_labels.append(f"{column_display_name}")
        header_row = [0] * len(combined_data)
        grid_data.append(header_row)
    
    # Create colorscale for combined view
    colorscale = [
        [0, '#e9ecef'],    # Gray for category headers
        [0.5, '#ffffff'],  # White for unselected
        [1, '#6c757d']     # Gray for selected (neutral color for combined view)
    ]
    
    # Create custom hover text with universe indices
    hover_text = []
    for i, row in enumerate(grid_data):
        hover_row = []
        for j, cell_value in enumerate(row):
            if j < len(combined_data):
                # combined_data contains pandas Series, so we can get the index
                universe_id = combined_data[j].name if hasattr(combined_data[j], 'name') else f"U{j}"
                hover_row.append(f"Universe: {universe_id}<br>Value: {cell_value}")
            else:
                hover_row.append(f"Value: {cell_value}")
        hover_text.append(hover_row)
    
    # Create the heatmap
    fig = go.Figure(data=go.Heatmap(
        z=grid_data,
        x=list(range(len(combined_data))),
        y=y_labels,
        colorscale=colorscale,
        showscale=False,
        hoverongaps=False,
        hoverinfo='text',
        text=hover_text,
        zmin=0,
        zmax=1
    ))
    
    # Add highlighting for the selected universe columns
    if highlighted_universe_1 is not None and highlighted_universe_1 < len(combined_data):
        # Only highlight if it's in the focal section of the combined grid
        # Count focal universes in the combined data
        focal_count = sum(1 for info in region_info if "Focal" in info)
        
        if highlighted_universe_1 < focal_count:
            # Add a vertical rectangle to highlight the selected column for profile 1
            fig.add_shape(
                type="rect",
                x0=highlighted_universe_1 - 0.5,
                x1=highlighted_universe_1 + 0.5,
                y0=-0.5,
                y1=len(y_labels) - 0.5,
                fillcolor="rgba(214, 51, 132, 0.4)",  # More visible pink for profile 1
                line=dict(color="#d63384", width=3),
                layer="above"  # Place above the heatmap
            )
    
    # Counterfactual profile highlighting disabled - no individual universe highlighting
    
    
    # Add region separators
    if len(combined_data) > 1:
        # Find region boundaries - add separators between different regions
        region_boundaries = []
        current_region = region_info[0] if region_info else ""
        
        for i in range(1, len(region_info)):
            if region_info[i] != current_region:
                region_boundaries.append(i)
                current_region = region_info[i]
        
        # Add vertical lines at region boundaries
        for boundary in region_boundaries:
            fig.add_shape(
                type="line",
                x0=boundary - 0.5, x1=boundary - 0.5,
                y0=-0.5, y1=len(y_labels) - 0.5,
                line=dict(color="white", width=10, dash="solid")
            )
        
        # Add region background colors
        current_region = region_info[0] if region_info else ""
        region_start = 0
        
        for i in range(len(region_info) + 1):
            if i == len(region_info) or region_info[i] != current_region:
                # End of current region, add background
                if i > region_start:
                    # Use blue background for CF regions, magenta for Focal regions
                    region_color = 'rgba(0, 123, 255, 0.1)' if 'CF' in current_region else 'rgba(255, 0, 255, 0.1)'
                    fig.add_shape(
                        type="rect",
                        x0=region_start - 0.5, x1=i - 0.5,
                        y0=-0.5, y1=len(y_labels) - 0.5,
                        fillcolor=region_color,
                        line=dict(width=0)
                    )
                
                if i < len(region_info):
                    current_region = region_info[i]
                    region_start = i
    
    # Add legend showing profile colors and region information
    legend_text = f"Showing {len(combined_data)} Universes"
    
    # Count regions for each profile
    focal_region_count = 0
    cf_region_count = 0
    
    if region_info:
        unique_regions = list(set(region_info))
        focal_region_count = len([r for r in unique_regions if 'Focal' in r])
        cf_region_count = len([r for r in unique_regions if 'CF' in r])
    
    if focal_region_count > 0:
        legend_text += f" (Focal: {focal_region_count} region{'s' if focal_region_count > 1 else ''})"
    if cf_region_count > 0:
        legend_text += f" (CF: {cf_region_count} region{'s' if cf_region_count > 1 else ''})"
    
    # Add region separator legend
    if len(region_info) > 1 and len(set(region_info)) > 1:
        legend_text += "<br>Orange lines separate regions"
    
    fig.add_annotation(
        x=0.02,
        y=0.98,
        xref='paper',
        yref='paper',
        text=legend_text,
        showarrow=False,
        font=dict(size=12, color='#333'),
        bgcolor='rgba(108, 117, 125, 0.1)',
        bordercolor='#6c757d',
        borderwidth=1
    )
    
    # Customize the layout
    fig.update_layout(
        title='Combined Specification Grid',
        xaxis_title='Analysis Pipeline Index',
        yaxis_title='Pipeline Choices',
        height=500,
        showlegend=False,
        margin=dict(l=150, r=50, t=80, b=50),  # Increased left margin for y-axis labels
        yaxis=dict(
            tickmode='array',
            tickvals=list(range(len(y_labels))),
            ticktext=y_labels,
            tickfont=dict(size=10)
        )
    )
    
    # Update x-axis to show universe labels
    if len(combined_data) <= 20:
        fig.update_xaxes(
            tickmode='array', 
            tickvals=list(range(len(combined_data))), 
            ticktext=universe_labels
        )
    else:
        step = max(1, len(combined_data) // 10)
        selected_positions = list(range(0, len(combined_data), step))
        fig.update_xaxes(
            tickmode='array', 
            tickvals=selected_positions, 
            ticktext=[universe_labels[i] for i in selected_positions]
        )
    
    return fig

def create_specification_grid(df, profile_num, selected_universes=None, highlighted_universe=None):
    """Create specification grid showing procedural choices vs universe indexes (sorted by probability)"""
    # Check if dataframe is empty (no analysis run yet)
    if df.empty:
        # Create empty figure with message
        fig = go.Figure()
        fig.add_annotation(
            text="No analysis data available.<br>Click 'Run Multiverse Analysis' to generate results.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, xanchor='center', yanchor='middle',
            showarrow=False, font=dict(size=16, color='gray')
        )
        fig.update_layout(
            title=f'{"Focal Profile" if profile_num == 1 else "Counterfactual Profile"} Specification Grid',
            xaxis_title='Analysis Pipeline Index',
            yaxis_title='Pipeline Choices',
            height=500,
            showlegend=False
        )
        return fig
    
    # Sort data by recidivism probability to match specification curve ordering
    df_sorted = df.sort_values('recidivism_prob')
    
    # Filter based on selected universes if provided
    if selected_universes is not None and len(selected_universes) > 0:
        # Convert selected universe indices to actual universe IDs
        try:
            selected_universe_ids = [df_sorted.index[i] for i in selected_universes if isinstance(i, int) and i < len(df_sorted)]
            if selected_universe_ids:
                # Sort selected universe IDs by their recidivism probability to maintain curve ordering
                selected_universe_ids.sort(key=lambda uid: df_sorted.loc[uid, 'recidivism_prob'])
                df_display = df_sorted.loc[selected_universe_ids]
                title_suffix = f" - {len(selected_universe_ids)} Selected Universes"
                # Identify regions within selections
                regions = identify_regions(selected_universes, df_sorted)
            else:
                df_display = df_sorted
                title_suffix = " - All Universes"
                regions = [list(range(len(df_sorted)))]
        except (IndexError, KeyError, TypeError):
            # Fallback to showing all universes if there's an error
            df_display = df_sorted
            title_suffix = " - All Universes"
            regions = [list(range(len(df_sorted)))]
    else:
        # Show all universes
        df_display = df_sorted
        title_suffix = " - All Universes"
        regions = [list(range(len(df_sorted)))]
    
    # Define procedural choices and their options, grouped as shown in the image
    procedural_choices = [
        # Recidivism time periods (grouped together)
        ('define_recid', ['5yr', '4yr', '3yr', '2yr', '1yr']),
        # Predictor methods (grouped together)
        ('predictor', ['fair', 'schmidt', 'full']),
        # Model types (grouped together)
        ('model', ['logistic', 'survival']),
        # Imbalancing methods (grouped together)
        ('imbalancing', ['weight', 'female only', 'male only', 'over', 'under']),
        # Age categories (grouped together)
        ('age_category', ['nij', 'compas', 'schmidt_age', 'year']),
        # Split ratios (grouped together)
        ('split', ['8:2', '7:3', '6:4', '1:2']),
        # Convert followup methods (grouped together)
        ('convert_followup', ['general_def', 'strict_def']),
        # Scaling methods (grouped together)
        ('scaling', ['no_scaling', 'schmidt_scale'])
    ]
    
    # Create the grid data matrix
    grid_data = []
    y_labels = []
    
    # Add category headers and their options
    for choice_name, options in procedural_choices:
        column_display_name = choice_name.replace('_', ' ').title()
        
        # Add the options FIRST (they will appear below the header in the visual)
        for option in options:
            y_labels.append(f"    {option}")  # More indented option for better hierarchy
            row = []
            
            for _, row_data in df_display.iterrows():
                # Determine if this option is selected for this universe
                if choice_name == 'define_recid':
                    selected = row_data['define_recid'] == option
                elif choice_name == 'predictor':
                    selected = row_data['predictor'] == option
                elif choice_name == 'imbalancing':
                    selected = row_data['imbalancing'] == option
                elif choice_name == 'age_category':
                    selected = row_data['age_category'] == option
                elif choice_name == 'split':
                    selected = row_data['split'] == option
                elif choice_name == 'scaling':
                    selected = row_data['scaling'] == option
                elif choice_name == 'convert_followup':
                    selected = row_data['convert_followup'] == option
                elif choice_name == 'model':
                    selected = row_data['model'] == option
                else:
                    selected = False
                
                # Simple binary selection: 1 if selected, 0.5 if not
                if selected:
                    row.append(1)
                else:
                    row.append(0.5)
            
            grid_data.append(row)
        
        # Add the category header AFTER the options (will appear above in visual due to y-axis ordering)
        y_labels.append(f"{column_display_name}")  # Category header
        # Create a header row (all 0 to indicate it's a header)
        header_row = [0] * len(df_display)
        grid_data.append(header_row)
    
    # Create colorscale: header rows, white for unselected, profile color for selected
    if profile_num == 1:
        colorscale = [
            [0, '#e9ecef'],    # Darker gray for category headers to make them stand out
            [0.5, '#ffffff'],  # White for unselected
            [1, '#d63384']     # Red for Profile 1 (matching specification curve)
        ]
    else:
        colorscale = [
            [0, '#e9ecef'],    # Darker gray for category headers to make them stand out
            [0.5, '#ffffff'],  # White for unselected
            [1, '#0d6efd']     # Blue for Profile 2
        ]
    
    # Create custom hover text with universe indices
    hover_text = []
    for i, row in enumerate(grid_data):
        hover_row = []
        for j, cell_value in enumerate(row):
            if j < len(df_display):
                universe_id = df_display.index[j]
                hover_row.append(f"Universe: {universe_id}<br>Value: {cell_value}")
            else:
                hover_row.append(f"Value: {cell_value}")
        hover_text.append(hover_row)
    
    # Create the heatmap
    fig = go.Figure(data=go.Heatmap(
        z=grid_data,
        x=list(range(len(df_display))),  # Use sequential positions (0, 1, 2, 3...) to match curve
        y=y_labels,
        colorscale=colorscale,
        showscale=False,  # Remove the colorbar to save space
        hoverongaps=False,
        hoverinfo='text',
        text=hover_text,
        zmin=0,
        zmax=1
    ))
    
    # Add highlighting for the selected universe column
    if highlighted_universe is not None and highlighted_universe < len(df_display):
        # Add a vertical rectangle to highlight the selected column
        fig.add_shape(
            type="rect",
            x0=highlighted_universe - 0.5,
            x1=highlighted_universe + 0.5,
            y0=-0.5,
            y1=len(y_labels) - 0.5,
            fillcolor="rgba(255, 255, 0, 0.4)",  # More visible yellow
            line=dict(color="orange", width=3),
            layer="above"  # Place above the heatmap
        )
    
    
    # Add region separators for single profile grid
    if len(regions) > 1 and len(df_display) > 1:
        # Find region boundaries in the displayed data
        region_boundaries = []
        current_region_idx = 0
        region_start = 0
        
        for i, (_, row) in enumerate(df_display.iterrows()):
            # Find which region this universe belongs to
            universe_idx = df_sorted.index.get_loc(row.name)
            region_idx = 0
            for j, region in enumerate(regions):
                if universe_idx in region:
                    region_idx = j
                    break
            
            if region_idx != current_region_idx:
                region_boundaries.append(i)
                current_region_idx = region_idx
        
        # Add vertical lines at region boundaries
        for boundary in region_boundaries:
            fig.add_shape(
                type="line",
                x0=boundary - 0.5, x1=boundary - 0.5,
                y0=-0.5, y1=len(y_labels) - 0.5,
                line=dict(color="white", width=10, dash="solid")
            )
        
        # Add region background colors
        current_region_idx = 0
        region_start = 0
        
        for i, (_, row) in enumerate(df_display.iterrows()):
            # Find which region this universe belongs to
            universe_idx = df_sorted.index.get_loc(row.name)
            region_idx = 0
            for j, region in enumerate(regions):
                if universe_idx in region:
                    region_idx = j
                    break
            
            if region_idx != current_region_idx or i == len(df_display) - 1:
                # End of current region, add background
                if i > region_start:
                    # Use blue background for CF regions (profile_num == 2), magenta for Focal regions
                    region_color = 'rgba(0, 123, 255, 0.1)' if profile_num == 2 else 'rgba(255, 0, 255, 0.1)'
                    fig.add_shape(
                        type="rect",
                        x0=region_start - 0.5, x1=i - 0.5,
                        y0=-0.5, y1=len(y_labels) - 0.5,
                        fillcolor=region_color,
                        line=dict(width=0)
                    )
                
                if i < len(df_display) - 1:
                    current_region_idx = region_idx
                    region_start = i
    
    # Add a legend for the current view
    if selected_universes is not None and len(selected_universes) > 0:
        region_count = len(regions)
        legend_text = f"Showing {len(selected_universes)} Selected Universes"
        if region_count > 1:
            legend_text += f" in {region_count} regions"
        legend_color = '#d63384' if profile_num == 1 else '#0d6efd'
    else:
        legend_text = "Showing All Universes"
        legend_color = '#666'
    
    # Add region separator legend for single profile
    if len(regions) > 1 and len(df_display) > 1:
        legend_text += "<br>Orange lines separate regions"
    
    # Parse color safely
    try:
        if legend_color.startswith('#') and len(legend_color) == 7:
            r = int(legend_color[1:3], 16)
            g = int(legend_color[3:5], 16)
            b = int(legend_color[5:7], 16)
            bgcolor = f'rgba({r}, {g}, {b}, 0.1)'
        else:
            bgcolor = 'rgba(102, 102, 102, 0.1)'  # Default gray
    except (ValueError, IndexError):
        bgcolor = 'rgba(102, 102, 102, 0.1)'  # Default gray
    
    fig.add_annotation(
        x=0.02,
        y=0.98,
        xref='paper',
        yref='paper',
        text=legend_text,
        showarrow=False,
        font=dict(size=12, color=legend_color),
        bgcolor=bgcolor,
        bordercolor=legend_color,
        borderwidth=1
    )
    
    # Customize the layout
    fig.update_layout(
        title=f'{"Focal Profile" if profile_num == 1 else "Counterfactual Profile"} Specification Grid{title_suffix}',
        xaxis_title=f'Analysis Pipeline Index',
        yaxis_title='Pipeline Choices',
        height=500,
        showlegend=False,
        margin=dict(l=150, r=50, t=80, b=50),  # Increased left margin for y-axis labels
        yaxis=dict(
            tickmode='array',
            tickvals=list(range(len(y_labels))),
            ticktext=y_labels,
            tickfont=dict(size=10)
        )
    )
    
    # Update x-axis to show appropriate tick labels
    if len(df_display.index) <= 20:  # Show all labels if 20 or fewer universes
        fig.update_xaxes(
            tickmode='array', 
            tickvals=list(range(len(df_display))), 
            ticktext=[f"{df_display.index[i]}" for i in range(len(df_display))]
        )
    else:  # Show every nth label for readability
        step = max(1, len(df_display.index) // 10)
        selected_positions = list(range(0, len(df_display), step))
        fig.update_xaxes(
            tickmode='array', 
            tickvals=selected_positions, 
            ticktext=[f"{df_display.index[i]}" for i in selected_positions]
        )
    
    return fig



def get_profile_summary(profile):
    """Generate a summary text for a profile"""
    return html.Div([
        html.P(f"Age: {profile['age']} years"),
        html.P(f"Gender: {profile['gender'].title()}"),
        html.P(f"Race: {profile['race'].replace('_', ' ').title()}")
    ])



def create_boxplot(df, profile, profile_num):
    """Create boxplot showing the distribution of recidivism probabilities for a specific profile"""
    fig = go.Figure()
    
    # Add boxplot for all data
    profile_color = '#d63384' if profile_num == 1 else '#0d6efd'
    fig.add_trace(go.Box(
        y=df['recidivism_prob'],
        name='All Specifications',
        boxpoints='outliers',
        marker_color=profile_color
    ))
    
    # Add vertical line for profile prediction
    profile_prob = df[
        (df['scaling'] == profile['scaling']) &
        (df['convert_followup'] == profile['convert_followup']) &
        (df['split'] == profile['split']) &
        (df['age_category'] == profile['age_category']) &
        (df['imbalancing'] == profile['imbalancing']) &
        (df['model'] == profile['model']) &
        (df['predictor'] == profile['predictor']) &
        (df['define_recid'] == profile['define_recid'])
    ]['recidivism_prob'].iloc[0] if len(df[
        (df['scaling'] == profile['scaling']) &
        (df['convert_followup'] == profile['convert_followup']) &
        (df['split'] == profile['split']) &
        (df['age_category'] == profile['age_category']) &
        (df['imbalancing'] == profile['imbalancing']) &
        (df['model'] == profile['model']) &
        (df['predictor'] == profile['predictor']) &
        (df['define_recid'] == profile['define_recid'])
    ]) > 0 else 0
    
    profile_color = '#d63384' if profile_num == 1 else '#0d6efd'
    profile_name = f'Profile {profile_num}'
    
    fig.add_hline(y=profile_prob, line_dash="dash", line_color=profile_color, 
                  annotation_text="")  # Remove profile annotation text
    
    fig.update_layout(
        #title=f'{profile_name} - Distribution of Recidivism Probabilities',
        yaxis_title='Predicted Risk Probability',
        height=400,
        showlegend=False,  # Remove legend
        yaxis=dict(
            range=[0, 1],  # Force y-axis range from 0 to 1 for consistent comparison
            tickmode='linear',
            dtick=0.2,  # Show ticks at 0, 0.2, 0.4, 0.6, 0.8, 1.0
            tickformat='.1f'  # Format ticks as 0.0, 0.2, 0.4, etc.
        )
    )
    
    return fig

def get_variable_importance_r(df):
    """Get variable importance using R decision tree analysis"""
    try:
        # Check if dataframe is empty or has insufficient data
        if df.empty or len(df) < 2:
            print("Warning: Insufficient data for variable importance analysis")
            return []
            
        with localconverter(robjects.default_converter + pandas2ri.converter):
            # Convert pandas DataFrame to R DataFrame
            # Reset index to avoid row names issues during conversion
            df_reset = df.reset_index(drop=True)
            r_df = robjects.conversion.py2rpy(df_reset)
            
            # Call R function to get variable importance
            r_get_var_importance = robjects.globalenv['get_variable_importance']
            try:
                var_importance = r_get_var_importance(r_df)
            except Exception as r_error:
                print(f"R error in variable importance: {r_error}")
                # Return default variable importance if R fails
                return [("scaling", 0), ("convert_followup", 0), ("split", 0), ("age_category", 0), 
                       ("imbalancing", 0), ("model", 0), ("predictor", 0), ("define_recid", 0)]
            
            # Convert R vector to Python dictionary
            var_importance_dict = {}
            
            # Debug: print what we get from R
            print(f"R var_importance type: {type(var_importance)}")
            print(f"R var_importance length: {len(var_importance)}")
            print(f"R var_importance names: {var_importance.names if hasattr(var_importance, 'names') else 'No names'}")
            print(f"R var_importance values: {list(var_importance)}")
            
            # Expected variable names in the order they should appear
            expected_vars = ["predictor", "define_recid", "imbalancing", "model",
                           "split", "scaling", "convert_followup", "age_category"]
            
            # Use the variable names from the R results (which should now be properly named)
            # The R function now returns named vectors, so we can extract the names
            if hasattr(var_importance, 'names') and var_importance.names is not None:
                print(f"Processing {len(var_importance.names)} names from R")
                # Use the names from R
                for i, name in enumerate(var_importance.names):
                    print(f"Processing name {i}: '{name}' (type: {type(name)})")
                    # Convert R string to Python string and clean it
                    if isinstance(name, str):
                        clean_name = name.strip()
                    else:
                        clean_name = str(name).strip() if name is not None else ""
                    
                    if clean_name and clean_name != "NA" and clean_name != "" and clean_name != "None":  # Skip null, NA, or empty names
                        var_importance_dict[clean_name] = float(var_importance[i])
                        print(f"Added {clean_name}: {float(var_importance[i])}")
                    else:
                        # If name is invalid, use expected variable name as fallback
                        if i < len(expected_vars):
                            fallback_name = expected_vars[i]
                        else:
                            fallback_name = f"Variable_{i+1}"
                        var_importance_dict[fallback_name] = float(var_importance[i])
                        print(f"Using fallback name {fallback_name}: {float(var_importance[i])}")
            else:
                print("No names from R, using expected variable names")
                # Fallback: use expected variable names
                for i in range(len(var_importance)):
                    if i < len(expected_vars):
                        var_importance_dict[expected_vars[i]] = float(var_importance[i])
                    else:
                        var_importance_dict[f"Variable_{i+1}"] = float(var_importance[i])
            
            # Sort by importance (descending) - ensure proper ordering
            sorted_importance = sorted(var_importance_dict.items(), key=lambda x: x[1], reverse=True)
            
            # Debug: print the sorted importance for verification
            print(f"Variable importance analysis completed:")
            print(f"  - Number of variables analyzed: {len(var_importance_dict)}")
            print(f"  - Variable names: {list(var_importance_dict.keys())}")
            print(f"  - Sorted importance: {sorted_importance}")
            
            return sorted_importance
    except Exception as e:
        print(f"Error in R variable importance calculation: {e}")
        import traceback
        traceback.print_exc()
        return []

def get_tree_nodes_r(df):
    """Get regression tree nodes and split rules using R analysis"""
    try:
        # Check if dataframe is empty or has insufficient data
        if df.empty or len(df) < 2:
            print("Warning: Insufficient data for tree nodes analysis")
            return []
            
        with localconverter(robjects.default_converter + pandas2ri.converter):
            # Convert pandas DataFrame to R DataFrame
            # Reset index to avoid row names issues during conversion
            df_reset = df.reset_index(drop=True)
            r_df = robjects.conversion.py2rpy(df_reset)
            
            # Call R function to get tree nodes
            r_get_tree_nodes = robjects.globalenv['get_tree_nodes']
            try:
                tree_nodes = r_get_tree_nodes(r_df)
            except Exception as r_error:
                print(f"R error in tree nodes: {r_error}")
                # Return empty list if R fails
                return []
            
            # Convert R list to Python list
            node_rules = []
            for i in range(len(tree_nodes)):
                node = tree_nodes[i]
                node_rules.append({
                    'rank': int(node[0]),  # rank
                    'node': int(node[1]),  # node number
                    'variable': str(node[2]),  # variable name
                    'rule': str(node[3])  # split rule
                })
            
            # Debug: print the tree nodes for verification
            print(f"Tree nodes: {node_rules}")
            
            return node_rules
    except Exception as e:
        print(f"Error in R tree nodes calculation: {e}")
        import traceback
        traceback.print_exc()
        return []

def identify_regions(selected_indices, df_sorted):
    """Identify separate regions within selected indices, preserving selection order"""
    if not selected_indices:
        return []
    
    # Don't sort the indices - preserve the original selection order
    # This ensures regions appear in the same order as selected in the curves
    regions = []
    current_region = [selected_indices[0]]
    
    for i in range(1, len(selected_indices)):
        # If the next index is adjacent to the current region, extend it
        if selected_indices[i] == selected_indices[i-1] + 1:
            current_region.append(selected_indices[i])
        else:
            # Gap found, start a new region
            regions.append(current_region)
            current_region = [selected_indices[i]]
    
    # Add the last region
    regions.append(current_region)
    
    return regions

def get_combined_regional_variable_importance_display(df, profile, profile_num, selected_universes=None):
    """Generate variable importance display treating all selected regions as one combined batch"""
    # Check if dataframe is empty
    if df.empty:
        return html.Div([
            html.P("No analysis data available. Click 'Run Multiverse Analysis' to generate results.", 
                   style={'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '50px'}),
            html.Hr(style={'margin': '15px 0'}),
            html.P(f"Profile: {'Focal' if profile_num == 1 else 'Counterfactual'}", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'}),
            html.P(f"Dataset Size: No data available", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'})
        ])
    
    # Sort data by recidivism probability to match specification curve ordering
    df_sorted = df.sort_values('recidivism_prob')
    
    # Identify regions within selections
    regions = identify_regions(selected_universes or [], df_sorted)
    
    if not regions:
        # No regions, use all data
        analysis_df = df_sorted
        region_title = "All Universes"
        region_count = 1
    else:
        # Combine all regions into one dataset for analysis
        all_selected_indices = []
        for region in regions:
            all_selected_indices.extend(region)
        
        # Get unique universe IDs for all selected regions
        region_universe_ids = [df_sorted.index[idx] for idx in all_selected_indices if idx < len(df_sorted)]
        analysis_df = df_sorted.loc[region_universe_ids] if region_universe_ids else df_sorted
        region_count = len(regions)
        region_title = f"Combined Analysis ({region_count} Region{'s' if region_count > 1 else ''})"
    
    # Create variable importance display
    var_importance_html = []
    
    # Add analysis status indicator
    if region_count > 1:
        status_indicator = html.Div([
            html.Span("", style={'fontSize': '14px', 'marginRight': '5px'}),
            html.Span(f"Analyzing {region_count} regions from {'focal' if profile_num == 1 else 'counterfactual'} profile", 
                     style={'fontSize': '11px', 'color': '#28a745', 'fontWeight': 'bold'})
        ], style={'marginBottom': '8px', 'padding': '4px 8px', 'backgroundColor': '#d4edda', 'borderRadius': '4px', 'border': '1px solid #c3e6cb'})
        var_importance_html.append(status_indicator)
    
    var_importance_html.append(html.H5(f"Decision Analysis ({region_title})", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '0', 'marginBottom': '6px', 'fontSize': '14px'}))
    var_importance_html.append(html.P("Most impactful methods on recidivism probability for the selected regions:", style={'fontSize': '12px', 'color': '#666', 'marginBottom': '5px'}))
    
    # Get variable importance for the combined dataset
    var_importance = get_variable_importance_r(analysis_df)
    
    # Get tree nodes for the combined dataset
    tree_nodes = get_tree_nodes_r(analysis_df)
    
    if var_importance:
        # Create bar plot for the combined analysis
        valid_vars = [(name, val) for name, val in var_importance if val is not None and val > 0]
        if valid_vars:
            top_vars_sorted = valid_vars[:5]
            var_names = [var[0] for var in top_vars_sorted]
            var_values = [var[1] for var in top_vars_sorted]
            
            # Create bar plot
            bar_fig = go.Figure(data=[
                go.Bar(
                    x=var_values,
                    y=var_names,
                    orientation='h',
                    marker_color='#d63384' if profile_num == 1 else '#0d6efd',
                    marker_line_color='#fff',
                    marker_line_width=1,
                    text=[f'{val:.3f}' for val in var_values],
                    textposition='auto',
                    textfont=dict(size=10, color='white')
                )
            ])
            
            bar_fig.update_layout(
                yaxis=dict(autorange='reversed')
            )
            
            bar_fig.update_layout(
                title=None,
                yaxis_title=None,
                height=120,
                margin=dict(l=8, r=8, t=5, b=5),
                showlegend=False,
                xaxis=dict(showgrid=False, zeroline=False, title_font=dict(size=9)),
                yaxis=dict(showgrid=False, zeroline=False),
                plot_bgcolor='rgba(0,0,0,0)',
                paper_bgcolor='rgba(0,0,0,0)'
            )
            
            # Wrap the chart in a horizontally scrollable container
            var_importance_html.append(html.Div([
                dcc.Graph(
                    figure=bar_fig,
                    config={'displayModeBar': False},
                    style={'height': '100px', 'marginBottom': '10px', 'minWidth': '300px'}
                )
            ], style={
                'marginBottom': '10px',
                'border': '1px solid #e9ecef',
                'borderRadius': '4px',
                'padding': '2px'
            }))
        
        # Add tree nodes
        if tree_nodes:
            tree_nodes_html = []
            for j, node in enumerate(tree_nodes[:5]):  # Show top 5 for combined analysis
                rank_text = "1st" if node['rank'] == 1 else "2nd" if node['rank'] == 2 else "3rd" if node['rank'] == 3 else f"#{node['rank']}"
                tree_nodes_html.append(
                    html.Div([
                        html.Span(f"{rank_text} ", style={'fontSize': '11px'}),
                        html.Span(f"Node {node['node']}: ", style={'fontSize': '10px', 'fontWeight': 'bold', 'color': '#333'}),
                        html.Span(node['rule'], style={'fontSize': '10px', 'color': '#666'})
                    ], style={
                        'marginBottom': '3px',
                        'padding': '2px 5px',
                        'backgroundColor': '#f8f9fa',
                        'borderRadius': '3px',
                        'borderLeft': f'2px solid {"#d63384" if profile_num == 1 else "#0d6efd"}'
                    })
                )
            
            var_importance_html.extend(tree_nodes_html)
    else:
        # Error for combined analysis
        error_msg = f"Decision tree analysis failed for combined regions. This might be due to insufficient data or identical values."
        var_importance_html.append(html.Div([
            html.Span("", style={'fontSize': '12px', 'marginRight': '5px'}),
            html.Span(error_msg, style={'fontSize': '10px', 'color': '#dc3545', 'fontStyle': 'italic'})
        ], style={'marginBottom': '8px', 'padding': '3px 6px', 'backgroundColor': '#f8d7da', 'borderRadius': '3px', 'border': '1px solid #f5c6cb'}))
    
    return html.Div([
        # Decision Importance section
        *var_importance_html,
        
        # Additional info
        html.Hr(style={'margin': '8px 0'})
    ])

def get_regional_variable_importance_display(df, profile, profile_num, selected_universes=None):
    """Generate variable importance display treating each selected region separately (legacy function)"""
    # Check if dataframe is empty
    if df.empty:
        return html.Div([
            html.P("No analysis data available. Click 'Run Multiverse Analysis' to generate results.", 
                   style={'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '50px'}),
            html.Hr(style={'margin': '15px 0'}),
            html.P(f"Profile: {'Focal' if profile_num == 1 else 'Counterfactual'}", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'}),
            html.P(f"Dataset Size: No data available", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'})
        ])
    
    # Sort data by recidivism probability to match specification curve ordering
    df_sorted = df.sort_values('recidivism_prob')
    
    # Identify regions within selections
    regions = identify_regions(selected_universes or [], df_sorted)
    
    if not regions:
        # No regions, use all data
        analysis_df = df_sorted
        regions = [list(range(len(df_sorted)))]
        region_title = "All Universes"
    else:
        region_title = f"{len(regions)} Region{'s' if len(regions) > 1 else ''}"
    
    # Create variable importance display
    var_importance_html = []
    
    # Add analysis status indicator
    if len(regions) > 1:
        status_indicator = html.Div([
            html.Span("", style={'fontSize': '14px', 'marginRight': '5px'}),
            html.Span(f"Analyzing {len(regions)} separate regions individually", 
                     style={'fontSize': '11px', 'color': '#28a745', 'fontWeight': 'bold'})
        ], style={'marginBottom': '8px', 'padding': '4px 8px', 'backgroundColor': '#d4edda', 'borderRadius': '4px', 'border': '1px solid #c3e6cb'})
        var_importance_html.append(status_indicator)
    
    var_importance_html.append(html.H5(f"Regional Decision Analysis ({region_title})", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '0', 'marginBottom': '6px', 'fontSize': '14px'}))
    var_importance_html.append(html.P("Most impactful methods on recidivism probability for each region:", style={'fontSize': '12px', 'color': '#666', 'marginBottom': '5px'}))
    
    # Analyze each region separately
    for i, region in enumerate(regions):
        if not region:  # Skip empty regions
            continue
            
        # Get data for this specific region
        region_universe_ids = [df_sorted.index[idx] for idx in region if idx < len(df_sorted)]
        region_df = df_sorted.loc[region_universe_ids] if region_universe_ids else df_sorted
        
        if region_df.empty:
            continue
        
        # Get variable importance for this region
        var_importance = get_variable_importance_r(region_df)
        
        # Get tree nodes for this region
        tree_nodes = get_tree_nodes_r(region_df)
        
        # Create region header
        region_header = html.Div([
            html.H6(f"Region {i+1} ({len(region_df)} universes)", 
                   style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '10px', 'marginBottom': '5px', 'fontSize': '13px', 'fontWeight': 'bold'}),
        ])
        var_importance_html.append(region_header)
        
        if var_importance:
            # Create bar plot for this region
            valid_vars = [(name, val) for name, val in var_importance if val is not None and val > 0]
            if valid_vars:
                top_vars_sorted = valid_vars[:5]
                var_names = [var[0] for var in top_vars_sorted]
                var_values = [var[1] for var in top_vars_sorted]
                
                # Create bar plot for this region
                bar_fig = go.Figure(data=[
                    go.Bar(
                        x=var_values,
                        y=var_names,
                        orientation='h',
                        marker_color='#d63384' if profile_num == 1 else '#0d6efd',
                        marker_line_color='#fff',
                        marker_line_width=1,
                        text=[f'{val:.3f}' for val in var_values],
                        textposition='auto',
                        textfont=dict(size=10, color='white')
                    )
                ])
                
                bar_fig.update_layout(
                    yaxis=dict(autorange='reversed')
                )
                
                bar_fig.update_layout(
                    title=None,
                    yaxis_title=None,
                    height=100,
                    margin=dict(l=8, r=8, t=5, b=5),
                    showlegend=False,
                    xaxis=dict(showgrid=False, zeroline=False, title_font=dict(size=9)),
                    yaxis=dict(showgrid=False, zeroline=False),
                    plot_bgcolor='rgba(0,0,0,0)',
                    paper_bgcolor='rgba(0,0,0,0)'
                )
                
                # Wrap the chart in a horizontally scrollable container
                var_importance_html.append(html.Div([
                    dcc.Graph(
                        figure=bar_fig,
                        config={'displayModeBar': False},
                        style={'height': '80px', 'marginBottom': '10px', 'minWidth': '300px'}
                    )
                ], style={
                    'marginBottom': '10px',
                    'border': '1px solid #e9ecef',
                    'borderRadius': '4px',
                    'padding': '2px'
                }))
            
            # Add tree nodes for this region
            if tree_nodes:
                tree_nodes_html = []
                for j, node in enumerate(tree_nodes[:3]):  # Show top 3 for each region
                    rank_text = "1st" if node['rank'] == 1 else "2nd" if node['rank'] == 2 else "3rd" if node['rank'] == 3 else f"#{node['rank']}"
                    tree_nodes_html.append(
                        html.Div([
                            html.Span(f"{rank_text} ", style={'fontSize': '11px'}),
                            html.Span(f"Node {node['node']}: ", style={'fontSize': '10px', 'fontWeight': 'bold', 'color': '#333'}),
                            html.Span(node['rule'], style={'fontSize': '10px', 'color': '#666'})
                        ], style={
                            'marginBottom': '3px',
                            'padding': '2px 5px',
                            'backgroundColor': '#f8f9fa',
                            'borderRadius': '3px',
                            'borderLeft': f'2px solid {"#d63384" if profile_num == 1 else "#0d6efd"}'
                        })
                    )
                
                var_importance_html.extend(tree_nodes_html)
        else:
            # Error for this region
            error_msg = f"Decision tree analysis failed for Region {i+1}. This might be due to insufficient data or identical values."
            var_importance_html.append(html.Div([
                html.Span("", style={'fontSize': '12px', 'marginRight': '5px'}),
                html.Span(error_msg, style={'fontSize': '10px', 'color': '#dc3545', 'fontStyle': 'italic'})
            ], style={'marginBottom': '8px', 'padding': '3px 6px', 'backgroundColor': '#f8d7da', 'borderRadius': '3px', 'border': '1px solid #f5c6cb'}))
        
        # Add separator between regions (except for the last one)
        if i < len(regions) - 1:
            var_importance_html.append(html.Hr(style={'margin': '8px 0', 'borderColor': '#eee'}))
    
    return html.Div([
        # Brief instruction
        html.P("This card shows separate decision tree analyses for each selected region (legacy mode)", 
               style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
        
        # Decision Importance section
        *var_importance_html,
        
        # Additional info
        html.Hr(style={'margin': '8px 0'}),
        html.P(f"Total Regions: {len(regions)}", style={'fontSize': '10px', 'color': '#666', 'fontStyle': 'italic'})
    ])

def get_combined_variable_importance_display(df1, df2, profile1, profile2, selected_universes_1=None, selected_universes_2=None):
    """Generate variable importance display for combined universe selections from focal profile only"""
    # Check if focal dataframe is empty
    if df1.empty:
        return html.Div([
            html.P("No analysis data available. Click 'Run Multiverse Analysis' to generate results.", 
                   style={'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '50px'}),
            html.Hr(style={'margin': '15px 0'}),
            html.P("Profile: Combined Analysis (Focal Only)", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'}),
            html.P("Dataset Size: No data available", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'})
        ])
    
    # Only use focal profile data (df1) - exclude counterfactual profile data
    combined_data = []
    
    # Process focal profile data (df1) only
    if not df1.empty:
        df1_sorted = df1.sort_values('recidivism_prob')
        if selected_universes_1 is not None and len(selected_universes_1) > 0:
            try:
                selected_universe_ids_1 = [df1_sorted.index[i] for i in selected_universes_1 if isinstance(i, int) and i < len(df1_sorted)]
                if selected_universe_ids_1:
                    df1_display = df1_sorted.loc[selected_universe_ids_1]
                    combined_data.extend(df1_display.to_dict('records'))
            except (IndexError, KeyError, TypeError):
                pass
    
    if not combined_data:
        # No selected universes, use all focal data
        if not df1.empty:
            combined_data.extend(df1.to_dict('records'))
    
    if not combined_data:
        return html.Div([
            html.P("No data available for analysis.", 
                   style={'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '50px'})
        ])
    
    # Convert to DataFrame for analysis
    import pandas as pd
    analysis_df = pd.DataFrame(combined_data)
    
    # Get variable importance using R decision tree
    var_importance = get_variable_importance_r(analysis_df)
    
    # Get tree nodes using R analysis
    tree_nodes = get_tree_nodes_r(analysis_df)
    
    # Create variable importance display
    var_importance_html = []
    if var_importance:
        # Add analysis status indicator
        focal_count = len([i for i in (selected_universes_1 or []) if isinstance(i, int) and i < len(df1)]) if not df1.empty else 0
        
        if focal_count > 0:
            status_indicator = html.Div([
                html.Span("", style={'fontSize': '14px', 'marginRight': '5px'}),
                html.Span(f"Analyzing {focal_count} selected universes from focal profile only", 
                         style={'fontSize': '11px', 'color': '#28a745', 'fontWeight': 'bold'})
            ], style={'marginBottom': '8px', 'padding': '4px 8px', 'backgroundColor': '#d4edda', 'borderRadius': '4px', 'border': '1px solid #c3e6cb'})
            var_importance_html.append(status_indicator)
        
        var_importance_html.append(html.H5("Key Decisions (Focal Profile Analysis)", style={'color': '#d63384', 'marginTop': '0', 'marginBottom': '6px', 'fontSize': '14px'}))
        var_importance_html.append(html.P("Most impactful methods on recidivism probability from focal profile regions:", style={'fontSize': '12px', 'color': '#666', 'marginBottom': '5px'}))
        
        # Create bar plot for top 5 variables
        if len(var_importance) >= 1:
            valid_vars = [(name, val) for name, val in var_importance if val is not None and val > 0]
            if valid_vars:
                top_vars_sorted = valid_vars[:5]
                var_names = [var[0] for var in top_vars_sorted]
                var_values = [var[1] for var in top_vars_sorted]
                
                # Create bar plot
                bar_fig = go.Figure(data=[
                    go.Bar(
                        x=var_values,
                        y=var_names,
                        orientation='h',
                        marker_color='#d63384',
                        marker_line_color='#fff',
                        marker_line_width=1,
                        text=[f'{val:.3f}' for val in var_values],
                        textposition='auto',
                        textfont=dict(size=10, color='white')
                    )
                ])
                
                bar_fig.update_layout(
                    yaxis=dict(autorange='reversed')
                )
                
                bar_fig.update_layout(
                    title=None,
                    yaxis_title=None,
                    height=120,
                    margin=dict(l=8, r=8, t=5, b=5),
                    showlegend=False,
                    xaxis=dict(showgrid=False, zeroline=False, title_font=dict(size=9)),
                    yaxis=dict(showgrid=False, zeroline=False),
                    plot_bgcolor='rgba(0,0,0,0)',
                    paper_bgcolor='rgba(0,0,0,0)'
                )
                
                # Wrap the chart in a horizontally scrollable container
                var_importance_html.append(html.Div([
                    dcc.Graph(
                        figure=bar_fig,
                        config={'displayModeBar': False},
                        style={'height': '100px', 'minWidth': '300px'}
                    )
                ], style={
                    'marginBottom': '10px',
                    'border': '1px solid #e9ecef',
                    'borderRadius': '4px',
                    'padding': '2px'
                }))
        
        # Add tree nodes section
        if tree_nodes:
            var_importance_html.append(html.Hr(style={'margin': '6px 0', 'borderColor': '#ddd'}))
            var_importance_html.append(html.H5("Tree Split Rules (Focal)", style={'color': '#d63384', 'marginTop': '6px', 'marginBottom': '5px', 'fontSize': '14px'}))
            var_importance_html.append(html.P("Most important decision splits from focal profile regions:", style={'fontSize': '12px', 'color': '#666', 'marginBottom': '5px'}))
            
            # Create list of tree nodes
            tree_nodes_html = []
            for i, node in enumerate(tree_nodes[:5]):
                rank_text = "1st" if node['rank'] == 1 else "2nd" if node['rank'] == 2 else "3rd" if node['rank'] == 3 else f"#{node['rank']}"
                tree_nodes_html.append(
                    html.Div([
                        html.Span(f"{rank_text} ", style={'fontSize': '12px'}),
                        html.Span(f"Node {node['node']}: ", style={'fontSize': '11px', 'fontWeight': 'bold', 'color': '#333'}),
                        html.Span(node['rule'], style={'fontSize': '11px', 'color': '#666'})
                    ], style={
                        'marginBottom': '4px',
                        'padding': '3px 6px',
                        'backgroundColor': '#f8f9fa',
                        'borderRadius': '4px',
                        'borderLeft': '3px solid #d63384'
                    })
                )
            
            var_importance_html.extend(tree_nodes_html)

    else:
        error_msg = "Decision tree analysis failed for combined universe selection. This might be due to insufficient data or identical values across all selected universes."
        var_importance_html.append(html.Div([
            html.Span("", style={'fontSize': '14px', 'marginRight': '5px'}),
            html.Span(error_msg, style={'fontSize': '11px', 'color': '#dc3545', 'fontStyle': 'italic'})
        ], style={'marginBottom': '8px', 'padding': '4px 8px', 'backgroundColor': '#f8d7da', 'borderRadius': '4px', 'border': '1px solid #f5c6cb'}))
    
    return html.Div([
        # Brief instruction
        
        # Decision Importance section
        *var_importance_html,
        
        # Additional info
        html.Hr(style={'margin': '8px 0'}),
        html.P(f"Focal: {profile1['age']} years, {profile1['gender'].title()}, {profile1['race'].replace('_', ' ').title()}", style={'fontSize': '10px', 'color': '#666', 'fontStyle': 'italic'})
    ])

def get_variable_importance_display(df, profile, profile_num, selected_universes=None):
    """Generate variable importance display for a specific profile"""
    # Check if dataframe is empty (no analysis run yet)
    if df.empty:
        return html.Div([
            html.P("No analysis data available. Click 'Run Multiverse Analysis' to generate results.", 
                   style={'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '50px'}),
            html.Hr(style={'margin': '15px 0'}),
            html.P(f"Dataset Size: No data available", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'})
        ])
    
    # Filter data based on selected universes if provided
    analysis_df = df.copy()
    analysis_title_suffix = ""
    
    if selected_universes is not None and len(selected_universes) > 0:
        # Convert selected universe indices to actual universe IDs
        try:
            # Sort data by recidivism probability to match specification curve ordering
            df_sorted = df.sort_values('recidivism_prob')
            selected_universe_ids = [df_sorted.index[i] for i in selected_universes if isinstance(i, int) and i < len(df_sorted)]
            if selected_universe_ids:
                analysis_df = df_sorted.loc[selected_universe_ids]
                analysis_title_suffix = f" (Selected {len(selected_universe_ids)} universes)"
            else:
                analysis_title_suffix = " (All universes)"
        except (IndexError, KeyError, TypeError):
            # Fallback to showing all universes if there's an error
            analysis_title_suffix = " (All universes)"
    else:
        analysis_title_suffix = " (All universes)"
    
    # Get variable importance using R decision tree
    var_importance = get_variable_importance_r(analysis_df)
    
    # Get tree nodes using R analysis
    tree_nodes = get_tree_nodes_r(analysis_df)
    
    # Debug: print the raw variable importance data
    print(f"Profile {profile_num} - Raw var_importance: {var_importance}")
    print(f"Profile {profile_num} - Tree nodes: {tree_nodes}")
    
    # Create variable importance display
    var_importance_html = []
    if var_importance:
        # Add analysis status indicator
        if selected_universes is not None and len(selected_universes) > 0:
            status_indicator = html.Div([
                html.Span("", style={'fontSize': '14px', 'marginRight': '5px'}),
                html.Span(f"Analyzing {len(selected_universes)} selected universes", 
                            style={'fontSize': '11px', 'color': '#28a745', 'fontWeight': 'bold'})
            ], style={'marginBottom': '8px', 'padding': '4px 8px', 'backgroundColor': '#d4edda', 'borderRadius': '4px', 'border': '1px solid #c3e6cb'})
            var_importance_html.append(status_indicator)
        
        var_importance_html.append(html.H5("Key Decisions (Regression Tree)", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '0', 'marginBottom': '6px', 'fontSize': '14px'}))
        var_importance_html.append(html.P("Most impactful methods on recidivism probability:", style={'fontSize': '12px', 'color': '#666', 'marginBottom': '5px'}))
        
        # Create bar plot for top 5 variables
        if len(var_importance) >= 1:
            # Filter out any None or invalid values and get top 5
            valid_vars = [(name, val) for name, val in var_importance if val is not None and val > 0]
            if valid_vars:
                # The data is already sorted by importance (descending) from the R function
                # Just take the top 5 to ensure proper display order
                top_vars_sorted = valid_vars[:5]
                var_names = [var[0] for var in top_vars_sorted]
                var_values = [var[1] for var in top_vars_sorted]
                
                # Create bar plot - ensure proper ordering from most to least important
                bar_fig = go.Figure(data=[
                    go.Bar(
                        x=var_values,
                        y=var_names,
                        orientation='h',  # Horizontal bars
                        marker_color='#d63384' if profile_num == 1 else '#0d6efd',
                        marker_line_color='#fff',
                        marker_line_width=1,
                        text=[f'{val:.3f}' for val in var_values],
                        textposition='auto',
                        textfont=dict(size=10, color='white')
                    )
                ])
                
                # Ensure y-axis is reversed so most important appears at the top
                bar_fig.update_layout(
                    yaxis=dict(autorange='reversed')
                )
                
                bar_fig.update_layout(
                    title=None,
                    # xaxis_title='Importance Score',
                    yaxis_title=None,
                    height=120,
                    margin=dict(l=8, r=8, t=5, b=5),
                    showlegend=False,
                    xaxis=dict(
                        showgrid=False, 
                        zeroline=False,
                        title_font=dict(size=9)  # Smaller font size
                    ),
                    yaxis=dict(showgrid=False, zeroline=False),
                    plot_bgcolor='rgba(0,0,0,0)',
                    paper_bgcolor='rgba(0,0,0,0)'
                )
                
                # Add the bar plot in a horizontally scrollable container
                var_importance_html.append(html.Div([
                    dcc.Graph(
                        figure=bar_fig,
                        config={'displayModeBar': False},
                        style={'height': '100px', 'minWidth': '300px'}
                    )
                ], style={
                    'marginBottom': '10px',
                    'border': '1px solid #e9ecef',
                    'borderRadius': '4px',
                    'padding': '2px'
                }))
        
        # Add tree nodes section
        if tree_nodes:
            var_importance_html.append(html.Hr(style={'margin': '6px 0', 'borderColor': '#ddd'}))
            var_importance_html.append(html.H5("Tree Split Rules", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '6px', 'marginBottom': '5px', 'fontSize': '14px'}))
            var_importance_html.append(html.P("Most important decision splits in the regression tree:", style={'fontSize': '12px', 'color': '#666', 'marginBottom': '5px'}))
            
            # Create list of tree nodes
            tree_nodes_html = []
            for i, node in enumerate(tree_nodes[:5]):  # Show top 5 nodes
                rank_text = "1st" if node['rank'] == 1 else "2nd" if node['rank'] == 2 else "3rd" if node['rank'] == 3 else f"#{node['rank']}"
                tree_nodes_html.append(
                    html.Div([
                        html.Span(f"{rank_text} ", style={'fontSize': '12px'}),
                        html.Span(f"Node {node['node']}: ", style={'fontSize': '11px', 'fontWeight': 'bold', 'color': '#333'}),
                        html.Span(node['rule'], style={'fontSize': '11px', 'color': '#666'})
                    ], style={
                        'marginBottom': '4px',
                        'padding': '3px 6px',
                        'backgroundColor': '#f8f9fa',
                        'borderRadius': '4px',
                        'borderLeft': f'3px solid {"#d63384" if profile_num == 1 else "#0d6efd"}'
                    })
                )
            
            var_importance_html.extend(tree_nodes_html)

    else:
        # Provide more helpful error message
        if selected_universes is not None and len(selected_universes) > 0:
            error_msg = f"Decision tree analysis failed for {len(selected_universes)} selected universes. This might be due to insufficient data or identical values across all selected universes."
        else:
            error_msg = "Variable importance analysis not available. Please ensure the multiverse analysis has been completed."
        
        var_importance_html.append(html.Div([
            html.Span("", style={'fontSize': '14px', 'marginRight': '5px'}),
            html.Span(error_msg, style={'fontSize': '11px', 'color': '#dc3545', 'fontStyle': 'italic'})
        ], style={'marginBottom': '8px', 'padding': '4px 8px', 'backgroundColor': '#f8d7da', 'borderRadius': '4px', 'border': '1px solid #f5c6cb'}))
    
    return html.Div([
        # Brief instruction in light gray
        
        # Decision Importance section
        *var_importance_html,
        
        # Additional info
        html.Hr(style={'margin': '8px 0'})
    ])

def get_dataset_overview_content(df_nc, df_low_risk):
    """Generate dataset overview content for the sidebar"""
    # Use either dataset since they have the same method combinations
    df = df_nc if not df_nc.empty else df_low_risk
    
    if df.empty:
        return html.Div([
            html.P("No analysis data available. Click 'Run Multiverse Analysis' to generate results.", 
                   style={'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '20px'})
        ])
    
    # Calculate total specifications across both profiles
    total_specs = len(df_nc) + len(df_low_risk) if not df_nc.empty and not df_low_risk.empty else len(df)
    
    return html.Div([
        html.P("Number of choices for each decision", style={'fontSize': '12px', 'color': '#333', 'marginBottom': '8px', 'fontWeight': 'bold', 'textAlign': 'center'}),
        html.P(f"Total Universes: {total_specs}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Scaling: {df['scaling'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Convert Followup: {df['convert_followup'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Split Ratios: {df['split'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Age Categories: {df['age_category'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Imbalancing: {df['imbalancing'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Model: {df['model'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Predictors: {df['predictor'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Recidivism Definitions: {df['define_recid'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.Hr(style={'margin': '8px 0', 'borderColor': '#ddd'}),
        html.P("Note: Both profiles use the same method combinations with different demographic parameters.", 
               style={'fontSize': '10px', 'color': '#666', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '5px'})
    ], style={'padding': '8px', 'backgroundColor': '#f8f9fa', 'borderRadius': '4px', 'border': '1px solid #dee2e6'})


# Selection handling is now integrated into the main callback



# Global variables to store accumulated selections
accumulated_selections_1 = []  # Store multiple regions for focal profile
accumulated_selections_2 = []  # Store multiple regions for counterfactual profile

# Combined callback to handle selection events
@app.callback(
    Output('highlighted-universe-1', 'children'),
    Output('highlighted-universe-2', 'children'),
    Input('spec-curve-1', 'clickData'),
    Input('spec-curve-2', 'clickData'),
    Input('combined-spec-grid', 'clickData'),
    Input('spec-curve-1', 'selectedData'),
    Input('spec-curve-2', 'selectedData'),
    State('highlighted-universe-1', 'children'),
    State('highlighted-universe-2', 'children'),
    prevent_initial_call=True
)
def update_highlighted_universe(click_data_1, click_data_2, click_data_grid, selected_data_1, selected_data_2, current_highlighted_1, current_highlighted_2):
    """Handle click events to highlight a single universe column"""
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    try:
        if trigger_id == 'spec-curve-1':
            if click_data_1 and 'points' in click_data_1 and click_data_1['points']:
                # Get the clicked point index
                clicked_index = click_data_1['points'][0]['pointIndex']
                # Clear counterfactual highlighting when focal is clicked
                return clicked_index, None
            else:
                return None, None
        
        elif trigger_id == 'spec-curve-2':
            # Disable counterfactual profile click functionality
            # Always clear counterfactual highlighting
            return dash.no_update, None
        
        elif trigger_id == 'combined-spec-grid':
            if click_data_grid and 'points' in click_data_grid and click_data_grid['points']:
                # For heatmaps, the click data structure is different
                point_data = click_data_grid['points'][0]
                
                if 'pointIndex' in point_data:
                    # Scatter plot style click data
                    clicked_index = point_data['pointIndex']
                elif 'x' in point_data:
                    # Heatmap style click data - x coordinate gives us the column index
                    clicked_index = int(point_data['x'])
                else:
                    print(f"Unexpected click data structure: {point_data}")
                    return dash.no_update, dash.no_update
                
                # For combined grid, only highlight the focal profile (profile 1)
                # Clear counterfactual highlighting when grid is clicked
                return clicked_index, None
            else:
                return dash.no_update, dash.no_update
        
        # Check if this is a selection event (not a click event)
        if ctx.triggered[0]['prop_id'].endswith('.selectedData'):
            print(f"Selection event detected: {ctx.triggered[0]['prop_id']}")
            # Clear both highlighting when any region is selected via dragging
            return None, None
        
    except (KeyError, IndexError, TypeError, ValueError) as e:
        print(f"Error handling click event: {e}")
        return dash.no_update, dash.no_update
    
    return dash.no_update, dash.no_update

@app.callback(
    Output('selected-universes-1', 'children'),
    Output('selected-universes-2', 'children'),
    Input('spec-curve-1', 'selectedData'),
    Input('spec-curve-2', 'selectedData'),
    State('selected-universes-1', 'children'),
    State('selected-universes-2', 'children'),
    prevent_initial_call=True
)
def update_selected_universes(selected_data_1, selected_data_2, current_selected_1, current_selected_2):
    global accumulated_selections_1, accumulated_selections_2
    
    # Get the context to determine which input triggered the callback
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if trigger_id == 'spec-curve-1':
        # Handle selection from specification curve 1 (Focal Profile)
        if selected_data_1 and 'points' in selected_data_1:
            selected_points = selected_data_1['points']
            new_selection = [point['pointIndex'] for point in selected_points if 'pointIndex' in point]
            
            # Add to accumulated selections (multiple regions)
            if new_selection:
                accumulated_selections_1.extend(new_selection)
                # Remove duplicates while preserving order
                accumulated_selections_1 = list(dict.fromkeys(accumulated_selections_1))
        else:
            # If no selection, clear accumulated selections for this profile
            accumulated_selections_1 = []
        
        return accumulated_selections_1, dash.no_update
    
    elif trigger_id == 'spec-curve-2':
        # Handle selection from specification curve 2 (Counterfactual Profile)
        if selected_data_2 and 'points' in selected_data_2:
            selected_points = selected_data_2['points']
            new_selection = [point['pointIndex'] for point in selected_points if 'pointIndex' in point]
            
            # Add to accumulated selections (multiple regions)
            if new_selection:
                accumulated_selections_2.extend(new_selection)
                # Remove duplicates while preserving order
                accumulated_selections_2 = list(dict.fromkeys(accumulated_selections_2))
        else:
            # If no selection, clear accumulated selections for this profile
            accumulated_selections_2 = []
        
        return dash.no_update, accumulated_selections_2
    
    return dash.no_update, dash.no_update

# Callback for clear selection buttons
@app.callback(
    Output('selected-universes-1', 'children', allow_duplicate=True),
    Output('selected-universes-2', 'children', allow_duplicate=True),
    Output('highlighted-universe-1', 'children', allow_duplicate=True),
    Output('highlighted-universe-2', 'children', allow_duplicate=True),
    Input('clear-selections-button', 'n_clicks'),
    Input('clear-focal-btn', 'n_clicks'),
    Input('clear-cf-btn', 'n_clicks'),
    prevent_initial_call=True
)
def clear_selections(clear_all_clicks, clear_focal_clicks, clear_cf_clicks):
    global accumulated_selections_1, accumulated_selections_2
    
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if trigger_id == 'clear-selections-button':
        # Clear all selections and highlighting
        accumulated_selections_1 = []
        accumulated_selections_2 = []
        return [], [], None, None
    elif trigger_id == 'clear-focal-btn':
        # Clear only focal selections and highlighting
        accumulated_selections_1 = []
        return [], dash.no_update, None, dash.no_update
    elif trigger_id == 'clear-cf-btn':
        # Clear only counterfactual selections and highlighting
        accumulated_selections_2 = []
        return dash.no_update, [], dash.no_update, None
    
    return dash.no_update, dash.no_update, dash.no_update, dash.no_update




# Global variables to store progress updates
progress_updates = {"focal": [], "counterfactual": []}
progress_counts = {"focal": 0, "counterfactual": 0}
total_universes = 9600

# Global variable to track if focal profile has been run
focal_profile_has_been_run = False

# Global variable to track the last focal profile parameters that were run
last_focal_profile_params = None

# Global variables to track selected universe analysis state
selected_universe_analysis_active = False
selected_universe_analysis_profile = None  # 1 for focal, 2 for counterfactual
selected_universe_analysis_data = None  # Store the selected universes for analysis

# Callback to show submit button feedback
@app.callback(
    Output('submit-profiles-button', 'children'),
    Output('submit-profiles-button', 'style'),
    Input('submit-profiles-button', 'n_clicks'),
    Input('progress-interval', 'n_intervals')
)
def update_submit_button(n_clicks, n_intervals):
    # Always show "Run Multiverse Analysis" button
    button_text = "Run Multiverse Analysis"
    button_style = {
        'fontSize': '16px',
        'padding': '12px 24px',
        'backgroundColor': '#6c757d',
        'color': 'white',
        'border': 'none',
        'borderRadius': '8px',
        'cursor': 'pointer',
        'fontWeight': 'bold',
        'width': '100%',
        'marginTop': '10px'
    }
    
    return button_text, button_style

# Callback to update progress bars
@app.callback(
    Output('focal-progress-bar', 'style'),
    Output('focal-progress-text', 'children'),
    Output('counterfactual-progress-bar', 'style'),
    Output('counterfactual-progress-text', 'children'),
    Input('progress-interval', 'n_intervals')
)
def update_progress_bars(n_intervals):
    global progress_counts, total_universes
    
    # Calculate progress percentages
    focal_progress = (progress_counts["focal"] / total_universes) * 100
    counterfactual_progress = (progress_counts["counterfactual"] / total_universes) * 100
    
    # Check if both profiles are running together (same progress)
    profiles_running_together = (progress_counts["focal"] == progress_counts["counterfactual"] and 
                                progress_counts["focal"] > 0)
    
    # Check if focal profile is completed and counterfactual is running
    focal_completed_counterfactual_running = (progress_counts["focal"] == total_universes and 
                                            progress_counts["counterfactual"] > 0 and 
                                            progress_counts["counterfactual"] < total_universes)
    
    # Update focal progress bar
    focal_bar_style = {
        'width': f'{focal_progress}%',
        'height': '8px',
        'backgroundColor': '#d63384',
        'borderRadius': '4px',
        'transition': 'width 0.3s ease'
    }
    
    # Update counterfactual progress bar
    counterfactual_bar_style = {
        'width': f'{counterfactual_progress}%',
        'height': '8px',
        'backgroundColor': '#0d6efd',
        'borderRadius': '4px',
        'transition': 'width 0.3s ease'
    }
    
    # Update text based on analysis status
    if profiles_running_together:
        focal_text = f"Running together: Universe {progress_counts['focal']}/{total_universes}"
        counterfactual_text = f"Running together: Universe {progress_counts['counterfactual']}/{total_universes}"
    elif focal_completed_counterfactual_running:
        focal_text = f"Completed (reused): Universe {progress_counts['focal']}/{total_universes}"
        counterfactual_text = f"Running: Universe {progress_counts['counterfactual']}/{total_universes}"
    else:
        focal_text = f"Universe {progress_counts['focal']}/{total_universes}"
        counterfactual_text = f"Universe {progress_counts['counterfactual']}/{total_universes}"
    
    return focal_bar_style, focal_text, counterfactual_bar_style, counterfactual_text

# Callback to update all outputs - now triggered by submit button
@app.callback(
    Output('spec-curve-1', 'figure'),
    Output('spec-curve-2', 'figure'),
    Output('combined-spec-grid', 'figure'),
    Output('variable-importance-content-1', 'children'),
    Output('variable-importance-content-2', 'children'),
    Output('dataset-overview-content', 'children'),
    Output('selection-status-indicator', 'children'),
    Input('submit-profiles-button', 'n_clicks'),
    Input('selected-universes-1', 'children'),
    Input('selected-universes-2', 'children'),
    Input('highlighted-universe-1', 'children'),
    Input('highlighted-universe-2', 'children'),
    State('age-slider-1', 'value'),
    State('gender-dropdown-1', 'value'),
    State('race-dropdown-1', 'value'),
    State('age-slider-2', 'value'),
    State('gender-dropdown-2', 'value'),
    State('race-dropdown-2', 'value')
)
def update_dashboard(submit_clicks, selected_1, selected_2, highlighted_1, highlighted_2, age_1, gender_1, race_1, age_2, gender_2, race_2):
    global current_grid_profile, df_nc, df_low_risk
    global selected_universe_analysis_active, selected_universe_analysis_profile, selected_universe_analysis_data
    global focal_profile_has_been_run, last_focal_profile_params
    
    try:
        # Handle None values for selected universes
        if selected_1 is None:
            selected_1 = []
        if selected_2 is None:
            selected_2 = []
        
        # Debug: print the selected universes
        print(f"Selected universes 1: {selected_1}")
        print(f"Selected universes 2: {selected_2}")
        
        # Only run analysis if submit button was clicked or if it's the initial load
        ctx = dash.callback_context
        if ctx.triggered:
            trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
            print(f"Callback triggered by: {trigger_id}")
            # Only proceed if submit button was clicked, or if it's selection events, highlighting events, or profile switcher
            if trigger_id not in ['submit-profiles-button', 'selected-universes-1', 'selected-universes-2', 'highlighted-universe-1', 'highlighted-universe-2']:
                # Return current state without updating if other inputs changed
                print(f"Preventing update for trigger: {trigger_id}")
                raise dash.exceptions.PreventUpdate
            
            # Handle highlighting updates without re-running analysis
            if trigger_id in ['highlighted-universe-1', 'highlighted-universe-2']:
                print(f"Highlighting update triggered by: {trigger_id}")
                # Skip the analysis and just update the plots with new highlighting
                # We'll fall through to the plot creation section below
            # If submit button was clicked, run the multiverse analysis
            elif trigger_id == 'submit-profiles-button':
                print("Submit button clicked - checking profile characteristics...")
                print(f"Focal profile params: age={age_1}, gender={gender_1}, race={race_1}")
                print(f"Counterfactual profile params: age={age_2}, gender={gender_2}, race={race_2}")
                print(f"Focal profile has been run: {focal_profile_has_been_run}")
                print(f"Last focal profile params: {last_focal_profile_params}")
                print(f"Current df_nc empty: {df_nc.empty}")
                
                # Check if profiles are identical
                profiles_identical = profiles_are_identical(age_1, gender_1, race_1, age_2, gender_2, race_2)
                
                if profiles_identical:
                    print("Profiles are identical - running shared analysis...")
                    # Run analysis once and use results for both profiles
                    shared_results = run_multiverse_analysis_with_profile(age_1, gender_1, race_1, is_shared_analysis=True)
                    
                    # Store the same results for both profiles
                    if shared_results is not None:
                        df_nc = shared_results
                        df_low_risk = shared_results
                        focal_profile_has_been_run = True  # Mark focal as run
                        last_focal_profile_params = (age_1, gender_1, race_1)  # Store parameters
                        print("Shared analysis completed - both profiles use same results!")
                    else:
                        print("Shared analysis failed!")
                else:
                    print("Profiles are different - checking if focal profile needs to be run...")
                    
                    # Check if focal profile has been run with the same parameters
                    current_focal_params = (age_1, gender_1, race_1)
                    
                    if (focal_profile_has_been_run and not df_nc.empty and 
                        last_focal_profile_params is not None and
                        current_focal_params == last_focal_profile_params):
                        print(f"Focal profile already run with same params {current_focal_params} - reusing results!")
                        print("Running only counterfactual analysis...")
                        # Set focal progress to completed since we're reusing results
                        progress_counts["focal"] = total_universes
                        progress_updates["focal"] = []  # Clear previous progress
                        progress_updates["focal"].append("Focal profile results reused from previous analysis")
                        progress_updates["focal"].append(f"Universe {total_universes}/{total_universes}")
                        
                        # Reset counterfactual progress before running
                        progress_counts["counterfactual"] = 0
                        progress_updates["counterfactual"] = []
                        
                        # Reuse existing focal results, only run counterfactual
                        # Pass is_shared_analysis=False explicitly to ensure only counterfactual updates
                        counterfactual_results = run_multiverse_analysis_with_profile(age_2, gender_2, race_2, is_shared_analysis=False)
                        
                        if counterfactual_results is not None:
                            df_low_risk = counterfactual_results
                            print("Counterfactual analysis completed - focal results reused!")
                        else:
                            print("Counterfactual analysis failed!")
                    else:
                        print(f"Focal profile needs to be run (has_been_run={focal_profile_has_been_run}, df_empty={df_nc.empty}, params_match={current_focal_params == last_focal_profile_params})")
                        print("Running sequential multiverse analysis...")
                        # Run analyses sequentially with proper progress tracking
                        focal_results = run_multiverse_analysis_with_profile(age_1, gender_1, race_1)
                        counterfactual_results = run_multiverse_analysis_with_profile(age_2, gender_2, race_2)
                        print("Sequential multiverse analysis completed!")
                        
                        # Store results globally for use in charts
                        if focal_results is not None:
                            df_nc = focal_results
                            # Mark focal as run and store parameters
                            focal_profile_has_been_run = True
                            last_focal_profile_params = current_focal_params
                        if counterfactual_results is not None:
                            df_low_risk = counterfactual_results
                
                # Reset selected universe analysis state when new analysis is run
                selected_universe_analysis_active = False
                selected_universe_analysis_profile = None
                selected_universe_analysis_data = None
            elif trigger_id in ['selected-universes-1', 'selected-universes-2']:
                # Handle selection changes - automatically activate analysis for new selections
                current_selection = selected_1 if trigger_id == 'selected-universes-1' else selected_2
                profile_num = 1 if trigger_id == 'selected-universes-1' else 2
                
                # Normalize selections for comparison (handle None vs empty list)
                current_selection_normalized = current_selection if current_selection is not None else []
                stored_selection_normalized = selected_universe_analysis_data if selected_universe_analysis_data is not None else []
                
                # If user has selected universes, automatically activate analysis
                if len(current_selection_normalized) > 0:
                    # Check if this is a new selection or change in selection
                    if (current_selection_normalized != stored_selection_normalized or 
                        not selected_universe_analysis_active or 
                        selected_universe_analysis_profile != profile_num):
                        
                        # Activate analysis for the selected universes
                        selected_universe_analysis_active = True
                        selected_universe_analysis_profile = profile_num
                        selected_universe_analysis_data = current_selection_normalized
                        print(f"Auto-activating decision tree analysis for profile {profile_num} with {len(current_selection_normalized)} selected universes")
                    else:
                        print(f"Selection unchanged, preserving analysis state (profile {profile_num}, {len(current_selection_normalized)} universes)")
                else:
                    # No universes selected, reset to default analysis
                    if selected_universe_analysis_active:
                        selected_universe_analysis_active = False
                        selected_universe_analysis_profile = None
                        selected_universe_analysis_data = None
                        print("No universes selected, resetting to default analysis")
                    else:
                        print("No universes selected, using default analysis")
        
        # Create profile dictionaries with demographic values
        profile1 = {
            'age': age_1,
            'gender': gender_1,
            'race': race_1
        }
        
        profile2 = {
            'age': age_2,
            'gender': gender_2,
            'race': race_2
        }
        
        # Create all the figures and content for both profiles
        # Profile 1 uses nc_multiverse.csv, Profile 2 uses low_risk_multiverse.csv
        # Always recreate curves to show previously selected regions
        # This ensures users can see their accumulated selections
        
        # Define selection state variables
        has_focal_selection = selected_1 is not None and len(selected_1) > 0
        has_cf_selection = selected_2 is not None and len(selected_2) > 0
        
        
        # Check if we have data to create plots
        if df_nc.empty and df_low_risk.empty:
            print("No data available for plotting - analysis may not have been run yet")
            # Create empty figures with message
            spec_curve_1 = create_specification_curve(pd.DataFrame(), profile1, 1, selected_1, highlighted_1)
            spec_curve_2 = create_specification_curve(pd.DataFrame(), profile2, 2, selected_2, None)
            combined_spec_grid = create_combined_specification_grid(pd.DataFrame(), pd.DataFrame(), selected_1, selected_2, highlighted_1, None)
        else:
            # Only show highlighting on the curve if that profile is the active one
            # Determine which profile is active based on highlighting or selections
            # Hide highlighting when there are selections (during drag operations)
            focal_is_active = highlighted_1 is not None and not has_focal_selection
            # Counterfactual highlighting is completely disabled
            cf_is_active = False
            
            # Show highlighting only for the active profile
            focal_highlighting = highlighted_1 if focal_is_active else None
            # Always disable counterfactual highlighting
            cf_highlighting = None
            
            spec_curve_1 = create_specification_curve(df_nc, profile1, 1, selected_1, focal_highlighting)
            spec_curve_2 = create_specification_curve(df_low_risk, profile2, 2, selected_2, None)
            
            # Create the combined specification grid
            if has_focal_selection and has_cf_selection:
                # Both profiles have selections - show combined view (no highlighting during drag)
                combined_spec_grid = create_combined_specification_grid(df_nc, df_low_risk, selected_1, selected_2, None, None)
            elif has_focal_selection:
                # Only focal profile has selection (no highlighting during drag)
                combined_spec_grid = create_specification_grid(df_nc, 1, selected_1, None)
            elif has_cf_selection:
                # Only counterfactual profile has selection
                combined_spec_grid = create_specification_grid(df_low_risk, 2, selected_2, None)
            elif highlighted_1 is not None:
                # Focal profile is highlighted but no selections - show all focal universes with highlighting
                combined_spec_grid = create_specification_grid(df_nc, 1, None, focal_highlighting)
            elif highlighted_2 is not None:
                # Counterfactual profile highlighting disabled - show all counterfactual universes without highlighting
                combined_spec_grid = create_specification_grid(df_low_risk, 2, None, None)
            else:
                # Default to focal profile if neither is explicitly selected
                selected_universes = selected_1 if selected_1 is not None else []
                combined_spec_grid = create_specification_grid(df_nc, 1, selected_universes, focal_highlighting)
        
        # Set variable importance content based on selection state
        if has_focal_selection and has_cf_selection:
            # Both profiles have selections - show focal analysis in first card and counterfactual in second card
            variable_importance_title_1 = "Focal Decision Importance"
            variable_importance_title_style_1 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
            # Wrap content with STEP 4 elements
            base_content_1 = get_combined_variable_importance_display(df_nc, df_low_risk, profile1, profile2, selected_1, selected_2)
            variable_importance_content_1 = html.Div([
                html.Div("STEP 4", style={
                    'fontSize': '12px', 
                    'fontWeight': 'bold', 
                    'color': 'white', 
                    'backgroundColor': '#6c757d', 
                    'padding': '4px 8px', 
                    'borderRadius': '4px',
                    'textAlign': 'center',
                    'marginBottom': '12px'
                }),
                html.P("View focal decision importance analysis for your selections. This shows which methodological choices most impact the focal profile results.", 
                       style={'fontSize': '11px', 'color': '#333', 'margin': '0 0 15px 0', 'lineHeight': '1.4', 'fontStyle': 'italic'}),
                html.H4("Focal Decision Importance", style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px', 'marginTop': '0px'}),
                base_content_1
            ])
            
            # Show counterfactual analysis in the second card
            variable_importance_title_2 = "Counterfactual Decision Importance"
            variable_importance_title_style_2 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px'}
            # Wrap content with STEP 4 elements
            base_content_2 = get_variable_importance_display(df_low_risk, profile2, 2, selected_2)
            variable_importance_content_2 = html.Div([
                html.Div("STEP 4", style={
                    'fontSize': '12px', 
                    'fontWeight': 'bold', 
                    'color': 'white', 
                    'backgroundColor': '#6c757d', 
                    'padding': '4px 8px', 
                    'borderRadius': '4px',
                    'textAlign': 'center',
                    'marginBottom': '12px'
                }),
                html.P("View counterfactual decision importance analysis for your selections. This shows which methodological choices most impact the counterfactual profile results.", 
                       style={'fontSize': '11px', 'color': '#333', 'margin': '0 0 15px 0', 'lineHeight': '1.4', 'fontStyle': 'italic'}),
                html.H4("CF Decision Importance", style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px', 'marginTop': '0px'}),
                base_content_2
            ])
        elif has_focal_selection:
            # Only focal profile has selection - use combined regional analysis
            variable_importance_title_1 = "Focal Decision Importance"
            variable_importance_title_style_1 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
            # Wrap content with STEP 4 elements
            base_content_1 = get_combined_regional_variable_importance_display(df_nc, profile1, 1, selected_1)
            variable_importance_content_1 = html.Div([
                html.Div("STEP 4", style={
                    'fontSize': '12px', 
                    'fontWeight': 'bold', 
                    'color': 'white', 
                    'backgroundColor': '#6c757d', 
                    'padding': '4px 8px', 
                    'borderRadius': '4px',
                    'textAlign': 'center',
                    'marginBottom': '12px'
                }),
                html.P("View focal decision importance analysis for your selections. This shows which methodological choices most impact the focal profile results.", 
                       style={'fontSize': '11px', 'color': '#333', 'margin': '0 0 15px 0', 'lineHeight': '1.4', 'fontStyle': 'italic'}),
                html.H4("Focal Decision Importance", style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px', 'marginTop': '0px'}),
                base_content_1
            ])
            
            variable_importance_title_2 = "CF Decision Importance"
            variable_importance_title_style_2 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px'}
            # Wrap content with STEP 4 elements
            base_content_2 = get_variable_importance_display(df_low_risk, profile2, 2, selected_2)
            variable_importance_content_2 = html.Div([
                html.Div("STEP 4", style={
                    'fontSize': '12px', 
                    'fontWeight': 'bold', 
                    'color': 'white', 
                    'backgroundColor': '#6c757d', 
                    'padding': '4px 8px', 
                    'borderRadius': '4px',
                    'textAlign': 'center',
                    'marginBottom': '12px'
                }),
                html.P("View counterfactual decision importance analysis for your selections. This shows which methodological choices most impact the counterfactual profile results.", 
                       style={'fontSize': '11px', 'color': '#333', 'margin': '0 0 15px 0', 'lineHeight': '1.4', 'fontStyle': 'italic'}),
                html.H4("CF Decision Importance", style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px', 'marginTop': '0px'}),
                base_content_2
            ])
        elif has_cf_selection:
            # Only counterfactual profile has selection - use combined regional analysis
            variable_importance_title_1 = "Focal Decision Importance"
            variable_importance_title_style_1 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
            # Wrap content with STEP 4 elements
            base_content_1 = get_variable_importance_display(df_nc, profile1, 1, selected_1)
            variable_importance_content_1 = html.Div([
                html.Div("STEP 4", style={
                    'fontSize': '12px', 
                    'fontWeight': 'bold', 
                    'color': 'white', 
                    'backgroundColor': '#6c757d', 
                    'padding': '4px 8px', 
                    'borderRadius': '4px',
                    'textAlign': 'center',
                    'marginBottom': '12px'
                }),
                html.P("View focal decision importance analysis for your selections. This shows which methodological choices most impact the focal profile results.", 
                       style={'fontSize': '11px', 'color': '#333', 'margin': '0 0 15px 0', 'lineHeight': '1.4', 'fontStyle': 'italic'}),
                html.H4("Focal Decision Importance", style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px', 'marginTop': '0px'}),
                base_content_1
            ])
            
            variable_importance_title_2 = "CF Decision Importance (Combined)"
            variable_importance_title_style_2 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px'}
            # Wrap content with STEP 4 elements
            base_content_2 = get_combined_regional_variable_importance_display(df_low_risk, profile2, 2, selected_2)
            variable_importance_content_2 = html.Div([
                html.Div("STEP 4", style={
                    'fontSize': '12px', 
                    'fontWeight': 'bold', 
                    'color': 'white', 
                    'backgroundColor': '#6c757d', 
                    'padding': '4px 8px', 
                    'borderRadius': '4px',
                    'textAlign': 'center',
                    'marginBottom': '12px'
                }),
                html.P("View counterfactual decision importance analysis for your selections. This shows which methodological choices most impact the counterfactual profile results.", 
                       style={'fontSize': '11px', 'color': '#333', 'margin': '0 0 15px 0', 'lineHeight': '1.4', 'fontStyle': 'italic'}),
                html.H4("CF Decision Importance", style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px', 'marginTop': '0px'}),
                base_content_2
            ])
        else:
            # Show individual profile analyses
            variable_importance_title_1 = "Focal Decision Importance"
            variable_importance_title_style_1 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
            # Wrap content with STEP 4 elements
            base_content_1 = get_variable_importance_display(df_nc, profile1, 1, selected_1)
            variable_importance_content_1 = html.Div([
                html.Div("STEP 4", style={
                    'fontSize': '12px', 
                    'fontWeight': 'bold', 
                    'color': 'white', 
                    'backgroundColor': '#6c757d', 
                    'padding': '4px 8px', 
                    'borderRadius': '4px',
                    'textAlign': 'center',
                    'marginBottom': '12px'
                }),
                html.P("View focal decision importance analysis for your selections. This shows which methodological choices most impact the focal profile results.", 
                       style={'fontSize': '11px', 'color': '#333', 'margin': '0 0 15px 0', 'lineHeight': '1.4', 'fontStyle': 'italic'}),
                html.H4("Focal Decision Importance", style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px', 'marginTop': '0px'}),
                base_content_1
            ])
            
            variable_importance_title_2 = "CF Decision Importance"
            variable_importance_title_style_2 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px'}
            # Wrap content with STEP 4 elements
            base_content_2 = get_variable_importance_display(df_low_risk, profile2, 2, selected_2)
            variable_importance_content_2 = html.Div([
                html.Div("STEP 4", style={
                    'fontSize': '12px', 
                    'fontWeight': 'bold', 
                    'color': 'white', 
                    'backgroundColor': '#6c757d', 
                    'padding': '4px 8px', 
                    'borderRadius': '4px',
                    'textAlign': 'center',
                    'marginBottom': '12px'
                }),
                html.P("View counterfactual decision importance analysis for your selections. This shows which methodological choices most impact the counterfactual profile results.", 
                       style={'fontSize': '11px', 'color': '#333', 'margin': '0 0 15px 0', 'lineHeight': '1.4', 'fontStyle': 'italic'}),
                html.H4("CF Decision Importance", style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px', 'marginTop': '0px'}),
                base_content_2
            ])
        
        
        
        # Generate dataset overview content
        dataset_overview_content = get_dataset_overview_content(df_nc, df_low_risk)
        
        # Generate selection status indicator
        has_focal_selection = selected_1 is not None and len(selected_1) > 0
        has_cf_selection = selected_2 is not None and len(selected_2) > 0
        
        if has_focal_selection and has_cf_selection:
            focal_count = len(selected_1)
            cf_count = len(selected_2)
            
            # Check for multiple regions in focal profile
            if not df_nc.empty:
                df_sorted = df_nc.sort_values('recidivism_prob')
                focal_regions = identify_regions(selected_1, df_sorted)
                region_info = f" ({len(focal_regions)} region{'s' if len(focal_regions) > 1 else ''})" if len(focal_regions) > 1 else ""
            else:
                region_info = ""
            
            # Get region count for focal profile
            if not df_nc.empty:
                df_sorted = df_nc.sort_values('recidivism_prob')
                focal_regions = identify_regions(selected_1, df_sorted)
                focal_region_count = len(focal_regions)
            else:
                focal_region_count = 1
            
            # Get region count for counterfactual profile
            if not df_low_risk.empty:
                df_sorted_cf = df_low_risk.sort_values('recidivism_prob')
                cf_regions = identify_regions(selected_2, df_sorted_cf)
                cf_region_count = len(cf_regions)
            else:
                cf_region_count = 1
            
            selection_status = html.Div([
                html.Span("", style={'fontSize': '14px'}),
                html.Span(f"Combined Selection Active: {focal_count} Focal ({focal_region_count} region{'s' if focal_region_count > 1 else ''}) + {cf_count} CF ({cf_region_count} region{'s' if cf_region_count > 1 else ''})", 
                         style={'fontWeight': 'bold', 'color': '#28a745'})
            ])
        elif has_focal_selection:
            focal_count = len(selected_1)
            
            # Check for multiple regions in focal profile
            if not df_nc.empty:
                df_sorted = df_nc.sort_values('recidivism_prob')
                focal_regions = identify_regions(selected_1, df_sorted)
                region_info = f" ({len(focal_regions)} region{'s' if len(focal_regions) > 1 else ''})" if len(focal_regions) > 1 else ""
            else:
                region_info = ""
            
            # Get region count for more detailed info
            if not df_nc.empty:
                df_sorted = df_nc.sort_values('recidivism_prob')
                focal_regions = identify_regions(selected_1, df_sorted)
                region_count = len(focal_regions)
            else:
                region_count = 1
            
            selection_status = html.Div([
                html.Span("", style={'fontSize': '14px'}),
                html.Span(f"Focal Selection Active: {focal_count} universes in {region_count} region{'s' if region_count > 1 else ''}", 
                         style={'fontWeight': 'bold', 'color': '#d63384'})
            ])
        elif has_cf_selection:
            cf_count = len(selected_2)
            
            # Check for multiple regions in counterfactual profile
            if not df_low_risk.empty:
                df_sorted = df_low_risk.sort_values('recidivism_prob')
                cf_regions = identify_regions(selected_2, df_sorted)
                cf_region_count = len(cf_regions)
            else:
                cf_region_count = 1
            
            selection_status = html.Div([
                html.Span("", style={'fontSize': '14px'}),
                html.Span(f"Counterfactual Selection Active: {cf_count} universes in {cf_region_count} region{'s' if cf_region_count > 1 else ''}", 
                         style={'fontWeight': 'bold', 'color': '#0d6efd'})
            ])
        else:
            selection_status = html.Div([
                html.Span("", style={'fontSize': '14px'}),
                html.Span("No selections - showing all universes", 
                         style={'color': '#6c757d'})
            ])
        
        return (spec_curve_1, spec_curve_2, combined_spec_grid, 
                variable_importance_content_1, variable_importance_content_2,
                dataset_overview_content, selection_status)
    
    except Exception as e:
        print(f"Error in update_dashboard: {e}")
        import traceback
        traceback.print_exc()
        # Return empty/default values to prevent the app from crashing
        empty_fig = go.Figure()
        empty_fig.add_annotation(text="Error loading data", xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        
        return (empty_fig, empty_fig, empty_fig, 
                html.Div("Error"), html.Div("Error"), 
                "Error", {}, html.Div("Error"),
                "Error", {}, html.Div("Error"),
                html.Div("Error"), html.Div("Error"))



app.layout = create_layout()

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8050)
