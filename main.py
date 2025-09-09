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
# Function to create decision tree and get variable importance
create_decision_tree <- function(data) {
    # Create the decision tree
    tree <- rpart(recidivism_prob ~ preprocessing + split + age_category + 
                  imbalancing_method + predictor_method + define_recid_method,
                  data = data,
                  method = "anova",
                  control = rpart.control(minsplit = 2, cp = 0))
    
    # Get variable importance
    var_importance <- tree$variable.importance
    
    # Get variable importance
    
    # Return both tree and variable importance
    return(list(tree = tree, var_importance = var_importance))
}

# Function to get variable importance for a dataset
get_variable_importance <- function(data) {
    result <- create_decision_tree(data)
    var_importance <- result$var_importance
    
    # Ensure we have the correct variable names in the right order
    # Based on your RStudio results, the order should be:
    # predictor_method, imbalancing_method, define_recid_method, split, preprocessing, age_category
    expected_names <- c("predictor_method", "imbalancing_method", "define_recid_method", 
                       "split", "preprocessing", "age_category")
    
    # Create a named vector with the correct names
    names(var_importance) <- expected_names
    
    # Sort by importance values (descending) to ensure proper ordering
    var_importance_sorted <- sort(var_importance, decreasing = TRUE)
    
    # Return sorted result
    return(var_importance_sorted)
}

# Function to get regression tree nodes and split rules
get_tree_nodes <- function(data) {
    result <- create_decision_tree(data)
    tree <- result$tree
    
    # 1. Get the tree frame and calculate the depth of each node
    tree_frame <- tree$frame
    node_depths <- floor(log2(as.numeric(row.names(tree_frame)))) + 1
    tree_frame$depth <- node_depths
    
    # 2. Filter for the top layers (e.g., top 5 layers)
    top_layers <- tree_frame[tree_frame$depth <= 5 & tree_frame$var != "<leaf>", ]
    
    # 3. Order the filtered nodes by depth
    ranked_nodes <- top_layers[order(top_layers$depth), ]
    
    # 4. Define the function to get the readable split rule
    get_node_rule <- function(tree, node_num) {
        tree_frame <- tree$frame
        splits_info <- tree$splits
        node_row <- which(rownames(tree_frame) == node_num)
        
        internal_nodes <- tree_frame[tree_frame$var != "<leaf>", ]
        split_index <- which(rownames(internal_nodes) == node_num)
        
        if (length(split_index) == 0) return("No split rule found.")
        
        split_var <- tree_frame$var[node_row]
        split_cutpoint <- splits_info[split_index, "index"]
        split_direction <- splits_info[split_index, "ncat"]
        
        if (split_direction > 1) {
            levels_all <- levels(data[[split_var]])
            left_levels <- levels_all[which(bitwAnd(2^(0:(length(levels_all) - 1)), split_cutpoint) != 0)]
            return(paste0(split_var, " is one of {", paste(left_levels, collapse = ", "), "}"))
        } else {
            return(paste0(split_var, " < ", split_cutpoint))
        }
    }
    
    # 5. Loop through the ranked nodes and collect the rules
    node_rules <- list()
    for (node_num in as.numeric(rownames(ranked_nodes))) {
        rule <- get_node_rule(tree, node_num)
        node_var <- tree_frame[rownames(tree_frame) == node_num, "var"]
        node_depth <- tree_frame[rownames(tree_frame) == node_num, "depth"]
        
        node_rules[[length(node_rules) + 1]] <- list(
            rank = node_depth,
            node = node_num,
            variable = node_var,
            rule = rule
        )
    }
    
    return(node_rules)
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
    
    try:
        # Reset progress for this profile
        progress_updates[profile_key] = []
        progress_counts[profile_key] = 0
        
        # If this is a shared analysis, also reset the other profile's progress
        if is_shared_analysis:
            other_profile_key = "counterfactual" if profile_key == "focal" else "focal"
            progress_updates[other_profile_key] = []
            progress_counts[other_profile_key] = 0
        
        # Update progress
        progress_updates[profile_key].append(f"Starting {profile_name} analysis...")
        progress_updates[profile_key].append(f"Generating 2700 universes...")
        
        # If shared analysis, update both profiles
        if is_shared_analysis:
            other_profile_name = "Counterfactual" if profile_name == "Focal" else "Focal"
            progress_updates[other_profile_key].append(f"Starting {other_profile_name} analysis...")
            progress_updates[other_profile_key].append(f"Generating 2700 universes...")
        
        with localconverter(robjects.default_converter + pandas2ri.converter):
            # Get the R function
            r_run_multiverse = robjects.globalenv['run_multiverse_analysis']
            
            # Define the procedural choice options
            preprocessing_options = robjects.StrVector(["Method_A", "Method_B", "Method_C"])
            split_options = robjects.StrVector(["1:2", "6:4", "7:3", "8:2"])
            age_cat_options = robjects.StrVector(["raw_age_year", "age_cat_compas", "age_cat_nij"])
            imbalancing_options = robjects.StrVector(["Undersampling", "Oversampling", "Male Only", "Female Only", "Weighting"])
            predictor_options = robjects.StrVector(["full", "schmidt", "protected"])
            define_recid_options = robjects.StrVector(["1yr", "2yr", "3yr", "4yr", "5yr"])
            
            # Convert profile parameters to R types
            profile_age_r = robjects.IntVector([profile_age])
            profile_gender_r = robjects.StrVector([profile_gender])
            profile_race_r = robjects.StrVector([profile_race])
            
            # Load the analysis data
            analysis_data = robjects.globalenv['analysis_1978']
            print(f"Running {profile_name} analysis: Age={profile_age}, Gender={profile_gender}, Race={profile_race}")
            
            # Start progress monitoring in a separate thread
            def monitor_progress():
                import os
                import time
                progress_file = f"progress_{profile_age}_{profile_gender}_{profile_race}.txt"
                
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
                                        total_universes = int(parts[3])
                                        progress_counts[profile_key] = universe_num
                                        
                                        # If shared analysis, update both profiles' progress
                                        if is_shared_analysis:
                                            other_profile_key = "counterfactual" if profile_key == "focal" else "focal"
                                            progress_counts[other_profile_key] = universe_num
                                        
                                        # Update progress text every 500 universes
                                        if universe_num % 500 == 0:
                                            progress_updates[profile_key].append(f"Processed {universe_num}/{total_universes} universes...")
                                            if is_shared_analysis:
                                                other_profile_name = "Counterfactual" if profile_name == "Focal" else "Focal"
                                                progress_updates[other_profile_key].append(f"Processed {universe_num}/{total_universes} universes...")
                                        
                                        # Check if analysis is complete
                                        if universe_num >= total_universes:
                                            break
                        except:
                            pass
                    time.sleep(0.5)  # Check every 0.5 seconds
            
            # Start progress monitoring
            import threading
            progress_thread = threading.Thread(target=monitor_progress)
            progress_thread.daemon = True
            progress_thread.start()
            
            # Call the R function with all parameters
            results = r_run_multiverse(analysis_data, preprocessing_options, split_options, 
                                     age_cat_options, imbalancing_options, predictor_options, 
                                     define_recid_options, profile_age_r, profile_gender_r, profile_race_r)
            
            # Convert R results to pandas DataFrame
            results_df = robjects.conversion.rpy2py(results)
            
            # Finalize progress
            progress_counts[profile_key] = 2700
            progress_updates[profile_key].append(f"Generated {len(results_df)} universes")
            progress_updates[profile_key].append(f"{profile_name} analysis completed!")
            
            # If shared analysis, finalize both profiles
            if is_shared_analysis:
                other_profile_key = "counterfactual" if profile_key == "focal" else "focal"
                other_profile_name = "Counterfactual" if profile_name == "Focal" else "Focal"
                progress_counts[other_profile_key] = 2700
                progress_updates[other_profile_key].append(f"Generated {len(results_df)} universes")
                progress_updates[other_profile_key].append(f"{other_profile_name} analysis completed!")
            
            # Clean up progress file
            import os
            progress_file = f"progress_{profile_age}_{profile_gender}_{profile_race}.txt"
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

# Constants for dropdown options - using the same values as defined in combinations.r
PREPROCESSING_OPTIONS = [{'label': method, 'value': method} for method in ["Method_A", "Method_B", "Method_C"]]
SPLIT_OPTIONS = [{'label': split, 'value': split} for split in ["1:2", "6:4", "7:3", "8:2"]]
AGE_CATEGORY_OPTIONS = [{'label': age_cat, 'value': age_cat} for age_cat in ["raw_age_year", "age_cat_compas", "age_cat_nij"]]
IMBALANCING_OPTIONS = [{'label': method, 'value': method} for method in ["Undersampling", "Oversampling", "Male Only", "Female Only", "Weighting"]]
PREDICTOR_OPTIONS = [{'label': method, 'value': method} for method in ["full", "schmidt", "protected"]]
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
        
        # Header
        html.H1("Multiverse Dashboard", style={'textAlign': 'center', 'marginBottom': '20px'}),
        
        # Main horizontal flex container
        html.Div([
            # Left Sidebar - Counterfactual Profiles
            html.Div([
                html.H3("Profiles", style={'textAlign': 'center', 'marginBottom': '12px', 'fontSize': '18px'}),
                html.P("Focal vs Counterfactual Demographic Profiles", style={'textAlign': 'center', 'marginBottom': '10px', 'fontSize': '12px', 'color': '#666'}),
                #html.P("(Read-only profiles - dropdowns disabled)", style={'textAlign': 'center', 'marginBottom': '15px', 'fontSize': '11px', 'color': '#999', 'fontStyle': 'italic'}),
                
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
                
                # Submit Button for Profile Analysis
                html.Div([
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
                            'marginTop': '10px'
                        }
                    ),
                    html.P("Click to analyze the selected profiles and generate multiverse results", 
                           style={'fontSize': '11px', 'color': '#666', 'textAlign': 'center', 'marginTop': '5px', 'fontStyle': 'italic'})
                ], style={'marginBottom': '15px'}),
                
                # Progress Indicator
                html.Div([
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
                            html.P(id='focal-progress-text', children="0/2700 universes", 
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
                            html.P(id='counterfactual-progress-text', children="0/2700 universes", 
                                   style={'fontSize': '10px', 'color': '#666', 'marginTop': '2px'})
                        ], style={'marginBottom': '10px'}),
                        
                        # Add a note about optimizations
                        html.Div([
                            html.P("Note: Identical profiles run once together. Focal profile results are reused when possible.", 
                                   style={'fontSize': '9px', 'color': '#888', 'textAlign': 'center', 'fontStyle': 'italic', 'marginTop': '5px'})
                        ], style={'marginBottom': '10px'})
                    ], style={'marginBottom': '10px'})
                ], style={'marginBottom': '15px'}),
                
                # Dataset Overview Section
                html.Div([
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
                ], style={'marginBottom': '15px'}),
                
            ], style={
                'flex': '0 0 400px',
                'padding': '15px',
                'borderRight': '1px solid #ccc',
                'marginRight': '20px',
                'backgroundColor': '#fefefe'
            }),
            
            # Central Charts & Right Cards Container
            html.Div([
                # Row 1: All Four Components in Horizontal Layout
                html.Div([
                    # 1. Focal Specification Curve
                    html.Div([
                        dcc.Graph(
                            id='spec-curve-1',
                            config={
                                'displayModeBar': True,
                                'modeBarButtonsToRemove': ['pan2d'],
                                'displaylogo': False
                            }
                        ),
                        # Clear Focal Selections button
                        html.Div([
                            html.Button(
                                'Clear Focal Selections',
                                id='clear-focal-btn',
                                n_clicks=0,
                                style={
                                    'backgroundColor': '#d63384',
                                    'color': 'white',
                                    'border': 'none',
                                    'padding': '8px 16px',
                                    'borderRadius': '4px',
                                    'cursor': 'pointer',
                                    'fontSize': '12px',
                                    'marginTop': '10px'
                                }
                            )
                        ], style={'textAlign': 'center'})
                    ], style={'flex': '1.2', 'minWidth': '0'}),
                    
                    # 2. Focal Variable Importance Card
                    html.Div([
                        html.H4(
                            id='variable-importance-title-1',
                            style={'textAlign': 'center', 'marginBottom': '0px', 'color': '#d63384', 'fontSize': '14px', 'marginTop': '60px'}
                        ),
                        html.Div(
                            id='variable-importance-content-1',
                            style={
                                'padding': '6px',
                                'border': '1px solid #ddd',
                                'borderRadius': '5px',
                                'minHeight': '0px',
                                'maxHeight': '300px',
                                'backgroundColor': '#f9f9f9',
                                'overflowY': 'auto'
                            }
                        ),
                        # Focal Profile Radio Button
                        html.Div([
                            html.Div([
                                dcc.RadioItems(
                                    id='grid-profile-switcher-1',
                                    options=[{'label': 'Show Focal Grid', 'value': 1}],
                                    value=1,  # Initially checked
                                    style={'marginTop': '8px'},
                                    labelStyle={
                                        'fontSize': '14px', 
                                        'color': '#d63384', 
                                        'fontWeight': 'bold',
                                        'padding': '10px 16px',
                                        'borderRadius': '25px',
                                        'backgroundColor': '#fff0f5',
                                        'border': '2px solid #d63384',
                                        'cursor': 'pointer',
                                        'display': 'inline-block',
                                        'margin': '0',
                                        'boxShadow': '0 2px 8px rgba(214, 51, 132, 0.2)',
                                        'textAlign': 'center',
                                        'minWidth': '160px'
                                    },
                                    inputStyle={
                                        'marginRight': '10px', 
                                        'transform': 'scale(1.3)',
                                        'accentColor': '#d63384'
                                    }
                                )
                            ], style={
                                'textAlign': 'center', 
                                'marginTop': '8px',
                                'padding': '6px'
                            })
                        ], style={'textAlign': 'center', 'marginTop': '8px'})
                    ], style={'flex': '0.6', 'minWidth': '0', 'alignSelf': 'flex-start'}),
                    
                    # 3. Counterfactual Specification Curve
                    html.Div([
                        dcc.Graph(
                            id='spec-curve-2',
                            config={
                                'displayModeBar': True,
                                'modeBarButtonsToRemove': ['pan2d'],
                                'displaylogo': False
                            }
                        ),
                        # Clear Counterfactual Selections button
                        html.Div([
                            html.Button(
                                'Clear CF Selections',
                                id='clear-cf-btn',
                                n_clicks=0,
                                style={
                                    'backgroundColor': '#0d6efd',
                                    'color': 'white',
                                    'border': 'none',
                                    'padding': '8px 16px',
                                    'borderRadius': '4px',
                                    'cursor': 'pointer',
                                    'fontSize': '12px',
                                    'marginTop': '10px'
                                }
                            )
                        ], style={'textAlign': 'center'})
                    ], style={'flex': '1.2', 'minWidth': '0'}),
                    
                    # 4. Counterfactual Variable Importance Card
                    html.Div([
                        html.H4(
                            id='variable-importance-title-2',
                            style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px', 'marginTop': '60px'}
                        ),
                        html.Div(
                            id='variable-importance-content-2',
                            style={
                                'padding': '6px',
                                'border': '1px solid #ddd',
                                'borderRadius': '5px',
                                'minHeight': '0px',
                                'maxHeight': '300px',
                                'backgroundColor': '#f9f9f9',
                                'overflowY': 'auto'
                            }
                        ),
                        # Counterfactual Profile Radio Button
                        html.Div([
                            html.Div([
                                dcc.RadioItems(
                                    id='grid-profile-switcher-2',
                                    options=[{'label': 'Show CF Grid', 'value': 2}],
                                    value=None,  # Initially unchecked
                                    style={'marginTop': '8px'},
                                    labelStyle={
                                        'fontSize': '14px', 
                                        'color': '#0d6efd', 
                                        'fontWeight': 'bold',
                                        'padding': '10px 16px',
                                        'borderRadius': '25px',
                                        'backgroundColor': '#f0f8ff',
                                        'border': '2px solid #0d6efd',
                                        'cursor': 'pointer',
                                        'display': 'inline-block',
                                        'margin': '0',
                                        'boxShadow': '0 2px 8px rgba(13, 110, 253, 0.2)',
                                        'textAlign': 'center',
                                        'minWidth': '160px'
                                    },
                                    inputStyle={
                                        'marginRight': '10px', 
                                        'transform': 'scale(1.3)',
                                        'accentColor': '#0d6efd'
                                    }
                                )
                            ], style={
                                'textAlign': 'center', 
                                'marginTop': '8px',
                                'padding': '6px'
                            })
                        ], style={'textAlign': 'center', 'marginTop': '8px'})
                    ], style={'flex': '0.6', 'minWidth': '0', 'alignSelf': 'flex-start'}),
                    
                ], style={'display': 'flex', 'flexDirection': 'row', 'alignItems': 'flex-start', 'gap': '15px', 'marginBottom': '20px'}),
                
                # Row 2: Combined Specification Grid
                html.Div([
                    # Combined Specification Grid
                    html.Div([
                        
                        # Instructions
                        html.P(
                            "The specification grid shows the procedural choices for all universes. Drag to select multiple regions in either or both specification curves above to filter the grid view. Multiple selections are accumulated - each new selection adds to previous ones. Each profile supports regional analysis with separate decision trees for each selected region.",
                            style={
                                'fontSize': '12px',
                                'color': '#666',
                                'textAlign': 'center',
                                'marginBottom': '10px',
                                'fontStyle': 'italic'
                            }
                        ),
                        
                        # Combined specification grid
                        html.Div([
                            dcc.Graph(
                                id='combined-spec-grid',
                                config={
                                    'displayModeBar': True,
                                    'modeBarButtonsToRemove': ['pan2d', 'lasso2d', 'select2d'],
                                    'displaylogo': False
                                }
                            )
                        ]),
                        
                        # Selection status indicator and controls (moved below grid)
                        html.Div([
                            html.Div(
                                id='selection-status-indicator',
                                style={
                                    'fontSize': '11px',
                                    'textAlign': 'center',
                                    'marginBottom': '5px'
                                }
                            ),
                            html.Div([
                                html.Button(
                                    "Clear All Selections",
                                    id='clear-selections-button',
                                    style={
                                        'fontSize': '10px',
                                        'padding': '4px 8px',
                                        'backgroundColor': '#dc3545',
                                        'color': 'white',
                                        'border': 'none',
                                        'borderRadius': '4px',
                                        'cursor': 'pointer'
                                    }
                                )
                            ], style={'textAlign': 'center', 'marginBottom': '10px'})
                        ], style={'marginTop': '15px'}),
                        
                    ], style={'flex': '1', 'minWidth': '0'}),
                    
                ], style={'display': 'flex', 'flexDirection': 'row', 'alignItems': 'flex-start', 'gap': '20px'}),
                
            ], style={'flexGrow': 1, 'display': 'flex', 'flexDirection': 'column', 'gap': '20px', 'maxWidth': '1400px', 'minWidth': '0'}),
            
        ], style={'display': 'flex', 'flexDirection': 'row', 'alignItems': 'flex-start', 'gap': '20px'}),
        
    ], style={'fontFamily': 'Arial, sans-serif', 'padding': '20px', 'maxWidth': '1800px', 'margin': 'auto'})

def create_specification_curve(df, profile, profile_num, previously_selected=None):
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
            xaxis_title='Universe Index',
            yaxis_title='Recidivism Probability',
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
                            line=dict(color="orange", width=2, dash="dash")
                        )
                        
                        # End boundary
                        fig.add_shape(
                            type="line",
                            x0=end_x, x1=end_x,
                            y0=y_min, y1=y_max,
                            line=dict(color="orange", width=2, dash="dash")
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
        xaxis_title='Universe Index',
        yaxis_title='Recidivism Probability',
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
            ticktext=[f"U{df_plot.index[i]}" for i in x_positions]
        )
    else:  # Show every nth position for readability
        step = max(1, len(x_positions) // 10)
        selected_positions = x_positions[::step]
        fig.update_xaxes(
            tickmode='array',
            tickvals=selected_positions,
            ticktext=[f"U{df_plot.index[i]}" for i in selected_positions]
        )
    

    
    return fig

def create_combined_specification_grid(df1, df2, selected_universes_1=None, selected_universes_2=None):
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
            xaxis_title='Universe Index',
            yaxis_title='Procedural Choices',
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
        for i, (_, row) in enumerate(df1_display.iterrows()):
            combined_data.append(row)
            universe_labels.append(f"F{i}")
            profile_colors.append('#d63384')  # Pink for focal
            
            # Determine which region this universe belongs to
            # Map the display index back to the original sorted index
            original_index = df1_sorted.index.get_loc(row.name)
            region_num = 1
            for j, region in enumerate(focal_regions):
                if original_index in region:
                    region_num = j + 1
                    break
            region_info.append(f"Focal Region {region_num}")
    
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
        for i, (_, row) in enumerate(df2_display.iterrows()):
            combined_data.append(row)
            universe_labels.append(f"CF{i}")
            profile_colors.append('#0d6efd')  # Blue for counterfactual
            
            # Determine which region this universe belongs to
            # Map the display index back to the original sorted index
            original_index = df2_sorted.index.get_loc(row.name)
            region_num = 1
            for j, region in enumerate(cf_regions):
                if original_index in region:
                    region_num = j + 1
                    break
            region_info.append(f"CF Region {region_num}")
    
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
            xaxis_title='Universe Index',
            yaxis_title='Procedural Choices',
            height=500,
            showlegend=False
        )
        return fig
    
    # Create the grid data matrix
    grid_data = []
    y_labels = []
    
    # Define procedural choices and their options
    procedural_choices = [
        ('define_recid_method', ['5yr', '4yr', '3yr', '2yr', '1yr']),
        ('predictor_method', ['protected', 'schmidt', 'full']),
        ('imbalancing_method', ['Female Only', 'Male Only', 'Oversampling', 'Undersampling', 'Weighting']),
        ('age_category', ['age_cat_nij', 'age_cat_compas', 'raw_age_year']),
        ('split', ['8:2', '7:3', '6:4', '1:2']),
        ('preprocessing', ['Method_C', 'Method_B', 'Method_A'])
    ]
    
    # Create the grid data matrix
    for choice_name, options in procedural_choices:
        column_display_name = choice_name.replace('_', ' ').title()
        
        # Add the options for this category
        for option in options:
            y_labels.append(f"    {option}")
            row = []
            
            for data_row in combined_data:
                # Determine if this option is selected for this universe
                if choice_name == 'define_recid_method':
                    selected = data_row['define_recid_method'] == option
                elif choice_name == 'predictor_method':
                    selected = data_row['predictor_method'] == option
                elif choice_name == 'imbalancing_method':
                    selected = data_row['imbalancing_method'] == option
                elif choice_name == 'age_category':
                    selected = data_row['age_category'] == option
                elif choice_name == 'split':
                    selected = data_row['split'] == option
                elif choice_name == 'preprocessing':
                    selected = data_row['preprocessing'] == option
                else:
                    selected = False
                
                # Binary selection: 1 if selected, 0.5 if not
                if selected:
                    row.append(1)
                else:
                    row.append(0.5)
            
            grid_data.append(row)
        
        # Add the category header
        y_labels.append(f"{column_display_name}")
        header_row = [0] * len(combined_data)
        grid_data.append(header_row)
    
    # Create colorscale for combined view
    colorscale = [
        [0, '#e9ecef'],    # Gray for category headers
        [0.5, '#ffffff'],  # White for unselected
        [1, '#6c757d']     # Gray for selected (neutral color for combined view)
    ]
    
    # Create the heatmap
    fig = go.Figure(data=go.Heatmap(
        z=grid_data,
        x=list(range(len(combined_data))),
        y=y_labels,
        colorscale=colorscale,
        showscale=False,
        hoverongaps=False,
        hoverinfo='z',
        zmin=0,
        zmax=1
    ))
    
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
                line=dict(color="orange", width=3, dash="solid")
            )
        
        # Add region background colors
        current_region = region_info[0] if region_info else ""
        region_start = 0
        
        for i in range(len(region_info) + 1):
            if i == len(region_info) or region_info[i] != current_region:
                # End of current region, add background
                if i > region_start:
                    region_color = 'rgba(220, 53, 69, 0.1)' if 'Focal' in current_region else 'rgba(0, 123, 255, 0.1)'
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
        legend_text += "<br>🟠 Orange lines separate regions"
    
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
        xaxis_title='Universe Index',
        yaxis_title='Procedural Choices',
        height=500,
        showlegend=False,
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

def create_specification_grid(df, profile_num, selected_universes=None):
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
            xaxis_title='Universe Index',
            yaxis_title='Procedural Choices',
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
        ('define_recid_method', ['5yr', '4yr', '3yr', '2yr', '1yr']),
        # Predictor methods (grouped together)
        ('predictor_method', ['protected', 'schmidt', 'full']),
        # Imbalancing methods (grouped together)
        ('imbalancing_method', ['Female Only', 'Male Only', 'Oversampling', 'Undersampling', 'Weighting']),
        # Age categories (grouped together)
        ('age_category', ['age_cat_nij', 'age_cat_compas', 'raw_age_year']),
        # Split ratios (grouped together)
        ('split', ['8:2', '7:3', '6:4', '1:2']),
        # Preprocessing methods (grouped together)
        ('preprocessing', ['Method_C', 'Method_B', 'Method_A'])
    ]
    
    # Create the grid data matrix
    grid_data = []
    y_labels = []
    
    # First, add category headers and their options
    for choice_name, options in procedural_choices:
        column_display_name = choice_name.replace('_', ' ').title()
        
        # Add the options for this category first (indented)
        for option in options:
            y_labels.append(f"    {option}")  # More indented option for better hierarchy
            row = []
            
            for _, row_data in df_display.iterrows():
                # Determine if this option is selected for this universe
                if choice_name == 'define_recid_method':
                    selected = row_data['define_recid_method'] == option
                elif choice_name == 'predictor_method':
                    selected = row_data['predictor_method'] == option
                elif choice_name == 'imbalancing_method':
                    selected = row_data['imbalancing_method'] == option
                elif choice_name == 'age_category':
                    selected = row_data['age_category'] == option
                elif choice_name == 'split':
                    selected = row_data['split'] == option
                elif choice_name == 'preprocessing':
                    selected = row_data['preprocessing'] == option
                else:
                    selected = False
                
                # Simple binary selection: 1 if selected, 0.5 if not
                if selected:
                    row.append(1)
                else:
                    row.append(0.5)
            
            grid_data.append(row)
        
        # Then add the category header (this will appear above the options due to Plotly's y-axis ordering)
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
    
    # Create the heatmap
    fig = go.Figure(data=go.Heatmap(
        z=grid_data,
        x=list(range(len(df_display))),  # Use sequential positions (0, 1, 2, 3...) to match curve
        y=y_labels,
        colorscale=colorscale,
        showscale=False,  # Remove the colorbar to save space
        hoverongaps=False,
        hoverinfo='z',
        zmin=0,
        zmax=1
    ))
    
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
                line=dict(color="orange", width=3, dash="solid")
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
                    region_color = 'rgba(220, 53, 69, 0.1)'
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
        legend_text += "<br>🟠 Orange lines separate regions"
    
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
        xaxis_title=f'Universe Index',
        yaxis_title='Procedural Choices',
        height=500,
        showlegend=False,
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
            ticktext=[f"U{df_display.index[i]}" for i in range(len(df_display))]
        )
    else:  # Show every nth label for readability
        step = max(1, len(df_display.index) // 10)
        selected_positions = list(range(0, len(df_display), step))
        fig.update_xaxes(
            tickmode='array', 
            tickvals=selected_positions, 
            ticktext=[f"U{df_display.index[i]}" for i in selected_positions]
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
        (df['preprocessing'] == profile['preprocessing']) &
        (df['split'] == profile['split']) &
        (df['age_category'] == profile['age_category']) &
        (df['imbalancing_method'] == profile['imbalancing_method']) &
        (df['predictor_method'] == profile['predictor_method']) &
        (df['define_recid_method'] == profile['define_recid_method'])
    ]['recidivism_prob'].iloc[0] if len(df[
        (df['preprocessing'] == profile['preprocessing']) &
        (df['split'] == profile['split']) &
        (df['age_category'] == profile['age_category']) &
        (df['imbalancing_method'] == profile['imbalancing_method']) &
        (df['predictor_method'] == profile['predictor_method']) &
        (df['define_recid_method'] == profile['define_recid_method'])
    ]) > 0 else 0
    
    profile_color = '#d63384' if profile_num == 1 else '#0d6efd'
    profile_name = f'Profile {profile_num}'
    
    fig.add_hline(y=profile_prob, line_dash="dash", line_color=profile_color, 
                  annotation_text="")  # Remove profile annotation text
    
    fig.update_layout(
        #title=f'{profile_name} - Distribution of Recidivism Probabilities',
        yaxis_title='Recidivism Probability',
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
        with localconverter(robjects.default_converter + pandas2ri.converter):
            # Convert pandas DataFrame to R DataFrame
            r_df = robjects.conversion.py2rpy(df)
            
            # Call R function to get variable importance
            r_get_var_importance = robjects.globalenv['get_variable_importance']
            var_importance = r_get_var_importance(r_df)
            
            # Convert R vector to Python dictionary
            var_importance_dict = {}
            
            # Use the variable names from the R results (which should now be properly named)
            # The R function now returns named vectors, so we can extract the names
            if hasattr(var_importance, 'names') and var_importance.names is not None:
                # Use the names from R
                for i, name in enumerate(var_importance.names):
                    var_importance_dict[name] = float(var_importance[i])
            else:
                # Fallback to expected names based on your RStudio results
                expected_names = ["predictor_method", "imbalancing_method", "define_recid_method", 
                                "split", "preprocessing", "age_category"]
                for i in range(len(var_importance)):
                    if i < len(expected_names):
                        var_importance_dict[expected_names[i]] = float(var_importance[i])
                    else:
                        var_importance_dict[f"Var_{i+1}"] = float(var_importance[i])
            
            # Sort by importance (descending) - ensure proper ordering
            sorted_importance = sorted(var_importance_dict.items(), key=lambda x: x[1], reverse=True)
            
            # Debug: print the sorted importance for verification
            print(f"Sorted variable importance: {sorted_importance}")
            
            return sorted_importance
    except Exception as e:
        print(f"Error in R variable importance calculation: {e}")
        import traceback
        traceback.print_exc()
        return []

def get_tree_nodes_r(df):
    """Get regression tree nodes and split rules using R analysis"""
    try:
        with localconverter(robjects.default_converter + pandas2ri.converter):
            # Convert pandas DataFrame to R DataFrame
            r_df = robjects.conversion.py2rpy(df)
            
            # Call R function to get tree nodes
            r_get_tree_nodes = robjects.globalenv['get_tree_nodes']
            tree_nodes = r_get_tree_nodes(r_df)
            
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
    """Identify separate regions within selected indices"""
    if not selected_indices:
        return []
    
    # Sort the indices to identify contiguous regions
    sorted_indices = sorted(selected_indices)
    regions = []
    current_region = [sorted_indices[0]]
    
    for i in range(1, len(sorted_indices)):
        # If the next index is adjacent to the current region, extend it
        if sorted_indices[i] == sorted_indices[i-1] + 1:
            current_region.append(sorted_indices[i])
        else:
            # Gap found, start a new region
            regions.append(current_region)
            current_region = [sorted_indices[i]]
    
    # Add the last region
    regions.append(current_region)
    
    return regions

def get_regional_variable_importance_display(df, profile, profile_num, selected_universes=None):
    """Generate variable importance display treating each selected region separately"""
    # Check if dataframe is empty
    if df.empty:
        return html.Div([
            html.P("This card shows the variable importance analysis", 
                   style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
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
                
                var_importance_html.append(dcc.Graph(
                    figure=bar_fig,
                    config={'displayModeBar': False},
                    style={'height': '80px', 'marginBottom': '10px'}
                ))
            
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
        html.P("This card shows separate decision tree analyses for each selected region", 
               style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
        
        # Variable Importance section
        *var_importance_html,
        
        # Additional info
        html.Hr(style={'margin': '8px 0'}),
        html.P(f"Total Regions: {len(regions)}", style={'fontSize': '10px', 'color': '#666', 'fontStyle': 'italic'})
    ])

def get_combined_variable_importance_display(df1, df2, profile1, profile2, selected_universes_1=None, selected_universes_2=None):
    """Generate variable importance display for combined universe selections from both profiles"""
    # Check if both dataframes are empty
    if df1.empty and df2.empty:
        return html.Div([
            html.P("This card shows the variable importance analysis", 
                   style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
            html.P("No analysis data available. Click 'Run Multiverse Analysis' to generate results.", 
                   style={'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '50px'}),
            html.Hr(style={'margin': '15px 0'}),
            html.P("Profile: Combined Analysis", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'}),
            html.P("Dataset Size: No data available", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'})
        ])
    
    # Combine data from both profiles
    combined_data = []
    
    # Process focal profile data (df1)
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
    
    # Process counterfactual profile data (df2)
    if not df2.empty:
        df2_sorted = df2.sort_values('recidivism_prob')
        if selected_universes_2 is not None and len(selected_universes_2) > 0:
            try:
                selected_universe_ids_2 = [df2_sorted.index[i] for i in selected_universes_2 if isinstance(i, int) and i < len(df2_sorted)]
                if selected_universe_ids_2:
                    df2_display = df2_sorted.loc[selected_universe_ids_2]
                    combined_data.extend(df2_display.to_dict('records'))
            except (IndexError, KeyError, TypeError):
                pass
    
    if not combined_data:
        # No selected universes, use all data
        if not df1.empty:
            combined_data.extend(df1.to_dict('records'))
        if not df2.empty:
            combined_data.extend(df2.to_dict('records'))
    
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
        cf_count = len([i for i in (selected_universes_2 or []) if isinstance(i, int) and i < len(df2)]) if not df2.empty else 0
        
        if focal_count > 0 or cf_count > 0:
            status_indicator = html.Div([
                html.Span("", style={'fontSize': '14px', 'marginRight': '5px'}),
                html.Span(f"Analyzing {focal_count + cf_count} selected universes (Focal: {focal_count}, CF: {cf_count})", 
                         style={'fontSize': '11px', 'color': '#28a745', 'fontWeight': 'bold'})
            ], style={'marginBottom': '8px', 'padding': '4px 8px', 'backgroundColor': '#d4edda', 'borderRadius': '4px', 'border': '1px solid #c3e6cb'})
            var_importance_html.append(status_indicator)
        
        var_importance_html.append(html.H5("Key Decisions (Combined Analysis)", style={'color': '#6c757d', 'marginTop': '0', 'marginBottom': '6px', 'fontSize': '14px'}))
        var_importance_html.append(html.P("Most impactful methods on recidivism probability across both profiles:", style={'fontSize': '12px', 'color': '#666', 'marginBottom': '5px'}))
        
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
                        marker_color='#6c757d',
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
                
                var_importance_html.append(dcc.Graph(
                    figure=bar_fig,
                    config={'displayModeBar': False},
                    style={'height': '100px'}
                ))
        
        # Add tree nodes section
        if tree_nodes:
            var_importance_html.append(html.Hr(style={'margin': '6px 0', 'borderColor': '#ddd'}))
            var_importance_html.append(html.H5("Tree Split Rules (Combined)", style={'color': '#6c757d', 'marginTop': '6px', 'marginBottom': '5px', 'fontSize': '14px'}))
            var_importance_html.append(html.P("Most important decision splits across both profiles:", style={'fontSize': '12px', 'color': '#666', 'marginBottom': '5px'}))
            
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
                        'borderLeft': '3px solid #6c757d'
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
        html.P("This card shows the variable importance analysis for combined universe selections", 
               style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
        
        # Variable Importance section
        *var_importance_html,
        
        # Additional info
        html.Hr(style={'margin': '8px 0'}),
        html.P("Profile: Combined Analysis", style={'fontSize': '10px', 'color': '#666', 'fontStyle': 'italic', 'marginBottom': '2px'}),
        html.P(f"Focal: {profile1['age']} years, {profile1['gender'].title()}, {profile1['race'].replace('_', ' ').title()}", style={'fontSize': '10px', 'color': '#666', 'fontStyle': 'italic', 'marginBottom': '2px'}),
        html.P(f"CF: {profile2['age']} years, {profile2['gender'].title()}, {profile2['race'].replace('_', ' ').title()}", style={'fontSize': '10px', 'color': '#666', 'fontStyle': 'italic', 'marginBottom': '2px'}),
        html.P(f"Dataset Size: {len(analysis_df)} specifications", style={'fontSize': '10px', 'color': '#666', 'fontStyle': 'italic'})
    ])

def get_variable_importance_display(df, profile, profile_num, selected_universes=None):
    """Generate variable importance display for a specific profile"""
    # Check if dataframe is empty (no analysis run yet)
    if df.empty:
        return html.Div([
            html.P("This card shows the variable importance analysis", 
                   style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
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
        
        var_importance_html.append(html.H5(f"Key Decisions (Regression Tree){analysis_title_suffix}", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '0', 'marginBottom': '6px', 'fontSize': '14px'}))
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
                
                # Add the bar plot
                var_importance_html.append(dcc.Graph(
                    figure=bar_fig,
                    config={'displayModeBar': False},
                    style={'height': '100px'}
                ))
        
        # Add tree nodes section
        if tree_nodes:
            var_importance_html.append(html.Hr(style={'margin': '6px 0', 'borderColor': '#ddd'}))
            var_importance_html.append(html.H5(f"Tree Split Rules{analysis_title_suffix}", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '6px', 'marginBottom': '5px', 'fontSize': '14px'}))
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
        html.P("This card shows the variable importance analysis", 
               style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
        
        # Variable Importance section
        *var_importance_html,
        
        # Additional info
        html.Hr(style={'margin': '8px 0'}),
        html.P(f"Dataset Size: {len(analysis_df)} specifications{analysis_title_suffix}", style={'fontSize': '10px', 'color': '#666', 'fontStyle': 'italic'})
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
        html.P(f"Preprocessing: {df['preprocessing'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Split Ratios: {df['split'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Age Categories: {df['age_category'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Imbalancing: {df['imbalancing_method'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Predictors: {df['predictor_method'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.P(f"Recidivism Definitions: {df['define_recid_method'].nunique()}", style={'fontSize': '11px', 'marginBottom': '3px'}),
        html.Hr(style={'margin': '8px 0', 'borderColor': '#ddd'}),
        html.P("Note: Both profiles use the same method combinations with different demographic parameters.", 
               style={'fontSize': '10px', 'color': '#666', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '5px'})
    ], style={'padding': '8px', 'backgroundColor': '#f8f9fa', 'borderRadius': '4px', 'border': '1px solid #dee2e6'})


# Selection handling is now integrated into the main callback



# Global variables to store accumulated selections
accumulated_selections_1 = []  # Store multiple regions for focal profile
accumulated_selections_2 = []  # Store multiple regions for counterfactual profile

# Combined callback to handle selection events and radio button clicks
@app.callback(
    Output('selected-universes-1', 'children'),
    Output('selected-universes-2', 'children'),
    Output('grid-profile-switcher-1', 'value'),
    Output('grid-profile-switcher-2', 'value'),
    Input('spec-curve-1', 'selectedData'),
    Input('spec-curve-2', 'selectedData'),
    Input('grid-profile-switcher-1', 'value'),
    Input('grid-profile-switcher-2', 'value'),
    State('selected-universes-1', 'children'),
    State('selected-universes-2', 'children'),
    prevent_initial_call=True
)
def update_selected_universes(selected_data_1, selected_data_2, switcher_1_value, switcher_2_value, 
                             current_selected_1, current_selected_2):
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
        
        return accumulated_selections_1, dash.no_update, dash.no_update, dash.no_update
    
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
        
        return dash.no_update, accumulated_selections_2, dash.no_update, dash.no_update
    
    elif trigger_id == 'grid-profile-switcher-1':
        # Handle manual radio button click for Focal Profile
        return dash.no_update, dash.no_update, 1, None  # Set focal to 1, clear counterfactual
    
    elif trigger_id == 'grid-profile-switcher-2':
        # Handle manual radio button click for Counterfactual Profile
        return dash.no_update, dash.no_update, None, 2  # Clear focal, set counterfactual to 2
    
    return dash.no_update, dash.no_update, dash.no_update, dash.no_update

# Callback for clear selection buttons
@app.callback(
    Output('selected-universes-1', 'children', allow_duplicate=True),
    Output('selected-universes-2', 'children', allow_duplicate=True),
    Input('clear-selections-button', 'n_clicks'),
    Input('clear-focal-btn', 'n_clicks'),
    Input('clear-cf-btn', 'n_clicks'),
    prevent_initial_call=True
)
def clear_selections(clear_all_clicks, clear_focal_clicks, clear_cf_clicks):
    global accumulated_selections_1, accumulated_selections_2
    
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if trigger_id == 'clear-selections-button':
        # Clear all selections
        accumulated_selections_1 = []
        accumulated_selections_2 = []
        return [], []
    elif trigger_id == 'clear-focal-btn':
        # Clear only focal selections
        accumulated_selections_1 = []
        return [], dash.no_update
    elif trigger_id == 'clear-cf-btn':
        # Clear only counterfactual selections
        accumulated_selections_2 = []
        return dash.no_update, []
    
    return dash.no_update, dash.no_update




# Global variables to store progress updates
progress_updates = {"focal": [], "counterfactual": []}
progress_counts = {"focal": 0, "counterfactual": 0}
total_universes = 2700

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
        focal_text = f"Running together: {progress_counts['focal']}/{total_universes} universes"
        counterfactual_text = f"Running together: {progress_counts['counterfactual']}/{total_universes} universes"
    elif focal_completed_counterfactual_running:
        focal_text = f"Completed (reused): {progress_counts['focal']}/{total_universes} universes"
        counterfactual_text = f"Running: {progress_counts['counterfactual']}/{total_universes} universes"
    else:
        focal_text = f"{progress_counts['focal']}/{total_universes} universes"
        counterfactual_text = f"{progress_counts['counterfactual']}/{total_universes} universes"
    
    return focal_bar_style, focal_text, counterfactual_bar_style, counterfactual_text

# Callback to update all outputs - now triggered by submit button
@app.callback(
    Output('spec-curve-1', 'figure'),
    Output('spec-curve-2', 'figure'),
    Output('combined-spec-grid', 'figure'),
    Output('variable-importance-title-1', 'children'),
    Output('variable-importance-title-1', 'style'),
    Output('variable-importance-content-1', 'children'),
    Output('variable-importance-title-2', 'children'),
    Output('variable-importance-title-2', 'style'),
    Output('variable-importance-content-2', 'children'),
    Output('dataset-overview-content', 'children'),
    Output('selection-status-indicator', 'children'),
    Input('submit-profiles-button', 'n_clicks'),
    Input('selected-universes-1', 'children'),
    Input('selected-universes-2', 'children'),
    Input('grid-profile-switcher-1', 'value'),
    Input('grid-profile-switcher-2', 'value'),
    State('age-slider-1', 'value'),
    State('gender-dropdown-1', 'value'),
    State('race-dropdown-1', 'value'),
    State('age-slider-2', 'value'),
    State('gender-dropdown-2', 'value'),
    State('race-dropdown-2', 'value')
)
def update_dashboard(submit_clicks, selected_1, selected_2, grid_profile_switcher_1, grid_profile_switcher_2,
                    age_1, gender_1, race_1, age_2, gender_2, race_2):
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
            # Only proceed if submit button was clicked, or if it's selection events, or profile switcher
            if trigger_id not in ['submit-profiles-button', 'selected-universes-1', 'selected-universes-2', 'grid-profile-switcher-1', 'grid-profile-switcher-2']:
                # Return current state without updating if other inputs changed
                print(f"Preventing update for trigger: {trigger_id}")
                raise dash.exceptions.PreventUpdate
            
            # If submit button was clicked, run the multiverse analysis
            if trigger_id == 'submit-profiles-button':
                print("Submit button clicked - checking profile characteristics...")
                
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
                    focal_profile_params = (age_1 == 25 and gender_1 == "male" and race_1 == "caucasian")
                    
                    if (focal_profile_has_been_run and not df_nc.empty and 
                        (focal_profile_params or current_focal_params == last_focal_profile_params)):
                        print("Focal profile already run - reusing results, running only counterfactual analysis...")
                        # Set focal progress to completed since we're reusing results
                        progress_counts["focal"] = total_universes
                        progress_updates["focal"].append("Focal profile results reused from previous analysis")
                        
                        # Reuse existing focal results, only run counterfactual
                        counterfactual_results = run_multiverse_analysis_with_profile(age_2, gender_2, race_2)
                        
                        if counterfactual_results is not None:
                            df_low_risk = counterfactual_results
                            print("Counterfactual analysis completed - focal results reused!")
                        else:
                            print("Counterfactual analysis failed!")
                    else:
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
            elif trigger_id in ['selected-universes-1', 'selected-universes-2', 'grid-profile-switcher-1', 'grid-profile-switcher-2']:
                # Handle selection changes or profile switcher changes
                if trigger_id in ['grid-profile-switcher-1', 'grid-profile-switcher-2']:
                    # Profile switcher changed - just update the grid display
                    if trigger_id == 'grid-profile-switcher-1':
                        print(f"Profile switcher changed to: Focal Profile")
                    else:
                        print(f"Profile switcher changed to: Counterfactual Profile")
                else:
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
        spec_curve_1 = create_specification_curve(df_nc, profile1, 1, selected_1)
        spec_curve_2 = create_specification_curve(df_low_risk, profile2, 2, selected_2)
        
        
        # Create the combined specification grid
        # Check if both profiles have selections
        has_focal_selection = selected_1 is not None and len(selected_1) > 0
        has_cf_selection = selected_2 is not None and len(selected_2) > 0
        
        if has_focal_selection and has_cf_selection:
            # Both profiles have selections - show combined view
            combined_spec_grid = create_combined_specification_grid(df_nc, df_low_risk, selected_1, selected_2)
        elif has_focal_selection:
            # Only focal profile has selection
            combined_spec_grid = create_specification_grid(df_nc, 1, selected_1)
        elif has_cf_selection:
            # Only counterfactual profile has selection
            combined_spec_grid = create_specification_grid(df_low_risk, 2, selected_2)
        elif grid_profile_switcher_1 == 1:
            # Focal profile is selected via radio button
            selected_universes = selected_1 if selected_1 is not None else []
            combined_spec_grid = create_specification_grid(df_nc, 1, selected_universes)
        elif grid_profile_switcher_2 == 2:
            # Counterfactual profile is selected via radio button
            selected_universes = selected_2 if selected_2 is not None else []
            combined_spec_grid = create_specification_grid(df_low_risk, 2, selected_universes)
        else:
            # Default to focal profile if neither is explicitly selected
            selected_universes = selected_1 if selected_1 is not None else []
            combined_spec_grid = create_specification_grid(df_nc, 1, selected_universes)
        
        # Set variable importance content based on selection state
        if has_focal_selection and has_cf_selection:
            # Both profiles have selections - show combined analysis in first card and counterfactual in second card
            variable_importance_title_1 = "Combined Variable Importance"
            variable_importance_title_style_1 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#6c757d', 'fontSize': '14px'}
            variable_importance_content_1 = get_combined_variable_importance_display(df_nc, df_low_risk, profile1, profile2, selected_1, selected_2)
            
            # Show counterfactual analysis in the second card
            variable_importance_title_2 = "Counterfactual Variable Importance"
            variable_importance_title_style_2 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px'}
            variable_importance_content_2 = get_variable_importance_display(df_low_risk, profile2, 2, selected_2)
        elif has_focal_selection:
            # Only focal profile has selection - use regional analysis
            variable_importance_title_1 = "Focal Variable Importance (Regional)"
            variable_importance_title_style_1 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
            variable_importance_content_1 = get_regional_variable_importance_display(df_nc, profile1, 1, selected_1)
            
            variable_importance_title_2 = "CF Variable Importance"
            variable_importance_title_style_2 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px'}
            variable_importance_content_2 = get_variable_importance_display(df_low_risk, profile2, 2, selected_2)
        elif has_cf_selection:
            # Only counterfactual profile has selection - use regional analysis
            variable_importance_title_1 = "Focal Variable Importance"
            variable_importance_title_style_1 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
            variable_importance_content_1 = get_variable_importance_display(df_nc, profile1, 1, selected_1)
            
            variable_importance_title_2 = "CF Variable Importance (Regional)"
            variable_importance_title_style_2 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px'}
            variable_importance_content_2 = get_regional_variable_importance_display(df_low_risk, profile2, 2, selected_2)
        else:
            # Show individual profile analyses
            variable_importance_title_1 = "Focal Variable Importance"
            variable_importance_title_style_1 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
            variable_importance_content_1 = get_variable_importance_display(df_nc, profile1, 1, selected_1)
            
            variable_importance_title_2 = "CF Variable Importance"
            variable_importance_title_style_2 = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px'}
            variable_importance_content_2 = get_variable_importance_display(df_low_risk, profile2, 2, selected_2)
        
        
        
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
                html.Span("🔗 ", style={'fontSize': '14px'}),
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
                html.Span("🎯 ", style={'fontSize': '14px'}),
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
                html.Span("🔄 ", style={'fontSize': '14px'}),
                html.Span(f"Counterfactual Selection Active: {cf_count} universes in {cf_region_count} region{'s' if cf_region_count > 1 else ''}", 
                         style={'fontWeight': 'bold', 'color': '#0d6efd'})
            ])
        else:
            selection_status = html.Div([
                html.Span("📊 ", style={'fontSize': '14px'}),
                html.Span("No selections - showing all universes", 
                         style={'color': '#6c757d'})
            ])
        
        return (spec_curve_1, spec_curve_2, combined_spec_grid, 
                variable_importance_title_1, variable_importance_title_style_1, variable_importance_content_1,
                variable_importance_title_2, variable_importance_title_style_2, variable_importance_content_2,
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
