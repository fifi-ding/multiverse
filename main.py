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

# Function to run multiverse analysis with profile parameters
def run_multiverse_analysis_with_profile(profile_age, profile_gender, profile_race):
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
        
        # Update progress
        progress_updates[profile_key].append(f"Starting {profile_name} analysis...")
        progress_updates[profile_key].append(f"Generating 2700 universes...")
        
        with localconverter(robjects.default_converter + pandas2ri.converter):
            # Get the R function
            r_run_multiverse = robjects.globalenv['run_multiverse_analysis']
            
            # Define the procedural choice options
            preprocessing_options = robjects.StrVector(["Method_A", "Method_B", "Method_C"])
            split_options = robjects.StrVector(["1:2", "6:4", "7:3", "8:2"])
            age_cat_options = robjects.StrVector(["raw_age_year", "age_cat_compas", "age_cat_nij"])
            imbalancing_options = robjects.StrVector(["Undersampling", "Oversampling", "Male Only", "Female Only", "Weighting"])
            predictor_options = robjects.StrVector(["full", "final", "protected"])
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
                                        
                                        # Update progress text every 500 universes
                                        if universe_num % 500 == 0:
                                            progress_updates[profile_key].append(f"Processed {universe_num}/{total_universes} universes...")
                                        
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
PREDICTOR_OPTIONS = [{'label': method, 'value': method} for method in ["full", "final", "protected"]]
RECID_METHOD_OPTIONS = [{'label': method, 'value': method} for method in ["1yr", "2yr", "3yr", "4yr", "5yr"]]

# Initialize the Dash app
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True)
app.title = "Multiverse Dashboard"


# Add state to track which profile's grid is currently displayed
current_grid_profile = 1  # 1 for Profile 1, 2 for Profile 2, 3 for Cross Comparison

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
                        html.Div(
                            id='selected-profile-summary-1',
                            style={'padding': '8px', 'border': '1px solid #eee', 'borderRadius': '5px', 'backgroundColor': '#fff0f5', 'marginBottom': '12px'}
                        ),
                        
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
                        html.P("🔒 This profile is locked for comparison", 
                               style={'fontSize': '10px', 'color': '#666', 'textAlign': 'center', 'marginTop': '8px', 'fontStyle': 'italic'})
                        
                    ], style={'flex': '1', 'padding': '12px', 'border': '2px solid #d63384', 'borderRadius': '5px', 'marginRight': '8px', 'backgroundColor': '#fff8fa'}),
                    
                    # Profile 2 - Counterfactual
                    html.Div([
                        html.Div("Counterfactual", style={'textAlign': 'center', 'marginBottom': '8px', 'fontSize': '12px', 'color': 'white', 'fontWeight': 'bold', 'backgroundColor': '#0d6efd', 'padding': '4px 8px', 'borderRadius': '12px'}),
                        html.Div(
                            id='selected-profile-summary-2',
                            style={'padding': '8px', 'border': '1px solid #eee', 'borderRadius': '5px', 'backgroundColor': '#f0f8ff', 'marginBottom': '12px'}
                        ),
                        
                        html.Div([
                            html.Label("Age:", style={'fontSize': '14px', 'marginBottom': '4px'}),
                            dcc.Slider(
                                id='age-slider-2',
                                min=0,
                                max=100,
                                step=1,
                                value=35,
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
                                value='female',
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
                                value='african_american',
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
                        ], style={'marginBottom': '10px'})
                    ], style={'marginBottom': '10px'})
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
                # Row 1: Profile 1 vs Profile 2 - Specification Curves and Boxplots
                html.Div([
                    # Profile 1 Specification Curve
                    html.Div([
                        # html.H3("Focal Profile", style={'textAlign': 'center', 'color': '#d63384'}),
                        dcc.Graph(
                            id='spec-curve-1',
                            config={
                                'displayModeBar': True,
                                'modeBarButtonsToRemove': ['pan2d'],
                                'displaylogo': False
                            }
                        ),

                    ], style={'flex': '1.2', 'minWidth': '0'}),
                    
                    # Profile 2 Specification Curve
                    html.Div([
                        # html.H3("Counterfactual Profile", style={'textAlign': 'center', 'color': '#0d6efd'}),
                        dcc.Graph(
                            id='spec-curve-2',
                            config={
                                'displayModeBar': True,
                                'modeBarButtonsToRemove': ['pan2d'],
                                'displaylogo': False
                            }
                        ),

                    ], style={'flex': '1.2', 'minWidth': '0'}),
                    
                    # Variable Importance Card
                    html.Div([
                        html.H4(
                            id='variable-importance-title',
                            style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
                        ),
                        html.Div(
                            id='variable-importance-content',
                            style={
                                'padding': '8px',
                                'border': '1px solid #ddd',
                                'borderRadius': '5px',
                                'minHeight': '400px',
                                'maxHeight': '400px',
                                'backgroundColor': '#f9f9f9',
                                'overflowY': 'auto'
                            }
                        )
                    ], style={'flex': '0.6', 'minWidth': '0', 'alignSelf': 'flex-start'}),
                    
                ], style={'display': 'flex', 'flexDirection': 'row', 'alignItems': 'flex-start', 'gap': '20px', 'marginBottom': '20px'}),
                
                # Row 2: Combined Specification Grid and Universe Cards
                html.Div([
                    # Combined Specification Grid with Toggle
                    html.Div([
                        # Toggle buttons for Profile selection
                        html.Div([
                            html.Button(
                                "Focal Profile Grid",
                                id='toggle-profile-1',
                                style={
                                    'fontSize': '14px',
                                    'padding': '8px 16px',
                                    'marginRight': '10px',
                                    'backgroundColor': '#d63384',
                                    'color': 'white',
                                    'border': 'none',
                                    'borderRadius': '5px',
                                    'cursor': 'pointer',
                                    'fontWeight': 'bold'
                                }
                            ),
                            html.Button(
                                "Counterfactual Profile Grid",
                                id='toggle-profile-2',
                                style={
                                    'fontSize': '14px',
                                    'padding': '8px 16px',
                                    'marginRight': '10px',
                                    'backgroundColor': '#6c757d',
                                    'color': 'white',
                                    'border': 'none',
                                    'borderRadius': '5px',
                                    'cursor': 'pointer'
                                }
                            ),
                            html.Button(
                                "Cross Comparison",
                                id='toggle-profile-3',
                                style={
                                    'fontSize': '14px',
                                    'padding': '8px 16px',
                                    'backgroundColor': '#6c757d',
                                    'color': 'white',
                                    'border': 'none',
                                    'borderRadius': '5px',
                                    'cursor': 'pointer'
                                }
                            )
                        ], style={'textAlign': 'center', 'marginBottom': '0px'}),
                        
                        # Instructions
                        html.P(
                            "Click on a profile button above to view its specification grid. Use 'Cross Comparison' to compare selected universes from both profiles side by side.",
                            style={
                                'fontSize': '12px',
                                'color': '#666',
                                'textAlign': 'center',
                                'marginBottom': '5px',
                                'fontStyle': 'italic'
                            }
                        ),
                        html.P(
                            "💡 Tip: Drag to select universes in the specification curves above to filter the grid view. The grid will automatically switch to show the selected profile. Use 'Cross Comparison' to see selected universes from both profiles together.",
                            style={
                                'fontSize': '11px',
                                'color': '#28a745',
                                'textAlign': 'center',
                                'marginBottom': '0px',
                                'fontStyle': 'italic',
                                'fontWeight': 'bold'
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
                        ])
                    ], style={'flex': '1.3', 'minWidth': '0'}),
                    
                    # Combined Universe Details Card
                    html.Div([
                        html.H4(
                            id='universe-card-title',
                            style={'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
                        ),
                        html.Div(
                            id='combined-universe-card-content',
                            style={
                                'padding': '8px',
                                'border': '1px solid #ddd',
                                'borderRadius': '5px',
                                'minHeight': '500px',
                                'maxHeight': '500px',
                                'backgroundColor': '#f9f9f9',
                                'overflowY': 'auto'
                            }
                        )
                    ], style={'flex': '0.33', 'minWidth': '0'}),
                    
                ], style={'display': 'flex', 'flexDirection': 'row', 'alignItems': 'flex-start', 'gap': '20px'}),
                
            ], style={'flexGrow': 1, 'display': 'flex', 'flexDirection': 'column', 'gap': '20px', 'maxWidth': '1400px', 'minWidth': '0'}),
            
        ], style={'display': 'flex', 'flexDirection': 'row', 'alignItems': 'flex-start', 'gap': '20px'}),
        
    ], style={'fontFamily': 'Arial, sans-serif', 'padding': '20px', 'maxWidth': '1800px', 'margin': 'auto'})

def create_specification_curve(df, profile, profile_num):
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
    

    
    # Profile highlighting removed since we now use dynamic data from R analysis
    
    fig.update_layout(
        title=f'{profile_name}',
        xaxis_title='Universe Index',
        yaxis_title='Recidivism Probability',
        height=400,
        showlegend=True,
        hovermode='x unified',
        dragmode='select',  # Enable selection mode by default
        yaxis=dict(
            range=[0, 1],  # Force y-axis range from 0 to 1 for consistent comparison
            tickmode='linear',
            dtick=0.2,  # Show ticks at 0, 0.2, 0.4, 0.6, 0.8, 1.0
            tickformat='.1f'  # Format ticks as 0.0, 0.2, 0.4, etc.
        )
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
            else:
                df_display = df_sorted
                title_suffix = " - All Universes"
        except (IndexError, KeyError, TypeError):
            # Fallback to showing all universes if there's an error
            df_display = df_sorted
            title_suffix = " - All Universes"
    else:
        # Show all universes
        df_display = df_sorted
        title_suffix = " - All Universes"
    
    # Define procedural choices and their options, grouped as shown in the image
    procedural_choices = [
        # Recidivism time periods (grouped together)
        ('define_recid_method', ['5yr', '4yr', '3yr', '2yr', '1yr']),
        # Predictor methods (grouped together)
        ('predictor_method', ['protected', 'final', 'full']),
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
        y_labels.append(f"📋 {column_display_name}")  # Category header with icon
        
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
    
    # Add a legend for the current view
    if selected_universes is not None and len(selected_universes) > 0:
        legend_text = f"🎯 Showing {len(selected_universes)} Selected Universes"
        legend_color = '#d63384' if profile_num == 1 else '#0d6efd'
    else:
        legend_text = "📊 Showing All Universes"
        legend_color = '#666'
    
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

def create_cross_comparison_grid(df1, df2, selected_1=None, selected_2=None):
    """Create a cross-comparison grid showing both profiles side by side with their selected universes"""
    # Check if dataframes are empty
    if df1.empty or df2.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No analysis data available.<br>Click 'Run Multiverse Analysis' to generate results.",
            xref="paper", yref="paper",
            x=0.5, y=0.5, xanchor='center', yanchor='middle',
            showarrow=False, font=dict(size=16, color='gray')
        )
        fig.update_layout(
            title='Cross Comparison - No Data Available',
            height=500,
            showlegend=False
        )
        return fig
    
    # Sort both dataframes by recidivism probability
    df1_sorted = df1.sort_values('recidivism_prob')
    df2_sorted = df2.sort_values('recidivism_prob')
    
    # Filter based on selected universes if provided
    if selected_1 is not None and len(selected_1) > 0:
        try:
            selected_universe_ids_1 = [df1_sorted.index[i] for i in selected_1 if isinstance(i, int) and i < len(df1_sorted)]
            if selected_universe_ids_1:
                df1_display = df1_sorted.loc[selected_universe_ids_1]
                title_suffix_1 = f" - {len(selected_universe_ids_1)} Selected"
            else:
                df1_display = df1_sorted
                title_suffix_1 = " - All Universes"
        except (IndexError, KeyError, TypeError):
            df1_display = df1_sorted
            title_suffix_1 = " - All Universes"
    else:
        df1_display = df1_sorted
        title_suffix_1 = " - All Universes"
    
    if selected_2 is not None and len(selected_2) > 0:
        try:
            selected_universe_ids_2 = [df2_sorted.index[i] for i in selected_2 if isinstance(i, int) and i < len(df2_sorted)]
            if selected_universe_ids_2:
                df2_display = df2_sorted.loc[selected_universe_ids_2]
                title_suffix_2 = f" - {len(selected_universe_ids_2)} Selected"
            else:
                df2_display = df2_sorted
                title_suffix_2 = " - All Universes"
        except (IndexError, KeyError, TypeError):
            df2_display = df2_sorted
            title_suffix_2 = " - All Universes"
    else:
        df2_display = df2_sorted
        title_suffix_2 = " - All Universes"
    
    # Define procedural choices and their options
    procedural_choices = [
        ('define_recid_method', ['5yr', '4yr', '3yr', '2yr', '1yr']),
        ('predictor_method', ['protected', 'final', 'full']),
        ('imbalancing_method', ['Female Only', 'Male Only', 'Oversampling', 'Undersampling', 'Weighting']),
        ('age_category', ['age_cat_nij', 'age_cat_compas', 'raw_age_year']),
        ('split', ['8:2', '7:3', '6:4', '1:2']),
        ('preprocessing', ['Method_C', 'Method_B', 'Method_A'])
    ]
    
    # Create the grid data matrix for both profiles
    grid_data_1 = []
    grid_data_2 = []
    y_labels = []
    
    # Process each procedural choice
    for choice_name, options in procedural_choices:
        column_display_name = choice_name.replace('_', ' ').title()
        
        # Add options for this category
        for option in options:
            y_labels.append(f"    {option}")
            
            # Profile 1 row
            row_1 = []
            for _, row_data in df1_display.iterrows():
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
                
                row_1.append(1 if selected else 0.5)
            
            # Profile 2 row
            row_2 = []
            for _, row_data in df2_display.iterrows():
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
                
                row_2.append(1 if selected else 0.5)
            
            grid_data_1.append(row_1)
            grid_data_2.append(row_2)
        
        # Add category header
        y_labels.append(f"📋 {column_display_name}")
        header_row_1 = [0] * len(df1_display)
        header_row_2 = [0] * len(df2_display)
        grid_data_1.append(header_row_1)
        grid_data_2.append(header_row_2)
    
    # Create subplot with two heatmaps side by side
    from plotly.subplots import make_subplots
    
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=(f"Focal Profile{title_suffix_1}", f"Counterfactual Profile{title_suffix_2}"),
        horizontal_spacing=0.1
    )
    
    # Add Profile 1 heatmap
    fig.add_trace(go.Heatmap(
        z=grid_data_1,
        x=list(range(len(df1_display))),
        y=y_labels,
        colorscale=[
            [0, '#e9ecef'],    # Gray for headers
            [0.5, '#ffffff'],  # White for unselected
            [1, '#d63384']     # Pink for Profile 1
        ],
        showscale=False,
        hoverongaps=False,
        hoverinfo='z',
        zmin=0,
        zmax=1
    ), row=1, col=1)
    
    # Add Profile 2 heatmap
    fig.add_trace(go.Heatmap(
        z=grid_data_2,
        x=list(range(len(df2_display))),
        y=y_labels,
        colorscale=[
            [0, '#e9ecef'],    # Gray for headers
            [0.5, '#ffffff'],  # White for unselected
            [1, '#0d6efd']     # Blue for Profile 2
        ],
        showscale=False,
        hoverongaps=False,
        hoverinfo='z',
        zmin=0,
        zmax=1
    ), row=1, col=2)
    
    # Update layout
    fig.update_layout(
        title='Cross Comparison: Selected Universes from Both Profiles',
        height=500,
        showlegend=False,
        yaxis=dict(
            tickmode='array',
            tickvals=list(range(len(y_labels))),
            ticktext=y_labels,
            tickfont=dict(size=10)
        ),
        yaxis2=dict(
            tickmode='array',
            tickvals=list(range(len(y_labels))),
            ticktext=y_labels,
            tickfont=dict(size=10)
        )
    )
    
    # Update x-axes
    if len(df1_display.index) <= 20:
        fig.update_xaxes(
            tickmode='array',
            tickvals=list(range(len(df1_display))),
            ticktext=[f"U{df1_display.index[i]}" for i in range(len(df1_display))],
            row=1, col=1
        )
    else:
        step = max(1, len(df1_display.index) // 10)
        selected_positions = list(range(0, len(df1_display), step))
        fig.update_xaxes(
            tickmode='array',
            tickvals=selected_positions,
            ticktext=[f"U{df1_display.index[i]}" for i in selected_positions],
            row=1, col=1
        )
    
    if len(df2_display.index) <= 20:
        fig.update_xaxes(
            tickmode='array',
            tickvals=list(range(len(df2_display))),
            ticktext=[f"U{df2_display.index[i]}" for i in range(len(df2_display))],
            row=1, col=2
        )
    else:
        step = max(1, len(df2_display.index) // 10)
        selected_positions = list(range(0, len(df2_display), step))
        fig.update_xaxes(
            tickmode='array',
            tickvals=selected_positions,
            ticktext=[f"U{df2_display.index[i]}" for i in selected_positions],
            row=1, col=2
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

def get_variable_importance_display(df, profile, profile_num):
    """Generate variable importance display for a specific profile"""
    # Check if dataframe is empty (no analysis run yet)
    if df.empty:
        return html.Div([
            html.P("This card shows the variable importance analysis", 
                   style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
            html.P("No analysis data available. Click 'Run Multiverse Analysis' to generate results.", 
                   style={'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '50px'}),
            html.Hr(style={'margin': '15px 0'}),
            html.P(f"Profile: {'Focal' if profile_num == 1 else 'Counterfactual'}", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'}),
            html.P(f"Demographics: {profile['age']} years, {profile['gender'].title()}, {profile['race'].replace('_', ' ').title()}", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'}),
            html.P(f"Dataset Size: No data available", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'})
        ])
    
    # Get variable importance using R decision tree
    var_importance = get_variable_importance_r(df)
    
    # Get tree nodes using R analysis
    tree_nodes = get_tree_nodes_r(df)
    
    # Debug: print the raw variable importance data
    print(f"Profile {profile_num} - Raw var_importance: {var_importance}")
    print(f"Profile {profile_num} - Tree nodes: {tree_nodes}")
    
    # Create variable importance display
    var_importance_html = []
    if var_importance:
        var_importance_html.append(html.H5("Key Decisions (Regression Tree)", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '0', 'marginBottom': '10px'}))
        var_importance_html.append(html.P("Most impactful methods on recidivism probability:", style={'fontSize': '12px', 'color': '#666', 'marginBottom': '8px'}))
        
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
                    height=200,
                    margin=dict(l=10, r=10, t=10, b=10),
                    showlegend=False,
                    xaxis=dict(
                        showgrid=False, 
                        zeroline=False,
                        title_font=dict(size=10)  # Smaller font size
                    ),
                    yaxis=dict(showgrid=False, zeroline=False),
                    plot_bgcolor='rgba(0,0,0,0)',
                    paper_bgcolor='rgba(0,0,0,0)'
                )
                
                # Add the bar plot
                var_importance_html.append(dcc.Graph(
                    figure=bar_fig,
                    config={'displayModeBar': False},
                    style={'height': '150px'}
                ))
        
        # Add tree nodes section
        if tree_nodes:
            var_importance_html.append(html.Hr(style={'margin': '10px 0', 'borderColor': '#ddd'}))
            var_importance_html.append(html.H5("Tree Split Rules", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '10px', 'marginBottom': '8px'}))
            var_importance_html.append(html.P("Most important decision splits in the regression tree:", style={'fontSize': '12px', 'color': '#666', 'marginBottom': '8px'}))
            
            # Create list of tree nodes
            tree_nodes_html = []
            for i, node in enumerate(tree_nodes[:5]):  # Show top 5 nodes
                rank_emoji = "🥇" if node['rank'] == 1 else "🥈" if node['rank'] == 2 else "🥉" if node['rank'] == 3 else f"#{node['rank']}"
                tree_nodes_html.append(
                    html.Div([
                        html.Span(f"{rank_emoji} ", style={'fontSize': '12px'}),
                        html.Span(f"Node {node['node']}: ", style={'fontSize': '11px', 'fontWeight': 'bold', 'color': '#333'}),
                        html.Span(node['rule'], style={'fontSize': '11px', 'color': '#666'})
                    ], style={
                        'marginBottom': '6px',
                        'padding': '4px 8px',
                        'backgroundColor': '#f8f9fa',
                        'borderRadius': '4px',
                        'borderLeft': f'3px solid {"#d63384" if profile_num == 1 else "#0d6efd"}'
                    })
                )
            
            var_importance_html.extend(tree_nodes_html)

    else:
        var_importance_html.append(html.P("Variable importance analysis not available", style={'color': '#999', 'fontStyle': 'italic'}))
    
    return html.Div([
        # Brief instruction in light gray
        html.P("This card shows the variable importance analysis", 
               style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
        
        # Variable Importance section
        *var_importance_html,
        
        # Additional info
        html.Hr(style={'margin': '15px 0'}),
        html.P(f"Profile: {'Focal' if profile_num == 1 else 'Counterfactual'}", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'}),
        html.P(f"Demographics: {profile['age']} years, {profile['gender'].title()}, {profile['race'].replace('_', ' ').title()}", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'}),
        html.P(f"Dataset Size: {len(df)} specifications", style={'fontSize': '11px', 'color': '#666', 'fontStyle': 'italic'})
    ])

def get_universe_details(df, profile, profile_num):
    """Generate universe details showing the demographic profile information"""
    # Check if dataframe is empty (no analysis run yet)
    if df.empty:
        return html.Div([
            html.P("This card shows the demographic profile information", 
                   style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
            
            # Demographic Profile section
            html.H5("Demographic Profile", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '0', 'marginBottom': '10px'}),
            html.P(f"Age: {profile['age']} years", style={'fontSize': '11px', 'marginBottom': '4px', 'fontWeight': 'bold'}),
            html.P(f"Gender: {profile['gender'].title()}", style={'fontSize': '11px', 'marginBottom': '4px'}),
            html.P(f"Race: {profile['race'].replace('_', ' ').title()}", style={'fontSize': '11px', 'marginBottom': '4px'}),
            
            # Dataset Overview section
            html.H5("Dataset Overview", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '15px', 'marginBottom': '10px'}),
            html.P("No analysis data available. Click 'Run Multiverse Analysis' to generate results.", 
                   style={'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginTop': '20px'})
        ])
    
    return html.Div([
        # Brief instruction in light gray
        html.P("This card shows the demographic profile information", 
               style={'fontSize': '12px', 'color': '#999', 'fontStyle': 'italic', 'textAlign': 'center', 'marginBottom': '15px', 'backgroundColor': '#f8f9fa', 'padding': '8px', 'borderRadius': '5px'}),
        
        # Demographic Profile section
        html.H5("Demographic Profile", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '0', 'marginBottom': '10px'}),
        html.P(f"Age: {profile['age']} years", style={'fontSize': '11px', 'marginBottom': '4px', 'fontWeight': 'bold'}),
        html.P(f"Gender: {profile['gender'].title()}", style={'fontSize': '11px', 'marginBottom': '4px'}),
        html.P(f"Race: {profile['race'].replace('_', ' ').title()}", style={'fontSize': '11px', 'marginBottom': '4px'}),
        
        # Dataset Overview section
        html.H5("Dataset Overview", style={'color': '#d63384' if profile_num == 1 else '#0d6efd', 'marginTop': '15px', 'marginBottom': '10px'}),
        html.P(f"Total Specifications: {len(df)}", style={'fontSize': '11px', 'marginBottom': '4px'}),
        html.P(f"Unique Preprocessing Methods: {df['preprocessing'].nunique()}", style={'fontSize': '11px', 'marginBottom': '4px'}),
        html.P(f"Unique Split Ratios: {df['split'].nunique()}", style={'fontSize': '11px', 'marginBottom': '4px'}),
        html.P(f"Unique Age Categories: {df['age_category'].nunique()}", style={'fontSize': '11px', 'marginBottom': '4px'}),
        html.P(f"Unique Imbalancing Methods: {df['imbalancing_method'].nunique()}", style={'fontSize': '11px', 'marginBottom': '4px'}),
        html.P(f"Unique Predictor Methods: {df['predictor_method'].nunique()}", style={'fontSize': '11px', 'marginBottom': '4px'}),
        html.P(f"Unique Recidivism Methods: {df['define_recid_method'].nunique()}", style={'fontSize': '11px', 'marginBottom': '4px'})
    ])

# Selection handling is now integrated into the main callback



# Callback to handle selection events from specification curve 1
@app.callback(
    Output('selected-universes-1', 'children'),
    Output('toggle-profile-1', 'style', allow_duplicate=True),
    Output('toggle-profile-2', 'style', allow_duplicate=True),
    Output('toggle-profile-3', 'style', allow_duplicate=True),
    Input('spec-curve-1', 'selectedData'),
    prevent_initial_call=True
)
def update_selected_universes_1(selected_data):
    global current_grid_profile
    
    # Update selection
    if selected_data and 'points' in selected_data:
        # Extract the point indices from the selection
        selected_points = selected_data['points']
        selected_indices = [point['pointIndex'] for point in selected_points if 'pointIndex' in point]
        
        # If there are selected universes, switch to profile 1
        if selected_indices:
            current_grid_profile = 1
    else:
        selected_indices = []
    
    # Update button styles based on current profile
    profile1_style = {
        'fontSize': '14px',
        'padding': '8px 16px',
        'marginRight': '10px',
        'backgroundColor': '#d63384' if current_grid_profile == 1 else '#6c757d',
        'color': 'white',
        'border': 'none',
        'borderRadius': '5px',
        'cursor': 'pointer',
        'fontWeight': 'bold' if current_grid_profile == 1 else 'normal'
    }
    
    profile2_style = {
        'fontSize': '14px',
        'padding': '8px 16px',
        'backgroundColor': '#0d6efd' if current_grid_profile == 2 else '#6c757d',
        'color': 'white',
        'border': 'none',
        'borderRadius': '5px',
        'cursor': 'pointer',
        'fontWeight': 'bold' if current_grid_profile == 2 else 'normal'
    }
    
    profile3_style = {
        'fontSize': '14px',
        'padding': '8px 16px',
        'backgroundColor': '#28a745' if current_grid_profile == 3 else '#6c757d',
        'color': 'white',
        'border': 'none',
        'borderRadius': '5px',
        'cursor': 'pointer',
        'fontWeight': 'bold' if current_grid_profile == 3 else 'normal'
    }
    
    return selected_indices, profile1_style, profile2_style, profile3_style

# Callback to handle selection events from specification curve 2
@app.callback(
    Output('selected-universes-2', 'children'),
    Output('toggle-profile-1', 'style', allow_duplicate=True),
    Output('toggle-profile-2', 'style', allow_duplicate=True),
    Output('toggle-profile-3', 'style', allow_duplicate=True),
    Input('spec-curve-2', 'selectedData'),
    prevent_initial_call=True
)
def update_selected_universes_2(selected_data):
    global current_grid_profile
    
    # Update selection
    if selected_data and 'points' in selected_data:
        # Extract the point indices from the selection
        selected_points = selected_data['points']
        selected_indices = [point['pointIndex'] for point in selected_points if 'pointIndex' in point]
        
        # If there are selected universes, switch to profile 2
        if selected_indices:
            current_grid_profile = 2
    else:
        selected_indices = []
    
    # Update button styles based on current profile
    profile1_style = {
        'fontSize': '14px',
        'padding': '8px 16px',
        'marginRight': '10px',
        'backgroundColor': '#d63384' if current_grid_profile == 1 else '#6c757d',
        'color': 'white',
        'border': 'none',
        'borderRadius': '5px',
        'cursor': 'pointer',
        'fontWeight': 'bold' if current_grid_profile == 1 else 'normal'
    }
    
    profile2_style = {
        'fontSize': '14px',
        'padding': '8px 16px',
        'backgroundColor': '#0d6efd' if current_grid_profile == 2 else '#6c757d',
        'color': 'white',
        'border': 'none',
        'borderRadius': '5px',
        'cursor': 'pointer',
        'fontWeight': 'bold' if current_grid_profile == 2 else 'normal'
    }
    
    profile3_style = {
        'fontSize': '14px',
        'padding': '8px 16px',
        'backgroundColor': '#28a745' if current_grid_profile == 3 else '#6c757d',
        'color': 'white',
        'border': 'none',
        'borderRadius': '5px',
        'cursor': 'pointer',
        'fontWeight': 'bold' if current_grid_profile == 3 else 'normal'
    }
    
    return selected_indices, profile1_style, profile2_style, profile3_style

# Combined callback to handle both initialization and toggle functionality
@app.callback(
    Output('toggle-profile-1', 'style'),
    Output('toggle-profile-2', 'style'),
    Output('toggle-profile-3', 'style'),
    Input('dummy-input', 'value'),
    Input('toggle-profile-1', 'n_clicks'),
    Input('toggle-profile-2', 'n_clicks'),
    Input('toggle-profile-3', 'n_clicks')
)
def handle_button_styles(dummy_value, n_clicks_1, n_clicks_2, n_clicks_3):
    global current_grid_profile
    
    # Determine which input triggered the callback
    ctx = dash.callback_context
    if not ctx.triggered:
        # Initial load - set Profile 1 as selected
        current_grid_profile = 1
    else:
        trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
        if trigger_id == 'toggle-profile-1':
            current_grid_profile = 1
        elif trigger_id == 'toggle-profile-2':
            current_grid_profile = 2
        elif trigger_id == 'toggle-profile-3':
            current_grid_profile = 3
    
    # Update button styles based on current selection
    profile1_style = {
        'fontSize': '14px',
        'padding': '8px 16px',
        'marginRight': '10px',
        'backgroundColor': '#d63384' if current_grid_profile == 1 else '#6c757d',
        'color': 'white',
        'border': 'none',
        'borderRadius': '5px',
        'cursor': 'pointer',
        'fontWeight': 'bold' if current_grid_profile == 1 else 'normal'
    }
    
    profile2_style = {
        'fontSize': '14px',
        'padding': '8px 16px',
        'backgroundColor': '#0d6efd' if current_grid_profile == 2 else '#6c757d',
        'color': 'white',
        'border': 'none',
        'borderRadius': '5px',
        'cursor': 'pointer',
        'fontWeight': 'bold' if current_grid_profile == 2 else 'normal'
    }
    
    profile3_style = {
        'fontSize': '14px',
        'padding': '8px 16px',
        'backgroundColor': '#28a745' if current_grid_profile == 3 else '#6c757d',
        'color': 'white',
        'border': 'none',
        'borderRadius': '5px',
        'cursor': 'pointer',
        'fontWeight': 'bold' if current_grid_profile == 3 else 'normal'
    }
    
    return profile1_style, profile2_style, profile3_style

# Global variables to store progress updates
progress_updates = {"focal": [], "counterfactual": []}
progress_counts = {"focal": 0, "counterfactual": 0}
total_universes = 2700
analysis_completed = False

# Callback to show submit button feedback
@app.callback(
    Output('submit-profiles-button', 'children'),
    Output('submit-profiles-button', 'style'),
    Input('submit-profiles-button', 'n_clicks'),
    Input('progress-interval', 'n_intervals')
)
def update_submit_button(n_clicks, n_intervals):
    global analysis_completed
    
    # Update button state
    if n_clicks and n_clicks > 0:
        if analysis_completed:
            button_text = "Complete"
            button_style = {
                'fontSize': '16px',
                'padding': '12px 24px',
                'backgroundColor': '#28a745',
                'color': 'white',
                'border': 'none',
                'borderRadius': '8px',
                'cursor': 'pointer',
                'fontWeight': 'bold',
                'width': '100%',
                'marginTop': '10px'
            }
        else:
            button_text = "Processing..."
            button_style = {
                'fontSize': '16px',
                'padding': '12px 24px',
                'backgroundColor': '#495057',
                'color': 'white',
                'border': 'none',
                'borderRadius': '8px',
                'cursor': 'not-allowed',
                'fontWeight': 'bold',
                'width': '100%',
                'marginTop': '10px'
            }
    else:
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
    
    # Update focal progress bar
    focal_bar_style = {
        'width': f'{focal_progress}%',
        'height': '8px',
        'backgroundColor': '#d63384',
        'borderRadius': '4px',
        'transition': 'width 0.3s ease'
    }
    focal_text = f"{progress_counts['focal']}/{total_universes} universes"
    
    # Update counterfactual progress bar
    counterfactual_bar_style = {
        'width': f'{counterfactual_progress}%',
        'height': '8px',
        'backgroundColor': '#0d6efd',
        'borderRadius': '4px',
        'transition': 'width 0.3s ease'
    }
    counterfactual_text = f"{progress_counts['counterfactual']}/{total_universes} universes"
    
    return focal_bar_style, focal_text, counterfactual_bar_style, counterfactual_text

# Callback to update all outputs - now triggered by submit button
@app.callback(
    Output('spec-curve-1', 'figure'),
    Output('spec-curve-2', 'figure'),
    Output('combined-spec-grid', 'figure'),
    Output('selected-profile-summary-1', 'children'),
    Output('selected-profile-summary-2', 'children'),
    Output('universe-card-title', 'children'),
    Output('universe-card-title', 'style'),
    Output('combined-universe-card-content', 'children'),
    Output('variable-importance-title', 'children'),
    Output('variable-importance-title', 'style'),
    Output('variable-importance-content', 'children'),
    Input('submit-profiles-button', 'n_clicks'),
    Input('toggle-profile-1', 'n_clicks'),
    Input('toggle-profile-2', 'n_clicks'),
    Input('toggle-profile-3', 'n_clicks'),
    Input('selected-universes-1', 'children'),
    Input('selected-universes-2', 'children'),
    State('age-slider-1', 'value'),
    State('gender-dropdown-1', 'value'),
    State('race-dropdown-1', 'value'),
    State('age-slider-2', 'value'),
    State('gender-dropdown-2', 'value'),
    State('race-dropdown-2', 'value')
)
def update_dashboard(submit_clicks, toggle_1, toggle_2, toggle_3, selected_1, selected_2,
                    age_1, gender_1, race_1, age_2, gender_2, race_2):
    global current_grid_profile
    
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
            # Only proceed if submit button was clicked, or if it's grid toggle events, or selection events
            if trigger_id not in ['submit-profiles-button', 'toggle-profile-1', 'toggle-profile-2', 'toggle-profile-3', 'selected-universes-1', 'selected-universes-2']:
                # Return current state without updating if other inputs changed
                raise dash.exceptions.PreventUpdate
            
            # If submit button was clicked, run the multiverse analysis
            if trigger_id == 'submit-profiles-button':
                global analysis_completed
                analysis_completed = False  # Reset completion state
                print("Submit button clicked - running sequential multiverse analysis...")
                # Run analyses sequentially with proper progress tracking
                focal_results = run_multiverse_analysis_with_profile(age_1, gender_1, race_1)
                counterfactual_results = run_multiverse_analysis_with_profile(age_2, gender_2, race_2)
                print("Sequential multiverse analysis completed!")
                
                # Store results globally for use in charts
                global df_nc, df_low_risk
                if focal_results is not None:
                    df_nc = focal_results
                if counterfactual_results is not None:
                    df_low_risk = counterfactual_results
                
                # Mark analysis as completed
                analysis_completed = True
        
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
        # Only recreate curves if not triggered by selection events (but allow toggle events)
        if ctx.triggered and ctx.triggered[0]['prop_id'].split('.')[0] in ['selected-universes-1', 'selected-universes-2']:
            # For selection events, return the current curves to preserve selection
            spec_curve_1 = dash.no_update
            spec_curve_2 = dash.no_update
        else:
            # For toggle events or other triggers, recreate curves
            spec_curve_1 = create_specification_curve(df_nc, profile1, 1)
            spec_curve_2 = create_specification_curve(df_low_risk, profile2, 2)
        
        # Create the combined specification grid based on current selection
        if current_grid_profile == 1:
            # Use selected universes from profile 1 if available
            selected_universes = selected_1 if selected_1 is not None else []
            combined_spec_grid = create_specification_grid(df_nc, 1, selected_universes)
            universe_card_title = "Focal Universe Details"
            universe_card_title_style = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
            combined_universe_content = get_universe_details(df_nc, profile1, 1)
            variable_importance_title = "Focal Variable Importance"
            variable_importance_title_style = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#d63384', 'fontSize': '14px'}
            variable_importance_content = get_variable_importance_display(df_nc, profile1, 1)
        elif current_grid_profile == 2:
            # Use selected universes from profile 2 if available
            selected_universes = selected_2 if selected_2 is not None else []
            combined_spec_grid = create_specification_grid(df_low_risk, 2, selected_universes)
            universe_card_title = "Counterfactual Universe Details"
            universe_card_title_style = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px'}
            combined_universe_content = get_universe_details(df_low_risk, profile2, 2)
            variable_importance_title = "Counterfactual Variable Importance"
            variable_importance_title_style = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#0d6efd', 'fontSize': '14px'}
            variable_importance_content = get_variable_importance_display(df_low_risk, profile2, 2)
        else:  # current_grid_profile == 3 (Cross Comparison)
            # Create cross comparison grid showing both profiles
            combined_spec_grid = create_cross_comparison_grid(df_nc, df_low_risk, selected_1, selected_2)
            universe_card_title = "Cross Comparison Details"
            universe_card_title_style = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#28a745', 'fontSize': '14px'}
            # Create combined universe content for cross comparison
            combined_universe_content = html.Div([
                html.H5("Focal Profile", style={'color': '#d63384', 'marginTop': '0', 'marginBottom': '10px'}),
                get_universe_details(df_nc, profile1, 1),
                html.Hr(style={'margin': '15px 0'}),
                html.H5("Counterfactual Profile", style={'color': '#0d6efd', 'marginTop': '15px', 'marginBottom': '10px'}),
                get_universe_details(df_low_risk, profile2, 2)
            ])
            variable_importance_title = "Cross Comparison Analysis"
            variable_importance_title_style = {'textAlign': 'center', 'marginBottom': '10px', 'color': '#28a745', 'fontSize': '14px'}
            # Create combined variable importance content for cross comparison
            variable_importance_content = html.Div([
                html.H5("Focal Profile", style={'color': '#d63384', 'marginTop': '0', 'marginBottom': '10px'}),
                get_variable_importance_display(df_nc, profile1, 1),
                html.Hr(style={'margin': '15px 0'}),
                html.H5("Counterfactual Profile", style={'color': '#0d6efd', 'marginTop': '15px', 'marginBottom': '10px'}),
                get_variable_importance_display(df_low_risk, profile2, 2)
            ])
    
        profile1_summary = get_profile_summary(profile1)
        profile2_summary = get_profile_summary(profile2)
        
        return spec_curve_1, spec_curve_2, combined_spec_grid, profile1_summary, profile2_summary, universe_card_title, universe_card_title_style, combined_universe_content, variable_importance_title, variable_importance_title_style, variable_importance_content
    
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
                "Error", {}, html.Div("Error"))



app.layout = create_layout()

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8050)
