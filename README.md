# NC Multiverse Dashboard

An interactive dashboard for exploring the NC Multiverse dataset with specification curves, specification grids, and counterfactual analysis.

## Features

### 🎯 **Counterfactual Profiles**
- **Profile 1**: Set as the baseline profile using **NC Multiverse** dataset (left side)
- **Profile 2**: Set as the counterfactual profile using **Low Risk Multiverse** dataset (right side)
- **Compact side-by-side layout** for easy comparison while maximizing chart space
- **Dual-dataset comparison**: Compare specifications across different risk populations
- **Read-only profiles**: Dropdowns are currently disabled for non-interactive display
- **Fixed configurations**:
  - Profile 1: Method_A, 7:3 split, raw_age_year, Undersampling, full, 1yr
  - Profile 2: Method_B, 8:2 split, age_cat_compas, Oversampling, final, 2yr
  - Preprocessing Method (Method_A, Method_B, Method_C)
  - Split Ratio (1:2, 6:4, 7:3, 8:2)
  - Age Category (raw_age_year, age_cat_compas, age_cat_nij)
  - Imbalancing Method (Undersampling, Oversampling, Male Only, Female Only, Weighting)
  - Predictor Method (full, final, protected)
  - Recidivism Method (1yr, 2yr, 3yr, 4yr, 5yr)

### 📊 **Visualizations**

#### **Profile 1 vs Profile 2 - Side by Side Comparison**

**Specification Curves**
- **Profile 1**: Pink-colored specification curve showing recidivism probabilities
- **Profile 2**: Blue-colored specification curve showing recidivism probabilities
- Each curve highlights the current profile's position with a diamond marker
- Interactive hover and click functionality
- Sorted by probability for better trend visualization

**Boxplots**
- **Profile 1**: Boxplot with pink horizontal line showing current prediction
- **Profile 2**: Boxplot with blue horizontal line showing current prediction
- Shows distribution of all recidivism probabilities
- Helps understand where each profile falls in the overall distribution

**Specification Grids**
- **Profile 1**: Grid showing methods vs universe numbers (NC Multiverse dataset)
- **Profile 2**: Grid showing methods vs universe numbers (Low Risk Multiverse dataset)
- **Y-axis**: Methods (preprocessing, split, age_category, imbalancing_method, predictor_method, define_recid_method)
- **X-axis**: Universe numbers ordered by recidivism probability (matching specification curve)
- **Color coding**: Distinct colors for each method value (e.g., Undersampling=pink, Oversampling=purple)
- **Data independence**: Grids are not affected by profile selections, showing actual method combinations
- **Interactive**: Hover over cells to see exact method values

### 🃏 **Information Cards**

#### **Profile 1 Universe Details Card**
- Profile 1 configuration summary
- Current recidivism probability
- Profile rank and percentile within the dataset
- Dataset overview statistics
- **Scrolling**: Vertical scrollbar for compact display (max height: 200px)

#### **Profile 2 Universe Details Card**
- Profile 2 configuration summary
- Current recidivism probability
- Profile rank and percentile within the dataset
- Dataset overview statistics
- **Scrolling**: Vertical scrollbar for compact display (max height: 200px)

**Note**: Each profile now has its own dedicated universe details card for easy comparison.

**Scrolling Feature**: Universe cards have a maximum height of 200px with vertical scrolling, ensuring compact layout while preserving all information accessibility.

## Layout Overview

The dashboard features a **side-by-side counterfactual profile layout**:

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           NC Multiverse Dashboard                          │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────────────────┐  ┌─────────────────────────────────┐  │
│  │      Counterfactual Profiles    │  │                                 │  │
│  │                                 │  │                                 │  │
│  │  ┌─────────────┐ ┌─────────────┐ │  │      Specification Curve      │  │
│  │  │  Profile 1  │ │  Profile 2  │ │  │      (Panel A)               │  │
│  │  │  (Baseline) │ │(Counterfact)│ │  │                                 │  │
│  │  │             │ │             │ │  │                                 │  │
│  │  │ [Dropdowns] │ │ [Dropdowns] │ │  │                                 │  │
│  │  │             │ │             │ │  │                                 │  │
│  │  └─────────────┘ └─────────────┘ │  │                                 │  │
│  └─────────────────────────────────┘  │                                 │  │
│                                       │                                 │  │
│  ┌─────────────────────────────────┐  │                                 │  │
│  │      Distribution Plot          │  │                                 │  │
│  │                                 │  │                                 │  │
│  └─────────────────────────────────┘  │                                 │  │
│                                       │                                 │  │
│  ┌─────────────────────────────────┐  │                                 │  │
│  │      Specification Grid         │  │                                 │  │
│  │      (Panel B)                  │  │                                 │  │
│  │                                 │  │                                 │  │
│  └─────────────────────────────────┘  │                                 │  │
│                                       │                                 │  │
│  ┌─────────────────────────────────┐  │                                 │  │
│  │      Universe Details           │  │                                 │  │
│  │                                 │  │                                 │  │
│  └─────────────────────────────────┘  │                                 │  │
│                                       │                                 │  │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Installation

1. **Clone or download the project files**
2. **Install dependencies** (if not already installed):
   ```bash
   pip install dash dash-bootstrap-components plotly pandas numpy
   ```
3. **Ensure you have the data file**: `nc_multiverse.csv` in the project directory

## Usage

### Starting the Dashboard

1. **Navigate to the project directory**:
   ```bash
   cd /path/to/nc_multiverse_dashboard
   ```

2. **Run the dashboard**:
   ```bash
   python main.py
   ```

3. **Open your web browser** and go to:
   ```
   http://localhost:8050
   ```

### Using the Dashboard

1. **View Profile 1** (Baseline):
   - Profile 1 is configured with fixed parameters (Method_A, 7:3 split, etc.)
   - Uses the NC Multiverse dataset
   - All parameters are visible but not editable (dropdowns disabled)

2. **View Profile 2** (Counterfactual):
   - Profile 2 is configured with fixed parameters (Method_B, 8:2 split, etc.)
   - Uses the Low Risk Multiverse dataset
   - All parameters are visible but not editable (dropdowns disabled)

3. **Explore the Visualizations**:
   - **Profile 1 vs Profile 2**: Side-by-side comparison of all visualizations
   - **Fixed comparison**: Compare two specific, predefined configurations
   - **Specification Curves**: Click on points to see specific specifications
   - **Specification Grids**: Hover over cells to see exact method values, methods vs universe numbers
   - **Boxplots**: See where each profile falls in the overall distribution

4. **Monitor the Information Cards**:
   - **Profile 1 Universe Details**: Configuration summary, current prediction, and ranking
   - **Profile 2 Universe Details**: Configuration summary, current prediction, and ranking

### Example Use Cases

#### **Current Fixed Configuration Comparison**
- **Profile 1 (NC)**: Method_A, 7:3 split, raw_age_year, Undersampling, full, 1yr
- **Profile 2 (Low Risk)**: Method_B, 8:2 split, age_cat_compas, Oversampling, final, 2yr
- **Analysis**: Compare how different preprocessing, split ratios, and methods affect different risk populations

#### **Dataset-Specific Insights**
- **NC Multiverse**: Baseline population with Method_A preprocessing
- **Low Risk Multiverse**: Low-risk population with Method_B preprocessing
- **Cross-population analysis**: Understand how specification choices impact different risk groups

#### **Age Category Analysis**
- Compare raw age years vs. categorized approaches
- Understand the trade-offs between different age representations across populations

#### **Imbalancing Method Comparison**
- Compare different approaches to handle class imbalance
- See which methods produce more stable predictions across different risk groups

## Data Structure

The dashboard works with two datasets:

**`nc_multiverse.csv`** (Profile 1):
- Contains the main multiverse specifications
- Used for baseline profile analysis

**`low_risk_multiverse.csv`** (Profile 2):
- Contains low-risk population specifications
- Used for counterfactual profile analysis

Both datasets contain the same structure:

- **preprocessing**: Method_A, Method_B, Method_C
- **split**: Train/test split ratios (1:2, 6:4, 7:3, 8:2)
- **age_category**: Age representation methods
- **imbalancing_method**: Class balancing techniques
- **predictor_method**: Feature selection approaches
- **define_recid_method**: Recidivism definition timeframes
- **recidivism_prob**: Predicted recidivism probability

## Technical Details

- **Framework**: Dash (Python web framework)
- **Styling**: Bootstrap components for responsive design
- **Charts**: Plotly for interactive visualizations
- **Data Processing**: Pandas for data manipulation
- **Port**: Default port 8050 (configurable in main.py)

## Customization

### Adding New Parameters
To add new specification parameters:

1. Add the new column to your CSV data
2. Update the dropdown options in `main.py`
3. Modify the profile creation logic
4. Update the visualization functions as needed

### Changing Visualizations
The dashboard uses Plotly for all charts. You can customize:
- Colors and styling
- Chart types and layouts
- Interactive features
- Hover information

### Modifying the Layout
The layout is defined in the `create_layout()` function and can be customized for:
- Different screen sizes
- Additional components
- Alternative arrangements

## Troubleshooting

### Common Issues

1. **Port already in use**:
   - Change the port in `main.py` (line 517)
   - Or kill the existing process: `lsof -ti:8050 | xargs kill`

2. **Data not loading**:
   - Ensure `nc_multiverse.csv` is in the project directory
   - Check file permissions and format

3. **Dependencies missing**:
   - Run: `pip install -r requirements.txt`
   - Or install manually: `pip install dash dash-bootstrap-components plotly pandas numpy`

4. **Dashboard not responding**:
   - Check the terminal for error messages
   - Ensure the Python process is running
   - Verify the correct URL in your browser

### Performance Tips

- For large datasets, consider sampling or aggregating data
- Use browser caching for static assets
- Monitor memory usage with large CSV files

## Contributing

Feel free to enhance the dashboard by:
- Adding new visualization types
- Improving the user interface
- Adding export functionality
- Implementing data filtering options
- Adding statistical analysis tools

## License

This project is open source and available under the MIT License.

---

**Enjoy exploring your multiverse of specifications! 🚀**
