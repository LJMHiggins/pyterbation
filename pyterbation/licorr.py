## Analysis functions ##
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
from adjustText import adjust_text

def list_datasets():
    """
    Returns:
        List of available datasets.
    """
    if not (Path.cwd().parent / "data").is_dir():
        data_dir = Path.cwd() / "data"
    else:
        data_dir = Path.cwd().parent / "data"
    files_in_directory = [f.name for f in data_dir.iterdir() if f.is_file() and ".pickle" in f.name]
    return files_in_directory

def retrieve_data():
    """Retrieve datasets.
    Args:
        Desired datasets
    Returns:
        Dictionary: Keys = dataset ids, values = dataframes
    """
    datasets = dict()
    data_dir = Path.cwd().parent / "data"
    files_in_directory = [file.name for file in data_dir.iterdir() if file.is_file() and ".pickle" in file.name]
    for f in files_in_directory:
        with open(data_dir / f, "rb") as file_p:
            datasets[f] = pickle.load(file_p)
    return datasets

def lm_fit_xdf(df, y):
    """Fit linear model to a dataframe of independent variables veruse a single vector (dependent variable)
    Args: 
        df: Pandas dataframe
        y: dataframe or series with ids to merge with df.
    Returns:
        pd.Dataframe of model parameters.
    """
    def lm_fit(x_col, y_var):
        mask = ~np.isnan(x_col) & ~np.isnan(y_var)
        try: 
            lm = stats.linregress(x_col[mask].astype(float), y_var[mask].values)
            x_var, slope, r, p, r2 = x_col.name, lm.slope, lm.rvalue, lm.pvalue,  lm.rvalue**2
            se_slope, intercept, se_intercept = lm.stderr, lm.intercept, lm.intercept_stderr
            new_row = pd.Series(
                {
                    "x_id": x_var,
                    "slope": slope,
                    "r_value": r,
                    "p_value": p,
                    "r2": r2,
                    "slope_se": se_slope,
                    "intercept": intercept,
                    "intercept_se": se_intercept
                }
            )
            return new_row
        except:
            next
    ## Match datasets
    cell_id_filter = df["CELL_LINE_NAME"].isin(y["CELL_LINE_NAME"])
    y = y.set_index("CELL_LINE_NAME")
    df_match = (df.loc[cell_id_filter, :]
                .set_index("CELL_LINE_NAME")
                .reindex(y.index)
                )
    ## Apply function
    df_lms = df_match.apply(lm_fit, args = (y.iloc[:,0],), axis = 0) # iloc seems to solve it
    return df_lms.T

def spearman_xdf(df, y):
    def spearman_calc(x_col, y_var):
        mask = ~np.isnan(x_col) & ~np.isnan(y_var)
        try:
            result = stats.spearmanr(x_col[mask].values, y_var[mask].values)
            x_id = x_col.name
            spr, p_value = result.correlation, result.pvalue
            new_row = pd.Series(
                {
                    "x_id": x_id,
                    "spearman_r": spr,
                    "p_value": p_value
                }
            )
            return new_row
        except:
            next
    ## Match datasets
    cell_id_filter = df["CELL_LINE_NAME"].isin(y["CELL_LINE_NAME"])
    y = y.set_index("CELL_LINE_NAME")
    df_match = (df.loc[cell_id_filter, :]
                .set_index("CELL_LINE_NAME")
                .reindex(y.index)
                )
    df_spear = df_match.apply(spearman_calc, args = (y.iloc[:,0],), axis = 0)
    return df_spear.T

# Modify above code to accept dataframe
# df x df
def lin_func_df(df1, df2):
    ## Inner function for applying across dataframe ##
    def lin_func(col):
        l = stats.linregress(col.astype(float), gene_i)
        phos_pro, slope, intercept, r, p, se_slope, se_intercept = (col.name, l.slope, l.intercept, l.rvalue, l.pvalue, l.stderr, l.intercept_stderr)
        new_row = {"phospho_protein":phos_pro, "slope":slope,"intercept":intercept, 
                   "r_value":r, "p_value":p, "std_error_slope":se_slope, "std_error_intercept":se_intercept}
        return pd.Series(new_row)
    gene_dict = {}
    for i, gene in enumerate(df1.columns):
        gene_i = df1.loc[:,gene]
        result_df = df2.apply(lin_func, axis = 0)
        gene_dict[gene] = result_df.T
        print("At RNAi {0}, number {1}".format(gene, i+1)) #Would be good to identify the core performing each calculation
    return gene_dict

#### Coeffienct volcano
def coefficient_volcano(df, 
                        x, 
                        y,
                        reg_or_corr = "reg",
                        title = None, 
                        x_label = None, 
                        label_points = True,
                        cmap = "viridis"):
    ## Plot construction ##
    df["p_value_transformed"] = -np.log10(df[y].astype(float))
    fig, ax = plt.subplots(figsize=(10, 8))
    if reg_or_corr == "reg":
        scatter = ax.scatter(x = x, 
                             y = "p_value_transformed",  # Updated to the transformed values
                             c="r2", cmap=cmap,
                             data=df,
                             alpha=0.7, 
                             edgecolors='w')
        cbar = plt.colorbar(scatter)
        cbar.set_label('Coefficient of determination (R^2)')
    elif reg_or_corr == "corr":
        if "r_value" in df.columns:
            x_column = 'r_value'
        elif "spearman_r" in df.columns:
            x_column = "spearman_r"
        else:
            raise ValueError("No suitable correlation coefficient column found in the dataframe.")
        scatter = ax.scatter(x = x_column, 
                             y = "p_value_transformed",  
                             data=df,
                             alpha=0.7, 
                             dgecolors='w')
    else: print("Clarify plot type.")
    ## Plot style formatting ##
    if x_label is not None:
        ax.set_xlabel(x_label)
    if title is not None:
        ax.set_title(title)
        
    ax.set_ylabel('-log10(p value)') 
    threshold = -np.log10(0.05)
    ax.axhline(threshold, color='red', linestyle='--', label='Threshold: -log10(0.05)')
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    border_color = '#808080'
    ax.spines['top'].set_color(border_color)
    ax.spines['right'].set_color(border_color)
    ax.spines['bottom'].set_color(border_color)
    ax.spines['left'].set_color(border_color)
    ax.grid(True, linestyle='-', linewidth=0.5, color='lightgrey')
    ## Point labelling ##
    if label_points:
        x_threshold = 0.3  # Change this to your desired x threshold
        y_threshold = -np.log10(0.01)  # Change this to your desired y threshold
        # Label the top points that exceed the given x and y value thresholds
        texts = []
        for index, row in df.iterrows():
            if (row[x] > x_threshold or row[x] < -x_threshold) and row['p_value_transformed'] > y_threshold:
                text = plt.text(row[x], row['p_value_transformed'], f"({row['x_id']})", ha='right')
                texts.append(text)
        # Use adjust_text to perform text repulsion and move only the labels up
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', lw=0.5), 
                    only_move={'points':'+y', 'text':'+y'},
                    force_text=(1, 2))
    # Increase axes span to accommodate text labels
    ax.set_xlim(ax.get_xlim()[0] - 0.5, ax.get_xlim()[1] + 0.5)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] + 0.5)

    plt.show()

#### Regression plot ####
def scatter_plot_with_regression(df, 
                                 x_label, 
                                 y_label, 
                                 ax = None,
                                 title = None):
    if ax is None:
        fig, ax = plt.subplots() 
    # Extracting data from the DataFrame
    x = df.loc[:, x_label]
    y = df.loc[:, y_label]
    # Calculating the least squares regression line
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line = slope * x + intercept
    spearman_corr, spearman_p_value = stats.spearmanr(x, y)
    # Plotting the scatter plot and the regression line
    ax.scatter(x, y, color="#000000", alpha=0.7, label =None)
    ax.plot(x, line, color='#FF7F50', alpha=0.7, label =None)
    # Adding labels and title
    ax.set_xlabel(x_label)
    if title is not None:
        ax.set_title(title)
    ax.set_xlabel(x_label, fontsize=14)
    ax.set_ylabel(y_label, fontsize=14)
    ax.set_title(x_label, fontsize=16, weight='bold')
    ax.set_ylim(0, 0.8)
    # Increasing border thickness
    ax.spines['top'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    # Making the borders darker
    border_color = '#808080'
    ax.spines['top'].set_color(border_color)
    ax.spines['right'].set_color(border_color)
    ax.spines['bottom'].set_color(border_color)
    ax.spines['left'].set_color(border_color)
    # Adjust tick label font sizes
    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelsize=12)
    # Displaying the grid
    ax.grid(True, linestyle='-', linewidth=0.5, color='lightgrey')
    # Displaying the legend
    legend_text = f'Lin. Reg.: Slope: {slope:.2f}, p-value: {p_value:.4f}\nSpearman: Coef: {spearman_corr:.2f}, p-value: {spearman_p_value:.4f}'
    ax.text(0.98, 0.95, legend_text, fontsize=13, ha='right', va='top', bbox=dict(facecolor='white', alpha=0.8), transform=ax.transAxes)
    
    if ax is None:
        plt.show()
## Regression panel ##
def scatter_panel(df,
                  protein_list,
                  y_label,
                  n_cols = 4,
                  figsize = (16, 4)):
    total_graphs = len(protein_list)
    n_rows = (total_graphs + n_cols - 1) // n_cols  # Calculate the number of rows needed

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten()  # Flatten

    for i, x_label in enumerate(protein_list):
        ax = axes[i]
        scatter_plot_with_regression(df=df, x_label=x_label, y_label=y_label, ax=ax)
        ax.set_title(f'{x_label}', fontsize=15)
    for j in range(total_graphs, n_rows * n_cols):
        fig.delaxes(axes[j])
        
    plt.tight_layout()
