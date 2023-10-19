import seaborn as sns
import matplotlib.pyplot as plt

def create_split_violin_plots(WT_df):
    # Filter DataFrame for variants with p-value < 0.3
    filtered_df = WT_df[WT_df['p_value'] < 0.2]
    
    # Calculate the number of rows and columns for the grid
    num_plots = len(filtered_df)
    num_cols = min(num_plots, 5)
    num_rows = math.ceil(num_plots / num_cols)
    
    print(f"{num_plots} variants to plot.")

    for i, row in filtered_df.iterrows():
        if type(row['case_values']) == str:
            case_values = [float(x) for x in row['case_values'].strip('[').strip(']').split(',')]
            control_values = [float(x) for x in row['control_values'].strip('[').strip(']').split(',')]
        else:
            case_values = row['case_values']
            control_values = row['control_values']
            
        print(row)
        avg_case = np.mean(case_values)
        avg_cont = np.mean(control_values)
        case_nozeros = [val for val in case_values if val>0]
        cont_nozeros = [val for val in control_values if val>0]

        print(f"case avg (nonzero): {avg_case} ({np.mean(case_nozeros)})")
        print(f"cont avg (nonzero): {avg_cont} ({np.mean(cont_nozeros)})")
        
        # Calculate normalized heights
        case_height = len(case_values) / len(row['case_values'])
        control_height = len(control_values) / len(row['control_values'])
        
        # Normalize values by sample size
        #case_values_normalized = [x / row['case_sample_size'] for x in case_values]
        #control_values_normalized = [x / row['control_sample_size'] for x in control_values]

        # Create a figure and axes for the plot
        fig, ax = plt.subplots()

        # Plot split violin plots for case and control values
        #sns.violinplot(data=[control_values, case_values], split=True, ax=ax, scale={"width": "count", "area": control_height + case_height})
        
        # Manually plot violins with adjusted heights
        #ax.violinplot([control_values], positions=[0], widths=0.8, vert=False, showmeans=False, showmedians=False)
        #ax.violinplot([case_values], positions=[1], widths=0.8, vert=False, showmeans=False, showmedians=False)
        
        # no scaling
        ax.violinplot([control_values], positions=[0], widths=0.8, vert=False, showmeans=True, showmedians=False)
        ax.violinplot([case_values], positions=[1], widths=0.8, vert=False, showmeans=True, showmedians=False)
        
        # Set labels and title for the plot
        #ax.set_xlabel('Control vs Case')
        #ax.set_ylabel('Normalized Allele 2 Size')
        #ax.set_title(f"Repeat Unit: {row['variant']} on {row['chrom']} ({row['effect']} effect)")
        
        # Set labels and title for the plot
        ax.set_yticks([0, 1])
        ax.set_yticklabels(['Control', 'Case'])
        ax.set_xlabel('Allele 2 Size')
        ax.set_title(f"Repeat Unit: {row['variant']} on {row['chrom']} ({row['statistic']} statistic)")

        # Adjust opacity of the violins
        for collection in ax.collections:
            collection.set_alpha(0.6)
        
        # Show the plot
        plt.show()
