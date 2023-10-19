import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def rev_complement(motif):
        rev_dict = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C"
        }
        return [rev_dict[base] for base in motif[::-1]]
    
def is_rotation(str1, str2):
    """
    Checks if string is a rotation of another
    """
    if len(str1) != len(str2):
        return False

    if str1 == str2:
        return True

    # to handle idexing in circular rotation
    concat = str1 + str1
    n = len(concat)

    # make prefix table as per Knuth-Morris-Pratt (KMP) algorithm
    prefix = [0] * n
    j = 0
    for i in range(1, n):
        if concat[i] == concat[j]:
            j += 1
            prefix[i] = j
        else:
            if j > 0:
                j = prefix[j-1]
                i -= 1
            else:
                prefix[i] = 0

    # Check if string 2 is a substring in the concatenated string with prefix table
    j = 0
    for i in range(n):
        if concat[i] == str2[j]:
            j += 1
            if j == len(str2):
                return True
        else:
            if j > 0:
                j = prefix[j-1]
                i -= 1

    # string 2 not found as substring
    return False

def create_3_plots(df1, df2, df3, slop=5000*10):
    # Filter DataFrame for variants with p-value < 0.3
    #filtered_df = WT_df1[WT_df1['variant'] == "ACGGGGAGAGGGAGAGGGAG"]
    #print(filtered_df)

     # Create a figure and axes for the plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 8))  # 12, 5

    fdf1 = df1[df1['p_corrected'] < 0.05]
    fdf1.sort_values(['chrom', 'variant', 'mean_left'], inplace=True)
    fdf1.reset_index(drop=True, inplace=True)

    ## First Dataframe
    for i, row1 in fdf1.iterrows():
        if i > 20:
            break
        print(row1)
        if type(row1['case_values']) == str:
            case_values1 = [float(x) for x in row1['case_values'].strip('[').strip(']').split(',')]
            control_values1 = [float(x) for x in row1['control_values'].strip('[').strip(']').split(',')]
        else:
            case_values1 = row1['case_values']
            control_values1 = row1['control_values']
            
        ###########################
        max_value1 = max(max(case_values1), max(control_values1))
        min_value1 = min(min(case_values1), min(control_values1))
        bins1 = np.linspace(min_value, max_value, num=20)
        
        # Calculate normalized heights by dividing counts by the total number of samples
        case_counts1, case_bins1 = np.histogram(case_values1, bins=bins1)
        control_counts1, control_bins1 = np.histogram(control_values1, bins=bins1)
        total_case_samples1 = len(case_values1)
        total_control_samples1 = len(control_values1)
        case_heights1 = case_counts1 / total_case_samples1
        control_heights1 = control_counts1 / total_control_samples1
        ##########################
        
        # Plot split violin plots for case and control values
        #sns.violinplot(data=[control_values, case_values], split=True, ax=ax, scale={"width": "count", "area": control_height + case_height})
        
        # no scaling
        ##ax.violinplot([control_values], positions=[0], widths=0.8, vert=False, showmeans=True, showmedians=False)
        ##ax.violinplot([case_values], positions=[1], widths=0.8, vert=False, showmeans=True, showmedians=False)
        axes[0].bar(bins1[:-1], case_heights1, width=np.diff(bins1), alpha=0.6, label='Case Strling')
        axes[0].bar(bins1[:-1], control_heights1, width=np.diff(bins1), alpha=0.6, label='Control Strling')
        cnt = 0

        for i2, row2 in df2.iterrows():
            if abs(row2['mean_left'] - row1['mean_left']) < 2*slop and (
                is_rotation(row1['variant'], row2['variant']) or is_rotation(row1['variant'], rev_complement(row2['variant']))
            ):
                cnt+=1
                
                if type(row['case_values']) == str:
                    case_values2 = [float(x) for x in row2['case_values'].strip('[').strip(']').split(',')]
                    control_values2 = [float(x) for x in row2['control_values'].strip('[').strip(']').split(',')]
                else:
                    case_values2 = row2['case_values']
                    control_values2 = row2['control_values']
                
                ################################
                max_value2 = max(max(case_values2), max(control_values2))
                min_value2 = min(min(case_values2), min(control_values2))
                bins2 = np.linspace(min_value2, max_value2, num=20)
                
                # Calculate normalized heights for comparison data
                case_counts2, case_bins2 = np.histogram(case_values2, bins=bins2)
                control_counts2, control_bins2 = np.histogram(control_values2, bins=bins2)
                total_case_samples2 = len(case_values2)
                total_control_samples2 = len(control_values2)
                case_heights2 = case_counts2 / total_case_samples2
                control_heights2 = control_counts2 / total_control_samples2

                ################################
                
                #ax.violinplot([comp_control_values], positions=[2*cnt], widths=0.8, vert=False, showmeans=True, showmedians=False)
                #ax.violinplot([comp_case_values], positions=[2*cnt+1], widths=0.8, vert=False, showmeans=True, showmedians=False)
                ##ax.violinplot([comp_control_values], positions=[2], widths=0.8, vert=False, showmeans=True, showmedians=False)
                ##ax.violinplot([comp_case_values], positions=[3], widths=0.8, vert=False, showmeans=True, showmedians=False)
                axes[1].bar(bins2[:-1], case_heights2, width=np.diff(bins2), alpha=0.6, label=f'Case EHDN-LRDN')
                axes[1].bar(bins2[:-1], control_heights2, width=np.diff(bins2), alpha=0.6, label=f'Control EHDN-LRDN')

        for i3, row3 in df3.iterrows():
            if abs(row3['mean_left'] - row2['mean_left']) < 2*slop and (
                is_rotation(row2['variant'], row3['variant']) or is_rotation(row2['variant'], rev_complement(row3['variant']))
            ):
                cnt+=1
                
                if type(row['case_values']) == str:
                    case_values3 = [float(x) for x in row3['case_values'].strip('[').strip(']').split(',')]
                    control_values3 = [float(x) for x in row3['control_values'].strip('[').strip(']').split(',')]
                else:
                    case_values2 = row3['case_values']
                    control_values2 = row3['control_values']
                
                ################################
                max_value3 = max(max(case_values3), max(control_values3))
                min_value3 = min(min(case_values3), min(control_values3))
                bins3 = np.linspace(min_value2, max_value3, num=20)
                
                # Calculate normalized heights for comparison data
                case_counts3, case_bins3 = np.histogram(case_values3, bins=bins3)
                control_counts3, control_bins3 = np.histogram(control_values3, bins=bins3)
                total_case_samples3 = len(case_values3)
                total_control_samples3 = len(control_values3)
                case_heights3 = case_counts3 / total_case_samples3
                control_heights3 = control_counts3 / total_control_samples3

                ################################
                
                #ax.violinplot([comp_control_values], positions=[2*cnt], widths=0.8, vert=False, showmeans=True, showmedians=False)
                #ax.violinplot([comp_case_values], positions=[2*cnt+1], widths=0.8, vert=False, showmeans=True, showmedians=False)
                ##ax.violinplot([comp_control_values], positions=[2], widths=0.8, vert=False, showmeans=True, showmedians=False)
                ##ax.violinplot([comp_case_values], positions=[3], widths=0.8, vert=False, showmeans=True, showmedians=False)
                axes[2].bar(bins3[:-1], case_heights3, width=np.diff(bins3), alpha=0.6, label=f'Case EHDN')
                axes[2].bar(bins3[:-1], control_heights3, width=np.diff(bins3), alpha=0.6, label=f'Control EHDN')


        # Set labels and title for the subplots
        for ax in axes:
            ax.set_xlabel('Large Allele Size')
            ax.set_ylabel('Frequency')
            ax.legend()

        axes[0].set_title(f"Repeat Unit: {row1['variant']} on {row1['chrom']} ({row1['statistic']} statistic)")
        #axes[1].set_title("Comparison")

        plt.tight_layout()
        plt.show()
        
        #image_format = 'svg' # e.g .png, .svg, etc.
        #image_name = f"Repeat Unit: {row1['variant']} on {row1['chrom']} ({row1['statistic']} statistic).svg"
        #fig.savefig(os.path.join("figures", image_name), format=image_format, dpi=1200)
        
