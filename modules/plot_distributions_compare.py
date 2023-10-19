import seaborn as sns
import matplotlib.pyplot as plt

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

def create_split_violin_plots_compare(WT_df1, WT_df2, slop=5000*10):
    # Filter DataFrame for variants with p-value < 0.3
    filtered_df = WT_df1[WT_df1['p_corrected'] < 0.05]
 
    for i, row in filtered_df.iterrows():
        if i > 20:
            break
        print(row)
        if type(row['case_values']) == str:
            case_values = [float(x) for x in row['case_values'].strip('[').strip(']').split(',')]
            control_values = [float(x) for x in row['control_values'].strip('[').strip(']').split(',')]
        else:
            case_values = row['case_values']
            control_values = row['control_values']
        print(f"CASE vs CONT LENGTH: {len(case_values)} vs {len(control_values)}")

        ###########################
        max_value = max(max(case_values), max(control_values))
        min_value = min(min(case_values), min(control_values))
        bins = np.linspace(min_value, max_value, num=20)
        
        case_counts, case_bins = np.histogram(case_values, bins=bins)
        control_counts, control_bins = np.histogram(control_values, bins=bins)
        total_case_samples = len(case_values)
        total_control_samples = len(control_values)
        case_heights = case_counts / total_case_samples
        control_heights = control_counts / total_control_samples
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        axes[0].bar(bins[:-1], case_heights, width=np.diff(bins), alpha=0.6, label='Normalized Case Strling')
        axes[0].bar(bins[:-1], control_heights, width=np.diff(bins), alpha=0.6, label='Normalized Control Strling')
        cnt = 0
        for comp_i, comp_row in WT_df2.iterrows():
            if abs(comp_row['mean_left'] - row['mean_left']) < 2*slop and (
                is_rotation(row['variant'], comp_row['variant']) or is_rotation(row['variant'], rev_complement(comp_row['variant']))
            ):
                cnt+=1
                
                if type(row['case_values']) == str:
                    comp_case_values = [float(x) for x in comp_row['case_values'].strip('[').strip(']').split(',')]
                    comp_control_values = [float(x) for x in comp_row['control_values'].strip('[').strip(']').split(',')]
                else:
                    comp_case_values = comp_row['case_values']
                    comp_control_values = comp_row['control_values']
                
                ################################
                max_value = max(max(comp_case_values), max(comp_control_values))
                min_value = min(min(comp_case_values), min(comp_control_values))
                comp_bins = np.linspace(min_value, max_value, num=20)
                
                # Calculate normalized heights for comparison data
                comp_case_counts, comp_case_bins = np.histogram(comp_case_values, bins=comp_bins)
                comp_control_counts, comp_control_bins = np.histogram(comp_control_values, bins=comp_bins)
                total_comp_case_samples = len(comp_case_values)
                total_comp_control_samples = len(comp_control_values)
                comp_case_heights = comp_case_counts / total_comp_case_samples
                comp_control_heights = comp_control_counts / total_comp_control_samples

                axes[1].bar(comp_bins[:-1], comp_case_heights, width=np.diff(comp_bins), alpha=0.6, label=f'Normalized Case EH')
                axes[1].bar(comp_bins[:-1], comp_control_heights, width=np.diff(comp_bins), alpha=0.6, label=f'Normalized Control EH')

        for ax in axes:
            ax.set_xlabel('Large Allele Size')
            ax.set_ylabel('Frequency')
            ax.legend()

        axes[0].set_title(f"Repeat Unit: {row['variant']} on {row['chrom']} ({row['statistic']} statistic)")
        axes[1].set_title("Comparison")
        
        plt.tight_layout()
        
        # Show the plot
        plt.show()
        
        
        image_format = 'svg' # e.g .png, .svg, etc.
        image_name = f"Repeat Unit: {row['variant']} on {row['chrom']} ({row['statistic']} statistic).svg"

        fig.savefig(os.path.join("figures", image_name), format=image_format, dpi=1200)