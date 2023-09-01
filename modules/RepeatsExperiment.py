import pandas as pd
import numpy as np
from itertools import zip_longest
import sys
import os
from scipy.stats import ranksums, ks_2samp
from statsmodels.stats.multitest import multipletests
import multiprocessing


class RepeatsExperiment:
    def __init__(self, tsv_dir, csv_metadata, chroms="All", sex=None, tissue=None,
                 dataset=None, cohort=None, race=None, ethnicity=None, apoe=None, 
                 slop=100, slop_modifier=1.5
    ):
        self.test_variable = None
        self.tsv_dir = tsv_dir
        self.csv_metadata = csv_metadata
        self.slop = slop
        self.slop_modifier = slop_modifier
        
        if chroms=="All":
            chroms =  ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
                       'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
                       'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
            
        if type(chroms) is not list:
            chroms = [chroms]
        
        self.chroms = chroms
            
        self.metadict = {
            "Sex": sex,
            "Tissue": tissue,
            "Dataset": dataset,
            "Cohort": cohort,
            "Race": race,
            "Ethnicity": ethnicity,
            "APOE": apoe,
        }

    @staticmethod
    def get_metadata_from_filename(filename):
        """
        Pull subject ID and tissue type from the filename
        """
        basename = os.path.basename(filename)
        subject = "-".join(basename.split("-")[0:3])
        tissue = basename.split("-")[3]
        
        return subject, tissue


    @staticmethod
    def rev_complement(motif):
        rev_dict = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "N": "N"
        }
        return [rev_dict[base] for base in motif[::-1]]

    
    @staticmethod
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
    

    def get_metadata_from_subject(self, subject):
        """
        Match subject name to metadata file and pull the rest of the useful metadata as a dictionary
        """
        meta_df = pd.read_csv(self.csv_metadata)
        subject_metadata = {
            'Dataset': meta_df.loc[meta_df['Subject'] == subject, 'Dataset'].values[0],
            'Disease': meta_df.loc[meta_df['Subject'] == subject, 'Disease'].values[0],
            'Cohort': meta_df.loc[meta_df['Subject'] == subject, 'Cohort'].values[0],
            'Sex': meta_df.loc[meta_df['Subject'] == subject, 'Sex'].values[0],
            'Race': meta_df.loc[meta_df['Subject'] == subject, 'Race'].values[0],
            'Ethnicity': meta_df.loc[meta_df['Subject'] == subject, 'Ethnicity'].values[0],
            'Diagnosis': meta_df.loc[meta_df['Subject'] == subject, 'Diagnosis'].values[0],
            'APOE': meta_df.loc[meta_df['Subject'] == subject, 'APOE'].values[0],
            'Assay': meta_df.loc[meta_df['Subject'] == subject, 'Assay'].values[0],
        }

        return subject_metadata


    def collapse_variants(self, df):
        """
        match adjacent variants from a subject and aggregate row based on chromosome, 
        proximity of the left coordinates and repeat unit and collapse the subject variant 
        into the aggregate variant by updating the mean and variance with the welford online 
        algorithm.
        """
        df.sort_values(['#chrom', 'repeatunit', 'left'], inplace=True)
        df.reset_index(drop=True, inplace=True)
        collapsed_variants = []

        prev_left = None
        prev_unit = None
        prev_counts = None
        prev_chrom = None
        prev_mean_left = None
        prev_m2 = None
        prev_variance_left = None

        for index, row in df.iterrows():
            
            # get new variant stats
            left = row['left']
            unit = row['repeatunit']
            count = row['counts']
            chrom = row['#chrom']
            mean_left = row['mean_left']
            variance_left = row['variance_left']
            m2 = row['m2']
    
            if (  # chromosomes match, repeat units match, within slop range of eachother
                prev_chrom is not None
                and prev_chrom == chrom
                and (self.is_rotation(unit, prev_unit) or self.is_rotation(self.rev_complement(unit), prev_unit))
                and abs(left - prev_left) <= self.slop
            ):
                # welford online algo
                n = prev_counts + count
                delta = mean_left - prev_mean_left
                delta2 = delta * delta

                new_mean = prev_mean_left + (delta * count) / n
                new_m2 = prev_m2 + prev_counts * count * (delta2 / n)
                new_variance_left = new_m2 / (n - 1) if n > 1 else 0

                # set next prev_ values for next iteration
                prev_left = new_mean
                prev_counts = n
                prev_mean_left = new_mean
                prev_m2 = new_m2
                prev_chrom = chrom
                prev_unit = unit
                prev_variance_left = new_variance_left
                
            # next variant not collapsible, add previous and set prev_ vals for next iter
            elif prev_chrom is not None:
                collapsed_variants.append({
                    '#chrom': prev_chrom,
                    'left': prev_left,
                    'repeatunit': prev_unit,
                    'counts': prev_counts,
                    'variance_left': prev_variance_left,
                    'mean_left': prev_mean_left,
                    'm2': prev_m2
                })
                prev_left = left
                prev_unit = unit
                prev_counts = count
                prev_chrom = chrom
                prev_mean_left = mean_left
                prev_m2 = m2
                prev_variance_left = variance_left
            # first iteration, set prev_ values for next iter
            else:
                prev_left = left
                prev_unit = unit
                prev_counts = count
                prev_chrom = chrom
                prev_mean_left = mean_left
                prev_m2 = m2
                prev_variance_left = variance_left
        # add last variant  
        collapsed_variants.append({
            '#chrom': prev_chrom,
            'left': prev_left,
            'repeatunit': prev_unit,
            'counts': prev_counts,
            'variance_left': prev_variance_left,
            'mean_left': prev_mean_left,
            'm2': prev_m2
        })
        
        new_df = pd.DataFrame(collapsed_variants)

        return new_df


    def summarize_experiment(self, summary_csv=None):
        """
        This creates a summarized df of variants accross all subjects by adding variants from a single subject to the aggregated summary df one at a time and collapsing them to the summary
        """
        if summary_csv:
            return pd.read_csv(summary_csv)
        # Initialize an empty DataFrame to store the results
        result_df = pd.DataFrame()
        iteration = 1

        for chrom in self.chroms:
            print(chrom)
            # Load, filter, combine and sort the first two TSV files into a DataFrame
            if iteration==1:
                df1 = self.filter_variants(pd.read_csv(self.tsvs[0], sep='\t'), chrom)
                df2 = self.filter_variants(pd.read_csv(self.tsvs[1], sep='\t'), chrom)
                dff = pd.concat([df1, df2])
                print(dff)
                dff["mean_left"] = dff["left"]
                dff["counts"] = np.ones(dff.shape[0])        
                dff["variance_left"] = np.zeros(dff.shape[0])
                dff["m2"] = np.zeros(dff.shape[0])
                #dff["mean_df"] = dff["left"]
                #print(dff)
                dff.drop(columns=self.cols_to_drop, inplace=True)
                #print(dff)
                #print("before: ", dff.shape)
                #print("--------------- INPUT -------------")
                #print(dff)
                dff = self.collapse_variants(dff)
                #print("--------------- OUTPUT -------------")
                #print(dff)
                dff['left'] = dff['mean_left']
                #print("after: ", dff.shape)
                #print(dff)

            if iteration == 1:
                start_idx = 2
            else:
                start_idx = 0
            files_processed = 2
            for file in self.tsvs[start_idx:]:
                files_processed+=1
                next_df = self.filter_variants(pd.read_csv(file, sep='\t'), chrom)
                next_df["counts"] = np.ones(next_df.shape[0])
                next_df["variance_left"] = np.zeros(next_df.shape[0])
                next_df["mean_left"] = next_df["left"]
                next_df["m2"] = np.zeros(next_df.shape[0])
                next_df.drop(columns=self.cols_to_drop, inplace=True)
                dff_before = dff
                dff = pd.concat([dff, next_df])
                
                #print("--------------- INPUT -------------")
                #print(dff)
                dff = self.collapse_variants(dff)
                # Check if counts are higher than the current iteration
                if any((dff['counts'] > files_processed) & (dff['#chrom'] == chrom)):
                    print(f"Error: Counts are higher than the current iteration ({iteration+1}).")
                    print("Offending rows in the aggregate DataFrame:")
                    print(dff[(dff['counts'] > files_processed) & (dff['#chrom'] == chrom)])
                    print('before')
                    print(dff_before[(dff_before['counts'] > iteration) & (dff_before['#chrom'] == chrom)])
                    print('next')
                    print(next_df[(next_df['repeatunit'] == "AT") & (next_df['#chrom'] == chrom)])
                    raise ValueError("Counts are higher than the current iteration.")

                    #print("--------------- OUTPUT -------------")
                    #print(dff)
                    dff['left'] = dff['mean_left']
            
            iteration+=1

        return dff
    

    def par_process_variant(self, variant_row):
        """
        Match each variant in a summary df to all corresponding variants in each subject df to extract case and control allele2 estimate sizes for each variant (defined by a row in the summary df). This can be done in parallel
        """
        variant = variant_row['repeatunit']
        mean_left = variant_row['mean_left']
        std_left = np.sqrt(variant_row['variance_left'])
        range_ = self.slop * self.slop_modifier
        chrom = variant_row['#chrom']
        matching_variants = []
        matching_filenames = []
        add_similar = True  # if true, when multiple expansions in same STR region, add together. Otherwise skip the sample
        
        # Get the total number of TSV files for this core
        total_files = len(self.tsvs)
        processed_files = 0
        warn_flag = False
        
        cnt_multi = 0
        for tsv_file in self.tsvs:
            tsv_df = pd.read_csv(tsv_file, sep='\t')
            variant_df = tsv_df[tsv_df['repeatunit'] == variant]
            add_variants = variant_df[
                (variant_df['left'] >= mean_left - range_ ) &
                (variant_df['left'] <= mean_left + range_ )
            ][self.test_variable].tolist()
            
            if len(add_variants) > 1:
                #print("WARNING: MORE THAN ONE MATCHING VARIANT FOUND IN SINGLE SUBJECT")
                warn_flag = True
                cnt_multi+=1
                if add_similar:
                    add_variants = np.sum(add_variants)
                    matching_variants.append(add_variants)
                else:
                    continue
                                                
            elif len(add_variants) > 0:
                matching_variants.append(*add_variants)
            else:
                matching_variants.append(0) # assume zero (same as ref) if no matching variants
                
            matching_filenames.append(tsv_file)
        
        case_vals = [matching_variants[idx] for idx, file in enumerate(matching_filenames) if file in self.case_tsvs]
        cont_vals = [matching_variants[idx] for idx, file in enumerate(matching_filenames) if file in self.cont_tsvs]

        total_found = len([val for val in case_vals if np.abs(val)>0] + [val for val in cont_vals if np.abs(val>0)])
        #print(f"counts: {variant_row['counts']} vs {total_found}")
        
        #effect, p_value = ranksums(case_vals, cont_vals)
        statistic, p_value = ks_2samp(case_vals, cont_vals)

        return {
            'chrom': chrom,
            'mean_left': mean_left,
            'std_dev': std_left,
            'variant': variant,
            #'effect': effect,
            'statistic': statistic,
            'p_value': p_value,
            'case_vals': case_vals,
            'cont_vals': cont_vals,
            'recovered_variants': total_found,
            'actual_variants': variant_row['counts'],
            'warning': warn_flag,
            'multi_expansions': cnt_multi
        }


    def perform_wilcoxon_test(self, summary_df, out_file=None):
        """
        Perform the Wilcoxon rank sum test on case and control allele2_size values to determine significant variants
        """
        significant_variants = []

        tsv_files = self.case_tsvs + self.cont_tsvs

        total_variants = len(summary_df)
        processed_variants = 0

        pool = multiprocessing.Pool()  # Create a multiprocessing Pool
        
        for result in pool.imap_unordered(self.par_process_variant, [row for _, row in summary_df.iterrows()]):
        # The rest of the code remains the same
            
            chrom = result['chrom']
            mean_left = result['mean_left']
            std_dev = result['std_dev']
            variant = result['variant']
            #effect = result['effect']
            statistic = result['statistic']
            p_value = result['p_value']
            case_values = result['case_vals']
            cont_values = result['cont_vals']
            act_vars = result['actual_variants']
            rec_vars = result['recovered_variants']
            warning = result['warning']
            multis = result['multi_expansions']

            # Perform the Wilcoxon rank sum test
            #effect, p_value = ranksums(case_values, control_values)

            significant_variants.append({
                'chrom': chrom,
                'mean_left': mean_left,
                'std_dev': std_dev,
                'variant': variant,
                #'effect': effect, 
                'statistic': statistic,
                'p_value': p_value,
                'case_values': case_values,
                'control_values': cont_values,
                'recovered_variants': rec_vars,
                'actual_variants': act_vars,
                'warning': warning,
                'multi_expansions': multis
            })

            processed_variants += 1
            progress = processed_variants / total_variants * 100
            #print(f"Progress: {progress:.2f}% [{processed_variants}/{total_variants}]\r", end='')

        pool.close()  # Close the multiprocessing Pool
        pool.join()   # Wait for all processes to finish

        # Perform multiple testing correction
        p_values = [variant['p_value'] for variant in significant_variants]
        reject, p_corrected, _, _ = multipletests(p_values, method='fdr_bh')

        for i, variant in enumerate(significant_variants):
            variant['p_corrected'] = p_corrected[i]
            variant['reject_null'] = reject[i]

        self.WT_df = pd.DataFrame(significant_variants)
        self.WT_df.sort_values('p_corrected', inplace=True, ascending=True)
        if out_file:
            self.WT_df.to_csv(out_file)

        return self.WT_df