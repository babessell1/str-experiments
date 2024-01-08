import pandas as pd
import numpy as np
from itertools import zip_longest
import sys
import os
from scipy.stats import ranksums, ks_2samp, anderson_ksamp
from statsmodels.stats.multitest import multipletests
import multiprocessing
import warnings

# Filter out the "p-value capped" warning
warnings.filterwarnings("ignore", category=UserWarning, message="p-value capped: true value larger than 0.25")

class RepeatsExperiment:
    def __init__(self, tsv_dir, csv_metadata, chroms="All", sex=None, tissue=None,
                 dataset=None, cohort=None, race=None, ethnicity=None, apoe=None, 
                 slop=100, slop_modifier=1.5, test="KS", motifs_to_drop=[]
    ):
        self.test_variable = None
        self.tsv_dir = tsv_dir
        self.csv_metadata = csv_metadata
        self.slop = slop
        self.slop_modifier = slop_modifier
        self.test = test
        self.caller = "None"
        self.motifs_to_drop = motifs_to_drop

        if cohort is None or cohort == "all_cohorts":
            cohort = "all_cohorts"
            self.cohort = cohort
        else:
            # make a manifest of only the specifid cohort
            manifest_df = pd.read_csv(self.csv_metadata)
            self.csv_metadata = os.path.join("manifests", os.path.basename(csv_metadata).split(".")[0] + "_" + cohort + ".csv")
            manifest_df[manifest_df['Cohort'] == cohort].to_csv(self.csv_metadata)

        if apoe is None or apoe == "all_apoe":
            apoe = "all_apoe"
        else:
            manifest_df = pd.read_csv(self.csv_metadata)
            # remove nans APOE
            manifest_df = manifest_df[manifest_df['APOE'].notna()]
            # APOE column should be an string
            manifest_df['APOE'] = manifest_df['APOE'].astype(int).astype(str)
            self.csv_metadata = os.path.join("manifests", os.path.basename(csv_metadata).split(".")[0] + "_" + apoe + ".csv")
            manifest_df[manifest_df['APOE'] == str(apoe)].to_csv(self.csv_metadata)

        if self.test == "AD":
            warnings.warn("")
        
        if chroms=="All":
            chroms =  ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
                       'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
                       'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
            
        if type(chroms) is not list:
            chroms = [chroms]
        
        self.chroms = chroms
        self.cohort = cohort

        # coerce apoe to string of int if possible
        if apoe is not None:
            if type(apoe) is not str:
                apoe = str(int(apoe))

        self.apoe = apoe

        print(self.csv_metadata)

        # if coha

            
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
        if basename.startswith("ADNI"):
            subject = "_".join(basename.split("_")[0:4])
            tissue = "unknown"
        else:
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
        return "".join([rev_dict[base] for base in motif[::-1]])

    
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
            'Assay': meta_df.loc[meta_df['Subject'] == subject, 'Assay'].values[0],
        }
        if self.apoe is not None and self.apoe != "all_apoe":
            subject_metadata['APOE'] = meta_df.loc[meta_df['Subject'] == subject, 'APOE'].values[0].astype(int).astype(str)
        else:
            subject_metadata['APOE'] = meta_df.loc[meta_df['Subject'] == subject, 'APOE'].values[0]

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
    
    @staticmethod
    def add_motif_dir_to_tsv_path(tsv, motif):
        """
        Add the motif directory to the second to last position in the tsv path
        """
        tsv_split = tsv.split("/")
        tsv_split.insert(-1, motif)
        return "/".join(tsv_split)


    def summarize(self, chrom, motif=None):
        # make log file
        #log_file = open(f"logs/{self.caller}_summary_{chrom}_{self.cohort}_{self.apoe}.log", "w")
    
        # get subset of expanded variants for the chromosome and motif
        tsv_to_read = self.add_motif_dir_to_tsv_path(self.tsvs[0], motif)
        df1 = pd.read_csv(tsv_to_read, sep='\t')
        df1 = df1[df1['#chrom'] == chrom]
        if motif is not None:
            df1 = df1[df1['repeatunit'] == motif]

        tsv_to_read = self.add_motif_dir_to_tsv_path(self.tsvs[1], motif)
        df2 = pd.read_csv(tsv_to_read, sep='\t')
        df2 = df2[df2['#chrom'] == chrom]
        if motif is not None:
            df2 = df2[df2['repeatunit'] == motif]

        # merge the two dataframes to create ongoing aggregate dataframe
        dff = pd.concat([df1, df2])
        dff["mean_left"] = dff["left"]
        dff["counts"] = np.ones(dff.shape[0])
        dff["variance_left"] = np.zeros(dff.shape[0])
        dff["m2"] = np.zeros(dff.shape[0])
        dff.drop(columns=self.cols_to_drop, inplace=True)
        dff = self.collapse_variants(dff)
        dff['left'] = dff['mean_left']

        start_idx = 2
        files_processed = 2
        iteration = 1
        for file in self.tsvs[start_idx:]:
            files_processed += 1
            tsv_to_read = self.add_motif_dir_to_tsv_path(file, motif)
            next_df = pd.read_csv(tsv_to_read, sep='\t')
            next_df = next_df[next_df['#chrom'] == chrom]
            next_df["counts"] = np.ones(next_df.shape[0])
            next_df["variance_left"] = np.zeros(next_df.shape[0])
            next_df["mean_left"] = next_df["left"]
            next_df["m2"] = np.zeros(next_df.shape[0])
            next_df.drop(columns=self.cols_to_drop, inplace=True)
            dff_before = dff
            dff = pd.concat([dff, next_df])
            dff = self.collapse_variants(dff)

            if any((dff['counts'] > files_processed) & (dff['#chrom'] == chrom)):
                print(f"Error: Counts are higher than the current iteration ({iteration+1}).")
                print("Offending rows in the aggregate DataFrame:")
                print(dff[(dff['counts'] > files_processed) & (dff['#chrom'] == chrom)])
                print('before')
                print(dff_before[(dff_before['counts'] > iteration) & (dff_before['#chrom'] == chrom)])
                print('next')
                print(next_df[(next_df['repeatunit'] == "AT") & (next_df['#chrom'] == chrom)])
                raise ValueError("Counts are higher than the current iteration.")
            
            iteration += 1

            progress = iteration / len(self.tsvs) * 100
            # write progress to logs
            #log_file.write(f"#### {chrom} Progress: {progress:.2f}% [{iteration}/{len(self.tsvs) }] ##\n")

            # print progress every 20th iteration
            if iteration % 20 == 0:
                print(f"## {chrom} Progress: {progress:.2f}% [{iteration}/{len(self.tsvs) }] ##", end='\r')

            # write to file
                
        if motif is None:
            motif = ""
        
        # don't write to file if no variants found
        if dff.shape[0] == 0:
            return
        
        dff.to_csv(f"results/{self.caller}_summary_{chrom}_{self.cohort}_{self.apoe}_{motif}.csv")
        
        return

    
    def summarize_experiment(self, summary_csv=None):
        print("Creating Summary file...")

        if summary_csv:
            return pd.read_csv(summary_csv)
        
        motif_df = pd.read_csv(os.path.join(f'{self.caller}_motifs.txt'))
        motifs = motif_df['motif'].tolist()
        chromosomes = self.chroms
    

        # Create a multiprocessing Pool
        pool = multiprocessing.Pool()

        # Parallelize the summarize function for each chromosome and motif combination
        pool.imap_unordered(self.summarize, [(chrom, motif) for chrom in chromosomes for motif in motifs])

        # Close the pool and wait for all processes to finish
        pool.close()
        pool.join()

        return None

    def par_process_variant(self, iter):
        """
        Match each variant in a summary df to all corresponding variants in each subject df to extract case and control allele2 estimate sizes for each variant (defined by a row in the summary df). This can be done in parallel
        """
        i = iter[0]
        variant_row = iter[1]
        variant = variant_row['repeatunit']
        mean_left = variant_row['mean_left']
        std_left = np.sqrt(variant_row['variance_left'])
        range_ = self.slop * self.slop_modifier
        chrom = variant_row['#chrom']
        matching_variants = []
        matching_filenames = []
        add_similar = True  # if true, when multiple expansions in same STR region, add together. Otherwise skip the sample
        
        warn_flag = False
        
        cnt_multi = 0
        for tsv_file in self.tsvs:
            tsv_df = pd.read_csv(tsv_file, sep='\t')
            variant_df = tsv_df[tsv_df['repeatunit'] == variant]
            add_variants = variant_df[
                (variant_df['left'] >= mean_left - range_ ) &
                (variant_df['left'] <= mean_left + range_ )
            ][self.test_variable].tolist()

            #if variant == "AT" and chrom == "chr1" and abs(mean_left - 52797421) < range_:
            if len(add_variants) > 1:
                #print(f"MULTI relative to {mean_left}-----------------------")
                for row in add_variants:
                    #print(row)
                    pass
                
                #print(f"IF SORTED")
                sorted_variant_df = variant_df.sort_values('left')
                add_variants_sorted = sorted_variant_df[
                    (sorted_variant_df['left'] >= mean_left - range_ ) &
                    (sorted_variant_df['left'] <= mean_left + range_ )
                ][self.test_variable].tolist()
                for row in add_variants_sorted:
                    #print(row)
                    pass
            
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
        if i % 20 == 0:
            print(f"## Progress: {i}/{self.total_to_process} | total_found: {total_found} ##", end='\r')
        # dont bother testing if the amount of expansions is low
        if sum(val > 0 for val in case_vals) + sum(val > 0 for val in cont_vals) < 10:
            return None

        if self.test == "KS":
            # perform Kolmogorov-Smirnov test
            statistic, p_value = ks_2samp(case_vals, cont_vals)
        elif self.test == "WT":
            # perform Wilcoxon rank sum test
            statistic, p_value = ranksums(case_vals, cont_vals)
        elif self.test == "AD":
            # perform Anderson-Darling test
            res = anderson_ksamp([case_vals, cont_vals])
            statistic  = res.statistic
            p_value = res.pvalue
        else:
            ValueError("Invalid test type. Must be one of 'KS', 'RS', or 'AD'")

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
    

    def perform_stat_test(self, summary_df):
        """
        Perform the Wilcoxon rank sum test on case and control allele2_size values to determine significant variants
        """
        significant_variants = []

        total_variants = len(summary_df)
        self.total_to_process = total_variants
        processed_variants = 0

        pool = multiprocessing.Pool()  # Create a multiprocessing Pool

        #summary_df = summary_df[summary_df['repeatunit'] == "AT"]
        
        for result in pool.imap_unordered(self.par_process_variant, [(i, row) for (i, row) in summary_df.iterrows()]):
        # The rest of the code remains the same
            if result is not None:
                significant_variants.append({
                    'chrom': result['chrom'],
                    'mean_left': result['mean_left'],
                    'std_dev': result['std_dev'],
                    'variant': result['variant'],
                    #'effect': effect, 
                    'statistic': result['statistic'],
                    'p_value': result['p_value'],
                    'case_values': result['case_vals'],
                    'control_values': result['cont_vals'],
                    'recovered_variants': result['recovered_variants'],
                    'actual_variants': result['actual_variants'],
                    'warning':  result['warning'],
                    'multi_expansions': result['multi_expansions']
                })

                processed_variants += 1
                progress = processed_variants / total_variants * 100
                #print(f"Progress: {progress:.2f}% [{processed_variants}/{total_variants}] ########################################\r", end='')

        pool.close()  # Close the multiprocessing Pool
        pool.join()   # Wait for all processes to finish

        # Perform multiple testing correction
        p_values = [variant['p_value'] for variant in significant_variants]
        reject, p_corrected, _, _ = multipletests(p_values, method='fdr_bh')

        for i, variant in enumerate(significant_variants):
            variant['p_corrected'] = p_corrected[i]
            variant['reject_null'] = reject[i]

        # convert list of chromosomes to string
        if len(self.chroms) == 1:
            chroms = self.chroms[0]
        else:
            chroms = "_".join(self.chroms)

        self.WT_df = pd.DataFrame(significant_variants)
        self.WT_df.sort_values('p_corrected', inplace=True, ascending=True)
        self.WT_df.to_csv(f"results/{self.caller}_{self.test}_{chroms}_{self.cohort}_{self.apoe}.csv", index=False)

        return