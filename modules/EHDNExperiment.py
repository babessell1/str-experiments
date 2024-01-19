import pandas as pd
import numpy as np
from itertools import zip_longest
import sys
import os
from scipy.stats import ranksums, ks_2samp
from statsmodels.stats.multitest import multipletests
import multiprocessing
from modules.RepeatsExperiment import RepeatsExperiment

class EHDNExperiment(RepeatsExperiment):
    def __init__(self, tsv_dir, csv_metadata, **kwargs):
        super().__init__(tsv_dir, csv_metadata, **kwargs)
        
        self.cols_to_drop = [
            'allele_est', 'right', 'merged_expansions'
        ]
        self.test_variable = 'allele_est'
        self.lrdn_files = None
        self.caller = "EHDN"
        self.locus_file_extension = 'locus.tsv'
        self.repeat_nomenclature = "motif"
        

    @staticmethod
    def filter_variants(tsv_df, chrom):
        """
        Remove single nucleotide expansions and select by chromosome
        """
        filtered_df = tsv_df[~tsv_df['motif'].isin(['A', 'G', 'C', 'T'])]
        if chrom != 'All':
            filtered_df = filtered_df[filtered_df['contig'] == chrom]

        return filtered_df


    def collapse_tsv_file(self, tsv_file, out_dir):
        """
        Process a single TSV file, collapsing variants.
        """

        if os.stat(tsv_file).st_size == 0:
            print(f"Skipping empty file: {tsv_file}")
            return
        
        subject, tissue = self.get_metadata_from_filename(tsv_file)
        df = pd.read_csv(tsv_file, sep='\t')
        df = self.filter_variants(df, "All")
        df.sort_values(['contig', 'motif', 'start'], inplace=True)
        df.reset_index(drop=True, inplace=True)
        df['counts'] = np.ones(df.shape[0])
        df['mean_left'] = df.start
        collapsed_variants = []

        # print tsv fie if it's empty
        if df.shape[0] == 0:
            print(f'{tsv_file} is empty')
            return

        prev_left = None
        prev_unit = None
        prev_counts = None
        prev_chrom = None
        prev_mean_left = None
        prev_sum_allele_est = None

        for index, row in df.iterrows():
            # get new variant stats
            left = row['start']
            unit = row['motif']
            count = row['counts']
            chrom = row['contig']
            mean_left = row['mean_left']
            new_allele_est = row['het_str_size']

            if (
                prev_chrom is not None
                and prev_chrom == chrom
                and (self.is_rotation(unit, prev_unit) or self.is_rotation(self.rev_complement(unit), prev_unit))
                and abs(left - prev_left) <= self.slop * self.slop_modifier
            ):
                # update counts, mean left, and add allele2 size estimate to the previous variant
                update_sum_allele_est = new_allele_est + prev_sum_allele_est
                n = prev_counts + count
                delta = mean_left - prev_mean_left
                new_mean = prev_mean_left + (delta * count) / n

                # set next prev_ values for the next iteration
                prev_left = new_mean
                prev_counts = n
                prev_mean_left = new_mean
                prev_chrom = chrom
                prev_unit = unit
                prev_sum_allele_est = update_sum_allele_est

            # next variant not collapsible, add previous and set prev_ vals for next iter
            elif prev_chrom is not None:
                collapsed_variants.append({
                    '#chrom': prev_chrom,
                    'left': int(prev_left),
                    'right': int(prev_left) + 1,
                    'repeatunit': prev_unit,
                    'allele_est': prev_sum_allele_est,
                    'counts': prev_counts,
                    'mean_left': prev_mean_left
                })
                prev_left = left
                prev_unit = unit
                prev_counts = count
                prev_chrom = chrom
                prev_mean_left = mean_left
                prev_sum_allele_est = new_allele_est

            # first iteration, set prev_ values for the next iteration
            else:
                prev_left = left
                prev_unit = unit
                prev_counts = count
                prev_chrom = chrom
                prev_mean_left = mean_left
                prev_sum_allele_est = new_allele_est
                
        # add last variant
        collapsed_variants.append({
            '#chrom': prev_chrom,
            'left': int(prev_left),
            'right': int(prev_left) + 1,
            'repeatunit': prev_unit,
            'allele_est': prev_sum_allele_est,
            'counts': prev_counts,
            'mean_left': prev_mean_left
        })

        new_df = pd.DataFrame(collapsed_variants)
        new_df['merged_expansions'] = new_df['counts'] > 1
        new_df['left'] = new_df['mean_left']
        new_df.drop(columns=['counts', 'mean_left'], inplace=True)
        
        # Write the result to an output file
        tsv_file_name = os.path.basename(tsv_file)
        out_file_path = os.path.join(out_dir, tsv_file_name)
        new_df.to_csv(out_file_path, sep='\t', index=False)

    
    def lrdn_process_chromosome(self, iter):
        """
        iteratate through rows of subset cdf from merge output of expansion hunter denovo-lrdn, get values for each sample in rows
        split into case/controls and perform statistical test
        row = looks like:
        HEADERS: contig  start   end counts
        SAMPLE_VAL: chr1    1000    1001    SAMPLE_NAME_1:{normalized counts},SAMPLE_NAME_2:{normalized counts},...
        """
        i = iter[0]
        #print(f"## {i} iter... ##", end='\r')
        row = iter[1]
        p_values = []
        processed_samples = []
        chrom = row[0]
        left = row['start']
        case_names = []
        case_counts = []
        cont_names = []
        cont_counts = []
        split = row["counts"].split(",")
        errors = 0
        for entry in split:
            s = entry.split(":")
            subject, _ = self.get_metadata_from_filename(s[0])
            subject_metadata = self.get_metadata_from_subject(subject)
            print(subject_metadata['Diagnosis'])
            if subject_metadata["Diagnosis"] == "Unknown":
                continue
            elif subject_metadata["Diagnosis"] == "Case":
                case_names.append(subject)
                case_counts.append(s[1])
            elif subject_metadata["Diagnosis"] == "Control":
                cont_names.append(subject)
                cont_counts.append(s[1])
            
            processed_samples.append(subject)

        if len(processed_samples) < 100:
            #print(f"## {i} not enough samples... ##", end='\r')
            return None

        print(f"## {i} adding uncalled samples... ##", end='\r')
        # add uncalled samples to the dataframe (assume they are unexpanded)
        print(f"number of lrdn files: {len(self.lrdn_files)}")
        for file in self.lrdn_files:
            subject, _ = self.get_metadata_from_filename(file)
            if i == 1:
                print(f"{i} {subject}")
            if subject not in processed_samples:
                subject_metadata = self.get_metadata_from_subject(subject)

                if subject_metadata["Diagnosis"] == "Unknown":
                    continue
                elif subject_metadata["Diagnosis"] == "Case":
                    case_names.append(subject)
                    case_counts.append(0)
                    if i == 1:
                        print(f"{i} {subject} case")
                elif subject_metadata["Diagnosis"] == "Control":
                    if i == 1:
                        print(f"{i} {subject} cont")
                    cont_names.append(0)
                    cont_counts.append(subject)
        
        print(f"## {i} performing statistical test... ##", end='\r')

        if self.test == "KS":
            # perform Kolmogorov-Smirnov test
            statistic, p_value = ks_2samp(case_counts, cont_counts)
        elif self.test == "WT":
            # perform Wilcoxon rank sum test
            statistic, p_value = ranksums(case_counts, cont_counts)
        elif self.test == "AD":
            # perform Anderson-Darling test
            res = anderson_ksamp([case_counts, cont_counts])
            statistic  = res.statistic
            p_value = res.pvalue
        else:
            ValueError("Invalid test type. Must be one of 'KS', 'RS', or 'AD'")
        
        p_values.append(p_value)

        progress = i / len(self.lrdn_files) * 100
        print(f"#### Progress: {progress:.2f}% [{i}/{len(self.lrdn_files)}] ########################################", end='\r')

        return {
            'chrom': chrom,
            'mean_left': left,
            'std_dev': "NA",
            'variant': "NA",
            #'effect': effect,
            'statistic': statistic,
            'p_value': p_value,
            'case_vals': case_counts,
            'cont_vals': cont_counts,
            'recovered_variants': len(case_counts) + len(cont_counts),
            'actual_variants': "NA",
            'warning': "NA",
            'multi_expansions': "NA"
        }


    def process_lrdn(self, merge_file):
        # Create a multiprocessing Pool

        significant_variants = []
        processed_variants = 0

        pool = multiprocessing.Pool()

        print("Processing LRDN variants...")

        for chrom in self.chroms:
            # Parallelize the summarize_chromosome function for each chromosome
            df = pd.read_csv(merge_file, sep='\t')
            df = df[df['contig'] == chrom]
            # sort the merge file by start position
            df.sort_values(['start'], inplace=True)
            for result in pool.imap_unordered(self.lrdn_process_chromosome, [(i, row) for (i, row) in df.iterrows()]):
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
                    
            # Close the pool and wait for all processes to finish
            pool.close()
            pool.join()

            # Perform multiple testing correction
            p_values = [variant['p_value'] for variant in significant_variants]
            reject, p_corrected, _, _ = multipletests(p_values, method='fdr_bh')

            for i, variant in enumerate(significant_variants):
                variant['p_corrected'] = p_corrected[i]
                variant['reject_null'] = reject[i]

            self.WT_df = pd.DataFrame(significant_variants)
            self.WT_df.sort_values('p_corrected', inplace=True, ascending=True)
            self.WT_df.to_csv(f"results/LRDN_{self.test}_{self.chroms}_{self.cohort}_{self.apoe}lrdn.csv", index=False)

        return  
