import pandas as pd
import numpy as np
from itertools import zip_longest
import sys
import os
from scipy.stats import ranksums, ks_2samp
from statsmodels.stats.multitest import multipletests
import multiprocessing
from modules.RepeatsExperiment import RepeatsExperiment

class STRlingExperiment(RepeatsExperiment):
    def __init__(self, tsv_dir, csv_metadata, chroms="All", sex=None, tissue=None,
                 dataset=None, cohort=None, race=None, ethnicity=None, apoe=None, 
                slop=100, slop_modifier=1.5
    ):
        super().__init__(tsv_dir, csv_metadata, chroms, sex, tissue, dataset, cohort, race, ethnicity, apoe, slop, slop_modifier)
        self.cols_to_drop = [
            'allele1_est', 'allele2_est', 'right', 'merged_expansions'
                ]
        self.test_variable = 'allele2_est'
   
    def filter_tsv_files(self):
        """
        Get case/control TSV file lists from the directory and filter files that adhere to desired covariates described by metadict
        """
        self.tsvs = []
        self.case_tsvs = []
        self.cont_tsvs = []
        
        for file in os.listdir(self.tsv_dir):
            if file.endswith('-genotype.txt'):
                subject, tissue = self.get_metadata_from_filename(file)
                if not self.metadict['Tissue'] or tissue == self.metadict['Tissue']:
                    subject_metadata = self.get_metadata_from_subject(subject)
                    add_flag = True
                    for key, val in self.metadict.items():
                        if val and val != subject_metadata[key]:
                            add_flag = False
                            break
                    if add_flag:
                        file_path = os.path.join(self.tsv_dir, file)
                        if subject_metadata["Diagnosis"] != "Unknown":
                            self.tsvs.append(file_path)
                            if subject_metadata["Diagnosis"] == "Case":
                                self.case_tsvs.append(file_path)
                            elif subject_metadata["Diagnosis"] == "Control":
                                self.cont_tsvs.append(file_path)


    def filter_variants(self, tsv_df, chrom):
        """
        Remove single nucleotide expansions and select by chromosome
        """
        filtered_df = tsv_df[~tsv_df['repeatunit'].isin(['A', 'G', 'C', 'T'])]
        filtered_df = filtered_df[filtered_df['#chrom'] == chrom]

        return filtered_df


    def collapse_sample_tsvs(self, in_dir, out_dir):
        """
        match adjacent variants from a single subject. Aggregate row based on chromosome,
        proximity of the left coordinates and repeat unit and collapse the subject variant
        into the aggregate variant by adding allele2 size estimate size, redefine left and
        right coordinates and average of aggregated variants.
        """
        for file in os.listdir(in_dir):
            if file.endswith('-genotype.txt'):
                tsv = os.path.join(in_dir, file)
                subject, tissue = self.get_metadata_from_filename(tsv)
                df = pd.read_csv(tsv, sep='\t')
                df.sort_values(['#chrom', 'repeatunit', 'left'], inplace=True)
                df.reset_index(drop=True, inplace=True)
                df['counts'] = np.ones(df.shape[0])
                df['mean_left'] = df.left
                collapsed_variants = []

                prev_left = None
                prev_unit = None
                prev_counts = None
                prev_chrom = None
                prev_mean_left = None
                prev_sum_allele2_est = None
                prev_sum_allele1_est = None
                

                for index, row in df.iterrows():

                    # get new variant stats
                    left = row['left']
                    unit = row['repeatunit']
                    count = row['counts']
                    chrom = row['#chrom']
                    mean_left = row['mean_left']
                    sum_allele1_est = row['allele1_est']
                    sum_allele2_est = row['allele2_est']

                    if (  # chromosomes match, repeat units match, within slop range of each other
                        prev_chrom is not None
                        and prev_chrom == chrom
                        and (self.is_rotation(unit, prev_unit) or self.is_rotation(self.rev_complement(unit), prev_unit))
                        and abs(left - prev_left) <= self.slop * self.slop_modifier
                    ):
                        #TODO: REFACTOR THESE REPETITIVE BLOCKS
                        # update counts, mean left and add allele2 size estimate to previous variant
                        new_sum_allele2_est = sum_allele2_est + prev_sum_allele2_est
                        new_sum_allele1_est = sum_allele1_est + prev_sum_allele1_est
                        n = prev_counts + count
                        delta = mean_left - prev_mean_left
                        new_mean = prev_mean_left + (delta * count) / n

                        # set next prev_ values for next iteration
                        prev_left = new_mean
                        prev_counts = n
                        prev_mean_left = new_mean
                        prev_chrom = chrom
                        prev_unit = unit
                        prev_sum_allele2_est = new_sum_allele2_est
                        prev_sum_allele1_est = new_sum_allele1_est

                    # next variant not collapsible, add previous and set prev_ vals for next iter
                    elif prev_chrom is not None:
                        collapsed_variants.append({
                            '#chrom': prev_chrom,
                            'left': int(prev_left),
                            'right': int(prev_left) + 1,
                            'repeatunit': prev_unit,
                            'allele1_est': prev_sum_allele1_est,
                            'allele2_est': prev_sum_allele2_est,
                            'counts': prev_counts,
                            'mean_left': prev_mean_left
                        })
                        prev_left = left
                        prev_unit = unit
                        prev_counts = count
                        prev_chrom = chrom
                        prev_mean_left = mean_left
                        prev_sum_allele2_est = sum_allele2_est
                        prev_sum_allele1_est = sum_allele1_est

                    # first iteration, set prev_ values for next iter
                    else:
                        prev_left = left
                        prev_unit = unit
                        prev_counts = count
                        prev_chrom = chrom
                        prev_mean_left = mean_left
                        prev_sum_allele2_est = sum_allele2_est
                        prev_sum_allele1_est = sum_allele1_est
                    
                # add last variant
                collapsed_variants.append({
                    '#chrom': prev_chrom,
                    'left': int(prev_left),
                    'right': int(prev_left) + 1,
                    'repeatunit': prev_unit,
                    'allele1_est': prev_sum_allele1_est,
                    'allele2_est': prev_sum_allele2_est,
                    'counts': prev_counts,
                    'mean_left': prev_mean_left
                })

                new_df = pd.DataFrame(collapsed_variants)
                new_df['merged_expansions'] = new_df['counts'] > 1
                new_df['left'] = new_df['mean_left']
                new_df.drop(columns=['counts', 'mean_left'], inplace=True)
            
                tsv_file = os.path.basename(tsv) 
                new_df.to_csv(os.path.join(out_dir, tsv_file), sep='\t', index=False)
    
