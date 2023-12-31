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
                slop=100, slop_modifier=1.5, test="KS", motifs_to_drop=[]
    ):
        super().__init__(tsv_dir, csv_metadata, chroms, sex, tissue, dataset, cohort, race, ethnicity, apoe, slop, slop_modifier, test, motifs_to_drop)
        self.cols_to_drop = [
            'allele1_est', 'allele2_est', 'right', 'merged_expansions'
                ]
        self.test_variable = 'allele2_est'
        self.caller = 'STRling'
   
    def filter_tsv_files(self):
        """
        Get case/control TSV file lists from the directory and filter files that adhere to desired covariates described by metadict
        """
        self.tsvs = []
        self.case_tsvs = []
        self.cont_tsvs = []

        cohort_subjects = pd.read_csv(self.csv_metadata).Subject.tolist()
        for file in os.listdir(self.tsv_dir):
            if file.endswith('-genotype.txt'):
                subject, tissue = self.get_metadata_from_filename(file)
                if (subject in cohort_subjects or cohort_subjects == "all_cohorts") and (self.metadict['Tissue'] == None or tissue == self.metadict['Tissue']):
                    subject_metadata = self.get_metadata_from_subject(subject)
                    add_flag = True
                    for key, val in self.metadict.items():
                        if val == "all_apoe": val = None
                        if val == "all_cohorts": val = None
                        if val and val != subject_metadata[key]:
                            print(val, subject_metadata[key])
                            add_flag = False
                            continue
                    if add_flag:
                        file_path = os.path.join(self.tsv_dir, file)
                        if subject_metadata["Diagnosis"] != "Unknown":
                            self.tsvs.append(file_path)
                            if subject_metadata["Diagnosis"] == "Case":
                                self.case_tsvs.append(file_path)
                            elif subject_metadata["Diagnosis"] == "Control":
                                self.cont_tsvs.append(file_path)


    @staticmethod
    def filter_variants(tsv_df, chrom):
        """
        Remove single nucleotide expansions and select by chromosome
        """
        filtered_df = tsv_df[~tsv_df['repeatunit'].isin(['A', 'G', 'C', 'T'])]
        if chrom != 'All':
            filtered_df = filtered_df[filtered_df['#chrom'] == chrom]

        return filtered_df


    def process_tsv_file(self, tsv_file, out_dir):
        """
        Process a single TSV file, collapsing variants.
        """
        subject, tissue = self.get_metadata_from_filename(tsv_file)
        df = pd.read_csv(tsv_file, sep='\t')
        df = self.filter_variants(df, "All")
        df.sort_values(['#chrom', 'repeatunit', 'left'], inplace=True)
        df.reset_index(drop=True, inplace=True)
        df['counts'] = np.ones(df.shape[0])
        df['mean_left'] = df.left
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
        prev_sum_allele2_est = None
        prev_sum_allele1_est = None

        try:
            for index, row in df.iterrows():
                # get new variant stats
                left = row['left']
                unit = row['repeatunit']
                count = row['counts']
                chrom = row['#chrom']
                mean_left = row['mean_left']
                new_allele1_est = row['allele1_est']
                new_allele2_est = row['allele2_est']

                #print(self.is_rotation(unit, prev_unit), self.is_rotation(self.rev_complement(unit), prev_unit))

                if (
                    prev_chrom is not None
                    and prev_chrom == chrom
                    and (self.is_rotation(unit, prev_unit) or self.is_rotation(self.rev_complement(unit), prev_unit))
                    and abs(left - prev_left) <= self.slop * self.slop_modifier
                ):
                    # update counts, mean left, and add allele2 size estimate to the previous variant
                    update_sum_allele2_est = new_allele2_est + prev_sum_allele2_est
                    update_sum_allele1_est = new_allele1_est + prev_sum_allele1_est
                    n = prev_counts + count
                    delta = mean_left - prev_mean_left
                    new_mean = prev_mean_left + (delta * count) / n

                    # set next prev_ values for the next iteration
                    prev_left = new_mean
                    prev_counts = n
                    prev_mean_left = new_mean
                    prev_chrom = chrom
                    prev_unit = unit
                    prev_sum_allele2_est = update_sum_allele2_est
                    prev_sum_allele1_est = update_sum_allele1_est

                # next variant not collapsible, add previous and set prev_ vals for the next iteration
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
                    prev_sum_allele2_est = new_allele2_est
                    prev_sum_allele1_est = new_allele1_est

                # first iteration, set prev_ values for the next iteration
                else:
                    # return if not a string
                    if not isinstance(unit, str):
                        print(unit)
                        return
                    prev_left = left
                    prev_unit = unit
                    prev_counts = count
                    prev_chrom = chrom
                    prev_mean_left = mean_left
                    prev_sum_allele2_est = new_allele2_est
                    prev_sum_allele1_est = new_allele1_est
                    
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

            # Write the result to an output file
            tsv_file_name = os.path.basename(tsv_file)
            # make sure the motif directory exists

            out_file_path = os.path.join(out_dir, tsv_file_name)
            new_df.to_csv(out_file_path, sep='\t', index=False)
        except:
            print("failed to process tsv file: ", tsv_file)
            return
        
    def seperate_tsvs_by_motif(self, tsvs):
        """
        For each tsv, get each unique motif, create a directory for
        that motif if it doesn't exist, filter the tsv by the motif,
        and write the filtered tsv to the motif directory.
        """
        all_motifs = set()
        for tsv in tsvs:
            df = pd.read_csv(tsv, sep='\t')
            motifs = df.repeatunit.unique()
            all_motifs.update(motifs)
            for motif in motifs:
                motif_dir = os.path.join(self.out_dir, motif)
                if not os.path.exists(motif_dir):
                    os.mkdir(motif_dir)
                motif_tsv = df[df.repeatunit == motif]
                motif_tsv.to_csv(os.path.join(motif_dir, os.path.basename(tsv)), sep='\t', index=False)

        # write all found motifs to a file
        motif_df = pd.DataFrame(list(all_motifs), columns=['motif'])
        motif_df.to_csv(os.path.join(f'{self.caller}_motifs.txt'), sep='\t', index=False)

    
    def collapse_sample_tsvs(self, in_dir, out_dir):
        """
        Parallelize the collapsing of sample TSV files.
        """
        print("Collapsing individual STRling TSV files...")
        for motif in os.listdir(in_dir):
            motif_dir = os.path.join(in_dir, motif)
            motif_out_dir = os.path.join(out_dir, motif)
            if not os.path.exists(motif_out_dir):
                os.mkdir(motif_out_dir)
            pool = multiprocessing.Pool()
            pool.starmap(self.process_tsv_file, [(os.path.join(motif_dir, file), motif_out_dir) for file in os.listdir(motif_dir) if file.endswith('-genotype.txt')])
            pool.close()
            pool.join()
            self.filter_tsv_files()
