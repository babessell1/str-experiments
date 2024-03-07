import pandas as pd
import numpy as np
from itertools import zip_longest, product
import sys
import os
from scipy.stats import ranksums, ks_2samp, anderson_ksamp
from statsmodels.stats.multitest import multipletests
import multiprocessing
import warnings
import pickle as pkl
import dask.dataframe as dd

# Filter out the "p-value capped" warning
warnings.filterwarnings("ignore", category=UserWarning, message="p-value capped: true value larger than 0.25")

#
class RepeatsExperiment:
    def __init__(self, tsv_dir, csv_metadata, chroms="All", sex=None, tissue=None,
                 dataset=None, cohort=None, race=None, ethnicity=None, apoe=None, 
                 slop=100, slop_modifier=1.5, test="KS", motifs_to_drop=[], count_cutoff=10,
                 assume_zero=True, repeat_nomenclature="", cores_per_node=36, mem_per_node="128GB",
                 partition="", account="", nodes=1, walltime="4:00:00", dask_log_directory="dask_logs",
                 corrupt_file_log=None
    ):
        self.test_variable = None
        self.tsv_dir = tsv_dir
        self.csv_metadata = csv_metadata
        self.slop = slop
        self.slop_modifier = slop_modifier
        self.test = test
        self.caller = "None"
        self.motifs_to_drop = motifs_to_drop
        self.count_cutoff = count_cutoff
        self.assume_zero = assume_zero
        self.repeat_nomenclature = repeat_nomenclature

        # cluster parameters
        self.cores_per_node = cores_per_node
        self.mem_per_node = mem_per_node
        self.partition = partition
        self.account = account
        self.nodes = nodes
        self.walltime = walltime
        self.dask_log_directory = dask_log_directory

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
        
        # Default to all chromosomes if none specified
        if chroms is None:
            self.chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        else:
            self.chroms = chroms if isinstance(chroms, list) else [chroms]
            
        if type(chroms) is not list:
            chroms = [chroms]
        
        self.chroms = chroms
        self.cohort = cohort

        # coerce apoe to string of int if possible
        if apoe is not None:
            if type(apoe) is not str:
                apoe = str(int(apoe))

        self.apoe = apoe
  
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
        Extract metadata from filename
        Arg: filename (str)
        Returns: subject (str), tissue (str)

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
        """
        Reverse complement a DNA sequence
        """
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

        for _, row in df.iterrows():
            
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
        # DEPRECATED
        tsv_split = tsv.split("/")
        tsv_split.insert(-1, motif)
        return "/".join(tsv_split)


    def summarize(self, chrom_motif):
        """
        Aggregate summary statisitcs for each variant in a chromosome and motif
        """
        # get subset of expanded variants for the chromosome and motif
        #tsv_to_read = self.add_motif_dir_to_tsv_path(self.tsvs[0], motif)
        chrom, motif = chrom_motif

        # do not summarize if too few variants found
        if len(self.tsvs[motif]) < self.count_cutoff:
            return
             
        # dont process nan motifs (from corrupted tsvs)
        if type(motif) is not str:
            return
        
        tsv_to_read = self.tsvs[motif][0]
        df1 = pd.read_csv(tsv_to_read, sep='\t')
        df1 = df1[df1['#chrom'] == chrom]
        if motif is not None:
            df1 = df1[df1['repeatunit'] == motif]


        #tsv_to_read = self.add_motif_dir_to_tsv_path(self.tsvs[1], motif)
        tsv_to_read = self.tsvs[motif][1]
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
        for file in self.tsvs[motif][start_idx:]:
            files_processed += 1
            #tsv_to_read = self.add_motif_dir_to_tsv_path(file, motif)
            tsv_to_read = file
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
                raise ValueError("Counts are higher than the current iteration.")
            
            iteration += 1

        # remove last row if last row is empty or nan
        if dff.iloc[-1].isnull().all() or dff.iloc[-1].isna().all():
            dff = dff.iloc[:-1]

        # remove rows with counts less than cutoff
        dff = dff[dff['counts'] >= self.count_cutoff]

        # if no variants found, return
        if dff.shape[0] == 0:
            return
        
        dff.to_csv(os.path.join(
            "results", "summaries", self.caller, chrom, self.cohort, self.apoe, 
            f"{self.caller}_summary_{chrom}_{self.cohort}_{self.apoe}_{motif}.csv"
        ), index=False)
        
        return

    
    def schedule_summarizing(self):
        """
        Handle the multiprocessing of the summarize function for each chromosome and motif combination
        """
        print("Creating Summary file...")
        if not hasattr(self, 'tsvs'): # if self.tsvs does not exist, get it from the pickled files
            with open(f"{self.caller}_tsvs.txt", 'rb') as handle:
                self.tsvs = pkl.load(handle)
            with open(f"{self.caller}_case_tsvs.txt", 'rb') as handle:
                self.case_tsvs = pkl.load(handle)
            with open(f"{self.caller}_cont_tsvs.txt", 'rb') as handle:
                self.cont_tsvs = pkl.load(handle)

        self.motifs = pd.read_csv(os.path.join(f'{self.caller}_motifs.txt'), sep='\t').iloc[:,0].dropna().tolist()
        chromosomes = self.chroms

        # create directories needed for summary funciton csv files
        for chrom in chromosomes:
            dir_path = os.path.join("results", "summaries", self.caller, chrom, self.cohort, self.apoe)
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)

        # Create a multiprocessing Pool
        pool = multiprocessing.Pool()

        # Parallelize the summarize function for each chromosome and motif combination
        #results = pool.imap_unordered(self.summarize, [(chrom, motif) for chrom, motif in product(chromosomes, motifs)])
        results = pool.imap_unordered(self.summarize, product(chromosomes, self.motifs))

        # Iterate over results to start the computation
        for _ in results:
            pass

        # Close the pool and wait for all processes to finish
        pool.close()
        pool.join()

        return None

    def par_process_variant(self, iter):
        """
        Match each variant in a summary df to all corresponding variants in each subject df to extract case and control allele2 estimate sizes for each variant (defined by a row in the summary df). This can be done in parallel
        """
        variant_row = iter[1]
        
        if variant_row['counts'] < self.count_cutoff:
            return None

        variant = variant_row['repeatunit']
        mean_left = variant_row['mean_left']
        std_left = np.sqrt(variant_row['variance_left'])
        range_ = self.slop * self.slop_modifier
        chrom = variant_row['#chrom']
        motif = variant_row['repeatunit']
        matching_variants = []
        matching_filenames = []
        add_similar = True  # if true, when multiple expansions in same STR region, add together. Otherwise skip the sample
        warn_flag = False 
        cnt_multi = 0

        # get all tsvs for the motif and chromosome
        # get files in tsv_dir/caller/motif
        tsvs_to_check = os.listdir(os.path.join(self.tsv_dir, motif))

        for tsv_file in tsvs_to_check:
            tsv_df = pd.read_csv(os.path.join(self.tsv_dir, motif, tsv_file), sep='\t')
            # subset to only the variants that match the current variant and same chromosome
            variant_df = tsv_df[(tsv_df['repeatunit'] == variant ) & (tsv_df['#chrom'] == chrom)]

            add_variants = variant_df[
                (variant_df['left'] >= mean_left - range_ ) &
                (variant_df['left'] <= mean_left + range_ )
            ][self.test_variable].tolist()
            
            if len(add_variants) > 1:
                warn_flag = True
                cnt_multi+=1
                if add_similar:
                    add_variants = np.sum(add_variants)
                    matching_variants.append(add_variants)                               
            elif len(add_variants) > 0:
                matching_variants.append(*add_variants)
                matching_filenames.append(os.path.join(self.tsv_dir, motif, tsv_file))
            elif self.assume_zero:
                matching_variants.append(0) # assume zero (same as ref) if no matching variants
                matching_filenames.append(os.path.join(self.tsv_dir, motif, tsv_file))
            else:
                matching_variants.append(np.nan)
                
        case_vals = [matching_variants[idx] for idx, file in enumerate(matching_filenames) if file in self.case_tsvs[motif]]
        cont_vals = [matching_variants[idx] for idx, file in enumerate(matching_filenames) if file in self.cont_tsvs[motif]]

        total_found = len([val for val in case_vals if np.abs(val)>0] + [val for val in cont_vals if np.abs(val>0)])

        # dont bother testing if the amount of expansions is low
        if sum(val > 0 for val in case_vals) + sum(val > 0 for val in cont_vals) < self.count_cutoff:
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
    
    def schedule_perform_stat_test(self, sum_file=None):
        """
        Handle the multiprocessing of the statistical test (par_process_variant) for each chromosome and motif combination
        """
        
        # if self.motifs does not exist, create it from the motifs.txt file
        if not hasattr(self, 'motifs'):
            self.motifs = pd.read_csv(os.path.join(f'{self.caller}_motifs.txt'), sep='\t').iloc[:,0].dropna().tolist()

        # get all combinations of motifs and chromosomes
        files_to_process = [
            os.path.join(
                "results", "summaries", self.caller, chrom, self.cohort, self.apoe,
                 f"{self.caller}_summary_{chrom}_{self.cohort}_{self.apoe}_{motif}.csv"
            ) for chrom, motif in product(self.chroms, self.motifs)
        ]

        # remove files that don't exist from the list
        files_to_process = [file for file in files_to_process if os.path.exists(file)]

        if sum_file is None:

            # create a dataframe from first file
            sum_df = pd.read_csv(files_to_process[0])
            # drop rows with counts less than cutoff
            sum_df = sum_df[sum_df['counts'] >= self.count_cutoff]
            
            # and then append the rest of the files
            for file in files_to_process[1:]:
                next_df = pd.read_csv(file)
                next_df = next_df[next_df['counts'] >= self.count_cutoff]
                sum_df = pd.concat([sum_df, next_df], ignore_index=True)

            # write the combined dataframe to a csv
            sum_df.to_csv(os.path.join(
                "results", f"full_summary_{self.caller}_{self.cohort}_{self.apoe}.csv"
            ), index=False)

        else:
            # read the file into df
            sum_df = pd.read_csv(sum_file)

        significant_variants = []

        rows2process = sum_df.shape[0]

        pool = multiprocessing.Pool()  # Create a multiprocessing Pool
        
        cnt = 0
        for result in pool.imap_unordered(self.par_process_variant, [(i, row) for (i, row) in sum_df.iterrows()]):
        # The rest of the code remains the same
            cnt += 1
            print(f"Progress: {cnt}/{rows2process}", end='\r')
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

        pool.close()  # Close the multiprocessing Pool
        pool.join()   # Wait for all processes to finish

        # Perform multiple testing correction
        p_values = [variant['p_value'] for variant in significant_variants]

        print("P values before correction:")
        print(p_values)

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
    
    def filter(self, motif):
        """
        Get case/control TSV file lists from the directory and filter files that adhere to desired covariates described by metadict
        """
        tsvs = []
        case_tsvs = []
        cont_tsvs = []
 
        cohort_subjects = pd.read_csv(self.csv_metadata).Subject.tolist()

        if type(motif) is not str:  # handle nans (may want to write csv to log file so we can rerun these samples)
            return tsvs, case_tsvs, cont_tsvs

        for file in os.listdir(os.path.join(self.tsv_dir, motif)):
            if file.endswith(self.locus_file_extension):
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
                        file_path = os.path.join(self.tsv_dir, motif, file)
                        if subject_metadata["Diagnosis"] != "Unknown":
                            tsvs.append(file_path)
                            if subject_metadata["Diagnosis"] == "Case":
                                case_tsvs.append(file_path)
                            elif subject_metadata["Diagnosis"] == "Control":
                                cont_tsvs.append(file_path)
            
        return motif, tsvs, case_tsvs, cont_tsvs
    
    def schedule_filter_tsv_files(self):
        """
        Filter TSV files for each motif in parallel using multiprocessing
        """
        print("Filtering TSV files...")

        self.motifs = pd.read_csv(os.path.join(f'{self.caller}_motifs.txt'), sep='\t').iloc[:,0].dropna().tolist()
        self.tsvs = {}
        self.case_tsvs = {}
        self.cont_tsvs = {}

        # Create a multiprocessing Pool
        pool = multiprocessing.Pool()

        # Parallelize the summarize function for each motif
        results = pool.imap_unordered(self.filter, self.motifs)

        # store results in class variables for summarizing function
        for result in results:
            self.tsvs[result[0]] = result[1]
            self.case_tsvs[result[0]] = result[2]
            self.cont_tsvs[result[0]] = result[3]
        
        pool.close()
        pool.join()

        # write self.tsvs, self.case_tsvs, and self.cont_tsvs dictionaries to pickles
        with open(f"{self.caller}_{self.apoe}_tsvs.txt", 'wb') as handle:
            pkl.dump(self.tsvs, handle)
        with open(f"{self.caller}_{self.apoe}_case_tsvs_{self.apoe}.txt", 'wb') as handle:
            pkl.dump(self.case_tsvs, handle)
        with open(f"{self.caller}_{self.apoe}_cont_tsvs.txt", 'wb') as handle:
            pkl.dump(self.cont_tsvs, handle)
    

    def seperate_tsv_by_motif(self, tsv):
        df = pd.read_csv(tsv, sep='\t')
        motifs = df[self.repeat_nomenclature].unique()
        motifs_to_remove = []
        for motif in motifs:
            if len(motif) > 9:
                motifs_to_remove.append(motif)
                continue
            try:
                motif_dir = os.path.join(self.seperated_dir, motif)
            except:
                warnings.warn(f"Motif {motif} is not a valid directory name. Corrupted file: {tsv}")
                continue

            if not os.path.exists(motif_dir):
                os.makedirs(motif_dir, exist_ok=True)

            motif_tsv = df[df[self.repeat_nomenclature] == motif]
            motif_tsv.to_csv(os.path.join(motif_dir, os.path.basename(tsv)), sep='\t', index=False)
        
        motifs = set(motifs) - set(motifs_to_remove)

        return list(motifs)


    def schedule_seperate_tsvs_by_motif(self, in_dir, seperated_dir):
        """
        For each tsv, get each unique motif, create a directory for
        that motif if it doesn't exist, filter the tsv by the motif,
        and write the filtered tsv to the motif directory.
        """
        print("Seperating TSVs by motif...")

        if not os.path.exists(seperated_dir):
            os.makedirs(seperated_dir, exist_ok=True)
        
        self.seperated_dir = seperated_dir

        all_motifs = set()
        tsvs = [os.path.join(in_dir, file) for file in os.listdir(in_dir) if file.endswith(self.locus_file_extension)]
        pool = multiprocessing.Pool()
        # run seperate_tsvs_by_motif for each tsv in parallel and return the set of all motifs
        for motifs in pool.imap_unordered(self.seperate_tsv_by_motif, tsvs):
            all_motifs.update(motifs)

        pool.close()
        pool.join()

        # write all found motifs to a file
        motif_df = pd.DataFrame(list(all_motifs), columns=['motif'])
        motif_df.to_csv(os.path.join(f'{self.caller}_motifs.txt'), sep='\t', index=False, header=False)

    
    def schedule_collapse_tsvs(self, seperated_dir, out_dir):
        """
        Parallelize the collapsing of sample TSV files.
        """
        print("Collapsing individual STRling TSV files...")

        if not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        pool = multiprocessing.Pool()  # Create a single pool here

        for motif in os.listdir(seperated_dir):
            motif_dir = os.path.join(seperated_dir, motif)
            motif_out_dir = os.path.join(out_dir, motif)

            if not os.path.exists(motif_out_dir):
                os.makedirs(motif_out_dir, exist_ok=True)

            pool.starmap(self.collapse_tsv_file, [(os.path.join(motif_dir, file), motif_out_dir) for file in os.listdir(motif_dir) if file.endswith(self.locus_file_extension)])

        pool.close()  # Close and join the pool after all tasks are submitted
        pool.join()


    def create_array(self, test_file):
        """
        Create a numpy array from the test file
        Converts categorical variables to integers
        """
        


