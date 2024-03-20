import pandas as pd
import numpy as np
from itertools import product
import sys
import os
from scipy.stats import ranksums, ks_2samp, anderson_ksamp
from statsmodels.stats.multitest import multipletests
import multiprocessing
import warnings
import pickle as pkl
import dask.dataframe as dd
from dask.distributed import Client
from modules.NIAGADSExperiment import NIAGADSExperiment
from tqdm.notebook import tqdm
from time import time

# Filter out the "p-value capped" warning
warnings.filterwarnings("ignore", category=UserWarning, message="p-value capped: true value larger than 0.25")

class RepeatsExperiment(NIAGADSExperiment):
    def __init__(
                self, 
                tsv_dir, 
                csv_metadata, 
                chroms="all_chroms", 
                slop=100, 
                slop_modifier=1.5,
                test="KS", 
                motifs_to_drop=None, 
                count_cutoff=10,
                assume_zero=True, 
                seperate_by_motif=False,
                **kwargs):
        
        super().__init__(csv_metadata, **kwargs)
    
        # data inputs
        self.tsv_dir = tsv_dir
        self.csv_metadata = csv_metadata

        # analysis parameters
        self.slop = slop
        self.slop_modifier = slop_modifier
        self.test = test
        self.motifs_to_drop = motifs_to_drop or []
        self.count_cutoff = count_cutoff
        self.assume_zero = assume_zero
        self.seperate_by_motif = seperate_by_motif
        self.aggregation_batch_size = 100

        self.meta1 = pd.DataFrame({
            '#chrom': pd.Series([], dtype='object'),
            'left': pd.Series([], dtype='float64'),
            'mean_left': pd.Series([], dtype='float64'),
            'repeatunit': pd.Series([], dtype='object'),
            'counts': pd.Series([], dtype='int64'),
            'variance_left': pd.Series([], dtype='float64'),
            'm2': pd.Series([], dtype='float64'),
        })

        self.meta2 = pd.DataFrame({
            '#chrom': pd.Series([], dtype='object'),
            'left': pd.Series([], dtype='float64'),
            'repeatunit': pd.Series([], dtype='object'),
            'counts': pd.Series([], dtype='int64'),
            'variance_left': pd.Series([], dtype='float64'),
            'mean_left': pd.Series([], dtype='float64'),
            'left_diff': pd.Series([], dtype='float64'),
        })



        # for subclasses not sure if these are needed
        #self.repeat_nomenclature = ""
        #self.test_variable = None
        #self.caller = "None"
        #self.locus_file_extension = ""
        #self.cols_to_drop = []
        
        chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'] if chroms == "all_chroms" else chroms
        self.chroms = [chroms] if type(chroms) is not list else chroms


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
        return len(str1) == len(str2) and str1 in str2 + str2
    

    def collapse_variants_sep(self, df):
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


    def collapse_variants(self, ddf):
        """
        Dask compatible method to match adjacent variants from a subject and aggregate 
        rows based on chromosome proximity of the left coordinates.
        """
        # Calculate the difference between successive 'left' values within groups
        ddf['left_diff'] = ddf.groupby(['#chrom', 'repeatunit'])['left'].diff().fillna(0)
        
        # Define a helper column to aid in identifying groups of rows to be collapsed
        ddf['group_id'] = (ddf['left_diff'] > self.slop).cumsum()
        
        # Aggregate the counts, mean and variance while also counting the number of samples
        agg_ddf = ddf.groupby(['#chrom', 'repeatunit', 'group_id']).agg({
            '#chrom': 'first',
            'left': 'mean',
            'repeatunit': 'first',
            'counts': 'sum',
            'variance_left': lambda x: np.var(x, ddof=1),
            'mean_left': 'mean',
            'left_diff': 'first',
        }).reset_index(drop=True)

        # Clean up helper columns
        #agg_ddf = agg_ddf.drop(columns=['left_diff', 'group_id'])
        print(agg_ddf.head())

        return agg_ddf

    def summarize_old(self, chrom_motif):
        """
        Aggregate summary statisitcs for each variant in a chromosome and motif
        """
        chrom, motif = chrom_motif

        if self.seperate_by_motif:
            # do not summarize if too few variants found
            if len(self.tsvs[motif]) < self.count_cutoff:
                return    
            # dont process nan motifs (from corrupted tsvs)
            if type(motif) is not str:
                return
            tsv_to_read = self.tsvs[motif][0]
        else:
            tsv_to_read = self.tsvs[0]

        df1 = pd.read_csv(tsv_to_read, sep='\t')
        df1 = df1[df1['#chrom'] == chrom]
        if motif is not None:
            df1 = df1[df1['repeatunit'] == motif]

        tsv_to_read = self.tsvs[motif][1] if self.seperate_by_motif else self.tsvs[1]
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
        tsv_list = self.tsvs[motif][start_idx:] if self.seperate_by_motif else self.tsvs[start_idx:]
        for file in tsv_list:
            files_processed += 1
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
        
        if self.seperate_by_motif:      
            csv_fname = os.path.join(
                "results", "summaries", self.caller, chrom, self.name_tag, 
                f"{self.caller}_summary_{chrom}_{self.name_tag}_{motif}.csv"
            )
        else:
            csv_fname = os.path.join(
                "results", "summaries", self.caller, chrom, self.name_tag, 
                f"{self.caller}_summary_{chrom}_{self.name_tag}.csv"
            )
        dff.to_csv(csv_fname, index=False)
        
        return

    
    def filter_df_by_chrom_and_motif(self, df, chrom, motif):
        """Utility method to filter dataframe by chromosome and motif."""
        return df[(df['#chrom'] == chrom) & 
                ((df['repeatunit'] == motif) if motif is not None else True)]


    def read_and_filter_tsv(self, filename, chrom, motif):
        """Read TSV file and apply filtering by chromosome and motif."""
        try:
            df = dd.read_csv(filename, sep='\t')
            return self.filter_df_by_chrom_and_motif(df, chrom, motif)
        except Exception as e:
            # Handle exceptions according to your logging and error handling strategy
            # Maybe log error message, filename that caused it, and continue or break.
            pass

    def summarize_hmmmm(self, chrom_motif):
        """
        Aggregate summary statistics for each variant in a chromosome and motif.
        """
        chrom, motif = chrom_motif

        # Keep the logic that checks for sufficient file count and correct motif format
        if self.seperate_by_motif:
            if motif not in self.tsvs or len(self.tsvs[motif]) < self.count_cutoff:
                return
            file_paths = self.tsvs[motif]
        else:
            # When not separating by motif, 'motif' will be None, process all files together
            motif = None
            file_paths = self.tsvs

        # Initialize `dff` with the content from the first file
        dff = self.read_and_filter_tsv(file_paths[0], chrom, motif)
        if dff is not None:
            # Add necessary columns and perform initial collapse
            dff["mean_left"] = dff["left"]
            dff["counts"] = np.ones(len(dff))
            dff["variance_left"] = np.zeros(len(dff))
            dff["m2"] = np.zeros(len(dff))
            dff.drop(columns=self.cols_to_drop, inplace=True, errors='ignore')
            dff = self.collapse_variants(dff)
            dff['left'] = dff['mean_left']
            
            # Process the rest of the files, starting from the second one
            for file_path in file_paths[1:]:
                next_df = self.read_and_filter_tsv(file_path, chrom, motif)
                if next_df is not None:
                    next_df["counts"] = np.ones(len(next_df))
                    next_df["variance_left"] = np.zeros(len(next_df))
                    next_df["mean_left"] = next_df["left"]
                    next_df["m2"] = np.zeros(len(next_df))
                    next_df.drop(columns=self.cols_to_drop, inplace=True, errors='ignore')
                    dff = dd.concat([dff, next_df])
                    dff = self.collapse_variants_sep(dff)

            # Conduct the final process on `dff`
            dff = dff[dff['counts'] >= self.count_cutoff].reset_index(drop=True) if 'counts' in dff.columns else dff

            # Save the result if not empty
            if not dff.empty:
                # Save results
                result_dir = os.path.join("results", "summaries", self.caller, chrom, self.name_tag)
                os.makedirs(result_dir, exist_ok=True)
                pattern = f"{self.caller}_summary_{chrom}_{self.name_tag}"
                csv_fname = f"{pattern}_{motif}.csv" if self.seperate_by_motif else f"{pattern}.csv"
                dff.to_csv(os.path.join(result_dir, csv_fname), index=False)


    @staticmethod
    def compute_divisions(ddf, ru_index, max_rows):
        """
        Compute the divisions for the dask dataframe such that divisions are near the max_rows
        and only occur where 'ru_index' changes.

        Parameters:
        - ddf: dask.DataFrame to repartition
        - max_rows: desired maximum number of rows per partition

        Returns: a list of divisions for repartitioning
        """
        # Compute diff to find where 'ru_index' changes
        # Fill NoneType values with a default value (0 in this case)
        ru_index_changes = (ru_index != ru_index.shift()).astype(int)

        # Now, get the change indices
        change_indices = ru_index_changes[ru_index_changes != 0].index
        # drop nonetype values
        change_indices = change_indices.dropna()
        change_indices = change_indices.compute().tolist()
        print("change idx: ", change_indices)

        # Get the actual divisions by finding indexes that are close to desired max_rows
        division_indices = [0]  # First division is always the first index
        current_max = max_rows
        for index in change_indices:
            if index > current_max:
                division_indices.append(index)
                current_max = index + max_rows
        division_indices.append(len(ddf))  # Last division is the last index

        # Convert to divisions which are the index values at division points
        divisions = ddf.index.compute()[division_indices].tolist()

        return divisions


    def summarize(self, chrom_motif):
        """
        Aggregate summary statistics for each variant in a chromosome and motif.
        Uses multiprocessing-safe progress reporting that works for Jupyter notebooks.
        """
        chrom, motif = chrom_motif

        if self.seperate_by_motif:
            if motif not in self.tsvs or len(self.tsvs[motif]) < self.count_cutoff:
                return
            file_paths = self.tsvs[motif]
        else:
            motif = None  # Explicitly handle the non-separated motif case
            file_paths = self.tsvs

        # Initialize a Dask DataFrame for the first file only
        dff = None
        progress_desc = f"Processing {chrom}-{motif if motif else 'All Motifs'}"

        # Iterate over file_paths and construct Dask DataFrame by concatenation
        for idx, file_path in enumerate(file_paths):
            next_df = self.read_and_filter_tsv(file_path, chrom, motif)

            if dff is None:
                dff = next_df
                dff["mean_left"] = dff["left"]
                dff["counts"] = 1
                dff["variance_left"] = 0
                dff["m2"] = 0
                dff = dff.map_partitions(lambda pdf: pdf.drop(columns=self.cols_to_drop, errors='ignore'))
            else:
                next_df["counts"] = 1
                next_df["variance_left"] = 0
                next_df["mean_left"] = next_df["left"]
                next_df["m2"] = 0
                next_df = next_df.map_partitions(lambda pdf: pdf.drop(columns=self.cols_to_drop, errors='ignore'))
                dff = dd.concat([dff, next_df])

            # Update tqdm progress manually for each file
            tqdm.write(f'{progress_desc}: {idx+1}/{len(file_paths)} files processed', end='\r')

        # to collapse the variants we first need to create an index string
        # the best way to do this will be to make the left coordinate a string and add 0's to make it a fixed length
        # then we can concatenate the epeatunit, and the left coordinate such that indexing will give the
        # same result as sorting by repeatunit and left coordinate

        # first extract the maximum length of the left coordinate
        max_len = dff['left'].apply(lambda x: len(str(int(x))), meta=('x', 'int')).max().compute()
        print("max_len: ", max_len)

        # max repeatunit length
        max_ru_len = dff['repeatunit'].apply(lambda x: len(x), meta=('x', 'object')).max().compute()
        print("max_ru_len: ", max_ru_len)

        ru_index = dff['repeatunit'].apply(lambda x: x.ljust(max_ru_len, '0'), meta=('x', 'object'))

        # then create the index string
        dff['index'] = ru_index + "_" + dff['left'].apply(lambda x: str(int(x)).zfill(max_len), meta=('x', 'object'))

        # then reindex the dataframe
        dff = dff.set_index('index')

        # repartition such the the breaks can only occur where the repeat unit part of the index changes
        divisions = self.compute_divisions(dff, ru_index, 20)
        print(divisions)

        print(dff.head())

        # repartitoin such the the breaks can only occur where the repeat unit part of the index changes
        # to do this we need a ballpark figure for the max size of a partition
        # which we will base off of the known memory usage of the dataframe
        # we will set a maximum parititon size of 2 million rows
        # this is a bit of a guess but it should be a good starting point
        max_rows = 2



        # Filter rows with counts below the cutoff (compute necessary if conditions are used)
        # dff = dff[dff['counts'] >= self.count_cutoff].compute()

        # Ensure dff is non-empty and compute to obtain the pandas DataFrame
        #if dff is not None and not dff.map_partitions(len).compute().sum() == 0:
            # Save results
        result_dir = os.path.join("results", "summaries", self.caller, chrom, self.name_tag)
        os.makedirs(result_dir, exist_ok=True)
        pattern = f"{self.caller}_summary_{chrom}_{self.name_tag}"
        csv_fname = f"{pattern}_{motif}.csv" if self.seperate_by_motif else f"{pattern}.csv"
        dff.to_csv(os.path.join(result_dir, csv_fname), single_file=True, index=False)

    
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
            dir_path = os.path.join("results", "summaries", self.caller, chrom, self.name_tag)
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)

        # Create a multiprocessing Pool
        pool = multiprocessing.Pool()

        # Parallelize the summarize function for each chromosome and motif combination
        #results = pool.imap_unordered(self.summarize, [(chrom, motif) for chrom, motif in product(chromosomes, motifs)])
        if self.seperate_by_motif:
            results = pool.imap_unordered(self.summarize, product(chromosomes, self.motifs))
        else:
            results = pool.imap_unordered(self.summarize, product(chromosomes, [None]))

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
                "results", "summaries", self.caller, chrom, self.name_tag,
                 f"{self.caller}_summary_{chrom}_{self.name_tag}_{motif}.csv"
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
                "results", f"full_summary_{self.caller}_{self.name_tag}.csv"
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
        self.WT_df.to_csv(f"results/{self.caller}_{self.test}_{chroms}_{self.name_tag}.csv", index=False)

        return

    
    def schedule_filter_tsv_files(self):
        """
        Filter TSV files for each motif in parallel using multiprocessing
        """
        print("Filtering TSV files...")

        self.motifs = pd.read_csv(os.path.join(f'{self.caller}_motifs.txt'), sep='\t').iloc[:,0].dropna().tolist()
        self.tsvs = {}
        self.case_tsvs = {}
        self.cont_tsvs = {}

        # Parallelize the summarize function for each motif
        if self.seperate_by_motif:
            # Create a multiprocessing Pool
            pool = multiprocessing.Pool()
            results = pool.imap_unordered(self.filter, self.motifs, self.tsv_dir)

            # store results in class variables for summarizing function
            for result in results:
                self.tsvs[result[0]] = result[1]
                self.case_tsvs[result[0]] = result[2]
                self.cont_tsvs[result[0]] = result[3]
        
            pool.close()
            pool.join()

        else: # no parallelization if not seperating by motif
            _, self.tsvs, self.case_tsvs, self.cont_tsvs = self.filter("no_seperation", self.tsv_dir)
    
        # write self.tsvs, self.case_tsvs, and self.cont_tsvs dictionaries to pickles
        with open(f"{self.caller}_{self.name_tag}_tsvs.txt", 'wb') as handle:
            pkl.dump(self.tsvs, handle)
        with open(f"{self.caller}_{self.name_tag}_case_tsvs_{self.name_tag}.txt", 'wb') as handle:
            pkl.dump(self.case_tsvs, handle)
        with open(f"{self.caller}_{self.name_tag}_cont_tsvs.txt", 'wb') as handle:
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
        if not self.seperate_by_motif:
            print("Seperate by motif is set to False. Skipping...")
            return

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

    
    def schedule_collapse_tsvs(self, dir, out_dir):
        """
        Parallelize the collapsing of sample TSV files.
        """
        print("Collapsing individual STRling TSV files...")

        if not os.path.exists(out_dir):
            os.makedirs(out_dir, exist_ok=True)

        pool = multiprocessing.Pool()  # Create a single pool here

        if self.seperate_by_motif:

            for motif in os.listdir(dir):
                motif_dir = os.path.join(dir, motif)
                motif_out_dir = os.path.join(out_dir, motif)

                if not os.path.exists(motif_out_dir):
                    os.makedirs(motif_out_dir, exist_ok=True)

                pool.starmap(self.collapse_tsv_file, [(os.path.join(motif_dir, file), motif_out_dir) for file in os.listdir(motif_dir) if file.endswith(self.locus_file_extension)])

        else:
            print(dir)
            pool.starmap(self.collapse_tsv_file, [(os.path.join(dir, file), out_dir) for file in os.listdir(dir) if file.endswith(self.locus_file_extension)])
        
        pool.close()  # Close and join the pool after all tasks are submitted
        pool.join()


    def create_array(self, test_file):
        """
        Create a numpy array from the test file
        Converts categorical variables to integers
        """
        pass
        


