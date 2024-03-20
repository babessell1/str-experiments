import pandas as pd
import numpy as np
from itertools import zip_longest, product
import os
from scipy.stats import ranksums, ks_2samp, anderson_ksamp
from statsmodels.stats.multitest import multipletests
import warnings
import dask.dataframe as dd
from dask.distributed import Client, progress, performance_report
from dask import delayed
from dask_jobqueue import SLURMCluster
import dask.bag as db
import dask.array as da
import pickle as pkl

# Filter out the "p-value capped" warning
warnings.filterwarnings("ignore", category=UserWarning, message="p-value capped: true value larger than 0.25")
# filter  UserWarning: Sending large graph of size 22.63 MiB. This may cause some slowdown. Consider scattering data ahead of time and using futures. warnings.warn(

warnings.filterwarnings("ignore", category=UserWarning, message="Sending large graph of size")


class RepeatsExperiment:
    def __init__(self, 
                tsv_dir, 
                csv_metadata, 
                chroms="all_chroms", 
                sex=None, 
                tissue=None,
                dataset=None, 
                cohort="all_cohorts", 
                race=None, 
                ethnicity=None, 
                apoe="all_apoe",
                slop=100, 
                slop_modifier=1.5,
                test="KS", 
                motifs_to_drop=None, 
                count_cutoff=10,
                assume_zero=True, 
                cores_per_node=36, 
                mem_per_node="128GB",
                partition="", 
                account="", 
                nodes=1, 
                walltime="12:00:00",
                dask_log_directory="/home/bbessell/str-analysis/dask_logs", 
                corrupt_file_log = "corrupt_files.log"
                ):
    
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

        # cluster parameters
        self.cores_per_node = cores_per_node
        self.mem_per_node = mem_per_node
        self.partition = partition
        self.account = account
        self.nodes = nodes
        self.walltime = walltime
        self.dask_log_directory = dask_log_directory

        # from subclasses
        self.repeat_nomenclature = ""
        self.test_variable = None
        self.caller = "None"
        self.locus_file_extension = ""
        self.cols_to_drop = []

        # cohort parameters
        self.metadict = {
            "Sex": sex,
            "Tissue": tissue,
            "Dataset": dataset,
            "Cohort": cohort,
            "Race": race,
            "Ethnicity": ethnicity,
            "APOE": apoe,
        }

        # file init
        self.corrupt_file_log = corrupt_file_log
        open(self.corrupt_file_log, 'w').close()
        
        if chroms=="all_chroms":
            self.chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'] if chroms in (None, "all_chroms") else chroms
           
        if type(chroms) is not list:
            self.chroms = [chroms]
        

    def _init_manifest(self):
        """
        Create a manifest file for the experiment from the NIAGADS metadata file
        """
        manifest_df = pd.read_csv(self.csv_metadata)
        if self.metadict["Cohort"] != "all_cohorts":
            manifest_df = manifest_df[manifest_df['Cohort'] == self.metadict["Cohort"]]
        if self.metadict["APOE"] != "all_apoe":
            manifest_df = manifest_df[manifest_df['APOE'] == self.metadict["APOE"]]
        if self.sex is not None:
            manifest_df = manifest_df[manifest_df['Sex'] == self.metadict['Sex']]
        if self.race is not None:
            manifest_df = manifest_df[manifest_df['Race'] == self.metadict['Race']]
        if self.ethnicity is not None:
            manifest_df = manifest_df[manifest_df['Ethnicity'] == self.metadict['Ethnicity']]
        if self.tissue is not None:
            # tissue is not in the actual manifest, so we need to search the filename
            # blood will have the string -BL-, brain will have -BR-
            if self.metadict["Tissue"] == "Blood" or self.metadict["Tissue"] == "BL":
                tissue = "BL"
            elif self.metadict["Tissue"] == "Brain" or self.metadict["Tissue"] == "BR":
                tissue = "BR"
            else:
                raise ValueError("Tissue not recognized. Must be 'Blood' or 'Brain'")
            manifest_df = manifest_df[manifest_df['location'].str.contains("-" + tissue + "-")]
        
        self.csv_metadata = os.path.join("manifests", os.path.basename(self.csv_metadata).split(".")[0] + "_" + self.metadict["Cohort"] + ".csv")
        manifest_df.to_csv(self.csv_metadata, index=False)
        self.meta_df = manifest_df


    @staticmethod
    def get_metadata_from_filename(filename):
        """
        Pull subject ID and tissue type from the filename (NIAGADS convention)
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
        return len(str1) == len(str2) and str1 in str2 + str2
    

    def start_slurm_cluster(self):
        """
        Start a SLURM cluster using the parameters specified in the constructor
        """
        return SLURMCluster(
            cores=self.cores_per_node,
            memory=self.mem_per_node,
            queue=self.partition,
            account=self.account,
            log_directory=self.dask_log_directory,
            walltime=self.walltime,
            interface="ib0",
            processes=self.cores_per_node,
            #job_extra_directives=[f'--nodes={self.nodes}'],
            death_timeout=1*24*60,  # 1 day
        )
    

    def get_metadata_from_subject(self, subject):
        """
        Match subject name to metadata file and pull the rest of the useful metadata as a dictionary.
        Assumes that self.meta_df is a DataFrame property of the class.
        """
        subject_row = self.meta_df[self.meta_df['Subject'] == subject]
        if subject_row.empty:
            # Handle case where subject is not found
            raise ValueError(f"Subject {subject} not found in metadata.")

        subject_data = subject_row.iloc[0]

        subject_metadata = {field: subject_data[field] for field in [
            'Dataset', 'Disease', 'Cohort', 'Sex', 'Race', 'Ethnicity', 'Diagnosis', 'Assay', 'APOE'
        ]}

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
    

    def summarize_chrom_motif(self, chrom_motif):
        """
        Aggregate summary statisitcs for each variant in a chromosome and motif
        """
        # get subset of expanded variants for the chromosome and motif
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
 
        # remove last row if last row is empty or nan
        if dff.iloc[-1].isnull().all() or dff.iloc[-1].isna().all():
            dff = dff.iloc[:-1]

        # remove rows with counts less than cutoff
        dff = dff[dff['counts'] >= self.count_cutoff]

        # if no variants found, return
        if dff.shape[0] == 0:
            return
        
        dff.to_csv(os.path.join(
            "results", "summaries", self.caller, chrom, self.metadict["Cohort"], self.metadict["APOE"], 
            f"{self.caller}_summary_{chrom}_{self.metadict["Cohort"]}_{self.metadict["APOE"]}_{motif}.csv"
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
            dir_path = os.path.join("results", "summaries", self.caller, chrom, self.metadict["Cohort"], self.metadict["APOE"])
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)

        # Create a SLURMCluster
        print("Waiting for workers...")
        cluster = self.start_slurm_cluster()
        cluster.scale(self.nodes)
        client = Client(cluster)
        #client.wait_for_workers(self.nodes)

        # Parallelize the summarize function for each chromosome and motif combination
        #results = pool.imap_unordered(self.summarize, [(chrom, motif) for chrom, motif in product(chromosomes, motifs)])
        #results = pool.imap_unordered(self.summarize, product(chromosomes, self.motifs))

        tasks = [delayed(self.summarize_chrom_motif)(chrom_motif) for chrom_motif in product(chromosomes, self.motifs)]

        futures = client.compute(tasks)
        progress(futures, notebook=False)

        with performance_report(filename=f"{self.dask_log_directory}/summarize_perf.html"):
            client.gather(futures)

        print("\nClosing cluster...")
        client.close()
        cluster.close()

        return

    def perform_stat_test(self, iter):
        """
        Match each variant in a summary df to all corresponding variants in each subject df to extract case and control allele2 estimate sizes for each variant (defined by a row in the summary df). This can be done in parallel
        """
        i = iter[0]
        variant_row = iter[1]
        variant = variant_row['repeatunit']
        mean_left = variant_row['mean_left']
        std_left = da.sqrt(variant_row['variance_left'])
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
            tsv_df = dd.read_csv(os.path.join(self.tsv_dir, motif, tsv_file), sep='\t')
            # subset to only the variants that match the current variant and same chromosome
            variant_df = tsv_df[(tsv_df['repeatunit'] == variant ) & (tsv_df['#chrom'] == chrom)]

            add_variants = variant_df[
                (variant_df['left'] >= mean_left - range_ ) &
                (variant_df['left'] <= mean_left + range_ )
            ][self.test_variable].compute().tolist()

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
                ][self.test_variable].compute().tolist()
                #for row in add_variants_sorted:
                #    #print(row)
                #    pass
            
            if len(add_variants) > 1:
                #print("WARNING: MORE THAN ONE MATCHING VARIANT FOUND IN SINGLE SUBJECT")
                warn_flag = True
                cnt_multi+=1
                if add_similar:
                    add_variants = da.sum(add_variants)
                    matching_variants.append(add_variants)
                else:
                    continue
                                                
            elif len(add_variants) > 0:
                matching_variants.append(*add_variants)
                matching_filenames.append(os.path.join(self.tsv_dir, motif, tsv_file))
            elif self.assume_zero:
                matching_variants.append(0) # assume zero (same as ref) if no matching variants
                matching_filenames.append(os.path.join(self.tsv_dir, motif, tsv_file))

        case_vals = [matching_variants[idx] for idx, file in enumerate(matching_filenames) if file in self.case_tsvs[motif]]
        cont_vals = [matching_variants[idx] for idx, file in enumerate(matching_filenames) if file in self.cont_tsvs[motif]]

        total_found = len([val for val in case_vals if da.abs(val)>0] + [val for val in cont_vals if da.abs(val>0)])

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
    

    def schedule_perform_stat_test(self):
        """
        Handle the multiprocessing of the statistical test (par_process_variant) for each chromosome and motif combination
        """

        # if self.motifs does not exist, create it from the motifs.txt file
        if not hasattr(self, 'motifs'):
            self.motifs = pd.read_csv(os.path.join(f'{self.caller}_motifs.txt'), sep='\t').iloc[:,0].dropna().tolist()

        if not hasattr(self, 'tsvs'): # if self.tsvs does not exist, get it from the pickled files
            with open(f"{self.caller}_tsvs.txt", 'rb') as handle:
                self.tsvs = pkl.load(handle)
            with open(f"{self.caller}_case_tsvs.txt", 'rb') as handle:
                self.case_tsvs = pkl.load(handle)
            with open(f"{self.caller}_cont_tsvs.txt", 'rb') as handle:
                self.cont_tsvs = pkl.load(handle)

        # get all combinations of motifs and chromosomes
        files_to_process = [
            os.path.join(
                "results", "summaries", self.caller, chrom, self.metadict["Cohort"], self.metadict["APOE"],
                 f"{self.caller}_summary_{chrom}_{self.metadict["Cohort"]}_{self.metadict["APOE"]}_{motif}.csv"
            ) for chrom, motif in product(self.chroms, self.motifs)
        ]

        # remove files that don't exist from the list
        files_to_process = [file for file in files_to_process if os.path.exists(file)]

        significant_variants = []

        # Create a SLURMCluster
        print("Waiting for workers...")
        cluster = self.start_slurm_cluster()
        cluster.scale(self.nodes)
        client = Client(cluster)
        #client.wait_for_workers(self.nodes)

        number_files = len(files_to_process)

        for file in files_to_process:

            # make a progress bar
            print(f"\nProcessing file {files_to_process.index(file)+1}/{number_files}", end='\r')
        
            summary_df = dd.read_csv(file)
            
            tasks = [delayed(self.perform_stat_test)(iter) for iter in summary_df.iterrows()]

            futures = client.compute(tasks)
            #progress(futures, notebook=False)

            with performance_report(filename=f"{self.dask_log_directory}/perform_stat_test_perf.html"):
                results = client.gather(futures)
            
            for result in results:
                if result is not None:
                    significant_variants.append(result)
            
        print("\nClosing cluster...")
        client.close()
        cluster.close()

        # Perform multiple testing correction
        p_values = [variant['p_value'] for variant in significant_variants]

        reject, p_corrected, _, _ = multipletests(p_values, method='fdr_bh')

        for i, variant in enumerate(significant_variants):
            variant['p_corrected'] = p_corrected[i]
            variant['reject_null'] = reject[i]

        if len(self.chroms) == 1:
            chroms = self.chroms[0]
        else:
            chroms = "_".join(self.chroms)

        self.WT_df = pd.DataFrame(significant_variants)
        self.WT_df.sort_values('p_corrected', ascending=True)
        self.WT_df.to_csv(f"results/{self.caller}_{self.test}_{chroms}_{self.metadict["Cohort"]}_{self.metadict["APOE"]}.csv", index=False)

        return
    

    def filter_tsvs(self, motif):
        """
        Get motif seperated case/control TSV file lists from the directory and filter files
        that adhere to desired covariates described by metadict
        """
        tsvs = []
        case_tsvs = []
        cont_tsvs = []
 
        cohort_subjects = dd.read_csv(self.csv_metadata).Subject.compute().tolist()

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
                            #print(val, subject_metadata[key])
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
        Schedule filtering of TSV files for each motif in parallel using multiprocessing
        """
        print("Filtering TSV files...")

        self.motifs = pd.read_csv(os.path.join(f'{self.caller}_motifs.txt'), sep='\t').iloc[:,0].dropna().tolist()
        self.tsvs = {}
        self.case_tsvs = {}
        self.cont_tsvs = {}

        # Create a SLURMCluster
        print("Waiting for workers...")
        cluster = self.start_slurm_cluster()
        cluster.scale(self.nodes)
        client = Client(cluster)
        #client.wait_for_workers(self.nodes)

        # create dask bag of motifs
        motifs_bag = db.from_sequence(self.motifs)

        # map filter function to each motif
        futures = motifs_bag.map(self.filter_tsvs).persist()

        # Print progress to the user
        progress(futures, notebook=False)

        # Generate performance report
        with performance_report(filename=f"{self.dask_log_directory}/filter_tsvs_perf.html"):
            results = client.compute(futures, sync=True)

        # TODO: find a more memory efficient way to do this
        # store results in class variables for summarizing function

        for result in results:
            self.tsvs[result[0]] = result[1]
            self.case_tsvs[result[0]] = result[2]
            self.cont_tsvs[result[0]] = result[3]

        # write self.tsvs, self.case_tsvs, and self.cont_tsvs dictionaries to pickles
        with open(f"{self.caller}_tsvs.txt", 'wb') as handle:
            pkl.dump(self.tsvs, handle)
        with open(f"{self.caller}_case_tsvs.txt", 'wb') as handle:
            pkl.dump(self.case_tsvs, handle)
        with open(f"{self.caller}_cont_tsvs.txt", 'wb') as handle:
            pkl.dump(self.cont_tsvs, handle)

        # close cluster
        print("\nClosing cluster...")
        client.close()
        cluster.close()


    def seperate_tsv_by_motif(self, tsv):
        """
        Seperate a TSV file by motif, writing each motif to a seperate file
        """

        df = dd.read_csv(tsv, sep='\t', assume_missing=True)

        # return empty list if df has missing values
        # remember to use compute() because lazy evaluation
        #if df.isnull().any().any().compute():
        #    return []

        motifs = df[self.repeat_nomenclature].unique().compute()

        for motif in motifs:
            try:
                motif_dir = os.path.join(self.seperated_dir, motif)
            except:
                warnings.warn(f"Motif {motif} is not a valid directory name. Corrupted file: {tsv}")
                with open(self.corrupt_file_log, 'a') as handle:
                    handle.write(f"{tsv}\n")
                continue

            if not os.path.exists(motif_dir):
                os.makedirs(motif_dir, exist_ok=True)

            motif_tsv = df[df[self.repeat_nomenclature] == motif]
            #print(f"Writing {os.path.join(motif_dir, os.path.basename(tsv))}")
            motif_tsv.to_csv(os.path.join(motif_dir, os.path.basename(tsv)), sep='\t', index=False, single_file=True)

        return motifs


    def schedule_seperate_tsvs_by_motif(self, in_dir, seperated_dir):
        """
        Create a run SLURM client t paralleize for each tsv, 
        getting each unique motif, creating a directory for
        that motif if it doesn't exist, filtering the tsv by the motif,
        and writing the filtered tsv to the motif directory.
        """
        print("Seperating TSVs by motif...")

        if not os.path.exists(seperated_dir):
            os.makedirs(seperated_dir, exist_ok=True)
        
        self.seperated_dir = seperated_dir

        all_motifs = set()
        tsvs = [os.path.join(in_dir, file) for file in os.listdir(in_dir) if file.endswith(self.locus_file_extension)]

        print("Waiting for workers...")
        cluster = self.start_slurm_cluster()
        cluster.scale(self.nodes)
        client = Client(cluster)
        #client.wait_for_workers(self.nodes)

        futures = client.map(self.seperate_tsv_by_motif, tsvs)

        # Print progress to the user
        progress(futures, notebook=False)

        # Generate performance report
        with performance_report(filename=f"{self.dask_log_directory}/seperate_tsvs_perf.html"):
            results = client.gather(futures)

        print("\nClosing cluster...")
        client.close()
        cluster.close()

        for motifs in results:
            all_motifs.update(motifs)

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

        print("Waiting for workers...")
        cluster = self.start_slurm_cluster()
        cluster.scale(self.nodes)
        client = Client(cluster)
        #client.wait_for_workers(self.nodes)

        tasks = []

        for motif in os.listdir(seperated_dir):
            motif_dir = os.path.join(seperated_dir, motif)
            motif_out_dir = os.path.join(out_dir, motif)

            if not os.path.exists(motif_out_dir):
                os.makedirs(motif_out_dir, exist_ok=True)
            
            for file in os.listdir(motif_dir):
                if file.endswith(self.locus_file_extension):
                    tasks.append(delayed(self.collapse_tsv_file)(
                        os.path.join(motif_dir, file),
                        motif_out_dir
                    ))

        futures = client.compute(tasks)
        progress(futures, notebook=False) # print progress
        
        with performance_report(filename=f"{self.dask_log_directory}/collapse_perf.html"):
            client.gather(futures)

        print("\nClosing cluster...")
        client.close()
        cluster.close()


def efficient_summarize(self, in_dir):
    """
    instead of seperating by motif and then summarizing with
    an online algorithm, merge tsv files with merge_asof for tolerance
    to generate the summary file
    """
    print("Creating Summary file...")

    tsvs = [os.path.join(in_dir, file) for file in os.listdir(in_dir) if file.endswith(self.locus_file_extension)]

    # merge all tsvs into one dataframe based on chromosome, repeat unit,
    # and left coordinate with a tolerance of slop
    df = pd.read_csv(tsvs[0], sep='\t')
    # keys must be sorted for merge_asof
    df.sort_values(['#chrom', 'repeatunit', 'left'], inplace=True)
    for tsv in tsvs[1:]:
        next_df = pd.read_csv(tsv, sep='\t')
        next_df.sort_values(['#chrom', 'repeatunit', 'left'], inplace=True)
        # need two operations to get full join (asof only works in one direction)
        df1 = pd.merge_asof(df, next_df, on='left', by=['#chrom', 'repeatunit'], tolerance=self.slop, direction='nearest')
        df2 = pd.merge_asof(next_df, df, on='left', by=['#chrom', 'repeatunit'], tolerance=self.slop, direction='nearest')
        df = pd.concat([df1, df2])
        df1.columns, df2.columns = df2.columns, df1.columns

        df.head()
        stop

        # in cases where there is a left coordinate match, but the repeat units are different,
        # take the mean of the left
        df['left'] = df[['left_x', 'left_y']].mean(axis=1)
        df.drop(columns=['left_x', 'left_y'], inplace=True)




        # in cases where the left coordinates are different, take the mean
        # of the two left coordina
    






