import os
import pandas as pd
from dask_jobqueue import SLURMCluster

class NIAGADSExperiment:
    def __init__(self, 
                 csv_metadata, 
                 sex=None, 
                 tissue=None,
                 dataset=None,
                 cohort="all_cohorts", 
                 race=None, 
                 ethnicity=None, 
                 apoe="all_apoe",
                 corrupt_file_log = "corrupt_files.log",
                 cores_per_node=36, 
                 mem_per_node="128GB",
                 partition="", 
                 account="", 
                 nodes=1, 
                 walltime="12:00:00",
                 dask_log_directory="/home/bbessell/str-analysis/dask_logs"):
        
        self.csv_metadata = csv_metadata

        self.metadict = {
            "Sex": sex,
            "Tissue": tissue,
            "Dataset": dataset,
            "Cohort": cohort,
            "Race": race,
            "Ethnicity": ethnicity,
            "APOE": apoe,
        }
        
        self.corrupt_file_log = corrupt_file_log
        open(self.corrupt_file_log, 'w').close()
        self._init_manifest()

        # SLURM cluster parameters
        self.cores_per_node = cores_per_node
        self.mem_per_node = mem_per_node
        self.partition = partition
        self.account = account
        self.nodes = nodes
        self.walltime = walltime
        self.dask_log_directory = dask_log_directory
        self.name_tag = f'{self.metadict["Cohort"]}_{self.metadict["APOE"]}'
        # add race, ethnicity, sex to name tag if not None
        if self.metadict["Sex"] is not None:
            self.name_tag += f"_{self.metadict['Sex']}"
        if self.metadict["Ethnicity"] is not None:
            self.name_tag += f"_{self.metadict['Ethnicity']}"
        if self.metadict["Race"] is not None:
            self.name_tag += f"_{self.metadict['Race']}"
        


    def _init_manifest(self):
            """
            Create a manifest file for the experiment from the NIAGADS metadata file
            """
            # read all columns as strings
            manifest_df = pd.read_csv(self.csv_metadata, dtype=str)

            if self.metadict["Cohort"] != "all_cohorts":
                manifest_df = manifest_df[manifest_df['Cohort'] == self.metadict["Cohort"]]
            if self.metadict["APOE"] != "all_apoe":
                manifest_df = manifest_df[manifest_df['APOE'] == self.metadict["APOE"]]
            if self.metadict["Sex"] is not None:
                manifest_df = manifest_df[manifest_df['Sex'] == self.metadict['Sex']]
            if self.metadict["Race"] is not None:
                manifest_df = manifest_df[manifest_df['Race'] == self.metadict['Race']]
            if self.metadict["Ethnicity"] is not None:
                manifest_df = manifest_df[manifest_df['Ethnicity'] == self.metadict['Ethnicity']]
            if self.metadict["Tissue"] is not None:
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
    
    def filter_raw_manifest(self, motif, file_dir):
        """
        DEPRACTED, FROM WHEN WE WERE USING THE ORIGINAL MANIFEST FILE
        Get case/control TSV file lists from the directory and filter files that adhere to desired covariates described by metadict
        """
        tsvs = []
        case_tsvs = []
        cont_tsvs = []
 
        cohort_subjects = pd.read_csv(self.csv_metadata).Subject.tolist()

        if motif != "no_seperation":

            if type(motif) is not str:  # handle nans (may want to write csv to log file so we can rerun these samples)
                return tsvs, case_tsvs, cont_tsvs
            
            dir_list = os.listdir(os.path.join(file_dir, motif))
        else:
            dir_list = os.listdir(file_dir)

        for file in dir_list:
            if file.endswith(self.locus_file_extension):
                subject, tissue = self.get_metadata_from_filename(file)
                print("SUBJECT: ", subject, tissue, self.metadict['Tissue'])
                if (subject in cohort_subjects or cohort_subjects == "all_cohorts") and (self.metadict['Tissue'] == None or tissue == self.metadict['Tissue']):
                    subject_metadata = self.get_metadata_from_subject(subject)
                    add_flag = True
                    for key, val in self.metadict.items():
                        print("KV: ", key, val, subject_metadata[key])
                        if val == "all_apoe": val = None
                        if val == "all_cohorts": val = None
                        if val and val != subject_metadata[key]:
                            print(val, subject_metadata[key])
                            add_flag = False
                            continue
                    if add_flag:
                        if motif == "no_seperation":
                            file_path = os.path.join(file_dir, file)
                        else:
                            file_path = os.path.join(file_dir, motif, file)

                        if subject_metadata["Diagnosis"] != "Unknown":
                            tsvs.append(file_path)
                            if subject_metadata["Diagnosis"] == "Case":
                                case_tsvs.append(file_path)
                            elif subject_metadata["Diagnosis"] == "Control":
                                cont_tsvs.append(file_path)

        print("FILTER TSVS: ", len(tsvs), len(case_tsvs), len(cont_tsvs))
            
        return motif, tsvs, case_tsvs, cont_tsvs
    

    def filter(self, motif, file_dir):
        """
        Get case/control TSV file lists from the directory, for a pre-filtered manifest file
        """
        tsvs = []
        case_tsvs = []
        cont_tsvs = []
        if motif != "no_seperation":
            if type(motif) is not str:
                return tsvs, case_tsvs, cont_tsvs
            dir_list = os.listdir(os.path.join(file_dir, motif))
        else:
            dir_list = os.listdir(file_dir)

        for file in dir_list:
            if file.endswith(self.locus_file_extension):
                subject, tissue = self.get_metadata_from_filename(file)
                if motif == "no_seperation":
                    file_path = os.path.join(file_dir, file)
                else:
                    file_path = os.path.join(file_dir, motif, file)
                subject_metadata = self.get_metadata_from_subject(subject)
                if subject_metadata["Diagnosis"] != "Unknown":
                    tsvs.append(file_path)
                    if subject_metadata["Diagnosis"] == "Case":
                        case_tsvs.append(file_path)
                    elif subject_metadata["Diagnosis"] == "Control":
                        cont_tsvs.append(file_path)

        return motif, tsvs, case_tsvs, cont_tsvs


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
