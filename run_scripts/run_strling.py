#!/usr/bin/env python

import os
import sys
import argparse

if __name__ == "__main__":

    import os
    from modules.STRlingExperiment import STRlingExperiment
    from modules.EHDNExperiment import EHDNExperiment
    from modules.plot_3_distributions import create_3_plots
    from modules.plot_distribution import create_1_plots
    import matplotlib
    import glob
    import shutil
    import pandas as pd
    import dask
    from s3_interface import download_tar_files_from_directory, unpack_tar_files

    cohort = "all_cohorts"
    test = "WT"
    apoe = "all_apoe"
    root_dir = '/scratch/remills_root/remills99/bbessell/str-data/'
    metadata_file = "manifests/ADNI_MIA_PR1066_VAN_StEP_manifest.csv"
    bucket_data_dir = '/mnt/data1/out'

    chromosomes = [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16',
        'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY'
    ]

    data_dir = os.path.join(root_dir, 'strling-calls')

    directory = os.path.join(data_dir, 'output')
    collapsed_directory = os.path.join(root_dir, 'collapsed', 'strling')
    seperated_dir = os.path.join(root_dir, 'seperated', 'strling')
    genotype_files = glob.glob(os.path.join(collapsed_directory, "*genotype*.txt"))

    # get strling summary_files
    # use chrom = "all" to run summaries in parallel by chrom

    strling_exp = STRlingExperiment(
        # data
        collapsed_directory,
        metadata_file,
        
        # cohort selection
        apoe=apoe,
        chroms=chromosomes,
        cohort=cohort,
        
        # analysis parameters
        slop=1000,
        slop_modifier=5,
        test=test,
        
        # cluster parameters
        cores_per_node=36,
        mem_per_node="180GB",
        partition="standard",
        account="remills99",
        nodes=2,
        walltime="8:00:00",
        dask_log_directory="/home/bbessell/str-analysis/dask_logs"
        
    )

    strling_exp.schedule_seperate_tsvs_by_motif(directory, seperated_dir)
    strling_exp.schedule_collapse_tsvs(seperated_dir, collapsed_directory)
    strling_exp.schedule_filter_tsv_files()
    strling_exp.schedule_summarizing()  # give csv to load instead of calculate
    strling_exp.schedule_perform_stat_test()
