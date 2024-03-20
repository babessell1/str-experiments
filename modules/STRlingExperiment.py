
import os
from scipy.stats import ranksums, ks_2samp
from statsmodels.stats.multitest import multipletests
import dask.dataframe as dd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import multiprocessing
from itertools import product
import pickle as pkl
from modules.RepeatsExperiment_pandas import RepeatsExperiment

class STRlingExperiment(RepeatsExperiment):
    def __init__(self, tsv_dir, csv_metadata, **kwargs):
        super().__init__(tsv_dir, csv_metadata, **kwargs)
        
        self.cols_to_drop = [
            'allele1_est',
            'allele2_est',
            'right',
            'merged_expansions'
        ]

        self.test_variable = 'allele2_est'
        self.caller = 'STRling'
        self.locus_file_extension = '-genotype.txt'
        self.repeat_nomenclature = 'repeatunit'
   
    @staticmethod
    def filter_variants(tsv_df, chrom):
        """
        Remove single nucleotide expansions and select by chromosome
        """
        filtered_df = tsv_df[~tsv_df['repeatunit'].isin(['A', 'G', 'C', 'T'])]
        if chrom != 'All':
            filtered_df = filtered_df[filtered_df['#chrom'] == chrom]

        return filtered_df

    def collapse_tsv_file_slow(self, tsv_file, out_dir):
        """
        Process a single TSV file, collapsing variants.
        """

        if os.stat(tsv_file).st_size == 0:
            print(f"Skipping empty file: {tsv_file}")
            return

        subject, tissue = self.get_metadata_from_filename(tsv_file)
        
        # NOTE: pandas is probably better than dask for this because the data is small and 
        # we have to sort and reset indexes 
        df = pd.read_csv(tsv_file, sep='\t')
        df = self.filter_variants(df, "All")   
        df = df.sort_values(['#chrom', 'repeatunit', 'left'])
        df = df.reset_index(drop=True)
        df['counts'] = np.ones(df.shape[0])
        df['mean_left'] = df.left

        # remove tsv fie if it's empty
        if df.shape[0] == 0:
            #print(f'{tsv_file} is empty')
            os.remove(tsv_file)
            return

        collapsed_variants = []

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
                        #print(unit)
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
            new_df = new_df.drop(columns=['counts', 'mean_left'])

            # Write the result to an output file
            tsv_file_name = os.path.basename(tsv_file)
            # make sure the motif directory exists

            out_file_path = os.path.join(out_dir, tsv_file_name)
            new_df.to_csv(out_file_path, sep='\t', index=False)
        except:
            print("failed to process tsv file: ", tsv_file)

            return

    
    def collapse_tsv_file(self, tsv_file, out_dir):
        if os.stat(tsv_file).st_size == 0:
            print(f"Skipping empty file: {tsv_file}")
            return
        
        df = pd.read_csv(tsv_file, sep='\t')
        df = self.filter_variants(df, "All")
        df = df.sort_values(['#chrom', 'repeatunit', 'left'])
        df = df.reset_index(drop=True)

        if df.empty:
            os.remove(tsv_file)
            return

        # Create a custom group key based on the slop_modifier
        df['group_key'] = df.groupby(['#chrom', 'repeatunit']).apply(
            lambda x: (x['left'] // (self.slop * self.slop_modifier)).astype(int)
        ).array

        # Define a function for collapsing within groups
        def collapse(group):
            sum_allele1_est = group['allele1_est'].sum()
            sum_allele2_est = group['allele2_est'].sum()
            return pd.Series({
                'left': group['left'].mean().astype(int),
                'right': group['left'].mean().astype(int) + 1,
                'repeatunit': group['repeatunit'].iloc[0],
                'allele1_est': sum_allele1_est,
                'allele2_est': sum_allele2_est,
                'counts': group.shape[0],
                'mean_left': group['left'].mean(),
                'merged_expansions': group.shape[0] > 1
            })

        # Apply the collapse function to each group and reset index
        collapsed_df = df.groupby("group_key").apply(collapse).reset_index(drop=True)

        out_file_path = os.path.join(out_dir, os.path.basename(tsv_file))
        collapsed_df.to_csv(out_file_path, sep='\t', index=False)

        print(f"Successfully collapsed file: {out_file_path}")


    def get_kmer_counts(self, iter):
        """
        get mappable normalized kmer count by getting values of 
        sum_str_count from -genotype file, dividing by the depth
        couln from the -genotype file. Also get the unplaced counts
        from the -unplaced file and divide by the average depth
        """
        motif = iter[0]
        tsv = iter[1]

        print("motif: ", motif)
        print("tsv: ", tsv)

        # get the depth from the -genotype file (is tsv)
        motif_file = tsv.replace('-genotype.txt', f'-unplaced.txt')
        motif_df = pd.read_csv(motif_file, sep='\t', header=None)
        genotype_df = pd.read_csv(tsv, sep='\t')
        avg_depth = genotype_df['depth'].mean()
        # get value of second column where first column is motif
        unplaced_counts = motif_df[motif_df[0] == motif][1].values[0] / avg_depth
        
        # get the placed counts from the -genotype file
        motif_match = genotype_df[genotype_df['repeatunit'] == motif]
        kmers = genotype_df[motif_match].sum_str_count.values[0] / genotype_df[motif_match].depth.values[0]

        return sum(kmers), unplaced_counts, motif, tsv
    

    def schedule_get_kmer_counts(self):
        """
        Schedule a get_kmer_counts job
        """
        if not hasattr(self, 'tsvs'): # if self.tsvs does not exist, get it from the pickled files
            with open(f"{self.caller}_tsvs.txt", 'rb') as handle:
                self.tsvs = pkl.load(handle)
            with open(f"{self.caller}_case_tsvs.txt", 'rb') as handle:
                self.case_tsvs = pkl.load(handle)
            with open(f"{self.caller}_cont_tsvs.txt", 'rb') as handle:
                self.cont_tsvs = pkl.load(handle)

        print(self.tsvs)
        #print(product(self.motifs, self.tsvs))


        self.motifs = pd.read_csv(os.path.join(f'{self.caller}_motifs.txt'), sep='\t').iloc[:,0].dropna().tolist()
        pool = multiprocessing.Pool()
        temp = product(self.motifs, self.tsvs)
        print([iter for iter in temp])
        results = pool.imap_unordered(self.get_kmer_counts, product(self.motifs, self.tsvs))

        case_kmer_counts = {}
        control_kmer_counts = {}

        unplaced_case_counts = {}
        unplaced_control_counts = {}

        sum_case_counts = {}
        sum_control_counts = {}

        for result in results:
            motif = result[0]
            tsv = result[1]

            if tsv in self.case_tsvs[motif]:
                # if motif not in case_kmer_counts, add it
                if motif not in case_kmer_counts:
                    case_kmer_counts[motif] = [result[2]]
                else:
                    case_kmer_counts[motif].append(result[2])
                if motif not in unplaced_case_counts:
                    unplaced_case_counts[motif] = [result[3]]
                else:
                    unplaced_case_counts[motif].append(result[3])
                if motif not in sum_case_counts:
                    sum_case_counts[motif] = [result[2]+result[3]]
                else:
                    sum_case_counts[motif].append(result[2]+result[3])

            elif tsv in self.control_tsvs[motif]:
                if motif not in control_kmer_counts:
                    control_kmer_counts[motif] = [result[2]]
                else:
                    control_kmer_counts[motif].append(result[2])
                if motif not in unplaced_control_counts:
                    unplaced_control_counts[motif] = [result[3]]
                else:
                    unplaced_control_counts[motif].append(result[3])
                if motif not in sum_control_counts:
                    sum_control_counts[motif] = [result[2]+result[3]]
                else:
                    sum_control_counts[motif].append(result[2]+result[3])
            
        # perform Wilcoxon rank sum test on case and control kmer counts
        case_control_pvals = {}
        for motif in self.motifs:
            case_kmers = case_kmer_counts[motif]
            control_kmers = control_kmer_counts[motif]
            case_control_pvals[motif] = ranksums(case_kmers, control_kmers)[1]
        
        # perform Wilcoxon rank sum test on unplacd case and control kmer counts
        unplaced_case_control_pvals = {}
        for motif in self.motifs:
            case_kmers = unplaced_case_counts[motif]
            control_kmers = unplaced_control_counts[motif]
            unplaced_case_control_pvals[motif] = ranksums(case_kmers, control_kmers)[1]
        
        # perform Wilcoxon rank sum test on sum case and control kmer counts
        sum_case_control_pvals = {}
        for motif in self.motifs:
            case_kmers = sum_case_counts[motif]
            control_kmers = sum_control_counts[motif]
            sum_case_control_pvals[motif] = ranksums(case_kmers, control_kmers)[1]
        
        # correct p-values
        case_control_pvals = multipletests(list(case_control_pvals.values()), method='fdr_bh')[1]
        unplaced_case_control_pvals = multipletests(list(unplaced_case_control_pvals.values()), method='fdr_bh')[1]
        sum_case_control_pvals = multipletests(list(sum_case_control_pvals.values()), method='fdr_bh')[1]

        # create dataframe of p-values
        case_control_df = pd.DataFrame(case_control_pvals, index=list(case_control_pvals.keys()), columns=['case_control_pval'])
        # add the case and control kmer counts to the dataframe
        case_control_df['case_kmer_counts'] = case_control_df.index.map(case_kmer_counts)
        case_control_df['control_kmer_counts'] = case_control_df.index.map(control_kmer_counts)

        # write dataframe to file
        case_control_df.to_csv(f'{self.caller}_case_control_pvals.txt', sep='\t')
        del case_control_df

        unplaced_case_control_df = pd.DataFrame(unplaced_case_control_pvals, index=list(unplaced_case_control_pvals.keys()), columns=['unplaced_case_control_pval'])
        # add the unplaced case and control kmer counts to the dataframe
        unplaced_case_control_df['unplaced_case_kmer_counts'] = unplaced_case_control_df.index.map(unplaced_case_counts)
        unplaced_case_control_df['unplaced_control_kmer_counts'] = unplaced_case_control_df.index.map(unplaced_control_counts)
        
        # write dataframe to file
        unplaced_case_control_df.to_csv(f'{self.caller}_unplaced_case_control_pvals.txt', sep='\t')
        del unplaced_case_control_pvals

        sum_case_control_df = pd.DataFrame(sum_case_control_pvals, index=list(sum_case_control_pvals.keys()), columns=['sum_case_control_pval'])
        # add the sum case and control kmer counts to the dataframe
        sum_case_control_df['sum_case_kmer_counts'] = sum_case_control_df.index.map(sum_case_counts)
        sum_case_control_df['sum_control_kmer_counts'] = sum_case_control_df.index.map(sum_control_counts)

        # write dataframe to file
        sum_case_control_df.to_csv(f'{self.caller}_sum_case_control_pvals.txt', sep='\t')
        del sum_case_control_df

    
    def plot_kmer_distributions(self):
        """
        Plot the kmer distributions for each significant motif
        """
        titles = [
            'Case vs Control Kmer Counts (placed)',
            'Case vs Control Kmer Counts (unplaced)',
            'Case vs Control Kmer Counts (sum)'
        ]

        def make_plot(title):
            for i, row in case_control_df.iterrows():
                motif = row.name
                case_values = row['case_kmer_counts']
                control_values = row['control_kmer_counts']
                # plot the distributions
                max_value = max(max(case_values), max(control_values))
                min_value = min(min(case_values), min(control_values))
                bins = np.linspace(min_value, max_value, num=20)
                
                # Calculate normalized heights by dividing counts by the total number of samples
                case_counts, case_bins = np.histogram(case_values, bins=bins)
                control_counts, control_bins = np.histogram(control_values, bins=bins)
                total_case_samples = len(case_values)
                total_control_samples = len(control_values)
                case_heights = case_counts / total_case_samples
                control_heights = control_counts / total_control_samples
                ##########################

                # plot the distributions
                fig, ax = plt.subplots()

                ax.bar(bins[:-1], case_heights, width=np.diff(bins), alpha=0.6, label='Case Kmers')
                ax.bar(bins[:-1], control_heights, width=np.diff(bins), alpha=0.6, label='Control Kmers')

                # title
                ax.set_title(title)

                # show
                plt.show()
            

        # read in the p-values
        case_control_df = pd.read_csv(f'{self.caller}_case_control_pvals.txt', sep='\t', index_col=0)

        # plot the distributions for each significant motif
        case_control_df = case_control_df[case_control_df['case_control_pval'] < 0.05]
        # sort by p-value
        case_control_df = case_control_df.sort_values(by=['case_control_pval'], ascending=True)
        
        # plot the distributions
        make_plot(titles[0])

        # next plot
        unplaced_case_control_df = pd.read_csv(f'{self.caller}_unplaced_case_control_pvals.txt', sep='\t', index_col=0)
        unplaced_case_control_df = unplaced_case_control_df[unplaced_case_control_df['unplaced_case_control_pval'] < 0.05]
        unplaced_case_control_df = unplaced_case_control_df.sort_values(by=['unplaced_case_control_pval'], ascending=True)
        make_plot(titles[1])

        # next plot
        sum_case_control_df = pd.read_csv(f'{self.caller}_sum_case_control_pvals.txt', sep='\t', index_col=0)
        sum_case_control_df = sum_case_control_df[sum_case_control_df['sum_case_control_pval'] < 0.05]
        sum_case_control_df = sum_case_control_df.sort_values(by=['sum_case_control_pval'], ascending=True)
        make_plot(titles[2])

            


