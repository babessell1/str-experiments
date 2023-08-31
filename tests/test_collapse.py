import pandas as pd
import numpy as np

def collapse_variants(df, slop):
    df.sort_values('left', inplace=True)
    df.reset_index(drop=True, inplace=True)
    collapsed_variants = []

    prev_left = None
    prev_unit = None
    prev_counts = None
    prev_chrom = None
    prev_mean_left = None
    prev_m2 = None

    for index, row in df.iterrows():
        left = row['left']
        unit = row['repeatunit']
        count = row['counts']
        chrom = row['#chrom']
        mean_left = row['mean_left']

        if (
            prev_chrom is not None
            and prev_chrom == chrom
            and is_rotation(unit, prev_unit)
            and left - prev_left <= slop
        ):
            # Pairwise algorithm
            n = prev_counts + count
            delta = mean_left - prev_mean_left
            delta2 = delta * delta

            new_mean = prev_mean_left + delta * count / n
            new_m2 = prev_m2 + prev_counts * count * (delta2 / n)

            prev_counts = n
            prev_mean_left = new_mean
            prev_m2 = new_m2
        else:
            if prev_chrom is not None:
                variance_left = prev_m2 / (prev_counts - 1) if prev_counts > 1 else 1e-6
                collapsed_variants.append({
                    '#chrom': prev_chrom,
                    'left': prev_left,
                    'repeatunit': prev_unit,
                    'counts': prev_counts,
                    'variance_left': variance_left,
                    'mean_left': prev_mean_left
                })

            prev_left = left
            prev_unit = unit
            prev_counts = count
            prev_chrom = chrom
            prev_mean_left = mean_left
            prev_m2 = 0.0

    if prev_chrom is not None:
        variance_left = prev_m2 / (prev_counts - 1) if prev_counts > 1 else 1e-6
        collapsed_variants.append({
            '#chrom': prev_chrom,
            'left': prev_left,
            'repeatunit': prev_unit,
            'counts': prev_counts,
            'variance_left': variance_left,
            'mean_left': prev_mean_left
        })

    new_df = pd.DataFrame(collapsed_variants)

    return new_df


# Input DataFrame
data = {
    '#chrom': ['chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1', 'chr1'],
    'left': [10179, 10179, 10180, 10182, 10185, 20000, 20001, 20002],
    'repeatunit': ['AACCCT', 'ACCCTA', 'AACCCT', 'AACCCT', 'AACCCT', 'AG', 'AG', 'AG'],
    'sum_left': [10179, 10179, 10180, 10182, 10185, 20000, 20001, 20002],
    'counts': [1, 1, 1, 1, 1, 1, 1, 1],
    'variance_left': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'mean_left': [10179.0, 10179.0, 10180.0, 10182.0, 10185.0, 20000.0, 20001.0, 20002.0]
}

input_df = pd.DataFrame(data)

# Expected Output DataFrame
expected_output_data = {
    '#chrom': ['chr1', 'chr1'],
    'left': [10179, 20000],
    'repeatunit': ['AACCCT', 'AG'],
    'counts': [5, 3],
    'variance_left': [6.5, 1.0],
    'mean_left': [10181., 20001.]
    #'sum_left': [50905, 60003],
}

expected_output_df = pd.DataFrame(expected_output_data)

# Set the slop value
slop = 150

# Apply the function
output_df = collapse_variants(input_df, slop)

print("--------------INPUT----------------")
print(input_df)
print("--------------OUTPUT----------------")
print(output_df)
print("---------EXPECTED OUTPUT------------")
print(expected_output_df)


# Compare the output with the expected output
pd.testing.assert_frame_equal(output_df, expected_output_df)

print("Test Passed!")
