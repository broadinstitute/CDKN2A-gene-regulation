import pandas as pd
import glob

# Load all the BEDGraph files
file_paths = glob.glob("sceptre_guides_nG50_new_ref_anril_guide*.bedGraph")

# Initialize an empty DataFrame for merging
merged_df = pd.DataFrame()

for file in file_paths:
    # Read each BEDGraph file into a DataFrame, skipping the header
    df = pd.read_csv(file, sep='\t', comment='t', header=None, names=['chr', 'start', 'end', 'value'])
    
    if merged_df.empty:
        merged_df = df
    else:
        # Merge the new DataFrame with the existing one, keeping the max value per bin
        merged_df = pd.merge(merged_df, df, on=['chr', 'start', 'end'], how='outer', suffixes=('', '_new'))
        merged_df['value'] = merged_df[['value', 'value_new']].max(axis=1)
        merged_df = merged_df.drop(columns=['value_new'])

# Write the final merged DataFrame back to a BEDGraph file
with open("combined_max_values_guide.bedGraph", "w") as f_out:
    f_out.write("track type=bedGraph\n")  # Writing the track line
    merged_df.to_csv(f_out, sep='\t', header=False, index=False)

