import pandas as pd
import matplotlib.pyplot as plt

# Manually add the file paths to a list
file_paths = ['rage.fps.0_42.0K.csv', 'rage.fps.1_42.0K.csv', 'rage.fps.2_42.0K.csv', 'rage.fps.3_42.0K.csv']  # Replace with actual file paths

# Set your desired `world_ix` value here
world_ix_value = 3

# Initialize empty lists to store combined values
fix_o_combined = []
fix_n_combined = []
fix_set_combined = []

# Loop through each file and process them
for file in file_paths:
    # Read the CSV file, ignoring the first line (header)
    df = pd.read_csv(file)
    
    # Filter rows where 'world_ix' matches the set value
    filtered_df = df[df['world_ix'] == world_ix_value]

    # Remove extremely large values
    filtered_df = filtered_df[(filtered_df['fix_o'] <= 1e6) & (filtered_df['fix_n'] <= 1e6) & (filtered_df['fix_set'] <= 1e6)]
    
    # Append 'fix_o', 'fix_n', and 'fix_set' values to the combined lists
    fix_o_combined.extend(filtered_df['fix_o'].values)
    fix_n_combined.extend(filtered_df['fix_n'].values)
    fix_set_combined.extend(filtered_df['fix_set'].values)

# Plot the combined values
plt.plot(fix_o_combined, label='fix_o')
plt.plot(fix_n_combined, label='fix_n')
plt.plot(fix_set_combined, label='fix_set')
plt.legend()
plt.xlabel('Index')
plt.ylabel('fix Values')
plt.title(f'Combined fix Values for world_ix = {world_ix_value}')
plt.show()
