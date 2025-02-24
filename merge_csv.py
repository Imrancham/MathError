import os
import pandas as pd

def merge_csv_files(data_folder):
    # List all CSV files in the data folder
    csv_files = [f for f in os.listdir(data_folder) if f.endswith('.csv')]
    
    # Read and concatenate all CSV files
    df_list = [pd.read_csv(os.path.join(data_folder, file)) for file in csv_files]
    merged_df = pd.concat(df_list, ignore_index=True)
    
    # Remove duplicates
    no_duplicates_df = merged_df.drop_duplicates()
    
    # Save the merged dataframes
    merged_df.to_csv(os.path.join(data_folder, 'merged_with_duplicates.csv'), index=False)
    no_duplicates_df.to_csv(os.path.join(data_folder, 'merged_no_duplicates.csv'), index=False)
    
    return merged_df, no_duplicates_df

if __name__ == "__main__":
    data_folder = 'data'  # Specify the path to your data folder
    merged_df, no_duplicates_df = merge_csv_files(data_folder)
    print("Merged CSV files with duplicates saved as 'merged_with_duplicates.csv'")
    print("Merged CSV files without duplicates saved as 'merged_no_duplicates.csv'")


    """import os
import pandas as pd
import json

def merge_json_files(data_folder):
    # List all JSON files in the data folder
    json_files = [f for f in os.listdir(data_folder) if f.endswith('.json')]
    
    # Read and concatenate all JSON files
    df_list = []
    for file in json_files:
        with open(os.path.join(data_folder, file), 'r', encoding='utf-8') as f:
            data = json.load(f)
            df_list.append(pd.DataFrame(data))
    
    merged_df = pd.concat(df_list, ignore_index=True)
    
    # Remove duplicates
    no_duplicates_df = merged_df.drop_duplicates()
    
    # Save the merged dataframes
    merged_df.to_json(os.path.join(data_folder, 'merged_with_duplicates.json'), orient='records', lines=True)
    no_duplicates_df.to_json(os.path.join(data_folder, 'merged_no_duplicates.json'), orient='records', lines=True)
    
    return merged_df, no_duplicates_df

if __name__ == "__main__":
    data_folder = 'data'  # Specify the path to your data folder
    merged_df, no_duplicates_df = merge_json_files(data_folder)
    print("Merged JSON files with duplicates saved as 'merged_with_duplicates.json'")
    print("Merged JSON files without duplicates saved as 'merged_no_duplicates.json'")"""