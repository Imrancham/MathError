import pandas as pd

# Read the JSON file into a Pandas DataFrame
json_file_path = 'error-log-2-999.json'
df = pd.read_json(json_file_path)

# Save the DataFrame as a CSV file
csv_file_path = 'output.csv'
df.to_csv(csv_file_path, index=False)

print("JSON file has been converted to CSV.")
