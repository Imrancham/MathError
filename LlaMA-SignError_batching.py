import transformers
import torch
import pandas as pd
from datasets import Dataset
from tqdm import tqdm
import os
import re
import logging

# Set up logging
logging.basicConfig(filename='sign_error_analysis.log', level=logging.INFO)

# Initialize tokenizer with proper settings for LLaMA
tokenizer = transformers.AutoTokenizer.from_pretrained(
    "meta-llama/Llama-3.3-70B-Instruct",
    padding_side="left"
)
tokenizer.pad_token = tokenizer.eos_token

# Initialize pipeline with proper configuration
pipeline = transformers.pipeline(
    "text-generation",
    model="meta-llama/Llama-3.3-70B-Instruct",
    tokenizer=tokenizer,
    torch_dtype=torch.bfloat16,
    device_map="auto",
    batch_size=2,
    return_full_text=False
)

def extract_sign_error(response_text: str) -> bool:
    """Extracts the sign error determination from the model response."""
    match = re.search(r"(sign\s*error|error)[:\-]?\s*(true|false)", response_text.strip(), re.IGNORECASE)
    if match:
        return match.group(2).lower() == "true"
    return False  # Default to False if unclear

def detect_sign_error(batch):
    """Process a batch of examples to detect sign errors"""
    prompts = [ f"""
    Analyze these two mathematical expressions and determine if the student made a sign error:
    Original: {original}
    Student's version: {student}

    Only respond with 'True' if there is definitely a sign error (changed + to - or - to +).
    Respond with 'False' if there is no sign error or if the changes are not sign-related.
    Response: """
                       for original, student in zip(batch['expression 1'], batch['expression 2'])
    ]
    #f"Original expression: {original}\nStudent's expression: {student}\nIs there a sign error? Answer only True or False."
#        f"""Analyze these mathematical expressions for sign errors:
#Original: {original}
#Student: {student}
#
#After careful analysis, respond ONLY with:
#"Sign error: True" if there's a sign inversion (+/-)
#"Sign error: False" otherwise"""


    
    messages = [[
        {"role": "system", "content": "You are a math tutor. You identify sign errors in student answers."},
        {"role": "user", "content": prompt}
    ] for prompt in prompts]

    try:
        responses = pipeline(
            messages,
            max_new_tokens=50,
            pad_token_id=tokenizer.eos_token_id,
            eos_token_id=tokenizer.eos_token_id
        )
        
        results = []
        for i, r in enumerate(responses):
            # Handle the response structure correctly
            if isinstance(r, list) and len(r) > 0:
                response_text = r[0].get("generated_text", "").strip().lower()
            else:
                response_text = str(r).strip().lower()
                
            logging.info(f"Original: {batch['expression 1'][i]}")
            logging.info(f"Student: {batch['expression 2'][i]}")
            logging.info(f"Response: {response_text}")
            
            error_status = extract_sign_error(response_text)
            results.append(error_status)
        
        return {'sign_error': results}
    
    except Exception as e:
        logging.error(f"Error processing batch: {str(e)}")
        return {'sign_error': [False] * len(batch)}  # Safe default

# Main processing function
def process_data(input_path, output_path, temp_path):
    # Load or initialize results
    if os.path.exists(temp_path):
        results_df = pd.read_pickle(temp_path)
        processed_count = results_df['sign_error'].notna().sum()
        print(f"Resuming from row {processed_count}")
    else:
        df = pd.read_csv(input_path)
        results_df = df.copy()
        results_df['sign_error'] = None
        processed_count = 0

    dataset = Dataset.from_pandas(results_df)
    batch_size = 2
    
    with tqdm(total=len(dataset), desc="Processing") as pbar:
        while processed_count < len(dataset):
            batch = dataset.select(range(processed_count, min(processed_count + batch_size, len(dataset))))
            results = detect_sign_error(batch)
            
            for i, error in enumerate(results['sign_error']):
                results_df.at[processed_count + i, 'sign_error'] = error
            
            processed_count += batch_size
            pbar.update(batch_size)
            
            # Save checkpoint every 10 batches
            if processed_count % (10 * batch_size) == 0:
                results_df.to_pickle(temp_path)
                pbar.set_postfix_str("Checkpoint saved")

    results_df.to_csv(output_path, index=False)
    return results_df

# Run the processing
if __name__ == "__main__":
    results = process_data(
        input_path='/data/ichamieh/mistral/similarity_results_full_Set_2.csv',
        output_path='analyzed_results_final_short_prompt_2.csv',
        temp_path='processing_checkpoint.pkl'
    )