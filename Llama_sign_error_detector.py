import re
from transformers import pipeline
import torch
import pandas as pd
from datasets import Dataset

# Initialize the pipeline
model_id = "meta-llama/Llama-3.3-70B-Instruct"  # Updated model ID

llm = pipeline(
    "text-generation",
    model=model_id,
    model_kwargs={"torch_dtype": torch.bfloat16},
    device_map="auto",
    batch_size=4,  # or higher, depending on VRAM
)

def extract_sign_error(response_text: str) -> bool | None:
    """
    Extracts the final 'Sign error: True/False' from the model response.
    Returns True, False, or None if unclear.
    """
    match = re.search(r"sign\s*error\s*[:\-]?\s*(true|false)", response_text.strip(), re.IGNORECASE)
    if match:
        return match.group(1).lower() == "true"
    return None  # Couldn't find a valid answer

def detect_sign_error(original_expr, student_expr):
    """Use LLM to detect if there's a sign error between two expressions"""


    prompt = f"""
You are an expert in reviewing math problem-solving. Your specialty is identifying sign errors in provided solution.
         You will receive the original math expression, and the student transformation.
1. Analyze the student's transformation step by step, and check if there are sig error.
2. be carfull about implicit sign error, like in the example: original expression: 2-4*s = -r+4*s Student's answer:  2 = -3. it contains a sign error, as the student moved the -4*s to the other side, but forgot to change its sign.
         
Original expression: {original_expr}
Student's answer: {student_expr}
        Only respond with 'True' if there is definitely a sign error (changed + to - or - to +).
        Respond with 'False' if there is no sign error or if the changes are not sign-related.
        Response:
"""

    
    messages = [
        {"role": "system", "content": "You are a math expert who detects sign errors in equations."},
        {"role": "user", "content": prompt},
    ]

    response = llm(
        messages,
        max_new_tokens=400
    )
    
    # Extract and clean the response
    answer = response[0]["generated_text"][-1]["content"].strip().lower()
    bool_answre = extract_sign_error(answer)
    if bool_answre is not None:
        return bool_answre
    return False

# Read CSV file
df = pd.read_csv('/data/ichamieh/mistral/similarity_results_full_Set_2.csv')  # Replace with your file path
# Process each row
results = []
for index, row in df.iterrows():
    original = row['expression 1']
    student = row['expression 2']
    has_error = detect_sign_error(original, student)
    results.append(has_error)
    print(f"Row {index+1}: Original: {original} | Student: {student} | Sign error: {has_error}")

# Add results to DataFrame and save
df['sign_error'] = results
df.to_csv('analyzed_results3.csv', index=False)
print("Analysis complete. Results saved to 'analyzed_results.csv'")