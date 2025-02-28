import re
import time
from mistralai import Mistral
from mistralai.models.sdkerror import SDKError
from sympy.parsing.latex import parse_latex

import os
from openai import OpenAI

import pandas as pd

api_key = os.getenv("api_key_deepSeek")
client = OpenAI(api_key, base_url="https://api.deepseek.com")

model = "mistral-large-latest"
matched_ids = []
df = pd.read_json("0_9999_v7.json")
output_file_path = "matched_ids2.txt"


#client = Mistral(api_key=api_key)

main_prompt = "You are an expert in analyzing mathematical errors. Below is an incorrect transformation in student solutions. answer by yes no. "

prompts = {'Arithmetic' : "Does the transformation contain a basic calculation mistake, such as incorrect addition, subtraction, multiplication, or division?"

,'Sign_Error' : "Does the transformation involve an incorrect handling of signs, such as missing or incorrect negative signs, or errors in inequality direction?"

,'Structure_Error' : "Does the transformation incorrectly rearrange terms or change the structure of the equation in a way that alters its meaning?"

,'Parenthesis_Expansion_Error' : "Does the transformation contain a mistake in factorization or distribution, such as incorrect expansion of parentheses or incorrect factoring?"

,'Exponent_Error' :  "Does the transformation involve an incorrect manipulation of exponents, such as misapplying exponent rules or incorrectly simplifying powers?"

,'Square_Root_Error' : "Does the transformation contain an incorrect handling of square roots, such as misapplying the square root property or incorrectly simplifying a radical expression?"

,'Careless_Error' : "Does the transformation contain a transcription error, such as missing symbols, miswritten numbers, or incorrect notation?"

,'Absolute_Value_Error' : "Does the transformation involve an incorrect handling of absolute values, such as missing absolute value symbols or misapplying absolute value properties?"}
for index, item in df.head(10).iterrows():
    t0 = parse_latex(item.t0)
    t1 = parse_latex(item.t1)

    for key in prompts.keys():
        prompt = f"""{main_prompt} \n {prompts[key]} \n - input: \n Original:  {t0} \n Transformation:  {t1}"""
       #prompt = f"""
       #- input: 
       #  Original:  {t0}
       #  Transformation:  {t1}
       #{prompts[key]}
       #"""
        retry_count = 0
        max_retries = 5
        print("--------------------------")
        print(f"t0: {t0}")
        print(f"t1: {t1}")
        while True:
            try:
                chat_response = client.chat.completions.create(
                    model="deepseek-reasoner",
                    messages=[
                        {
                            "role": "user",
                            "content": prompt,
                        },
                    ]
                )
                break
            except SDKError as e:
                if e.status_code == 429:
                    print("Rate limit exceeded. Waiting for 1 minute...")
                    wait_time = 5
                    time.sleep(wait_time)  # Wait for 5 minutes
                else:
                    raise e

        print(index,key,  chat_response.choices[0].message.content)
    

    # Regular expression to find "True" or "False" and capture the ID
    # pattern that find True or ture or TRUE in any where
    pattern = r"(?i)(True|true|TRUE)"
    

    # if the above pattern is found return the ID
    match = re.search(pattern, chat_response.choices[0].message.content)
    print(match)
    if match:
        matched_ids.append(item["id"])
        print("Match found")
        print(item["id"])
        # Append the matched ID to the JSON file
        with open(output_file_path, "a") as file:
            file.write(f"{item['id']}\n")
    else:
        print("No match found")
    print("--------------------")
# Save the matched IDs to a JSON file
print(matched_ids)



print(f"Matched IDs saved to {output_file_path}")

