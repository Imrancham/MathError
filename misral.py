import json
import os
import re
import time
from mistralai import Mistral
from mistralai.models.sdkerror import SDKError

import pandas as pd

api_key = "dXMrrJyyV2G2RjabskIGmylVcCuap55E"
model = "open-mistral-nemo"  # "pixtral-12b-2409"
matched_ids = {'id': [1,2]}
df = pd.read_json("/data/ichamieh/mistral/1_200_annotated.csv")

client = Mistral(api_key=api_key)

for index, item in df[1:2].iterrows():
    t0 = item["t0"]
    t1 = item.t1

    examale = f"""
Given an expression and its corresponding transformation in latex, check wheather the transformation contain a typo errorlike using similar letters like q and g, p and P, or o instead of 0:


Example:
- Input:  
  Original: \( x + 4 = -x - 50 \)  
  Transformation: \( 2X = -54 \)  
- output:  
  Error:True. using X instead og x.

  - Input: 
    Original:   4=3d+4v
    Transformation:  o=3d+4v-4
  -  Output:
    Error: True. Use 0 instead of o 

   - Input: 
    Original:   \(2x+4\)\(2x-4)
    Transformation:  4x-16
  -  Output:
    Error: False. Other error


    Check the following experssoin, return only true or false:
    - input: 
      Original:  {t0}
      Transformation:  {t1}
    - output:
    
    """

    while True:
        try:
            chat_response = client.chat.complete(
                model=model,
                messages=[
                    {
                        "role": "user",
                        "content": examale,
                    },
                ]
            )
            break
        except SDKError as e:
            if e.status_code == 429:
                print("Rate limit exceeded. Waiting for 1 minute...")
                wait_time = 60
                time.sleep(wait_time)  # Wait for 5 minutes
            else:
                raise e

    print(index, chat_response.choices[0].message.content)
    print("--------------------")

    # Regular expression to find "True" or "False" and capture the ID
    pattern = r'\s*Error:\s*(True)'
    # if the above pattern is found return the ID
    match = re.search(pattern, chat_response.choices[0].message.content)
    print(match)
    if match:
        matched_ids.append(item["id"])
    else:
        print("No match found")

# Save the matched IDs to a JSON file
output_file_path = "matched_ids.json"

with open(output_file_path, "w") as file:
    json.dump(matched_ids, file)

print(f"Matched IDs saved to {output_file_path}")

"""model = "mistral-embed"

client = Mistral(api_key=api_key)

embeddings_response = client.embeddings.create(
    model=model,
    inputs=["Embed this sentence.", "As well as this one."]
)

print(embeddings_response)"""