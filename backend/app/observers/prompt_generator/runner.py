# setup to initialize app using path
import os
import sys

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../../..')))

from app.observers.reaction_iter import provide_reaction_details
from app.observers.prompt_generator.base_prompt import BASE_PROMPT
import base64
from openai import OpenAI
import pandas
import tqdm

def encode_image(image_path):
    """
    Reads an image file and encodes it in Base64 format.
    """
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode('utf-8')

def run_prompt_generator(source_df_path: str, column_names=["Serial No.", "Reaction Type", "Source"]):
    print("Running prompt generator...")
    # create an empty df assign the serial no. and the generated prompt that will be again merged with the source df later
    dummy_df = pandas.DataFrame(columns=["Serial No.", "Generated Prompt"])
    for value_dict in provide_reaction_details(source_df_path, column_to_access=column_names):
        reaction_type = value_dict["Reaction Type"]
        source = value_dict["Source"]
        sl_no = value_dict["Serial No."]
        image_path = f"/home/kalki/src/valency/reactor/analysis/data/rxn_data/{sl_no}/sources/image.png"
        prompt = BASE_PROMPT.format(reaction_type=reaction_type)

        image_str = encode_image(image_path) if os.path.exists(image_path) else None
        # print("len of image_str >> ", len(image_str)) if image_str is not None else print("image_str is None")

        client = OpenAI()

        response = client.chat.completions.create(
            model="gpt-4o", #  Or another multimodal model like gpt-4o-mini
            messages=[
                {
                    "role": "user",
                    "content": [
                        {
                            "type": "text",
                            "text": prompt
                        },
                        {
                            "type": "image_url",
                            "image_url": {
                                "url": f"data:image/jpeg;base64,{image_str}" # Specify the image type (jpeg, png, etc.)
                            },
                        },
                    ],
                },
            ],
        ) 

        # print(response.choices[0].message.content)
        dummy_df.loc[len(dummy_df)] = {"Serial No.": sl_no, "Generated Prompt": response.choices[0].message.content}
        print('completed [sl. no. >> {}]'.format(sl_no), "/67")

    # read the source df path and merge the dummy df based on the Serial No. and then export to the same location with a modified.csv tag
    source_df = pandas.read_csv(source_df_path)
    merged_df = source_df.merge(dummy_df, on="Serial No.", how="left")
    merged_df.to_csv(source_df_path.replace(".csv", "_modified.csv"), index=False)
    print("Modified csv saved at >>", source_df_path.replace(".csv", "_modified.csv"))

if __name__ == "__main__":
    import dotenv
    dotenv.load_dotenv()
    print("Prompt generator runner started")
    run_prompt_generator(
        source_df_path="/home/kalki/src/valency/reactor/analysis/data/chemical_reaction_types_completed_urls.csv"
    )
