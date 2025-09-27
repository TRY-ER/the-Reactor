#!/usr/bin/env python3
"""
Python script to make a request to Groq API
Converts the curl command to Python using requests library
"""
import os
import requests
import openai
import json


def make_groq_request(model_name: str,
                      api_key: str,
                      query: str):
    """Make a request to Groq API"""
    # API endpoint
    url = "https://api.groq.com/openai/v1/responses"  # for response

    # Headers
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}"
    }

    # Request payload
    data = {
        "model": model_name,
        "input": query,
    }

    try:
        # Make the POST request
        response = requests.post(url, headers=headers, json=data)

        # Check if request was successful
        response.raise_for_status()

        # Parse and return JSON response
        return response.json()

    except requests.exceptions.RequestException as e:
        print(f"Error making request: {e}")
        if hasattr(e, 'response') and e.response is not None:
            print(f"Response status code: {e.response.status_code}")
            print(f"Response content: {e.response.text}")
        return None


def make_openai_groq_request(model_name: str,
                             api_key: str,
                             query: str):
    client = openai.OpenAI(
        api_key=api_key,
        base_url="https://api.groq.com/openai/v1"
    )

    response = client.responses.create(
        model=model_name,
        input=query
    )

    return response

# if __name__ == "__main__":
#     import dotenv
#     dotenv.load_dotenv()

#     import os
#     api_key = os.getenv('GROQ_API_KEY', "<key>")

#     # # Make the request
#     # result = make_groq_request(model_name="openai/gpt-oss-120b",
#     #                            api_key=api_key,
#     #                            query="how are you?")

#     result = make_openai_groq_request(model_name="openai/gpt-oss-120b",
#                                        api_key=api_key,
#                                        query="how are you?")

#     if result:
#         # Pretty print the response
#         print("Response from Groq API:")
#         try:
#             print(json.dumps(result, indent=2))
#         except:
#             print("Raw response >>", result)
#             print("type of Raw response >>", type(result))

#     else:
#         print("Failed to get response from Groq API")
