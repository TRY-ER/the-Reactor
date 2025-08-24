import pandas as pd
from os import path


def provide_reaction_details(reaction_df_path,
                             column_to_access: list[str] | None = None,
                             custom_range: list[int] | None = None):
    if isinstance(reaction_df_path, str) and path.exists(reaction_df_path) and reaction_df_path.endswith('.csv'):
        df = pd.read_csv(reaction_df_path)
        # yield an iterator with a dictionary containing column name and the value same as column_to_access list
        for index, row in df.iterrows():
            if custom_range and index not in custom_range:
                continue
            yield {col: row[col] for col in column_to_access} if column_to_access else row.to_dict()

# if __name__ == "__main__":
#     for values in provide_reaction_details(
#         "/home/kalki/src/valency/reactor/analysis/data/chemical_reaction_types_completed_urls.csv",
#         column_to_access=["Reaction Type", "Source"],
#         custom_range=[i for i in range(63, 67)]):
#         print("data received", values)