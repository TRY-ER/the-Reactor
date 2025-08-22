import pandas as pd
import time
import os

class PersistentDataLogger:
    def __init__(self, log_dir: str, filename: str, headers: list):
        self.log_dir = log_dir
        self.filename = filename
        self.headers = headers
        self.check_create_csv()


    def check_create_csv(self):
        file_path = os.path.join(self.log_dir, self.filename)
        if not os.path.exists(self.log_dir):
            os.makedirs(self.log_dir)
        if os.path.isfile(file_path):
            exist_df = pd.read_csv(file_path)
            if all([x in self.headers for x in exist_df.columns]):
                self.target_csv = exist_df
            else:
                empty_df = pd.DataFrame(columns=self.headers)
                self.target_csv = empty_df
                self.filename = self.filename[:-3]+"_"+str(time.time())+".csv"
                empty_df.to_csv(os.path.join(self.log_dir, self.filename),index=False)
        else:
            empty_df = pd.DataFrame(columns=self.headers)
            self.target_csv = empty_df
            empty_df.to_csv(os.path.join(self.log_dir, self.filename),index=False)

    def log_one(self, data: dict):
        new_row = pd.DataFrame([data])
        self.target_csv = pd.concat([self.target_csv, new_row], ignore_index=True)
        self.target_csv.to_csv(os.path.join(self.log_dir, self.filename), index=False)

    def log_df(self, df: pd.DataFrame):
        self.target_csv = pd.concat([self.target_csv, df], ignore_index=True)
        self.target_csv.to_csv(os.path.join(self.log_dir, self.filename), index=False)


if __name__ == "__main__":
    logger = PersistentDataLogger(
        log_dir="/home/kalki/src/valency/reactor/backend/app/observers/one_shot/outputs/observation_0",
        filename="findings.csv",
        headers=["column1", "column2", "column3"]
    )

    # testing logging one
    logger.log_one({"column1": "value1", "column2": "value2", "column3": "value88"})

    # testing logging dataframe
    # logger.log_df(pd.DataFrame([{"column1": "value4", "column2": "value5", "column3": "value6"}]))