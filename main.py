import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from read_data import hist

def process_data(filename):
    path = "DataAvg/" + filename
    df = pd.read_csv(path)
    df_noise = pd.read_csv("DataAvg/xenonbaseline.csv")
    # subtract noise
    df["p"] = df["p"] - df_noise["p"]
    df["err"] = df["err"] - df_noise["err"]
    hist(df, filename, 1, plot=True)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    process_data("Xenon.csv")