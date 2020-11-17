import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from read_data import hist

def process_data(filename, baseline_filename):
    path = "DataAvg/" + filename
    baseline_path = "DataAvg/" + baseline_filename
    df = pd.read_csv(path)
    df_noise = pd.read_csv(baseline_path)
    # subtract noise
    df["p"] = df["p"] - df_noise["p"]
    df["err"] = np.sqrt(df["err"]**2 + df_noise["err"]**2)
    hist(df, filename, combine=1, plot=True)
    

if __name__ == '__main__':
    process_data("Xenon.csv", "xenonbaseline.csv")