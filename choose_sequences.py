import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from read_data import hist


def sequences(filename, baseline_filename, combine=5, start=0, end=140, relative=True, plot=True):
    path = "DataAvg/" + filename
    df = pd.read_csv(path)
    if baseline_filename == False:
        df_noise = df * 0
    else:
        baseline_path = "DataAvg/" + baseline_filename
        df_noise = pd.read_csv(baseline_path)
        # subtract noise
    df["p"] = df["p"] - df_noise["p"]
    df["err"] = np.sqrt(df["err"] ** 2 + df_noise["err"] ** 2)
    amu, p, err = hist(df, filename, combine, plot=True)
    if relative:
        err = err / sum(p) * 100
        p = p / sum(p) * 100
    end = int(end * 5 / combine)
    start = int(start * 5 / combine)
    if end > len(amu):
        end = len(amu)
    amu = amu[start:end + 1]
    p = p[start:end + 1]
    err = err[start:end + 1]

    if plot == True:
        plt.title(filename.replace(".csv", ""))
        plt.bar(amu, p, width=0.2 * combine, color="red", label=str(df.type[0]))
        plt.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
        plt.legend()
        if relative:
            plt.ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
        else:
            plt.ylabel(r"$p$ [Torr]")
        plt.xlabel("amu")
        plt.savefig("Report/DataResultsPlots/" + filename.replace(".csv", ".pdf"))
        plt.show()
    else:
        return amu, p, err


if __name__ == '__main__':
    sequences("Xenon.csv", "xenonbaseline.csv", 5, 0, 150, False, True)
