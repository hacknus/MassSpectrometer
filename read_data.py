import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def hist(df, file, combine=1, plot=False):
    # could have used np.hist function...
    n = len(df.amu) % combine
    if n == 0:
        n = combine
    amu = np.array(df.amu)[:-n:combine]
    p = np.array(df.p)[:-n:combine]
    err_quat = np.array(df.err)[:-n:combine] ** 2
    for i in range(1, combine):
        amu += np.array(df.amu)[i:-n:combine]
        p += np.array(df.p)[i:-n:combine]
        err_quat += np.array(df.err)[i:-n:combine] ** 2
    amu /= combine
    err = np.sqrt(err_quat)

    plt.title(file.replace(".csv", ""))
    plt.bar(amu, p, width=0.2 * combine, color="red", label=str(df.type[0]))
    plt.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')

    plt.legend()
    #    plt.semilogy()
    plt.ylabel(r"$p$ [Torr]")
    plt.xlabel("amu")
    if not plot:
        plt.savefig("Report/DataAvgPlots/" + file.replace(".csv", ".pdf"))
        plt.clf()
    else:
        plt.show()


def avg_all():
    files = os.listdir("Data")
    for file in files:
        filename = "Data/" + file
        df = convert_file(filename, True)
        df.to_csv("DataAvg/" + file, index=0)
        print(f"completed {file}")
        # make a histogram of 50 bars
        # hist(df, file, len(df)//50)
        # do not create a histogram (bin combining = 1), just a bar plot
        hist(df, file, 1)


def convert_file(filename, average=False):
    df = pd.read_csv(filename, delimiter=";", header=21)
    cycles = max(df.Cycle)
    amu = np.array(df[df.Cycle == 1]["mass amu"])
    if not average:
        try:
            p1 = np.array(df[df.Cycle == 1]["Faraday torr"])
            d = {"amu": amu,
                 "p1": p1,
                 "type": ["Faraday"] * len(amu)}
            for c in range(2, cycles + 1):
                d[f"p{c}"] = np.array(df[df.Cycle == c]["Faraday torr"])
        except KeyError:
            print("found SEM")
            p1 = np.array(df[df.Cycle == 1]["SEM torr"])
            d = {"amu": amu,
                 "p1": p1,
                 "type": ["SEM"] * len(amu)}
            for c in range(2, cycles + 1):
                d[f"p{c}"] = np.array(df[df.Cycle == c]["SEM torr"])
    else:
        df_mean = df[df.Cycle > 2].groupby("mass amu").mean()
        df_std = df[df.Cycle > 2].groupby("mass amu").std()
        try:
            d = {"amu": amu,
                 "p": np.array(df_mean["Faraday torr"]),
                 "err": np.array(df_std["Faraday torr"]),
                 "type": ["Faraday"] * len(amu)}
        except KeyError:
            print("found SEM")
            d = {"amu": amu,
                 "p": np.array(df_mean["SEM torr"]),
                 "err": np.array(df_std["SEM torr"]),
                 "type": ["SEM"] * len(amu)}
    df = pd.DataFrame(data=d)
    return df


def plot(df):
    """
    for i in range(1, cycles + 1):
        plt.scatter(df.amu, df[f"p{i}"])
    """
    plt.semilogy()
    plt.bar(df.amu, df.p, yerr=df.err, color="red")
    # plt.errorbar(df.amu, df.p, yerr=df.err, color="red", label="Data", fmt='o', markeredgecolor="black",
    #             ecolor='black', capthick=2, capsize=2, elinewidth=1, markeredgewidth=0.5, ms=3)
    plt.ylim(1e-10, 0.6e-6)
    plt.xlabel("amu")
    plt.ylabel("p [Torr]")
    plt.show()


if __name__ == "__main__":
    avg_all()
