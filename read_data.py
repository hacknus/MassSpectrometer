import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def hist(df, file, combine=1, plot=False, start=0, end=-1):
    # could have used np.hist function...    
    amu0 = np.array(df.amu)
    # amu0 = np.concatenate(
    #     (np.array([0.6, 0.8]), amu0, np.array([len(amu0) * 0.2 + 1, (len(amu0) + 1) * 0.2 + 1])))  # complete vector
    n = len(amu0) % combine
    if n == 0:
        n = combine
    amu = np.array(amu0)[:-n:combine]
    p0 = np.array(df.p)
    # p0 = np.concatenate(([p0[2]], p0, [p0[-2],p0[-3],p0[-4]]))  # startet nicht bei p0[2] da daten satz um eins verschoben ist
    p = np.array(p0)[:-n:combine]
    err_quat0 = np.array(df.err) ** 2
    # err_quat0 = np.concatenate(([err_quat0[2]], err_quat0, [err_quat0[-2], err_quat0[-3],err_quat0[-4]]))
    err_quat = np.array(err_quat0)[:-n:combine]
    for i in range(1, combine):
        amu += amu0[i:-n:combine]
        p += p0[i:-n:combine]
        err_quat += err_quat0[i:-n:combine] ** 2
    amu /= combine
    err = np.sqrt(err_quat)
    if not plot:
        plt.title(file.replace(".csv", ""))
        plt.bar(amu, p, width=0.2 * combine, color="red", label=str(df.type[0]))
        plt.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
        plt.legend()
        plt.ylabel(r"$p$ [Torr]")
        plt.xlabel("amu")
        plt.savefig("Report/DataAvgPlots/" + file.replace(".csv", ".pdf"))
        plt.clf()
    else:
        return amu, p, err


def avg_all(new=False):
    if new:
        path = "Data2"
    else:
        path = "Data"
    files = os.listdir(path)
    for file in files:
        filename = f"{path}/" + file
        df = convert_file(filename, True)
        df.to_csv(f"{path}Avg/" + file, index=0)
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
            p1 = np.array(df[df.Cycle == 1]["SEM torr"])
            d = {"amu": amu,
                 "p1": p1,
                 "type": ["SEM"] * len(amu)}
            for c in range(2, cycles + 1):
                d[f"p{c}"] = np.array(df[df.Cycle == c]["SEM torr"])
    else:
        i=2
        if filename[0:42]=='Data/restgasspektrum_FARx7_20prozent_5res_' :
            print('true')
            i=3
        if filename == 'Data2/xenonbaseline_highres.csv':
            i=3
        if filename == 'Data2/mix_baseline2':
            i = 3
        if filename == 'co2_premature_balloon_loss':
            i = 0
        df_mean = df[df.Cycle > i].groupby("mass amu").mean()
        df_std = df[df.Cycle > i].groupby("mass amu").std()
        try:
            d = {"amu": amu,
                 "p": np.array(df_mean["Faraday torr"]),
                 "err": np.array(df_std["Faraday torr"]),
                 "type": ["Faraday"] * len(amu)}
        except KeyError:
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
    avg_all(new=True)
