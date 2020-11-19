from jcamp import JCAMP_reader
import pandas as pd


def get_nist_peaks(name, p_number=1):
    path = "NIST/" + name + ".jdx"
    d = JCAMP_reader(path)
    m = d["x"]
    y = d["y"]
    nist = {"m": m,
            "y": y}
    n = pd.DataFrame(data=nist)
    n.sort_values("y", inplace=True, ascending=False)
    n.reset_index()
    n = n[:p_number]
    return n


if __name__ == "__main__":
    print(get_nist_peaks("oxygen", 2).head())
