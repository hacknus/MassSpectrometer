import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


filename='Data2/air2.csv'
df = pd.read_csv(filename, delimiter=";", header=21)
cycles = max(df.Cycle)
for i in np.arange(1,8):
    df_mean = df[df.Cycle == i]['Faraday torr']
    plt.plot(df_mean)
    plt.title('cycle {}'.format(i))
    plt.show()
