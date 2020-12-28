import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from read_data import hist
from choose_sequences import sequences
from voltage_variation import double_plot, sequence_plot
from utils import fit_peak
from nist import get_nist_peaks, nist_approx
from gas_analysis import gas_analysis

nist_approx(*gas_analysis('xenon_highres.csv', 'xenonbaseline_highres.csv'), 'xenon')
nist_approx(*gas_analysis('xenonbaseline_highres.csv', False), 'resiudal gas')
nist_approx(*gas_analysis('argon2.csv', 'argonbaseline.csv'), 'argon')
nist_approx(*gas_analysis('mix2.csv', 'mix_baseline2.csv'), 'mix')
nist_approx(*gas_analysis('co2_premature_balloon_loss.csv', 'co2baseline.csv'), r'CO_2', True)
nist_approx(*gas_analysis('deo.csv', 'airbaseline.csv'), 'Deo')
nist_approx(*gas_analysis('ethanol.csv', 'co2baseline.csv'), 'Ethanol', True)

nist_approx(*gas_analysis('air00.csv', 'airbaseline.csv'), 'air00')
nist_approx(*gas_analysis('air1.csv', 'airbaseline.csv'), 'air1')
nist_approx(*gas_analysis('air2.csv', 'airbaseline.csv'), 'air2')
nist_approx(*gas_analysis('air3.csv', 'airbaseline.csv'), 'air3')
nist_approx(*gas_analysis('air4.csv', 'airbaseline.csv'), 'air4')
nist_approx(*gas_analysis('air5.csv', 'airbaseline.csv'), 'air5')

nist_approx(*gas_analysis('vodka_pure.csv', 'air_background.csv', False), 'vodka pure', True)

nist_approx(*gas_analysis('vodka.csv', 'air_background.csv', False), 'vodka', True)
