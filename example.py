from nist import nist_approx
from gas_analysis import gas_analysis





if __name__ == "__main__":
    params = gas_analysis('deo.csv', 'airbaseline.csv')
    air = nist_approx(*params, 'deo')
