from smere.smere import SMERE
from utils.logger import Logger

if __name__ == '__main__':
    sm = SMERE(num_observed=91)
    sm.load_records(file = 'data/RL1.csv', fields=[str, str, int, int])
    sm.fit()


