from smere.smere import SMERE
from utils.logger import Logger

if __name__ == '__main__':
    sm = SMERE(num_observed=91)
    sm.load_records(file = 'data/RLdata500.csv', fields=[str, str, int, int, int])
    sm.fit(33)
    print(sm.lambda_)

