import random


def sample_index_from_multinomial_list(multinomial):
    mark = random.random()
    cummulative = 0.0
    for i in range(len(multinomial)):
        cummulative += multinomial[i]
        if cummulative > mark:
            return i
    assert (False)


def sample_key_from_multinomial_dict(multinomial):
    mark = random.random()
    cummulative = 0.0
    for k in multinomial.keys():
        cummulative += multinomial[k]
        if cummulative > mark:
            return k
