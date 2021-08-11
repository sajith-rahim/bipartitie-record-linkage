def list_to_unit_vector(plist: list):
    _sum = sum(plist)
    if _sum == 0:
        print_linkage_structure()
    assert _sum != 0
    for i in range(len(plist)):
        plist[i] /= _sum


def dict_to_unit_vector(pdict: dict):
    _sum = sum(pdict.values())
    for k in pdict.keys():
        pdict[k] /= _sum


def print_linkage_structure(X, y, lambda_inverse, logger_=None):
    for j in range(len(lambda_inverse)):
        if len(lambda_inverse[j]) > 0:
            obs = '\n'.join([str(X[i]) for i in lambda_inverse[j]])
            if logger_ is not None:
                logger_.info('latent_index={}\nobserved_indexes={}\nlatent_record={}\nobserved_records={}\n\n'
                        .format(j,
                                lambda_inverse[j],
                                y[j],
                                obs
                                ))
            else:
                print(
                    'latent_index={}\nobserved_indexes={}\nlatent_record={}\nobserved_records={}\n\n'
                        .format(j,
                                lambda_inverse[j],
                                y[j],
                                obs
                                ))
