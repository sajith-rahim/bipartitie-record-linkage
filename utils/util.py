def list_to_unit_vector(alist):
    sum = sum(alist)
    if sum == 0:
        print_linkage_structure()
    assert sum != 0
    for i in range(len(alist)):
        alist[i] /= sum


def dict_to_unit_vector(adict):
    sum = sum(adict.values())
    for k in adict.keys():
        adict[k] /= sum


def print_linkage_structure(X, y, lambda_inverse, logger_):
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
