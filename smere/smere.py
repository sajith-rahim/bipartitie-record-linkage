from collections import defaultdict
import io
import csv

import scipy
import scipy.stats
import scipy.special

from metrics.levenshtein import levenshtein_ratio
from sampling.sampler import sample_index_from_multinomial_list, sample_key_from_multinomial_dict
from utils.logger import Logger
from utils.util import list_to_unit_vector, dict_to_unit_vector, print_linkage_structure


class SMERE:

    def __init__(self, num_observed, al=1, bl=1, max_iter=100):
        # params of beta prior
        self.al = 1
        self.bl = 1

        self.c = [-1.0, -1.0]

        if num_observed is not None:
            self.num_observed = num_observed
        else:
            raise ValueError("N: number of observed individuals not defined!")

        self.max_iter = 33

        self.z = []
        '''
        self.z : list of list 
        z[j][l] = true => x[j][l] is a distorted form of y[Lambda[j]][l]
        z[j][l] = false implies that x[j][l] = y[Lambda[j]][l], but
        z[j][l] = true doesn't imply that x[j][l] != y[Lambda[j]][l]
        '''
        self.y = []
        '''
        self.y : list
        list of latent vectors/records
        '''
        self.lambda_ = []
        '''
        self.lambda_ : list
        list of "pointers" to the latent record that corresponds to an observed record
        lambda_[i] = j <=> observed record x[i] was generated from latent record y[j]
        '''

        self.lambda_inverse = []
        '''
         self.lambda_inverse :list of sets.
         for each latent record y[i], lambda_inverse[i] gives the indexes {j} of aligned
         observed records x[j].
         lambda_inverse[j] = set([i1, i2, i3]) <=> lambda_[i1] = lambda_[i2] = lambda_[i3] = j, lambda_[else] != j
        '''

        self.empirical = []
        '''
        the empirical distribution for each field in the observed records.
        empirical[l] is a dictionary that maps a field value to the percentage of records with that value
        '''

        self.distortion = {}
        '''
        distortion[(l,y')] => the distortion distribution of x' | y' for field l which is a defaultdict
        '''
        self.M_l = []
        '''
        list of sets. M_l[l] is the set of observed values at field l.
        '''

        self.X = []  # records

        self.max_iter = max_iter

        self.logger = Logger().get_logger(file_name='smere')
        self.prob_logger = Logger().get_logger(file_name='sampling')
        self.gibbs_logger = Logger().get_logger(file_name='gibbs')

    def load_records(self, file, fields):
        """ fields: data type of fields"""
        self.fields = fields
        try:
            csv_file = io.open(file, encoding='utf8', mode='r')
            csv_reader = csv.reader(csv_file, delimiter=',')
            headers = next(csv_reader)
            assert len(headers) == len(fields) + 1

            for record in csv_reader:
                for i in range(1, len(record)):
                    _field = record[i]
                    if self.fields[i - 1] == int:
                        _field = int(_field)
                    elif self.fields[i - 1] == float:
                        _field = float(_field)
                    elif self.fields[i - 1] == str:
                        pass
                    else:
                        self.logger.error("Unsupported datatype in record: {}".format(record));
                        raise TypeError("Unsupported datatype in record!")
                    record[i] = _field
                self.X.append(record[1:])
                # print record[1:]
        except IOError as e:
            self.logger.error(e)
        except AttributeError as e:
            self.logger.error(e)

        self.logger.info("Number of records read = {}".format(len(self.X)))

        self.M_l = [set([record[l] for record in self.X]) for l in range(len(self.fields))]

    def initialize_latent_variables(self):
        """ initialize y, lambda and z"""
        # copy the latent records from the read records
        self.y = [list(self.X[record]) for record in range(self.num_observed)]
        _num_of_records = len(self.X)
        # each observed points at the corresponding latent;
        # if the number of unique latents (num_observed) < len(x), then use num_observed-1
        self.lambda_ = [min(i, self.num_observed - 1) for i in range(_num_of_records)]
        self.lambda_inverse = [set() for i in range(self.num_observed)]
        for i in range(_num_of_records):
            self.lambda_inverse[self.lambda_[i]].add(i)

            # initially, no distortion happens
            self.z.append([])
            for l in range(len(self.fields)):
                if i < self.num_observed:
                    self.z[i].append(False)
                else:
                    # except when Lambda[i] != i
                    self.z[i].append(True)

    def compute_lambda_posteriors(self, posteriors, flags, i):
        for j in range(self.num_observed):
            # compute the posterior of p(lambda_i = j | everything else)
            posterior = 1.0
            for l in range(len(self.fields)):
                if not self.z[i][l]:
                    delta = 1.0 if self.X[i][l] == self.y[j][l] else 0.0
                    self.prob_logger.info(
                        "z_{}_{} is False ==> X[i][l] = {} y[j][l] = {} delta = {}".
                            format(i, l, self.X[i][l], self.y[j][l], delta))
                    posterior *= scipy.special.beta(self.al, self.bl + 1) * delta
                else:
                    _distortion = self.distortion[(l, self.y[j][l])][self.X[i][l]]
                    posterior *= scipy.special.beta(self.al + 1, self.bl) * _distortion
            posteriors[j] = posterior
            flags[j] = True

    def resample_lambda(self, i):

        # first, compute the posterior distribution of lambda_i
        posteriors = [0.0 for index in range(self.num_observed)]
        flags = [False for index in range(self.num_observed)]
        self.compute_lambda_posteriors(posteriors, flags, i)

        if sum(posteriors) == 0:
            raise ValueError("posterior --> 0")
        # normalize the posteriors
        list_to_unit_vector(posteriors)

        # now that we computed the posteriors, sample a value for lambda_i
        old_lambda_i = self.lambda_[i]

        lambda_i = sample_index_from_multinomial_list(posteriors)
        self.lambda_[i] = lambda_i

        # update inverse lambda
        self.lambda_inverse[old_lambda_i] -= set([i])
        self.lambda_inverse[self.lambda_[i]] |= set([i])

    def resample_z(self, i, l):
        # first, compute the posterior distribution of z_{i,l}
        delta = 1.0 if self.X[i][l] == self.y[self.lambda_[i]][l] else 0.0
        posteriors = []
        posterior_prob_of_z_eq_zero = scipy.special.beta(self.al, self.bl + 1) * delta
        posteriors.append(posterior_prob_of_z_eq_zero)
        _distortion = self.distortion[(l, self.y[self.lambda_[i]][l])][self.X[i][l]]
        posterior_prob_of_z_eq_one = scipy.special.beta(self.al + 1, self.bl) * _distortion
        posteriors.append(posterior_prob_of_z_eq_one)

        # normalize the posteriors
        list_to_unit_vector(posteriors)

        # now sample a value for z_{i,l}
        if sample_index_from_multinomial_list(posteriors) == 0:
            self.z[i][l] = False
        else:
            self.z[i][l] = True

    def resample_y(self, i, l):
        # find the indexes of observed records which are currently aligned to this latent record
        aligned_observed_record_indexes = self.lambda_inverse[i]

        # first, compute posteriors
        posteriors = defaultdict(float)
        for v in self.M_l[l]:
            posteriors[v] = 1.0
            for j in aligned_observed_record_indexes:
                if not self.z[j][l]:
                    delta = 1.0 if self.X[j][l] == v else 0.0
                    posteriors[v] *= float(scipy.special.beta(self.al, self.bl + 1)) * delta
                else:
                    posteriors[v] *= float(scipy.special.beta(self.al + 1, self.bl)) * self.distortion[(l, v)][
                        self.X[j][l]]

        # normalize posteriors
        dict_to_unit_vector(posteriors)

        # now, sample a value from the posterior
        self.y[i][l] = sample_key_from_multinomial_dict(posteriors)

    def score(self, l, x_value, y_value):
        if self.fields[l] == str:
            return self.c[0] * levenshtein_ratio(x_value, y_value)
        else:
            print("{}:{}-{}:{}".format(x_value,type(x_value),y_value,type(y_value)))
            return self.c[1] * (x_value - y_value) * (x_value - y_value)

    def precompute_distortions(self):
        for l in range(len(self.fields)):
            for y_value in self.M_l[l]:
                self.distortion[(l, y_value)] = defaultdict(float)
                for x_value in self.M_l[l]:
                    self.distortion[(l, y_value)][x_value] = self.score(l, x_value, y_value)
                sum_ly = sum(self.distortion[(l, y_value)].values())
                for x_value in self.M_l[l]:
                    self.distortion[(l, y_value)][x_value] /= sum_ly

    def precompute_empiricals(self):
        for l in range(len(self.fields)):
            self.empirical.append(defaultdict(float))
        for i in range(len(self.X)):
            for l in range(len(self.fields)):
                self.empirical[l][self.X[i][l]] += 1.0
        for l in range(len(self.fields)):
            sum_l = sum(self.empirical[l].values())
            for value in self.empirical[l].keys():
                self.empirical[l][value] /= sum_l

    def fit(self, max_iter=100):

        self.max_iter = max_iter

        self.precompute_distortions()
        self.precompute_empiricals()
        self.initialize_latent_variables()

        iter = 0
        while True:

            for i in range(len(self.X)):
                self.resample_lambda(i)

                for l in range(len(self.fields)):
                    self.resample_z(i, l)

            for j in range(self.num_observed):
                for l in range(len(self.fields)):
                    self.resample_y(j, l)

            if iter > self.max_iter:
                break
            else:
                print("{}/{}".format(iter, self.max_iter))
                print_linkage_structure(self.X, self.y, self.lambda_inverse, self.gibbs_logger)

        print("Converged")
        print_linkage_structure(self.X, self.y, self.lambda_inverse)
