from collections import defaultdict
import io
import csv

from utils.logger import Logger


class SMERE:

    def __init__(self, num_observed, al=1, bl=1, max_iter=100):
        # params of beta prior
        self.al = 1
        self.bl = 1

        c = defaultdict(float)

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

        self.logger = Logger().get_logger(file_name='smere')

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
                self.X.append(record[1:])
                # print record[1:]
        except IOError as e:
            self.logger.error(e)
        except AttributeError as e:
            self.logger.error(e)


        self.logger.info("Number of records read = {}".format(len(self.X)))

        self.M_l = [set([record[l] for record in self.X]) for l in range(len(self.fields))]
