import getopt
import re
import math
import numpy as np
import sys

# thresholds K to use for pre-processing; values for coverage 10, 20, 30, ... 200
# always the largest values for any tumour proportion, rounded up to one decimal place
# Ks = [1.42, 3.38, 4.61, 5.42, 5.24, 5.04, 5.24, 5.11, 3.77, 4.06, 4.33, 4.55, 4.69, 4.89, 5.06, 5.22, 5.31, 2.78, 3.17, 3.66]
# thresholds for tumour proportion 0.5
Ks = [0.872, 0.872, 0.872, -1.238, -1.358, -5.749, -10.139, -13.912, -25.998, -33.936, -34.469, -54.084, -44.453,
      -56.842, -63.352, -88.297, -90.711, -96.841, -115.730, -112.601]


# log factorial(n) using Stirling formula
# !!! DOES NOT CHECK IF N IS INTEGER !!!
def log_fact(n):
    if n > 170:
        return 0.5 * math.log(2 * math.pi * n) + n * math.log(n / math.e)
    else:
        return math.log(math.factorial(n))


# function to do pre-processing; Bayesian approach, homoOnly
# returns True -- the position is kept, or False -- the position is discarded
# data is a vector of length 4 with counts of A,C,G, and T in the pooled data at a fixed position
def preprocessing(data, theta, K):
    logTheta = math.log(theta / 3)
    logOneMinusTheta = math.log(1 - theta)

    # sort data; sorts in ascending order
    data = np.sort(data)
    # total coverage
    cov = np.sum(data)
    # if there are no reads or if we have just one base, do not keep
    if cov > 0 and np.count_nonzero(data) > 1:
        # the multinomial coefficient
        logMultinom_coef = log_fact(cov) - log_fact(data[0]) - log_fact(data[1]) - log_fact(data[2]) - log_fact(data[3])
        # log P(data|wt) for wt=homo
        logProb_homo = logMultinom_coef + data[3] * logOneMinusTheta + (cov - data[3]) * logTheta
        # add prior on the most probable genotype (1/4, because all four homozygous genotypes are equally likely)
        logProb_homo += math.log(1 / 4)
        # add prior on null hypothesis (0.998)
        logProb_homo += math.log(0.998)
        # the normalizing coefficient (normalizing for read depth)
        logNorm_coef = log_fact(cov + 3) - math.log(6) - log_fact(cov)
        # the result
        print('log logMultinom_coef ', logMultinom_coef)
        print('logNorm_coef ', logNorm_coef)
        print('log prob homo ', logProb_homo, 'logNorm_coef ', logNorm_coef)
        res = logNorm_coef + logProb_homo
        if res < K:
            return True
        else:
            return False
    else:
        return False


def main(argv):
    print(preprocessing([51, 1, 0, 0], 0.01, Ks[4]))
    # parse the arguments
    try:
        opts, args = getopt.getopt(argv, "n:t:")
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-n':
            ncells = int(arg)
        elif opt == '-t':
            theta = float(arg)

    # the input is coming on stdin
    for line in sys.stdin:
        # if line is empty, end of file is reached
        if not line:
            break
        # process the line
        else:
            line = line.rstrip("\n")
            # split by tabs
            splitLine = line.split("\t")
            # first column: chromosome
            chrom = splitLine[0]
            # second column: position
            position = splitLine[1]

            bases = ""
            readNames = ""
            cells = ""
            for i in range(ncells):
                cov = int(splitLine[3 + 4 * i])
                if cov > 0:
                    bases_i = splitLine[4 + 4 * i]
                    readNames_i = splitLine[6 + 4 * i].split(",")

                    # in "bases", there are some characters we don't want
                    # starts of reads are denoted by ^ and the next character denotes the mapping quality -- delete those
                    bases_i = re.sub('\^.', '', bases_i)

                    # ends of reads denoted by $
                    bases_i = re.sub('\$', '', bases_i)

                    # delete the insertions - + followed by length of the insertion and the inserted bases_i, e.g. +1a
                    p = re.compile('\+[0-9]+')
                    while True:
                        m = re.search(p, bases_i)
                        # no more match
                        if m == None:
                            break
                        else:
                            # position in the string
                            pos = m.start()
                            # length of the insertion
                            ins_len = int(m.group()[1:])
                            # delete the +, the number and the inserted bases
                            bases_i = bases_i[:pos] + bases_i[pos + ins_len + len(m.group()):]

                    # delete the deletions: - followed by length of the insertion and n's
                    p = re.compile('-[0-9]+')
                    while True:
                        m = re.search(p, bases_i)
                        # no more match
                        if m == None:
                            break
                        else:
                            # position in the string
                            pos = m.start()
                            # length of the deletion
                            del_len = int(m.group()[1:])
                            # delete the +, the number and the inserted bases_i
                            bases_i = bases_i[:pos] + bases_i[pos + del_len + len(m.group()):]

                    # if there is "*", the base was deleted; delete the corresponding read names, and then also delete the stars
                    star_pos = [x for x in range(len(bases_i)) if bases_i[x] == "*"]
                    tmp = [readNames_i[x] for x in range(len(readNames_i)) if x not in star_pos]
                    readNames_i = tmp
                    bases_i = re.sub('\*', '', bases_i)

                    # if we have some properly sequenced bases, append to the existing strings
                    if len(bases_i) > 0:
                        if bases != "":
                            bases = bases + bases_i
                            readNames = ','.join([readNames, ','.join(readNames_i)])
                            cells = ','.join([cells, ','.join([str(i)] * len(readNames_i))])
                        else:
                            bases = bases_i
                            readNames = ','.join(readNames_i)
                            cells = ','.join([str(i)] * len(readNames_i))

            # ---- pre-processing ----
            # count number of As, Cs, Gs, Ts
            nAs = bases.count('a') + bases.count('A')
            nCs = bases.count('c') + bases.count('C')
            nGs = bases.count('g') + bases.count('G')
            nTs = bases.count('t') + bases.count('T')
            # print([nAs, nCs, nGs, nTs])
            # coverage
            cov = nAs + nCs + nGs + nTs
            # print(cov)
            # choose K for the closest coverage
            i = round(cov / 10) - 1
            if i == -1:
                i = 0
            elif i > 19:
                i = 19
            K = Ks[i]
            # do pre-processing
            if preprocessing([nAs, nCs, nGs, nTs], theta, K):
                # write the line to the pre-processed file
                line = '\t'.join([chrom, position, str(cov), bases, cells, readNames])
                print(line)
        # print("kept")
    # else:
    # print("discarded")


if __name__ == "__main__":
    main(sys.argv[1:])
