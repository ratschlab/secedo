# Computes how many positions in a VCF file match the significant positions as output by SVC in the following format:
# CHROMOSOME\tPOS

import argparse


def parse_file(name, out=None):
    file1 = open(name)
    list1 = []
    while True:
        line = file1.readline()
        if not line:
            break
        if line.startswith('#'):
            if out:
                out.write(line)
            continue
        tokens = line.split('\t')
        if not (tokens[0] == 'X' or tokens[0] == 'Y' or int(tokens[0]) < 25):
            print('Invalid chromosome name: ', tokens[0])
            exit(1)
        if tokens[0] == '23':
            tokens[0] = 'X'
        elif tokens[0] == '24':
            tokens[0] = 'Y'
        if tokens[1][-1] == '\n':
            tokens[1] = tokens[1][:-1]
        list1.append((tokens[0], int(tokens[1])))
    return list1


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf1', help='First VCF file', default=None)
    parser.add_argument('--vcf2', help='Second VCF file', default=None)
    parser.add_argument('--out', help='Where to write the resulting intersection and diff', default="./res")

    args = parser.parse_args()

    intersect_file = open(args.out + '.intersect', 'w')
    diff_file = open(args.out + '.diff', 'w')
    print(f'Reading {args.vcf1}')
    list1 = parse_file(args.vcf1)
    print(f'Reading {args.vcf2}')
    list2 = parse_file(args.vcf2, intersect_file)
    print('Sorting...')
    list1.sort()
    list2.sort()
    i = 0
    j = 0
    matches = 0
    diffs = 0
    while i < len(list1) and j < len(list2):
        while i < len(list1) and list1[i] < list2[j]:
            diff_file.write(list1[i][0] + '\t' + str(list1[i][1]) + '\n')
            diffs += 1
            i += 1
        if i < len(list1) and list1[i] == list2[j]:
            intersect_file.write(list1[i][0] + '\t' + str(list1[i][1]) + '\n')
            matches += 1
            i += 1

        j += 1

    print('Total matches: ', matches, ' out of ', len(list1), '/', len(list2))
    print(matches/len(list2), matches/len(list1), ' ', matches, ' SNVs\n', len(list1), ' ', len(list2))
    print('Total positions in first file that are not in the 2nd: ', diffs)
    print(f'Intersection written to {args.out}.intersect')
    print(f'Diff written to {args.out}.diff')
