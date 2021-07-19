# Computes how many positions in a VCF file match the significant positions as output by SVC in the following format:
# CHROMOSOME\tPOS

import argparse


def uniq(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not ((x[0], x[1]) in seen or seen_add((x[0], x[1])))]

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
        if tokens[0] != 'X' and tokens[0] != 'Y' and (not tokens[0].isdigit() or int(tokens[0]) > 24):
            # print('Invalid chromosome name: ', tokens[0])
            # exit(1)
            continue
        if tokens[0] == '23':
            tokens[0] = 'X'
        elif tokens[0] == '24':
            tokens[0] = 'Y'
        if tokens[1][-1] == '\n':
            tokens[1] = tokens[1][:-1]
        list1.append((tokens[0], int(tokens[1]), '\t'.join(tokens[2:])[:-1]))
    print(len(list1), ' positions.')
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
    list1.sort(key=lambda c: (c[0], c[1]))
    list2.sort(key=lambda c: (c[0], c[1]))

    print('De-duping')
    list1 = uniq(list1)
    list2 = uniq(list2)
    print(f'De-duped lists have {len(list1)} and {len(list2)} items respectively')
    i = 0
    j = 0
    matches = 0
    diffs = 0
    while i < len(list1) and j < len(list2):
        while i < len(list1) and (
                list1[i][0] < list2[j][0] or (list1[i][0] == list2[j][0] and list1[i][1] < list2[j][1])):
            diff_file.write(list1[i][0] + '\t' + str(list1[i][1]) + '\t' + list1[i][2] + '\t' + list2[j][2] + '\n')
            diffs += 1
            i += 1
        if i < len(list1) and list1[i][0] == list2[j][0] and list1[i][1] == list2[j][1]:
            intersect_file.write(
                list1[i][0] + '\t' + str(list1[i][1]) + '\t' + list1[i][2] + '\t' + list2[j][2] + '\n')
            matches += 1
            i += 1

        j += 1

print('Total matches: ', matches, ' out of ', len(list1), '/', len(list2))
print(matches / len(list2), matches / len(list1), ' ', matches, ' SNVs\n', len(list1), ' ', len(list2))
print('Total positions in first file that are not in the 2nd: ', diffs)
print(f'Intersection written to {args.out}.intersect')
print(f'Diff written to {args.out}.diff')
