# Sorts a VCF file (for human) by position
# It is assumed that the VCF file is sorted lexicographically by an inane tool such as bcftools, i.e. the order of
# the chromosomes is [1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 3, 4, 5, 6, 7, 8, 9, X, Y] and this
# script wil rearrange the chromosomes to [1-22, X, Y]

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--i', help='Input VCF file, sorted lexicographically by position', default=None)
    parser.add_argument('--o', help='Output VCF file, sorted numerically by position', default="./sorted.vcf")

    args = parser.parse_args()

    if not args.i:
        print('Please specify an input file via --i <input_file>')
        exit(1)

    f = open(args.i)
    out = open(args.o, 'w')

    print("First pass, writing chromosomes 1-10...")

    for line in f.readlines():
        if line[0] == '#':
            out.write(line)
        if line[1] == '\t' and line[0] != 'X' and line[0] != 'Y':
            out.write(line)

    print('Second pass, writing chromosomes 11-22,X,Y')
    f = open(args.i)

    for line in f.readlines():
        if line[0] != '#' and (line[2] == '\t' or line[0] == 'X' or line[0] == 'Y'):
            out.write(line)

    out.close()

    print('Done.')
