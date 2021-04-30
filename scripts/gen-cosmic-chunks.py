# split the cosmic vcf file (containing somatic variations) into several pieces
# this allows us to generate somatic (cancer) simulated reads corresponding to different variants

num_variants = 7  # how many variants to generate
num_lines = 23  # number of lines in each variations

fcoding = open('/Users/dd/Downloads/CosmicCodingMuts.vcf')
fnoncoding = open('/Users/dd/Downloads/CosmicNonCodingVariants.vcf')

chromosomes = ['1'] + list(range(10, 20)) + ['2', '21', '22', 'X']

outf = [open(f'/Users/dd/work/svc/data/cosmic_{i}.vcf', 'w') for i in range(num_variants)]

line1 = fcoding.readline()
while line1.startswith('#'):
    for f in outf:
            f.write(line1)
    line1 = fcoding.readline()

for chr in chromosomes:
    while not line1.startswith(str(chr)):
        line1 = fcoding.readline()

    for j in range(int(num_lines/len(chromosomes))):
        for i in range(num_variants):
            outf[i].write(line1)
            line1 = fcoding.readline()

line2 = fnoncoding.readline()
while line2.startswith('#'):
    for f in outf:
            f.write(line2)
    line2 = fnoncoding.readline()

for chr in chromosomes:
    while not line2.startswith(str(chr)):
        line2 = fnoncoding.readline()

    for j in range(int(num_lines/len(chromosomes))):
        for i in range(num_variants):
            outf[i].write(line2)
            line2 = fnoncoding.readline()

for f in outf:
    f.close()