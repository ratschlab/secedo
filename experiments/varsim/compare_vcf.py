# Computes how many positions in a VCF file match the significant positions as output by SVC in the following format:
# CHROMOSOME\tPOS

def parse_file(name):
    file1 = open(name)
    list1 = []
    while True:
        line = file1.readline()
        if not line:
            break
        if line.startswith('#'):
            continue
        tokens = line.split('\t')
        if not (tokens[0] == 'X' or tokens[0] == 'Y' or int(tokens[0]) < 23):
            print(tokens[0])
            exit(1)
        if tokens[1][-1] == '\n':
            tokens[1] = tokens[1][:-1]
        list1.append((tokens[0], tokens[1]))
    return list1


print('Reading...')
list1 = parse_file('/Users/dd/work/svc/data/tumor-40K-1_somatic.vcf')  # cosmic.vcf')
print('Reading...')
list2 = parse_file('/Users/dd/work/svc/data/svc_05_01/significant_positions')
print('Sorting...')
list1.sort()
list2.sort()
i = 0
j = 0
matches = 0
while i < len(list1) and j < len(list2):
    while i < len(list1) and list1[i] < list2[j]:
        #print(list1[i])
        #exit(0)
        i += 1
    if i < len(list1) and list1[i] == list2[j]:
        matches += 1
        i += 1

    j += 1

print('Total matches: ', matches, 'out of ', len(list1), '/', len(list2))
