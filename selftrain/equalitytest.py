import codecs


count = 0
file1 = codecs.open('datasets/dataB/dataset-B.fasta', 'r')
name1 = ''
for line1 in file1:
    if '>' in line1:
        name1 = line1#.strip('>0123456789\n')
        continue
    file2 = codecs.open('datasets/dataB/dataset-B.fasta', 'r')
    name2 = ''
    for line2 in file2:
        if '>' in line2:
            # print('tick')
            name2 = line2#.strip('>0123456789\n')
            continue
        if line1 == line2 and (not name1 == name2):
            count += 1
            print(line1)
            print(name1)
            print(name2)
print(count)
