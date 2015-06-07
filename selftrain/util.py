import codecs
import prefs


def numeratesequences():
    filtered = codecs.open(
        'datasets/a_domains_filtered.fa', 'r', encoding=prefs.ENCODING)
    out = codecs.open(
        'datasets/a_domains_enumerated.fa', 'w', encoding=prefs.ENCODING)
    sequences = []
    for line in filtered:
        if not line.startswith('>'):
            sequences.append(line)
    sequencesset = set(sequences)
    print(len(sequences))
    print(len(sequencesset))

    sequencenumber = 0
    for sequence in sequencesset:
        out.write('>seq%i\n%s' % (sequencenumber, sequence))
        # out.write(line)
        sequencenumber += 1


if __name__ == '__main__':
    numeratesequences()
