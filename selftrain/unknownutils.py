import codecs
import prefs
import makefolds
import copy
from collections import OrderedDict
import re


def sequencesfromfasta():
    fastafile = codecs.open(
        'datasets/a_domains_enumerated.fa', 'r', encoding=prefs.ENCODING)
    sequences = []
    for line in fastafile:
        if not line.startswith('>'):
            sequences.append(line)
    return sequences


def makecombinedfastafile(dataset, unknownsequences, resfasta):
    datasetfile = codecs.open(
        dataset, 'r', encoding=prefs.ENCODING)
    unknownsequences = codecs.open(
        unknownsequences, 'r', encoding=prefs.ENCODING)
    resfile = codecs.open(
        resfasta, 'w', encoding=prefs.ENCODING)

    seqnumber = -0

    for line in datasetfile:
        resfile.write(line)
        if line.startswith('>'):
            seqnumber += 1

    for line in unknownsequences:
        # To have a limited initial set for testing.
        # if seqnumber > 2000:
        #     return
        if not line.startswith('>'):
            resfile.write(line)
        else:
            resfile.write('>seq%i\n' % seqnumber)
            seqnumber += 1


def makecombinedselftrainingset(fold, datasetfile):
    if not prefs.limitinitialtrainingsetforunknown:
        fold['trainingset'] = copy.deepcopy(fold['referenceset'])
    domainnames = copy.deepcopy(prefs.domainnames)
    domainnames.append('seq')
    dictionary = unknowndictfromfasta(datasetfile, domainnames)
    selftrainingset = fold['selftrainingset']
    for substrate in domainnames:
        if substrate == 'seq':
            selftrainingset[substrate] = dictionary['seq']
        elif not prefs.includeknownsequencesinselftrainingset:
            selftrainingset[substrate] = []
        else:
            pass


def unknowndictfromfasta(dataset, domainnames):
    data = codecs.open(dataset, encoding=prefs.ENCODING)

    substrates = OrderedDict()
    for substrate in domainnames:
        substrates[str(substrate)] = []

    '''
    Objective: Make a list of numbers for each substrate.
    The number goes into what?
    '''

    for line in data:
        if line.startswith('>'):
            match = re.match(r'([a-z-]+)([0-9]+)', line.strip(' >\n'), re.I)
            if match:
                items = match.groups()
            substrates[items[0]].append(int(items[1]))
    return substrates


def addunknowntotrainingset(bestsequences, index, trainingset):
    substrates = prefs.domainnames
    for sequence in bestsequences:
        addtosubstrate = substrates[int(sequence[0])]
        trainingset[addtosubstrate].append(index[int(sequence[1])])


if __name__ == '__main__':
    makecombinedfastafile(
        prefs.datasetmakefolds,
        'datasets/a_domains_enumerated.fa',
        'datasets/test.fa')
