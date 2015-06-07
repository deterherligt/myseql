
import prefs
import codecs
import re
import random
# import pprint
from collections import OrderedDict
# import sys


def makemap(domainnames):
    substrates = dictfromfasta(prefs.dataset, domainnames)
    resmap = {}
    for substratenumber in range(len(substrates)):
        for sequencenumber in substrates[domainnames[substratenumber]]:
            resmap[sequencenumber] = substratenumber
    return resmap


# Append to sequencefiles
def appendtoseqfile(
    bestsequences,
    index,
    sequencearray,
    substratearray,
    trainingdir,
    log
):
    log.debug('appending to sequncefiles')
    substrates = prefs.domainnames

    for substrate in substrates:
        seqfilename = '%s%s.seq' % (trainingdir, substrate)
        seqfile = codecs.open(seqfilename, 'a+', encoding=prefs.ENCODING)

        for sequence in bestsequences:
            sequencesubstrate = substratearray[sequence[0]]
            seqline = '+1 ' if sequencesubstrate == substrate else '-1 '
            seqline += sequencearray[index[sequence[1]]] + '\n'
            seqfile.write(seqline)
        seqfile.flush()
        seqfile.close()


# Move from selftrainingset to trainingset
def addtotrainingset(
    bestsequences,
    index,
    substratearray,
    trainingset,
    selftrainingset,
    addascorrectsubstrate,
    log
):
    substrates = prefs.domainnames
    for sequence in bestsequences:
        removefromsubstrate = substratearray[index[int(sequence[1])]]
        addtosubstrate = substrates[int(sequence[0])]
        if addascorrectsubstrate:
            addtosubstrate = removefromsubstrate
        # remove
        # print(sequence)
        selftrainingset[removefromsubstrate].remove(index[int(sequence[1])])
        # add
        trainingset[addtosubstrate].append(index[int(sequence[1])])


def dictfromfasta(dataset, domainnames):
    data = codecs.open(dataset, encoding=prefs.ENCODING)
    # substratesfile = codecs.open(
    #     prefs.substratenames,
    #     encoding=prefs.ENCODING)
    substrates = OrderedDict()
    for substrate in domainnames:
        # substrate = line.strip()
        substrates[str(substrate)] = []
        # substrates.append(currentsubstrate)
    # print(substrates)

    # substrates = prefs.domainnames
    # substrates.append('fig')

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

    # pprint.pprint(substrates)
    # substrates now contains every number of each substrate
    return substrates


def getchunks(list, folds):
    return [list[i:: folds] for i in range(folds)]


def makenfolds(substratedict, logger, randseed=prefs.randomseed):
    '''
Should return a dict for each fold. N i given in prefs file.
Random seed is given as argument
    '''
    randseed += 'makefolds'
    logger.info('making folds with random seed: ' + randseed)
    makefoldsrandom = random.Random()
    makefoldsrandom.seed(randseed)
    # print(makefoldsrandom.getstate())
    # pp = pprint.PrettyPrinter(width=1000)

    for substrate in substratedict:
        makefoldsrandom.shuffle(substratedict[substrate])

    folds = OrderedDict()
    for substrate in substratedict:
        chunks = getchunks(substratedict[substrate], prefs.numberoffolds)
        folds[substrate] = chunks

    resfolds = []
    for foldnumber in range(prefs.numberoffolds):

        traindict = OrderedDict()
        selftraindict = OrderedDict()
        predictdict = OrderedDict()
        referencedict = OrderedDict()

        # Calculate witch parts goes to this fold:
        predictpart = foldnumber % prefs.numberoffolds
        selftrainparts = []
        trainparts = []
        for y in range((prefs.numberoffolds - 1) // 2):
            numbertoaddtoselftrain = (foldnumber + y + 1) % prefs.numberoffolds
            numbertoaddtotrain = (
                foldnumber + y + 1 +
                (prefs.numberoffolds - 1) // 2) % prefs.numberoffolds
            selftrainparts.append(numbertoaddtoselftrain)
            trainparts.append(numbertoaddtotrain)
        # If even number of folds. just in case
        if prefs.numberoffolds % 2 == 0:
            trainparts.append(
                (foldnumber +
                 (prefs.numberoffolds - 1)) % prefs.numberoffolds)

        # Add parts for this fold
        for substrate in folds:
            seqs = folds[substrate]
            # substrate is substrate (key) and seqs is sequences.

            # Predict part is always only one ot the n parts. We can therefore
            # just set it
            predictdict[substrate] = seqs[predictpart]

            # Init substrate arrays
            selftraindict[substrate] = []
            traindict[substrate] = []
            referencedict[substrate] = []

            for selftrainpart in selftrainparts:
                if prefs.limitinitialtrainingset:
                    for number in seqs[selftrainpart]:
                        # This is incredebly ugly, but the linter says its the
                        # prettiest way...
                        if len(
                                traindict[substrate]
                        ) < prefs.initaltrainingsetlimit:
                            traindict[substrate].append(number)
                        else:
                            selftraindict[substrate].append(number)
                else:
                    traindict[substrate] += seqs[selftrainpart]
                # Always add to reference
                referencedict[substrate] += seqs[selftrainpart]

            for trainpart in trainparts:
                selftraindict[substrate] += seqs[trainpart]
                # Always add to reference
                referencedict[substrate] += seqs[trainpart]

        resfolds.append(dict(
            trainingset=traindict,
            selftrainingset=selftraindict,
            predictset=predictdict,
            referenceset=referencedict))
    # pp.pprint(resfolds)
    logger.debug('Makefolds done')
    return resfolds


def makeclassifyfile(sequences, directory, log):
    classifyfilename = directory + 'classify.fasta'
    classifyfile = codecs.open(classifyfilename, 'w', encoding=prefs.ENCODING)
    for sequence in sequences:
        classifyfile.write(sequence + '\n')
    classifyfile.close()


def makenumbertosequence(logger):
    logger.debug('making dict for conversion between numbers and seqences')
    data = codecs.open(prefs.dataset, encoding=prefs.ENCODING)
    sequences = []
    substrates = []
    currentnumber = 0
    currentsubstrate = ''
    for line in data:
        if line.startswith('>'):
            # TODO: RE should proberly be compiled
            # as the same RE is used multiply times
            # re.I means ignore case.
            match = re.match(r'([a-z-]+)([0-9]+)', line.strip(' >\n'), re.I)
            if match:
                items = match.groups()
            currentnumber = int(items[1])
            currentsubstrate = str(items[0])
        else:
            sequences.insert(currentnumber, str(line.strip()))
            substrates.insert(currentnumber, currentsubstrate)
            # sequences[currentnumber] = [currentsubstrate, line.strip()]
            # substrates[items[0]].append(int(items[1]))
    # pprint.pprint(sequences)
    return sequences, substrates


def makeseqeuncefiles(fold, foldnumber, directory, sequencearray, logger):
    logger.debug('making sequencefiles for fold' + str(foldnumber))

    sequencedir = directory

    for substrate in fold:

        sequencefilename = sequencedir + substrate + '.seq'
        sequencefile = codecs.open(
            sequencefilename, 'w', encoding=prefs.ENCODING)

        # print(sequencefilename)
        for substrates in fold:
            for sequence in fold[substrates]:
                seqstring = '+1 ' if substrates == substrate else '-1 '
                sequencefile.write(
                    seqstring + sequencearray[int(sequence)] + '\n')
        sequencefile.flush()
        sequencefile.close()
        # currentsequencenumber = 0


if __name__ == '__main__':
    # makenfolds()
    pass
