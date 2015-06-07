
import prefs
import numpy
import random


# Check correctness:
def checkcorrectness(
    bestsequences,
    index,
    resultscores,
    substratearray,
    log
):
    resstring = ''
    log.debug('checking correctness')
    substrates = prefs.domainnames
    for sequence in bestsequences:
        # print(sequence)
        addedassubstrate = substrates[sequence[0]]
        if substrates[int(sequence[0])] == substratearray[
                index[int(sequence[1])]]:
            # print('correct')
            resstring += addedassubstrate + ' '
        else:
            # print('incorrect')
            resstring += addedassubstrate.upper() + ' '
    return resstring.strip()


def checkscores(scores, classifyindex, substratearray, fold, log):

    substrates = prefs.domainnames
    numberofsubstrates = len(substrates)
    bestclassification = []

    bestscore = -1000
    bestcounter = 0
    counter = 0
    for score in numpy.nditer(scores.T.copy(order='C')):
        if score > bestscore:
            bestscore = score
            bestcounter = counter
        counter += 1
        if counter % numberofsubstrates is 0:
            bestclassification.append(substrates[bestcounter])
            bestcounter = counter = 0
            bestscore = -1000
    correct = 0
    counter = 0
    # print(substratearray)
    for sequencenumber in range(len(bestclassification)):
        if bestclassification[sequencenumber] \
                == substratearray[classifyindex[sequencenumber]]:
            correct += 1
        counter += 1
    correctpercent = 0.0 + float(correct)/float(counter) * 100
    return bestclassification, correctpercent


def getbestscores(scores, randomseed, log):
    results = []
    resultscores = []

    if prefs.randomselection:
        log.debug('finding %i random sequences to add' % prefs.numbertoadd)
        rand = random.Random()
        # randseed = prefs.randomseed + 'getbestscores'
        rand.seed(randomseed)
        log.debug('seeding random function with seed: "%s"' % randomseed)
        numbertoadd = prefs.numbertoadd if scores.shape[
            1] >= prefs.numbertoadd else scores.shape[1]
        sequencestoadd = rand.sample(
            range(scores.shape[1]), numbertoadd)
        for seqnumber in sequencestoadd:
            seqscore = numpy.max(scores[:, seqnumber])
            substratenumber = numpy.argmax(scores[:, seqnumber])
            results.append([substratenumber, seqnumber])
            resultscores.append(seqscore)
    else:
        for n in range(prefs.numbertoadd):
            log.debug('finding ' + str(n + 1) + 'th best score')
            substratenumber, sequencenumber = numpy.unravel_index(
                scores.argmax(),
                scores.shape
            )
            score = scores[substratenumber, sequencenumber]
            if score <= prefs.seqlsminscore:
                log.debug(
                    'best score below %f only found %i results.' % (
                        prefs.seqlsminscore,
                        n)
                )
                return results, resultscores
            log.debug(
                'best score found with value ' +
                str(score) +
                ' in position (' +
                str(substratenumber) +
                ',' +
                str(sequencenumber) +
                ')'
            )

            # print(score)
            row = numpy.zeros(30)
            row.fill(-1000)
            scores[:, sequencenumber] = row
            if score < prefs.seqlmaxscore:
                results.append([substratenumber, sequencenumber])
                resultscores.append(score)

    return results, resultscores
