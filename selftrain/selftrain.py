
import logging
import prefs
import makefolds
import copy
import pprint
import seqlcaller
from multiprocessing import Process
import scorechecker
import subprocess
import sys
# import numpy
import datetime
import time
import outputscripts
from collections import deque
import unknownutils
import metadata
import randomseeds


def replaylastrun(fromiteration):
    # initloggers()
    cleanup()
    # print('replayprocces running')
    foldresults = []

    for foldnumber in range(prefs.numberoffolds):
        foldfilename = '%s%i/selftraingresults.log' % (
            prefs.logdir, foldnumber)
        foldresults.append(outputscripts.parseoutput(foldfilename))

    sequencearray, substratearray = makefolds.makenumbertosequence(
        logging.getLogger('numbertoseqeuncesdict'))

    dict = makefolds.dictfromfasta(prefs.dataset, prefs.domainnames)
    log = logging.getLogger('makefolds')
    randomizeddict = copy.deepcopy(dict)
    folds = makefolds.makenfolds(randomizeddict, log)
    # print(folds[0])
    for foldnumber in range(prefs.numberoffolds):
        # for foldnumber in range(1):
        selftrainingset = folds[foldnumber]['selftrainingset']
        trainingset = folds[foldnumber]['trainingset']

        # print(selftrainingset)
        # print(trainingset)
        # print(folds[foldnumber]['referenceset'])
        # print(folds[foldnumber]['predictset'])

        for iterationnumber in range(fromiteration):

            sequencestoadd = foldresults[
                foldnumber].sequences[iterationnumber].split(' ')

            addsequencesass = foldresults[
                foldnumber].substrates[iterationnumber].lower().split(' ')

            for sequence in range(len(sequencestoadd)):
                sequencetoadd = sequencestoadd[sequence]
                addassubstrate = addsequencesass[sequence]

                # remove from selftrainingset
                removed = False
                for substrate in prefs.domainnames:
                    intsequence = int(sequencetoadd)
                    if intsequence in selftrainingset[substrate]:
                        selftrainingset[substrate].remove(intsequence)
                        removed = True
                if not removed:
                    print('FAULT! %s' % sequencetoadd)
                    sys.exit()
                else:
                    pass
                    # print('added to selftrainingset %s' % sequencetoadd)
                # Add to trainingset
                trainingset[addassubstrate].append(sequencetoadd)
    # folds now contains the same as the folds after the iteration specified.
    # lets get the sequences we need to replay the run:
    sequencestoadd = []
    for foldnumber in range(prefs.numberoffolds):
        foldsequences = []
        for iteration in range(len(foldresults[
                foldnumber].iterationnumber) - fromiteration-1):

            foldsequences.append(foldresults[foldnumber].sequences[
                iteration+fromiteration].split(' '))
        sequencestoadd.append(foldsequences)
    # cleanup()
    replayprocs = []
    for foldnumber in range(prefs.numberoffolds):
        foldsets = folds[foldnumber]
        replayprocs.append(
            Process(
                target=replayprocess,
                args=(
                    foldsets,
                    foldnumber,
                    sequencearray,
                    substratearray,
                    sequencestoadd[foldnumber],
                    fromiteration)
                )
            )

    for proc in replayprocs:
        proc.start()
    for proc in replayprocs:
        proc.join()


def replayprocess(
    fold,
    foldnumber,
    sequencearray,
    substratearray,
    sequences,
    fromiteration
):

    resultlog = logging.getLogger('selftraingresults%i' % foldnumber)
    resultlog.debug(
        '%%==========================================================================')
    resultlog.debug('%%=== Replaying from iteration %i' % fromiteration)
    processlog = logging.getLogger('fold%i' % foldnumber)
    processlog.info('started selflearning')

    # Get directory strings
    foldstring = str(foldnumber) + '/'
    currenttraindir = prefs.sequencedirs + foldstring
    currentmodeldir = prefs.modelsdirs + foldstring
    currentbindir = prefs.bindirs + foldstring
    currentpredictorsdir = prefs.predictorsdirs + foldstring
    classifydir = prefs.classifydirs + foldstring

    iteration = 0
    iterationsleft = 0

    while True:
        processlog.info('iteration %i\t iterations left %i' % (
            iteration, iterationsleft))

        # Make sequencefiles from trainingset
        log = logging.getLogger('makesseqfiles' + str(foldnumber))
        makefolds.makeseqeuncefiles(
            fold['trainingset'],
            foldnumber,
            currenttraindir,
            sequencearray,
            log)

        # Train on sequencefiles
        log = logging.getLogger('train' + str(foldnumber))
        processlog.debug('training')
        seqlcaller.train(currenttraindir, currentmodeldir, log)

        # Make binary model
        log = logging.getLogger('binarymodel' + str(foldnumber))
        processlog.debug('binarymodel')
        seqlcaller.makebinarymodel(
            currentmodeldir,
            currentbindir,
            currentpredictorsdir,
            log)

        # Classify selftrainingset
        selftrainseqs = []
        selftrainindex = []
        selftrainingset = fold['selftrainingset']
        # print(selftrainingset)

        for substrate in selftrainingset:
            for sequence in selftrainingset[substrate]:
                selftrainseqs.append(sequencearray[sequence])
                selftrainindex.append(int(sequence))

        iterationsleft = len(selftrainindex) // prefs.numbertoadd

        if len(selftrainseqs) > 0:
            processlog.debug('Making classifyfiles')
            log = logging.getLogger('makeclassifyfiles' + str(foldnumber))
            makefolds.makeclassifyfile(selftrainseqs, classifydir, log)

            processlog.debug('classifying sequencefiles')
            scores = seqlcaller.classifysequencefile(
                classifydir + 'classify.fasta',
                len(selftrainseqs),
                currentbindir
            )

            processlog.debug('Get best')
            log = logging.getLogger('getbestscores' + str(foldnumber))
            randomseed = prefs.randomseed + str(foldnumber) + str(iteration)
            bestsequences, resultscores = scorechecker.getbestscores(
                scores, randomseed,  log)
            # if len(resultscores) < prefs.numbertoadd:
            #     print('fold%i last iteration' % foldnumber)
            #     lastiteration = True

            # bestsequences contains pairs of [substratenumber, sequncenumber]
            # Substratearray is the substrate of sequence number.
            # ei substratearray[sequencenumber] -> the substratename
            # sequencearray[sequence] contains the actual sequence.
            # the bestseqeunces should be read throught the index.

            # Check correctness:
            # print(selftrainindex)
            # print(sequences[iteration])
            sequencesthisiteration = []
            for sequence in sequences[iteration]:
                sequencenumber = selftrainindex.index(int(sequence))
                sequencesubstrate = substratearray[int(sequence)]
                sequencesubstratenumber = prefs.domainnames.index(
                    sequencesubstrate)
                # print(sequencesubstrate)
                sequencesthisiteration.append([
                    int(sequencesubstratenumber), int(sequencenumber)])
            # print(sequencesthisiteration)
            # print(sequencesthisiteration)
            # print(selftrainindex)
            processlog.debug('check correctness')
            correctness = scorechecker.checkcorrectness(
                sequencesthisiteration,
                selftrainindex,
                resultscores,
                substratearray,
                logging.getLogger('checkcorrect%i' % foldnumber)
            )
            processlog.debug('moving from selftrainingset')
            makefolds.addtotrainingset(
                sequencesthisiteration,
                selftrainindex,
                substratearray,
                fold['trainingset'],
                fold['selftrainingset'],
                True,
                logging.getLogger('addtoselftrain%i' % foldnumber)
            )
            sequencesadded = ''
            for sequence in sequencesthisiteration:
                sequencesadded += str(selftrainindex[int(sequence[1])]) + ' '
            sequencesadded = sequencesadded.strip()
        # end if len(selftrainseqs) > 0
        # run classifier on referenceset
        classifysequences = []
        classifyindex = []
        for substrate in fold['predictset']:
            for sequence in fold['predictset'][substrate]:
                classifysequences.append(sequencearray[sequence])
                classifyindex.append(sequence)

        makefolds.makeclassifyfile(classifysequences, classifydir, log)
        processlog.debug('classifying reference')
        scores = seqlcaller.classifysequencefile(
            classifydir + 'classify.fasta',
            len(classifysequences),
            currentbindir
        )

        classification, correctpercent = scorechecker.checkscores(
            scores,
            classifyindex,
            substratearray,
            fold,
            log)

        if len(selftrainseqs) > 0:
            resultlog.debug('%i\t%f\t%s\t%s\t%f\t%f' % (
                iteration,
                correctpercent,
                correctness,
                sequencesadded,
                min(resultscores),
                max(resultscores)
                )
            )
        else:
            # Exiting. Logging last classification.
            resultlog.debug('%i\t%f' % (
                iteration,
                correctpercent
                )
            )
            processlog.info('fold%i, terminated' % foldnumber)
            sys.exit()
        iteration += 1



def main(randseed=prefs.randomseed):
    """ Main method
Should setup parser and loggers befor making datasets.
Returns?
    """
    # Setup logging:
    # initloggers()
    # Cleanup:
    cleanup()

    # Each method that needs logging should set up their own
    # logger calling setup
    # Make datasets
    sequencearray, substratearray = makefolds.makenumbertosequence(
        logging.getLogger('numbertoseqeuncesdict'))

    dict = makefolds.dictfromfasta(prefs.datasetmakefolds, prefs.domainnames)
    log = logging.getLogger('makefolds')
    randomizeddict = copy.deepcopy(dict)
    folds = makefolds.makenfolds(randomizeddict, log, randseed=randseed)

    # make reference for every fold:
    referencetrainproc = []
    for fold in range(prefs.numberoffolds):
        referenceset = folds[fold]['referenceset']
        predictionset = folds[fold]['predictset']
        referencetrainproc.append(
            Process(
                target=referenceprocess,
                args=(
                    fold,
                    referenceset,
                    sequencearray,
                    substratearray,
                    predictionset)
                )
            )
    for proc in referencetrainproc:
        proc.start()
    for proc in referencetrainproc:
        proc.join()

    cleanup()
    # print(folds[0])
    selftrainproc = []
    for fold in range(prefs.numberoffolds):
        foldsets = folds[fold]
        selftrainproc.append(
            Process(
                target=selftrainprocess,
                args=(
                    foldsets,
                    fold,
                    sequencearray,
                    substratearray)
                )
            )
    for proc in selftrainproc:
        proc.start()
    for proc in selftrainproc:
        proc.join()
    cleanup()


def rununknown(randseed=prefs.randomseed):

    cleanup()
    # print(prefs.domainnames)
    sequencearray, substratearray = makefolds.makenumbertosequence(
        logging.getLogger('numbertoseqeuncesdict'))

    # make on dataset without the other sequences appended on.
    dict = makefolds.dictfromfasta(prefs.datasetmakefolds, prefs.domainnames)
    log = logging.getLogger('makefolds')
    randomizeddict = copy.deepcopy(dict)
    folds = makefolds.makenfolds(randomizeddict, log, randseed=randseed)
    for foldnumber in range(prefs.numberoffolds):
        unknownutils.makecombinedselftrainingset(
            folds[foldnumber], 'datasets/test.fa'
        )
    # print(prefs.domainnames)
    referencetrainproc = []
    for fold in range(prefs.numberoffolds):
        referenceset = folds[fold]['referenceset']
        predictionset = folds[fold]['predictset']
        referencetrainproc.append(
            Process(
                target=referenceprocess,
                args=(
                    fold,
                    referenceset,
                    sequencearray,
                    substratearray,
                    predictionset)
                )
            )
    for proc in referencetrainproc:
        proc.start()
    for proc in referencetrainproc:
        proc.join()

    cleanup()

    selftrainprocs = []
    for foldnumber in range(prefs.numberoffolds):
        fold = folds[foldnumber]
        selftrainprocs.append(
            Process(
                target=unknowprocess,
                args=(
                    foldnumber,
                    fold,
                    sequencearray,
                    substratearray
                    )
                ))
    for proc in selftrainprocs:
        proc.start()
    for proc in selftrainprocs:
        proc.join()


def unknowprocess(
    foldnumber,
    fold,
    sequencearray,
    substratearray
):
    resultlog = logging.getLogger('selftraingresults%i' % foldnumber)
    processlog = logging.getLogger('fold%i' % foldnumber)
    processlog.info('started selftraining')

    foldstring = str(foldnumber) + '/'
    currenttraindir = prefs.sequencedirs + foldstring
    currentmodeldir = prefs.modelsdirs + foldstring
    currentbindir = prefs.bindirs + foldstring
    currentpredictorsdir = prefs.predictorsdirs + foldstring
    classifydir = prefs.classifydirs + foldstring

    iteration = 0
    totaliterationsleft = prefs.unknowniterations
    iterationsleft = 0
    iterationtime = time.time()

    iterationtimes = deque(2*[0], 2)
    noresultscores = False

    while True:
        # pprint.pprint(len(fold['trainingset']['a']))
        # iterationsleft = totaliterationsleft - iteration
        thisiterationtime = time.time()
        timesincelastiteration = thisiterationtime - iterationtime
        iterationtimes.appendleft(timesincelastiteration)
        if iteration > 1:
            averagetime = sum(iterationtimes)/len(iterationtimes)
            dt = datetime.datetime.utcfromtimestamp(averagetime)
            dts = dt.time()
        else:
            dts = 'estimating. Wait for iteration 2'
        processlog.info(
            'iteration %i of %i\t approxtime per iteration: %s' % (
                iteration,
                totaliterationsleft,
                dts
                )
            )
        iterationtime = thisiterationtime
        # Make sequencefiles from trainingset
        log = logging.getLogger('makesseqfiles' + str(foldnumber))
        makefolds.makeseqeuncefiles(
            fold['trainingset'],
            foldnumber,
            currenttraindir,
            sequencearray,
            log)

        # Train on sequencefiles
        log = logging.getLogger('train' + str(foldnumber))
        processlog.debug('training')
        seqlcaller.train(currenttraindir, currentmodeldir, log)

        # Make binary model
        log = logging.getLogger('binarymodel' + str(foldnumber))
        processlog.debug('binarymodel')
        seqlcaller.makebinarymodel(
            currentmodeldir,
            currentbindir,
            currentpredictorsdir,
            log)

        # Classify selftrainingset
        selftrainseqs = []
        selftrainindex = []
        selftrainingset = fold['selftrainingset']

        for substrate in selftrainingset:
            for sequence in selftrainingset[substrate]:
                selftrainseqs.append(sequencearray[sequence])
                selftrainindex.append(int(sequence))

        if iteration == 0:
            iterationsuntilnomore = (len(selftrainindex)//prefs.numbertoadd)+1
            totaliterationsleft = min(
                iterationsuntilnomore, totaliterationsleft)
        if iteration <= prefs.unknowniterations and len(selftrainseqs) > 0:
            processlog.debug('Making classifyfiles')
            log = logging.getLogger('makeclassifyfiles' + str(foldnumber))
            makefolds.makeclassifyfile(selftrainseqs, classifydir, log)

            # processlog.info('classifying sequencefiles')

            processlog.info(
                'classifying %i sequencefiles' % len(selftrainseqs))
            scores = seqlcaller.classifysequencefile(
                classifydir + 'classify.fasta',
                len(selftrainseqs),
                currentbindir
            )
            # pprint.pprint(scores)
            # processlog.info('get best')
            processlog.debug('Get best')
            log = logging.getLogger('getbestscores' + str(foldnumber))
            randomseed = prefs.randomseed + str(foldnumber) + str(iteration)
            bestsequences, resultscores = scorechecker.getbestscores(
                scores, randomseed, log)

            if len(resultscores) > 0:

                # print(bestsequences)
                # print(resultscores)
                # sys.exit()
                # if len(resultscores) < prefs.numbertoadd:
                #     print('fold%i last iteration' % foldnumber)
                #     lastiteration = True
                processlog.debug('check correctness')
                correctness = scorechecker.checkcorrectness(
                    bestsequences,
                    selftrainindex,
                    resultscores,
                    substratearray,
                    logging.getLogger('checkcorrect%i' % foldnumber)
                )
                # bestsequences contains pairs of
                # [substratenumber, sequncenumber]
                # Substratearray is the substrate of sequence number.
                # ei substratearray[sequencenumber] -> the substratename
                # sequencearray[sequence] contains the actual sequence.
                # the bestseqeunces should be read throught the index.

                # Move from selftrainingset to trainingset
                processlog.debug('moving from selftrainingset')
                makefolds.addtotrainingset(
                    bestsequences,
                    selftrainindex,
                    substratearray,
                    fold['trainingset'],
                    fold['selftrainingset'],
                    prefs.addsequenceascorrectsubstrate,
                    logging.getLogger('addtoselftrain%i' % foldnumber)
                )
                sequencesadded = ''
                for sequence in bestsequences:
                    sequencesadded += str(selftrainindex[sequence[1]]) + ' '
                sequencesadded = sequencesadded.strip()
                noresultscores = False
            else:
                noresultscores = True
        # end if len(selftrainseqs) > 0
        else:
            # noresultscores = True
            processlog.info('last iteration, classifying and terminating')

        # run classifier on referenceset
        classifysequences = []
        classifyindex = []
        for substrate in fold['predictset']:
            for sequence in fold['predictset'][substrate]:
                classifysequences.append(sequencearray[sequence])
                classifyindex.append(sequence)

        makefolds.makeclassifyfile(classifysequences, classifydir, log)
        processlog.debug('classifying reference')
        scores = seqlcaller.classifysequencefile(
            classifydir + 'classify.fasta',
            len(classifysequences),
            currentbindir
        )

        classification, correctpercent = scorechecker.checkscores(
            scores,
            classifyindex,
            substratearray,
            fold,
            log)
        processlog.info('\tacc:%f' % (correctpercent))

        if iteration <= prefs.unknowniterations and len(
                resultscores) > 0 and len(selftrainseqs) > 0:
            # processlog.info('fold%i\tacc:%f' % (foldnumber, correctpercent))
            resultlog.debug('%i\t%f\t%s\t%s\t%f\t%f' % (
                iteration,
                correctpercent,
                correctness,
                sequencesadded,
                min(resultscores),
                max(resultscores)
                )
            )
            noresultscores = False
        else:
            # Exiting. Logging last classification.
            resultlog.debug('%i\t%f' % (
                iteration,
                correctpercent
                )
            )
            # processlog.info('fold%i, terminated' % foldnumber)
            if len(selftrainseqs) == 0 or iteration > prefs.unknowniterations:
                processlog.info('fold%i, terminated' % foldnumber)
                sys.exit()
            elif not noresultscores:
                noresultscores = True
            else:
                processlog.info('fold%i, terminated' % foldnumber)
                sys.exit()
        iteration += 1


def selftrainprocess(
    fold,
    foldnumber,
    sequencearray,
    substratearray
):
    resultlog = logging.getLogger('selftraingresults%i' % foldnumber)
    processlog = logging.getLogger('fold%i' % foldnumber)
    processlog.info('started selftraining')

    # Get directory strings
    foldstring = str(foldnumber) + '/'
    currenttraindir = prefs.sequencedirs + foldstring
    currentmodeldir = prefs.modelsdirs + foldstring
    currentbindir = prefs.bindirs + foldstring
    currentpredictorsdir = prefs.predictorsdirs + foldstring
    classifydir = prefs.classifydirs + foldstring

    iteration = 0
    totaliterationsleft = 0
    iterationsleft = 0
    iterationtime = time.time()
    iterationtimes = deque(5*[0], 5)
    noresultscores = False
    while True:
        iterationsleft = totaliterationsleft - iteration
        thisiterationtime = time.time()
        timesincelastiteration = thisiterationtime - iterationtime
        iterationtimes.appendleft(timesincelastiteration)
        if iteration > 4:
            averagetime = sum(iterationtimes)/len(iterationtimes)
            # print(iterationtimes)
            dt = datetime.datetime.utcfromtimestamp(
                averagetime*(iterationsleft+1))
            dts = dt.time()
        else:
            dts = 'estimating. Waiting for iteration 5'

        processlog.info(
            'iteration %i of %i\t approxtimeleft: %s' % (
                iteration,
                totaliterationsleft,
                dts
                )
            )
        iterationtime = thisiterationtime

        # Make sequencefiles from trainingset
        log = logging.getLogger('makesseqfiles' + str(foldnumber))
        makefolds.makeseqeuncefiles(
            fold['trainingset'],
            foldnumber,
            currenttraindir,
            sequencearray,
            log)

        # Train on sequencefiles
        log = logging.getLogger('train' + str(foldnumber))
        processlog.debug('training')
        seqlcaller.train(currenttraindir, currentmodeldir, log)

        # Make binary model
        log = logging.getLogger('binarymodel' + str(foldnumber))
        processlog.debug('binarymodel')
        seqlcaller.makebinarymodel(
            currentmodeldir,
            currentbindir,
            currentpredictorsdir,
            log)

        # Classify selftrainingset
        selftrainseqs = []
        selftrainindex = []
        selftrainingset = fold['selftrainingset']

        for substrate in selftrainingset:
            for sequence in selftrainingset[substrate]:
                selftrainseqs.append(sequencearray[sequence])
                selftrainindex.append(sequence)
        if iteration == 0:
            totaliterationsleft = len(selftrainindex) // prefs.numbertoadd

        if len(selftrainseqs) > 0:
            processlog.debug('Making classifyfiles')
            log = logging.getLogger('makeclassifyfiles' + str(foldnumber))
            makefolds.makeclassifyfile(selftrainseqs, classifydir, log)

            processlog.debug('classifying sequencefiles')
            scores = seqlcaller.classifysequencefile(
                classifydir + 'classify.fasta',
                len(selftrainseqs),
                currentbindir
            )

            processlog.debug('Get best')
            log = logging.getLogger('getbestscores' + str(foldnumber))
            randomseed = prefs.randomseed + str(foldnumber) + str(iteration)
            bestsequences, resultscores = scorechecker.getbestscores(
                scores, randomseed, log)
            # if len(resultscores) < prefs.numbertoadd:
            #     print('fold%i last iteration' % foldnumber)
            #     lastiteration = True

            # bestsequences contains pairs of [substratenumber, sequncenumber]
            # Substratearray is the substrate of sequence number.
            # ei substratearray[sequencenumber] -> the substratename
            # sequencearray[sequence] contains the actual sequence.
            # the bestseqeunces should be read throught the index.
            if len(resultscores) > 0:
                # Check correctness:
                processlog.debug('check correctness')
                correctness = scorechecker.checkcorrectness(
                    bestsequences,
                    selftrainindex,
                    resultscores,
                    substratearray,
                    logging.getLogger('checkcorrect%i' % foldnumber)
                )
                # Append to sequencefiles
                # Does not work yet, as we rewrite the sequencefiles
                # every iteration. That works but this would be nicer.
                # This code, for some reason makes seql training slower
                # by a lot. This makes no sense at all.
                # processlog.debug('appending to sequencefiles')
                # makefolds.appendtoseqfile(
                #     bestsequences,
                #     selftrainindex,
                #     sequencearray,
                #     substratearray,
                #     currenttraindir,
                #     logging.getLogger('appendtoseqfile%i' % foldnumber)
                # )

                # Move from selftrainingset to trainingset
                processlog.debug('moving from selftrainingset')
                makefolds.addtotrainingset(
                    bestsequences,
                    selftrainindex,
                    substratearray,
                    fold['trainingset'],
                    fold['selftrainingset'],
                    prefs.addsequenceascorrectsubstrate,
                    logging.getLogger('addtoselftrain%i' % foldnumber)
                )
                sequencesadded = ''
                for sequence in bestsequences:
                    sequencesadded += str(selftrainindex[sequence[1]]) + ' '
                sequencesadded = sequencesadded.strip()
            else:
                noresultscores = True
        # end if len(selftrainseqs) > 0

        # run classifier on referenceset
        classifysequences = []
        classifyindex = []
        for substrate in fold['predictset']:
            for sequence in fold['predictset'][substrate]:
                classifysequences.append(sequencearray[sequence])
                classifyindex.append(sequence)

        makefolds.makeclassifyfile(classifysequences, classifydir, log)
        processlog.debug('classifying reference')
        scores = seqlcaller.classifysequencefile(
            classifydir + 'classify.fasta',
            len(classifysequences),
            currentbindir
        )

        classification, correctpercent = scorechecker.checkscores(
            scores,
            classifyindex,
            substratearray,
            fold,
            log)
        processlog.info('\tacc:%f' % (correctpercent))

        if len(selftrainseqs) > 0 and len(resultscores) > 0:
            # processlog.info('fold%i\tacc:%f' % (foldnumber, correctpercent))
            resultlog.debug('%i\t%f\t%s\t%s\t%f\t%f' % (
                iteration,
                correctpercent,
                correctness,
                sequencesadded,
                min(resultscores),
                max(resultscores)
                )
            )
            noresultscores = False
        else:
            # Exiting. Logging last classification.
            resultlog.debug('%i\t%f' % (
                iteration,
                correctpercent
                )
            )
            # processlog.info('fold%i, terminated' % foldnumber)
            if len(selftrainseqs) == 0:
                processlog.info('fold%i, terminated' % foldnumber)
                sys.exit()
            elif not noresultscores:
                noresultscores = True
            else:
                processlog.info('fold%i, terminated' % foldnumber)
                sys.exit()
        iteration += 1


def referenceprocess(
    fold,
    referenceset,
    sequencearray,
    substratearray,
    predictionset
):
    resultlog = logging.getLogger('selftraingresults%i' % fold)

    foldstring = str(fold) + '/'
    currentrefdir = prefs.referencedirs + foldstring

    log = logging.getLogger('makesrefseqfiles' + str(fold))
    makefolds.makeseqeuncefiles(
        referenceset,
        fold,
        currentrefdir,
        sequencearray,
        log)

    currentmodeldir = prefs.modelsdirs + foldstring
    log = logging.getLogger('reftrain' + str(fold))

    seqlcaller.train(currentrefdir, currentmodeldir, log)

    currentbindir = prefs.bindirs + foldstring
    log = logging.getLogger('binarymodels' + str(fold))
    currentpredictorsdir = prefs.predictorsdirs + foldstring

    seqlcaller.makebinarymodel(
        currentmodeldir,
        currentbindir,
        currentpredictorsdir,
        log)

    classifysequences = []
    classifyindex = []
    for substrate in predictionset:
        for sequence in predictionset[substrate]:
            classifysequences.append(sequencearray[sequence])
            classifyindex.append(sequence)

    classifydir = prefs.classifydirs + foldstring
    makefolds.makeclassifyfile(classifysequences, classifydir, log)

    scores = seqlcaller.classifysequencefile(
        classifydir + 'classify.fasta',
        len(classifysequences),
        currentbindir
    )

    # scorechecker must log the results..
    classification, correctpercent = scorechecker.checkscores(
        scores,
        classifyindex,
        substratearray,
        fold,
        log)
    processlog = logging.getLogger('fold%i' % fold)
    processlog.info('referencescore: %f' % correctpercent)
    resultlog.debug(metadata.getmetadatastring(predictionset, correctpercent))


def initloggers():
    """ Method for setting up logger for a method.
    This log writes to a file in logdir
    """
    # Loggging init taken from
    # https://docs.python.org/3/howto/logging-cookbook.html#logging-cookbook
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
        datefmt='%m-%d %H:%M',
        filename=prefs.logdir + prefs.mainlog,
        filemode='a+')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(prefs.verbositylevel)
    # set a format which is simpler for console use
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)
    logging.info('Logging initialised')

    for number in range(prefs.numberoffolds):
        setupresultsloggers(
            'selftraingresults%i' % number,
            '../out/%i/selftraingresults.log' % number
        )


def setupresultsloggers(logname, filename, level=logging.DEBUG):
    l = logging.getLogger(logname)
    formatter = logging.Formatter('%(message)s')
    fileHandler = logging.FileHandler(filename, mode='a+')
    fileHandler.setFormatter(formatter)
    l.setLevel(level)
    l.addHandler(fileHandler)


def cleanup():
    log = logging.getLogger('cleanup')
    log.info('Removing previous experiments')
    rmcall = [
        'rm',
        '-r',
        prefs.sequencedirs,
        prefs.modelsdirs,
        prefs.bindirs,
        prefs.predictorsdirs,
        prefs.referencedirs,
        prefs.classifydirs
    ]
    subprocess.call(rmcall)

    mkdircall = ['mkdir', '-p']
    for number in range(prefs.numberoffolds):
        mkdircall += [prefs.sequencedirs + '/' + str(number)]
        mkdircall += [prefs.predictorsdirs + '/' + str(number)]
        mkdircall += [prefs.modelsdirs + '/' + str(number)]
        mkdircall += [prefs.bindirs + '/' + str(number)]
        mkdircall += [prefs.referencedirs + '/' + str(number)]
        mkdircall += [prefs.classifydirs + '/' + str(number)]

    log.info(
        'Making new directories for a ' + str(prefs.numberoffolds) + '-fold')
    subprocess.call(mkdircall)


def runmainandunknown(times=10):
    if times > 100:
        print('not enoght ramdomseeds. Exiting')
        sys.exit()
    for iteration in range(times):
        main(randseed=randomseeds.randomseeds[iteration])
        outputscripts.parseoutputs(
            randseed=randomseeds.randomseeds[iteration],
            writestdres=True)
        rununknown(randseed=randomseeds.randomseeds[iteration])
        outputscripts.parseoutputs(
            randseed=randomseeds.randomseeds[iteration],
            writestdres=True)


if __name__ == '__main__':
    initloggers()
    main()
    # outputscripts.parseoutputs()
    # rununknown()
    outputscripts.parseoutputs()
    # replaylastrun(0)
    # runmainandunknown(6)

