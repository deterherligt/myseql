import codecs
import prefs
import re
import pprint
from collections import namedtuple
from operator import attrgetter
import statistics
import makefolds
import metadata
import copy
from scipy import stats


def printtest():
    foldresults = []

    for foldnumber in range(prefs.numberoffolds):
        foldfilename = '%s%i/selftraingresults.log' % (
            prefs.logdir, foldnumber)
        foldresults.append(parseoutput(foldfilename))
    # print(foldresults)


def parseoutput(filename):
    outputfile = codecs.open(filename, 'r', encoding=prefs.ENCODING)

    lines = []
    for line in outputfile:
        lines.append(line.strip())
    lengthofoutput = len(lines)
    lastline = lines[-1]
    numberofiterations = int(lastline.split('\t')[0])
    firstlineindex = lengthofoutput - numberofiterations - 1

    iterationnumber = []
    scores = []
    substrates = []
    substratenumbers = []
    sequences = []
    minscores = []
    maxscores = []

    for line in lines[firstlineindex:]:
        rematch = re.match(
            r'([0-9]+\t)'
            r'([\d.]+\t*)'
            r'([\D-]+\t)?'
            r'([0-9 ]+\t)?'
            r'([\d.]+\t)?'
            r'([\d.]+)?',
            line.strip(),
            re.I)
        if rematch:
            items = rematch.groups()

        iterationnumber.append(int(items[0].strip()))
        scores.append(items[1].strip())
        if items[2] is not None:
            substrates.append(items[2].strip())
            substratenumbers.append(parsesubstrates(items[2].strip()))
            sequences.append(items[3].strip())
            minscores.append(items[4].strip())
            maxscores.append(items[5].strip())
    try:
        predictionsetline = [int(line.strip()) for line in lines[
            firstlineindex-2].split('[')[1].split(']')[0].split(',')]
    except:
        predictionsetline = ''

    referencescore = float(lines[firstlineindex-1].split(' ')[-1])

    # print(predictionsetline)

    Results = namedtuple('Results', [
        'iterationnumber',
        'scores',
        'substrates',
        'substratenumbers',
        'sequences',
        'minscores',
        'maxscores',
        'predictionset',
        'referencescore'])

    results = Results(
        iterationnumber,
        scores,
        substrates,
        substratenumbers,
        sequences,
        minscores,
        maxscores,
        predictionsetline,
        referencescore)
    return results


def parseoutputs(randseed=prefs.randomseed, writestdres=False):
    domainnames = copy.deepcopy(prefs.domainnames)
    domainnames.append('seq')
    # print(domainnames)
    foldmap = makefolds.makemap(domainnames)

    foldresults = []
    for foldnumber in range(prefs.numberoffolds):
        foldfilename = '%s%i/selftraingresults.log' % (
            prefs.logdir, foldnumber)
        foldresults.append(parseoutput(foldfilename))
    # sorted(seq, key=lambda x: x.age)
    maxit = max(max(foldresults, key=attrgetter('iterationnumber'))[0])
    results = []
    seqlminscores = []
    seqlmaxscores = []
    referencescores = []
    for foldnumber in range(prefs.numberoffolds):
        results.append(foldresults[foldnumber].scores)
        seqlminscores.append(foldresults[foldnumber].minscores)
        seqlmaxscores.append(foldresults[foldnumber].maxscores)
        referencescores.append(foldresults[foldnumber].referencescore)
    # print(seqlminscores)
    # print(seqlmaxscores)
    avgres = []
    seqlminres = []
    seqlmaxres = []

    typesofsequencesadded = []
    addedassubstrates = []
    progressaddedas = []
    progresscorrect = []
    progressfalse = []
    for substrate in domainnames:
        typesofsequencesadded.append([])
        addedassubstrates.append([])
        progressaddedas.append([])
        progresscorrect.append([])
        # progressfalse.append([])

    for iteration in range(maxit+1):
        iterationres = []
        iterationseqlmin = []
        iterationseqlmax = []
        correct = 0
        total = 0

        iterationrescorrect = []
        iterationresaddedas = []

        for substrate in domainnames:
            iterationrescorrect.append(0)
            iterationresaddedas.append(0)

        for foldnumber in range(prefs.numberoffolds):
            if len(foldresults[foldnumber].sequences) <= iteration:
                continue
            else:
                for sequencenumber in foldresults[
                        foldnumber].sequences[iteration].split(' '):
                    iterationrescorrect[foldmap[int(sequencenumber)]] += 1
                for substratenumber in foldresults[
                        foldnumber].substratenumbers[iteration]:
                    iterationresaddedas[substratenumber] += 1
        wrongforsubstrate = []
        for substrate in range(len(prefs.domainnames)):
            wrongforsubstrate.append(0)
        for substratenumber in range(len(iterationrescorrect)):
            progresscorrect[
                substratenumber].append(iterationrescorrect[
                    substratenumber])
            progressaddedas[substratenumber].append(
                iterationresaddedas[substratenumber])
            # print(iterationrescorrect)
            # if not iterationrescorrect[
            #         substratenumber] == iterationresaddedas[substratenumber]:
            #     # print('%i==%i' % (
            #     #     iterationrescorrect[substratenumber],
            #     #     iterationresaddedas[substratenumber]))
            #     # wrongforsubstrate += 1
            #     # print(substrate)
            #     # print(sum(iterationrescorrect))
            #     # print(sum(iterationresaddedas))
            #     wrongforsubstrate[substratenumber] = iterationresaddedas[substratenumber] - iterationrescorrect[substratenumber]
                # print('wrong')
            # print(substrate)
            # progressfalse[substrate].append(wrong)
            # print(substratenumber)
            # print(wrong)
            # I cannot use wrong as this is for all substrates.
            #  i need the number of wrongs for only this substrate
            # progressfalse[substratenumber].append(wrong)
            if iteration == 0:
                typesofsequencesadded[
                    substratenumber].append(iterationrescorrect[
                        substratenumber])
                addedassubstrates[substratenumber].append(
                    iterationresaddedas[substratenumber])
            else:
                typesofsequencesadded[
                    substratenumber].append(typesofsequencesadded[
                        substratenumber][-1] + iterationrescorrect[
                            substratenumber])
                addedassubstrates[substratenumber].append(
                    addedassubstrates[
                        substratenumber][-1] + iterationresaddedas[
                            substratenumber])
        # print(wrongforsubstrate)
        # wrong = 0
        # wrongforsubstrate = []
        # for substrate in range(len(prefs.domainnames)):
        #     wrongforsubstrate.append(0)
        # for foldnumber in range(prefs.numberoffolds):
        #     if len(foldresults[foldnumber].substrates) > iteration:
        #         substrates = foldresults[foldnumber].substrates[iteration]
        #         substratelist = substrates.split()
        #         for substrate in substratelist:
        #             # if substrate.lower() == 'seq':
        #             #     continue
        #             if substrate.isupper():
        #                 # wrong += 1
        #                 wrongforsubstrate[
        #                     prefs.domainnames.index(substrate.lower())] += 1
        # # wrong = sum(wrongforsubstrate)
        wrongforsubstrate = []
        for substrate in range(len(prefs.domainnames)):
            wrongforsubstrate.append(0)

        for foldnumber in range(prefs.numberoffolds):
            while len(results[foldnumber]) < maxit+1:
                results[foldnumber].append(results[foldnumber][-1])
            while len(seqlmaxscores[foldnumber]) < maxit+1:
                seqlmaxscores[foldnumber].append(seqlmaxscores[foldnumber][-1])
            while len(seqlminscores[foldnumber]) < maxit+1:
                seqlminscores[foldnumber].append(seqlminscores[foldnumber][-1])
            iterationres.append(float(results[foldnumber][iteration]))
            iterationseqlmin.append(
                float(seqlminscores[foldnumber][iteration]))
            iterationseqlmax.append(
                float(seqlmaxscores[foldnumber][iteration]))
            if len(foldresults[foldnumber].substrates) > iteration:
                substrates = foldresults[foldnumber].substrates[iteration]
                substratelist = substrates.split()
                for substrate in substratelist:
                    total += 1
                    if substrate.islower():
                        correct += 1
                    else:  # substrate.isupper():
                        # wrong += 1
                        wrongforsubstrate[
                            prefs.domainnames.index(substrate.lower())] += 1
            else:
                continue
        # TODO: Maybe i should also log median for each iteration?
        stdev = statistics.stdev(iterationres)
        iterationres.append(statistics.mean(iterationres))
        iterationres.append(stdev)
        iterationres.append(total)
        iterationres.append(total-correct)
        iterationseqlmin.append(statistics.mean(iterationseqlmin))
        iterationseqlmax.append(statistics.mean(iterationseqlmax))
        seqlminres.append(iterationseqlmin)
        seqlmaxres.append(iterationseqlmax)
        avgres.append(iterationres)
        progressfalse.append(wrongforsubstrate)
    # pprint.pprint(len(avgres))
    # print(len(progressaddedas))
    # pprint.pprint((progressfalse[-4]))
    # pprint.pprint(sum(progressfalse[-4]))

    # Subtract unknown from wrong, and add number at the end of iteration
    # print(avgres[-1])
    # print(avgres[-1][:5])
    # print(statistics.stdev(avgres[-1][:5]))
    stdev = statistics.stdev(avgres[-1][:5])
    for iteration in range(len(avgres)):
        avgres[iteration][-1] -= progresscorrect[-1][iteration]
        avgres[iteration].append(progresscorrect[-1][iteration])
    makeavgdat(
        avgres,
        stdev,
        referencescores,
        seqlminres,
        seqlmaxres,
        'average.dat')
    makeheatmapdat(
        typesofsequencesadded,
        stdev,
        referencescores,
        'heatmaptypes.dat')
    makeheatmapdat(
        addedassubstrates,
        stdev,
        referencescores,
        'heatmapaddedas.dat')
    # # print(len(progresscorrect))
    # # print(progresscorrect)
    # # print(avgres)
    makeprogressheatmap(
        progressaddedas,
        avgres,
        progressfalse,
        stdev,
        referencescores,
        'heatmapprogressaddedas.dat')
    makeprogressheatmap(
        progresscorrect,
        avgres,
        progressfalse,
        stdev,
        referencescores,
        'heatmapprogresstypes.dat')
    if writestdres:
        appendtoresfile(avgres, stdev, 'stdevres.dat', randseed)


def appendtoresfile(avgres, stdev, filename, randseed):
    output = codecs.open(
        'texdatout/%s' % filename, 'a', encoding=prefs.ENCODING)
    lastline = avgres[-1]
    line = ''
    for res in lastline[:6]:
        line += str(res) + '\t'
    line += '%f\t%s\n' % (stdev, randseed)
    output.write(line)


def makeprogressheatmap(
    heatmap,
    avgres,
    progressfalse,
    stdev,
    referencescores,
    filename
):
    output = codecs.open(
        'texdatout/%s' % filename, 'a', encoding=prefs.ENCODING)
    # output.write('\n')
    titleline = metadata.getoutputmetadatastring(
        statistics.mean(referencescores),
        statistics.stdev(referencescores),
        stdev)
    output.write(titleline + '\n')
    totaladdedperiteration = []
    for it in range(len(heatmap[0])):
        totaladdedperiteration.append(0)
        for substrate in prefs.domainnames:
            value = heatmap[prefs.domainnames.index(substrate)][it]
            totaladdedperiteration[it] += value
    maxwrong = 0
    for iteration in progressfalse:
        maxiteration = max(iteration)
        if maxiteration > maxwrong:
            maxwrong = maxiteration
    for substrate in prefs.domainnames:
        line = ' {:<8} &'.format(substrate)
        for iteration in range(len(heatmap[0])):
            value = heatmap[prefs.domainnames.index(substrate)][iteration]
            if not totaladdedperiteration[iteration] == 0:
                line += '  {:<8.4} &'.format(
                    value/totaladdedperiteration[iteration]*100.0)
            # else:
            #     line += ''
        line = '{:s}{:>9}'.format(line.strip('&'), '\\\\[-2ex]\n')
        # line = line.strip('& ') + '\\\\[-2ex]\n'
        output.write(line)

        # Now for the wrongly added substrates
        line = '          &'
        # print(len(progressfalse[0]))
        for iteration in range(len(progressfalse)-1):
            # print(prefs.domainnames.index(substrate))
            value = progressfalse[iteration][
                prefs.domainnames.index(substrate)] / maxwrong
            line += ' -{:<8.4} &'.format(value*100.0)
        line = line.rstrip('&') + '\\\\[-3ex]\n'
        output.write(line)

    avgresarray = []
    for it in range(len(avgres)):
        itres = avgres[it][prefs.numberoffolds]
        avgresarray.append(itres)
    minavg = min(element for element in avgresarray)
    maxavg = max(element for element in avgresarray)
    resarraynormalized = []
    for element in avgresarray:
        resarraynormalized.append((element-minavg)/(maxavg-minavg))
    # print(resarraynormalized)
    line = 'average   &'
    for itres in resarraynormalized:
        line += '  {:<8.4} &'.format(itres*100.0)
    output.write(line.strip('&') + '\\\\\n')
    output.close()


def makeavgdat(average, stdev, referencescores, seqlmin, seqlmax, filename):
    output = codecs.open(
        'texdatout/%s' % filename, 'a', encoding=prefs.ENCODING)
    # TODO. Need referencescore in output
    titleline = metadata.getoutputmetadatastring(
        statistics.mean(referencescores),
        statistics.stdev(referencescores),
        stdev) + '\nite  '
    for fold in range(prefs.numberoffolds):
        titleline += '{1}{0:<3d}{2}{0:<2d}{3}{0:<3d}'.format(
            fold, 'fold', 'min', 'max')
    titleline += ' avgscor stddev avgsmin avgsmax added wrong unknowns correctpct\n'
    output.write(titleline)
    for iteration in range(len(average)):
        iterationres = average[iteration]
        iterationseqlmin = seqlmin[iteration]
        iterationseqlmax = seqlmax[iteration]

        line = '{:<4}'.format(iteration)
        for fold in range(prefs.numberoffolds):
            line += '{:^7.4}{:^5.2}{:^6.2}'.format(
                iterationres[fold],
                iterationseqlmin[fold],
                iterationseqlmax[fold])
        line += ' {:^7.4} '.format(iterationres[prefs.numberoffolds])
        line += ' {:^7.4} '.format(iterationres[prefs.numberoffolds + 1])
        line += ' {:<7.3}'.format(iterationseqlmin[prefs.numberoffolds])
        line += ' {:<7.3}'.format(iterationseqlmax[prefs.numberoffolds])

        for element in range(len(iterationres)-prefs.numberoffolds-2):
            line += ' {:^6}'.format(
                iterationres[element+prefs.numberoffolds+2])
        addedandwrong = iterationres[-3:][:2]
        if addedandwrong[0] > 0:
            wrongpct = addedandwrong[1]/addedandwrong[0]*100
        else:
            wrongpct = 100.0
        # print(iterationres[-3:][:2])
        line += ' {:<7.4}'.format(wrongpct)
        output.write(line.strip() + '\n')
    output.close()


def makeheatmapdat(heatmap, stdev, referencescores, filename):
    output = codecs.open(
        'texdatout/%s' % filename, 'a', encoding=prefs.ENCODING)
    titleline = metadata.getoutputmetadatastring(
        statistics.mean(referencescores),
        statistics.stdev(referencescores),
        stdev)
    output.write(titleline + '\n')
    for substrate in prefs.domainnames:
        line = '{:<9}& '.format(substrate)
        # line = ''
        for iteration in range(len(heatmap[0])):
            value = heatmap[prefs.domainnames.index(substrate)][iteration]
            if heatmap[prefs.domainnames.index(substrate)][-1] > 0:
                percent = (value / heatmap[
                    prefs.domainnames.index(substrate)][-1])*100
            else:
                percent = 0.0
            line += '{:<7.4} & '.format(percent)
        line = line.strip('& ') + '\\\\\n'
        output.write(line)
    output.close()


def calcstddev(filename):
    file = codecs.open(
        'texdatout/%s' % filename, 'r', encoding=prefs.ENCODING)
    outputfile = codecs.open(
        'texdatout/stedev.tex', 'w', encoding=prefs.ENCODING)
    foldresknown = []
    resavgknown = []

    foldresunknown = []
    resavgunknown = []

    outputline = 'seed                  '
    for foldnumber in range(prefs.numberoffolds):
        outputline += 'known%i ' % foldnumber
    for foldnumber in range(prefs.numberoffolds):
        outputline += 'unkno%i ' % foldnumber

    outputline += ' knownavg unknownavg knownstdev unknowstdev\n'

    outputfile.write(outputline)
    linenumber = 0
    knownstdev = 0.0
    for line in file:
        line = line.split('\t')
        if linenumber % 2 == 0:
            outputline = line[-1].strip() + '  '
            # even which means its with knowns
            for i in range(prefs.numberoffolds):
                foldresknown.append(float(line[i]))
                outputline += '{:<7.5}'.format(float(line[i]))
            resavgknown.append(float(line[i+1]))
            knownstdev = float(line[-2])
        else:
            # unknows
            for i in range(prefs.numberoffolds):
                foldresunknown.append(float(line[i]))
                outputline += '{:<7.5}'.format(float(line[i]))
            resavgunknown.append(float(line[i+1]))
            outputline += ' {:<8.5} {:<8.5}'.format(
                resavgknown[-1], resavgunknown[-1])
            outputline += '   {:<8.5}   {:<8.5}\n'.format(
                knownstdev, float(line[-2]))
            outputfile.write(outputline)
        linenumber += 1
    # print(foldresknown)
    # print(resavgknown)

    # print(foldresunknown)
    # print(resavgunknown)

    # baseline = [67.2, 67.4, 71.5, 77.6, 86.0, 89.1, 59.5, 81.9, 105.5]
    # follow_up = [62.4, 64.6, 70.4, 62.6, 80.1, 73.2, 58.2, 71.0, 101.0]

    # paired_sample = stats.ttest_rel(baseline, follow_up)

    res = 'known:   Average of every fold: {:5.4}\n'.format(
        statistics.mean(foldresknown))
    res += 'known:   Standard deviation of every fold: {:5.4}\n'.format(
        statistics.stdev(foldresknown))
    res += 'known:   Average of avereges: {:5.4}\n'.format(
        statistics.mean(resavgknown))
    res += 'known:   Standard deviation of averages: {:5.4}\n'.format(
        statistics.stdev(resavgknown))
    res += 'unknown: Average of every fold: {:5.4}\n'.format(
        statistics.mean(foldresunknown))
    res += 'unknown: Standard deviation of every fold: {:5.4}\n'.format(
        statistics.stdev(foldresunknown))
    res += 'unknown: Average of avereges: {:5.4}\n'.format(
        statistics.mean(resavgunknown))
    res += 'unknown: Standard deviation of averages: {:5.4}\n'.format(
        statistics.stdev(resavgunknown))
    pairedres = stats.ttest_rel(resavgknown, resavgunknown)
    res += 'Average: The t-statistic is {:5.4} with p-value {:5.4}\n'.format(
        pairedres[0], pairedres[1])
    pairedres = stats.ttest_rel(foldresknown, foldresunknown)
    res += 'Every fold: The t-statistic is {:5.4} with p-value {:5.4}'.format(
        pairedres[0], pairedres[1])
    # print(pairedres)
    print(len(resavgknown))
    print(len(resavgunknown))
    print(res)


def parsesubstrates(substratestring):
    substrates = str(substratestring).split()
    res = []
    domains = prefs.domainnames
    for substrate in substrates:
        res.append(domains.index(substrate.lower()))
    return res

if __name__ == '__main__':
    # parseoutputs()
    # printtest()
    calcstddev('stdevres.dat')
