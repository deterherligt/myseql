import codecs
import subprocess
import prefs
import numpy


def train(inputdir, models, log):
    substrates = codecs.open(prefs.substratenames, encoding=prefs.ENCODING)
    verbosity = 0
    wildcards = 0
    traversalStrategy = 0
    regularizer = 0.1
    '''
    seql_learn help:
    [-o objective_function][-m minsup] [-l minpat]
    [-L maxpat] [-g gap] = wildcards [-r traversal_strategy ] [-T #round]
    [-n token_type] [-c convergence_threshold] [-C regularizer_value]
    [-a l1_vs_l2_regularizer_weight] [-v verbosity]
    train model
    '''
    for substrate in substrates:
        seqlCall = ['../seql_learn']
        seqlCall += ['-v', str(verbosity)]
        seqlCall += ['-g', str(wildcards)]
        seqlCall += ['-r', str(traversalStrategy)]
        seqlCall += ['-C', str(regularizer)]
        # seqlCall += ['-a', str(0.5)]
        seqlCall += [inputdir + substrate.strip() + '.seq']
        seqlCall += [models + substrate.strip() + '.model']
        out = subprocess.Popen(seqlCall, stdout=subprocess.PIPE)
        log.debug(out.communicate()[0].decode(prefs.ENCODING))


def makebinarymodel(models, bin, predictors, log):
    substrates = codecs.open(prefs.substratenames, encoding=prefs.ENCODING)
    '''
    seql_mkmodel help:
    [-i model_file] [-o binary_model_file] [-O predictors_file]
    '''
    for substrate in substrates:
        mkModelCall = ['../seql_mkmodel']
        mkModelCall += ['-i', models + (substrate.strip()) + '.model']
        mkModelCall += ['-o', bin + (substrate.strip()) + '.model.bin']
        mkModelCall += ['-O', predictors + substrate.strip() +
                        '.model.predictors']
        out = subprocess.Popen(mkModelCall, stdout=subprocess.PIPE)
        log.debug(out.communicate()[0].decode(prefs.ENCODING))


def classifysequence(sequences, substrate, binD):
    '''
    [-n token_type: 0 word tokens, 1 char tokens]
    [-t classif_threshold]
    [-v verbose] 0 will print scores, 1 will print posibilities
    test_file binary_model_file
    '''
    # 0 will print scores, 1 will also print proberbilities
    verbose = 0
    results = ""
    classifyCall = ['../seql_classify']
    classifyCall += ['-v', str(verbose)]
    classifyCall += ['-t', str(-1)]
    classifyCall += [sequences, binD + substrate.strip() + '.model.bin']
    p = subprocess.Popen(classifyCall, stdout=subprocess.PIPE)
    output = p.stdout.read()
    results += output.decode(prefs.ENCODING)
    return results


def classifysequencefile(classifyfile, numberofseqs, bindir):

    substrates = prefs.domainnames

    numberofsubstrates = len(substrates)
    scores = numpy.zeros((numberofsubstrates, numberofseqs))

    substrateNumber = 0
    for substrate in substrates:
        score = classifysequence(classifyfile, substrate, bindir)
        sequenceNumber = 0
        for line in score.split('\n'):
            if len(line) > 0:
                scores[substrateNumber][sequenceNumber] = line
            sequenceNumber += 1
        substrateNumber += 1

    # This is done because seql classifies avery sequence,
    # and it cant read fasta format.
    # Instead we only give it sequences from the test file.
    # print(scores.shape)
    return scores
