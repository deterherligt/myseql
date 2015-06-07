import datetime
import prefs


def getmetadatastring(predictionset, referencescore):
    time = datetime.datetime.now().time()
    predictionsize = []
    for substrate in predictionset:
        predictionsize.append(len(predictionset[substrate]))
    logstring = '\
%%==========================================================================\n\
%%=== New run time: %s\n\
%%=== Dataset: %s \n\
%%=== Randomseed: \'%s\'\n\
%%=== Limittrainingset: %r\n\
%%=== Limit: %i\n\
%%=== number added per iteration: %i\n\
%%=== random selection: %s\n\
%%=== add as correct substrate: %s\n\
%%=== seql threshold: %f\n\
%%=== predictionset %s\n\
%%=== Unknown: limit intial trainingset: %s\n\
%%=== Unknown: include known sequences in selftraingset: %s\n\
%%=== Unknown: Max iterations: %s\n\
%%=== Referencescore: %f' % (
        time,
        prefs.dataset,
        prefs.randomseed,
        prefs.limitinitialtrainingset,
        prefs.initaltrainingsetlimit,
        prefs.numbertoadd,
        prefs.randomselection,
        prefs.addsequenceascorrectsubstrate,
        prefs.seqlsminscore,
        predictionsize,
        prefs.limitinitialtrainingsetforunknown,
        prefs.includeknownsequencesinselftrainingset,
        prefs.unknowniterations,
        referencescore)
    return logstring


def getoutputmetadatastring(referencescore, refstdev, stdev):
    time = datetime.datetime.now().time()
    logstring = '\
%%==========================================================================\n\
%%=== New run time: %s\n\
%%=== Dataset: %s \n\
%%=== Randomseed: \'%s\'\n\
%%=== Limittrainingset: %r\n\
%%=== Limit: %i\n\
%%=== number added per iteration: %i\n\
%%=== random selection: %s\n\
%%=== add as correct substrate: %s\n\
%%=== seql threshold: %f\n\
%%=== Unknown: limit intial trainingset: %s\n\
%%=== Unknown: include known sequences in selftraingset: %s\n\
%%=== Unknown: Max iterations: %s\n\
%%=== Standart deviation last iteration: %f\n\
%%=== Referencescore standart deviation: %f\n\
%%=== Referencescore: %f' % (
        time,
        prefs.dataset,
        prefs.randomseed,
        prefs.limitinitialtrainingset,
        prefs.initaltrainingsetlimit,
        prefs.numbertoadd,
        prefs.randomselection,
        prefs.addsequenceascorrectsubstrate,
        prefs.seqlsminscore,
        prefs.limitinitialtrainingsetforunknown,
        prefs.includeknownsequencesinselftrainingset,
        prefs.unknowniterations,
        stdev,
        refstdev,
        referencescore)
    return logstring
