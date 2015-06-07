import logging

randomseed = '8b4Dxf5IDkV9A7zSOawyeNJ5GR9npeXZbcYofSv355HmqAXvDAD12wsxCpuyb5la'

numberoffolds = 5

limitinitialtrainingset = True
initaltrainingsetlimit = 5

limitinitialtrainingsetforunknown = True
includeknownsequencesinselftrainingset = True

numbertoadd = 25

addsequenceascorrectsubstrate = False
randomselection = False

unknowniterations = 20

# For score based termination. 1.1 seems like the magic number
seqlsminscore = -1000
seqlmaxscore = 100000

# datasets
datasetmakefolds = 'datasets/dataset-B.fasta'

dataset = 'datasets/test.fa'

datasetwithunknows = 'dataset/test.fa'


# LOG setup
logdir = '../out/'
mainlog = 'main.log'
verbositylevel = logging.INFO

ENCODING = 'utf-8'

# Directories
substratenames = '../domainA.names'
sequencedirs = '../training/'
modelsdirs = '../models/'
bindirs = '../bin/'
predictorsdirs = '../predictors/'
referencedirs = '../reference/'
classifydirs = '../classifytemp/'
selftraindirs = '../selftraining/'


domainnames = [
    'a',
    'aad',
    'beta-ala',
    'bht',
    'c',
    'd',
    'dab',
    'dhb',
    'dhpg',
    'dht',
    'e',
    'f',
    'g',
    'horn',
    'hpg',
    'hyv-d',
    'i',
    'k',
    'l',
    'n',
    'orn',
    'p',
    'pip',
    'q',
    'r',
    's',
    't',
    'v',
    'w',
    'y']


# if __name__ == '__main__':
#     init()
