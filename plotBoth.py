import numpy as np
import matplotlib.pyplot as plt
import stats

# Import all the parameters from input.dat 
file = open("input.dat")
inputDict={}
for line in file:
    tokens=line.split("=")
    try:
        inputDict[tokens[0]]=float(tokens[1])
    except:
        inputDict[tokens[0]]=tokens[1]
file.close()
markovSteps = int(inputDict["markovSteps"])
TimeStepsTaken = int(inputDict["TimeStepsTaken"])
StepSize=inputDict["StepSize"]
beta = 2*TimeStepsTaken*StepSize
NumWalkers=int( inputDict["NumWalkers"])
numThreads=1
burnIn = 20
print "Beta =",beta

observableList = [0 for i in range(1000000)]

for i in range(numThreads):
    for j in range(NumWalkers):
        observableFile=open("MyRunobservables."+str(i)+"_"+str(j))
        index = 0
        for line in observableFile:
            observableList[index]+=float(line)/(numThreads*NumWalkers)
            index+=1
        observableFile.close()

observableList = observableList[0:index]

print stats.Stats(np.array(observableList[int(.2*index):]))
