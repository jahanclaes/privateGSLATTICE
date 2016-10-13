import numpy as np
import stats
import matplotlib.pyplot as plt

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
print "Beta =",beta

observableList = [0 for i in range(1000000)]
for i in range(numThreads):
    for j in range(NumWalkers):
        observableFile=open("MyRunobservables."+str(i)+"_"+str(j))
        index = 0
        for line in observableFile:
            observableList[index]+=float(line)/(numThreads*NumWalkers)
            index+=1
        print index
        observableFile.close()

burnIn=min(20,index/5)
observableList = observableList[burnIn:index]

print stats.Stats(np.array(observableList))
