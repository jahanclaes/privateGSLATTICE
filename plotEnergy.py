import matplotlib.pyplot as plt
import matplotlib

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
StepSize = float(inputDict["StepSize"])

energyFile = open("Energy.dat")
energyList = []
varianceList = []
acceptanceRatioList = []
for line in energyFile:
    pieces = line.split(" ")
 #   energyList.append(float(pieces[1]))
    energyList.append(float(pieces[3]))
    if -float(pieces[3])<5.:
        varianceList.append(-float(pieces[4]))
    else:
        varianceList.append(5.)
    acceptanceRatioList.append(float(pieces[5]))
index= len(energyList)
print index

numCycles = 15
rng1 = [i for i in range(index)]
rng2 = [i*TimeStepsTaken-1 for i in range(1,markovSteps+1) if i*TimeStepsTaken-1 <index]
rng3 = [i*TimeStepsTaken for i in range(markovSteps) if i*TimeStepsTaken<index]

average2 = sum([energyList[i] for i in rng2[len(rng2)/5:]])/len(rng2[len(rng2)/5:])
average3 = sum([energyList[i] for i in rng3[len(rng3)/5:]])/len(rng3[len(rng3)/5:])
print average2, average3

lineweight=2
markersize=7
fontsize=20
matplotlib.rcParams.update({'text.usetex':'True','font.family':'Computer Modern Roman','font.size': fontsize})


plt.figure(num=None, figsize=(8, 6))
#for i in range((len(rng3))/5+5,len(rng3)/5+9):
#    if energyList[rng3[i+1]]>energyList[rng3[i]]:
#        plt.annotate('', xytext=(StepSize*rng3[i],energyList[rng3[i]]),xy=(StepSize*rng3[i+1],energyList[rng3[i+1]]),arrowprops=dict(color='g',ls='dotted',lw=lineweight,arrowstyle='-|>',connectionstyle='angle3,angleA=80,angleB=-10',shrinkA=10,shrinkB=5))
#    else:
#        plt.annotate('', xytext=(StepSize*rng3[i],energyList[rng3[i]]),xy=(StepSize*rng3[i+1],energyList[rng3[i+1]]),arrowprops=dict(color='g',ls='dotted',lw=lineweight,arrowstyle='-|>',connectionstyle='angle3,angleA=-10,angleB=100',shrinkA=10,shrinkB=5))
plt.subplot(2,1,1)
plt.plot([StepSize*i for i in rng1], [energyList[i] for i in rng1],lw=lineweight,label="Energy")
plt.plot([StepSize*i for i in rng3], [energyList[i] for i in rng3], 'go',lw=lineweight,ms=markersize,label="Start of Markov step")
plt.plot([StepSize*i for i in rng2], [energyList[i] for i in rng2], 'rs',lw=lineweight,ms=markersize,label="End of Markov step")
plt.plot([StepSize*i for i in rng2], [average2 for i in rng2],'r--',lw=2,label="Average Energy" )
plt.subplot(2,1,2)
plt.plot([StepSize*i for i in rng1], [varianceList[i] for i in rng1],lw=lineweight,label="Variance")
plt.plot([StepSize*i for i in rng1], [acceptanceRatioList[i] for i in rng1],lw=lineweight,label="AcceptanceRatio")
plt.ylim([-.5,3])
plt.legend()
plt.show()
