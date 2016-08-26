import numpy
import numpy.linalg

def write_rPoints(a,b,Lx,Ly,rPoints):
    rPoints_file=open("rPoints.txt","w")
    rPoints_file.write(str(a[0])+" "+str(a[1])+'\n')
    rPoints_file.write(str(b[0])+" "+str(b[1])+'\n')
    rPoints_file.write(str(Lx)+"\n")
    rPoints_file.write(str(Ly)+"\n")
    for i in range(0,len(rPoints)):
        rPoints_file.write(str(rPoints[i][0])+" "+str(rPoints[i][1])+"\n")
    rPoints_file.close()


def write_eigs(s):
    t=-1.0
    mySize=len(s.rPoints)
    H=numpy.zeros((mySize,mySize))
    ll=len(s.rPoints)
    for i in range(0,len(s.rPoints)):
        for j in range(0,len(s.neighbors[i])):
            myNeighbor=s.neighbors[i][j]
            H[i,myNeighbor]=t
    (e,v)=numpy.linalg.eigh(H[0:ll,0:ll])
    eigs_file=open("eigs.txt",'w')
    A=range(0,s.numUpElectrons)
    for i in A:
        for j in v[:,i]:
            eigs_file.write(str(j)+" ")
        eigs_file.write("\n")


class SquareHubbardModel:
    Lx=4
    Ly=4
    rPoints_dict=dict()
    rPoints=[]
    doping=0.875
    def Init(self):
        self.numSites=self.Lx*self.Ly
        self.numUpElectrons=int(round(self.doping*(self.Lx*self.Ly)/2.)) # actually num up electrons in spin layer
    def build_rPoints(self):
        countMe=0
        for ix in range(0,self.Lx):
            for iy in range(0,self.Ly):
                self.rPoints_dict[(ix,iy)]=countMe
                self.rPoints.append((ix,iy))
                countMe=countMe+1
        return 
    def buildHopping(self):
        self.hops=[]
        for n in range(0,len(self.neighbors)):
            for j in self.neighbors[n]:
                self.hops.append((n,j))
    def buildNeighbors(self,xBoundary="periodic",yBoundary="periodic"):
        self.neighbors=[]
        for idx,(ix,iy) in enumerate(self.rPoints):
            self.neighbors.append([])
            if (ix!=self.Lx or xBoundary=="periodic"):
                myNeighbor=((ix+1) % self.Lx,iy)
                self.neighbors[-1].append(self.rPoints_dict[myNeighbor])   
            if (ix!=0 or xBoundary=="periodic"):
                myNeighbor=(( (ix-1)+self.Lx) % self.Lx,iy)
                self.neighbors[-1].append(self.rPoints_dict[myNeighbor])   
            if (iy!=self.Ly or yBoundary=="periodic"):
                myNeighbor=(ix,(iy+1) % self.Ly)
                self.neighbors[-1].append(self.rPoints_dict[myNeighbor])   
            if (iy!=0 or yBoundary=="periodic"):
                myNeighbor=(ix,(iy-1+self.Ly) % self.Ly)
                self.neighbors[-1].append(self.rPoints_dict[myNeighbor])   

def write_wavefunction(outFile,wf_type,wf_dict):
    outFile.write("WaveFunction:\n")
    if wf_type=="RVB":
        outFile.write("\t name=RVB\n")
    if wf_type=="PEPS":
        outFile.write("\t name=PEPS\n")
        outFile.write("\t Dx="+wf_dict["Dx"]+"\n")
        outFile.write("\t Dy="+wf_dict["Dy"]+"\n")
        outFile.write("\t chi="+wf_dict["chi"]+"\n")
        outFile.write("\t L="+wf_dict["L"]+"\n")
        outFile.write("\t W="+wf_dict["W"]+"\n")
        outFile.write("\t tol="+wf_dict["tol"]+"\n")
    outFile.write("END\n")

def write_Hamiltonian(outFile,h_type,ham_dict):
    outFile.write("Hamiltonian:\n")
    if h_type=="Hubbard":
        outFile.write("\t name=Hubbard\n")
        outFile.write("\t bondFile="+ham_dict["bondFile"]+'\n')
        outFile.write("\t U="+ham_dict["U"]+'\n')
        outFile.write("\t t="+ham_dict["t"]+'\n')
    outFile.write("END"+'\n')



def write_input(fileName,paramsDict):
    outFile=open(fileName,'w')
    outFile.write("RunType="+paramsDict["RunType"]+'\n')
    if paramsDict["RunType"]=="OPT":
        outFile.write("OptType="+paramsDict["OptType"]+'\n')
    outFile.write("StepSize="+paramsDict["StepSize"]+'\n')
    outFile.write("NumSteps="+paramsDict["NumSteps"]+'\n')
    outFile.write("EquilSweeps="+paramsDict["EquilSweeps"]+'\n')
    outFile.write("SampleSweeps="+paramsDict["SampleSweeps"]+'\n')
    outFile.write("NumWalkers="+paramsDict["NumWalkers"]+'\n')
    outFile.write("doping="+paramsDict["doping"]+'\n')
    outFile.write("MoveType="+paramsDict["MoveType"]+'\n')
    write_wavefunction(outFile,"PEPS",p["wavefunction"])
    
    write_Hamiltonian(outFile,"Hubbard",paramsDict["hamiltonian"])
    outFile.write("ReadParams="+paramsDict["ReadParams"]+'\n')
    if paramsDict["ReadParams"]=="True":
        outFile.write("ParamFileName="+paramsDict["ParamFileName"]+'\n')
p=dict()
p["RunType"]="OPT"
p["OptType"]="GRADIENT"
p["StepSize"]=str(0.01)
p["NumSteps"]=str(1000)
p["EquilSweeps"]=str(100)
p["SampleSweeps"]=str(1000)
p["NumWalkers"]=str(1)
p["doping"]=str(0.875)
p["MoveType"]="HOP"
p["ReadParams"]="False"

p["hamiltonian"]=dict()
p["hamiltonian"]["U"]="8.0"
p["hamiltonian"]["t"]="1.0"
p["hamiltonian"]["bondFile"]="bond.dat"

p["wavefunction"]=dict()
p["wavefunction"]["Dx"]=str(2)
p["wavefunction"]["Dy"]=str(2)
p["wavefunction"]["chi"]=str(20)
p["wavefunction"]["L"]=str(4)
p["wavefunction"]["W"]=str(4)
p["wavefunction"]["tol"]=str(1e-4)

write_input("input_hubbard.dat",p)


s=SquareHubbardModel()
s.Init()
s.build_rPoints()
write_rPoints([1,0],[0,1],s.Lx,s.Ly,s.rPoints)
s.buildNeighbors(xBoundary="open")
s.buildHopping()
write_eigs(s)

a=open("bond.dat",'w')
for i in s.hops:
    a.write(str(i[0])+" "+str(i[1])+"\n")
a.close()
