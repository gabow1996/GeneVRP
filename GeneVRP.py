import random, math, sys
from copy import deepcopy
from tqdm import *  
import numpy as np
import pandas as pd

DEBUG = False

geneNum = 100 # 種族數量
generationNum = 2000  # 迭代次數
CENTER = 0  # 配送中心
HUGE = 9999999
VARY = 0.05  # 變異機率
n = 90  # 客戶點數量
k = 8   # 車輛數量
Q = 264 # 額定载重量, t
carC=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
carx=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21]
Tkmax=480
M=999999
h=[]
h.append(0)
for i in range(90):
    h.append(20)


def getdata(filename):
    data=pd.read_csv(filename)
    dArray=data.values[:,1:]
     
    return dArray


def getalldata():
    # a_ik=getdata('a_ik.csv') #車輛能到達之站點
    a_k=getdata('a_k.csv')   #車輛之燃油效率
    c_ij=getdata('c_ij.csv') #配送路徑
    d_i=getdata('d_i.csv')   #站點需求量
    e_i=getdata('e_i.csv')   #各站點最早到達時間
    l_i=getdata('l_i.csv')   #各站點最晚到達時間
    Q_k=getdata('Q_k.csv')   #車容量
    t_ij=getdata('t_ij.csv') #行駛路徑之運輸時間

    a_k=a_k.flatten()
    d_i=d_i.flatten()
    e_i=e_i.flatten()
    l_i=l_i.flatten()
    Q_k=Q_k.flatten()
    
    c_ij=c_ij/1000
    t_ij=t_ij/60
    
    return a_k ,c_ij ,d_i ,e_i ,l_i ,Q_k ,t_ij
a_k ,distB ,t ,eh ,lh ,q ,t_ij =getalldata()

class Gene:
    def __init__(self, name='Gene', data=None):
        self.name = name
        self.length = n + 1
        self.usedcar=[]
        self.subpathT=[]
        self.subpathDist=[]
        self.subpathPos=[]
        self.pathtime=[]
        self.twoTimeIndex=[]
        self.carSubpath=[]
        self.TimeIneed=[]
        self.subPTNeed=[]
        self.genefit=0
        self.pathtimeadd20=[]
        self.costlist=[]
        self.possible=[]
        if data is None:
            self.data = self._getGene(self.length)
        else:
            self.data = data
        self.fit = self.getFit()
        self.chooseProb = 0  
        
        
    def _generate(self, length):
        data = [i for i in range(1, length)]
        random.shuffle(data)
        data.insert(0, CENTER)
        data.append(CENTER)
        return data

    def _insertZeros(self, data):
        sum = 0
        newData = []
        usedcar=[]
        global distB
        Q=270
        
        for index, pos in enumerate(data):
            sum += t[pos]
            if sum > Q:
                newData.append(CENTER)
                sum = t[pos]
            newData.append(pos)
        
        return newData

    def _getGene(self, length):
        data = self._generate(length)
        data = self._insertZeros(data)
        
        return data
    
    def getFit(self):
        fit = distCost = timeCost = overloadCost = overTimeCost=pathcost = 0
        culposT = culdist = x =0
        sumtime =sumtimeadd20 =0
        usedcar=[]
        subpathT=[]
        pathDist=[]
        pathPos=[]
        subpathPos=[]
        pathtime=[]
        pathtimeadd20=[]
        global distB
        for i,pos in  enumerate(self.data):
            if i==0:
                pathPos.append(pos)
                continue
            else:
                culposT+=t[pos]
                culdist+=distB[self.data[i-1]][self.data[i]]
                pathPos.append(pos)
                sumtime+=(t_ij[self.data[i-1]][self.data[i]])
                sumtimeadd20+=(t_ij[self.data[i-1]][self.data[i]])+20
                if(self.data[i]==0):
                    subpathT.append(culposT)
                    culposT=t[pos]
                    pathDist.append(culdist)
                    culdist=0
                    subpathPos.append(pathPos)
                    pathPos=[]
                    pathPos.append(0)
                    pathtime.append(sumtime)
                    sumtime=0
                    pathtimeadd20.append(sumtimeadd20)
                    sumtimeadd20=0
        self.subpathT=subpathT
        self.subpathDist=pathDist
        self.subpathPos=subpathPos
        self.pathtime=pathtime 
        self.pathtimeadd20=pathtimeadd20
        
        timeSpent = 0
        #時間成本
        for i, pos in enumerate(self.data):
            if i == 0:
                continue
            elif pos == CENTER:
                timeSpent = 0
            timeSpent += (t_ij[self.data[i-1]][self.data[i]])
            if timeSpent > lh[pos]:
                timeCost += (timeSpent - lh[pos])
            if timeSpent > 540:
                timeCost+=(timeSpent*HUGE)
            timeSpent += h[pos]
        
        choosetime=[]
        sumtwoT=[]
        twoTimeI=[]
        choosetime=np.array(pathtimeadd20)
        choosetime.sort
        choosetimeBig=choosetime[choosetime>270]
        choosetimeLess=choosetime[choosetime<=270]
        culchTime=0
        pathIndex=[]
        pathTimeindex=[]
        need=[]
        Tneed=[]
        maxT=[]
        subPTindex=[]
        if(len(choosetimeLess))==1:
            twoTimeI.append(choosetimeLess[-1])
            find=pathtimeadd20.index(choosetimeLess[-1])
            pathIndex.append(find)
            pathNeed=subpathT[find]
            need.append(pathNeed)
            pathTimeindex.append(pathIndex)
            maxT.append(pathNeed)
        else:   
            for i in range(len(choosetimeLess)):
                culchTime+=choosetimeLess[i]
                if culchTime<540:
                    sumtwoT.append(choosetimeLess[i])
                    find=pathtimeadd20.index(choosetimeLess[i])
                    pathIndex.append(find)
                    pathNeed=subpathT[find]
                    need.append(pathNeed)
                    if i==len(choosetimeLess)-1:
                        twoTimeI.append(sumtwoT)
                        if len(pathIndex)!=1:
                            pathTimeindex.append(pathIndex)
                            Tneed.append(need)
                            maxT.append(max(need))
                        else:
                            subPTindex.append(pathIndex)
                            pathTimeindex.append(subPTindex)
                            subPTindex=[]
                            Tneed.append(need)
                            maxT.append(max(need))
                elif culchTime>540:
                    twoTimeI.append(sumtwoT)
                    sumtwoT=[]
                    culchTime=choosetimeLess[i]
                    sumtwoT.append(choosetimeLess[i])
                    pathTimeindex.append(pathIndex)
                    Tneed.append(need)
                    maxT.append(max(need))
                    need=[]
                    pathIndex=[]
                    find=pathtimeadd20.index(choosetimeLess[i])
                    pathIndex.append(find)
                    pathNeed=subpathT[find]
                    need.append(pathNeed)
                    if i==len(choosetimeLess)-1:
                        twoTimeI.append(choosetimeLess[i])
                        subPTindex.append(find)
                        pathTimeindex.append(subPTindex)
                        subPTindex=[]
                        Tneed.append(pathNeed)
                        maxT.append(pathNeed)
        for i in range(len(choosetimeBig)):
            twoTimeI.append(choosetimeBig[i])
            find=pathtimeadd20.index(choosetimeBig[i])
            subPTindex.append(find)
            pathTimeindex.append(subPTindex)
            pathNeed=subpathT[find]
            Tneed.append(pathNeed)
            maxT.append(pathNeed)
            subPTindex=[]
        self.twoTimeIndex=twoTimeI
        self.carSubpath=pathTimeindex
        self.TimeIneed=Tneed
        self.subPTNeed=maxT
        
        fitcar=[]
        for i in range(len(maxT)):
            near=(np.abs(q -maxT[i] )).argmin()
            if q[near]<maxT[i]:
                near+=1
                if carC[near]==carC[-1]:
                    fitcar.append(carC[near])
                else:
                    while carC[near] in fitcar:
                            near+=1
                            if carC[near]==carC[-1]:
                                fitcar.append(carC[near])
                                break
                    fitcar.append(carC[near])
                    
            else:
                while carC[near] in fitcar:
                            near+=1
                            if carC[near]==carC[-1]:
                                fitcar.append(carC[near])
                                break
                fitcar.append(carC[near])
        self.usedcar=fitcar
        
        pathcost=0

        for i in range(len(pathTimeindex)):
            carFuel=a_k[fitcar[i]]            
            for j in range(len(pathTimeindex[i])):
                subDist=pathDist[pathTimeindex[i][j]]
                pathcost+=subDist*carFuel*19.2

        load = 0
        distAfterCharge = 0
        Q_i=0
        for i, pos in enumerate(self.data):
            if i == 0:
                continue
            if pos > n:
                distAfterCharge = 0
            elif pos == CENTER:
                load = 0
                distAfterCharge = 0
            else:
                load += t[pos]
                if Q_i<2:
                    Q=416
                    overloadCost += (HUGE * (load > Q))
                    Q_i+=1
                else:
                    Q=264
                    overloadCost += (HUGE * (load > Q))
        
        carcost=len(fitcar)*966.667

        sumTotaltime=sum(pathtimeadd20)
        fit = pathcost + timeCost + overloadCost + carcost + sumTotaltime
        
        costlist=[]
        costlist.append(fit)
        costlist.append(pathcost)
        costlist.append(sumTotaltime)
        costlist.append(carcost)
        costlist.append(timeCost)
        costlist.append(overloadCost)
        self.genefit=costlist
        
        return 1/fit

    def updateChooseProb(self, sumFit):
        self.chooseProb = self.fit / sumFit

    def moveRandSubPathLeft(self):
        path = random.randrange(k-1)
        index = self.data.index(CENTER, path+1)
        locToInsert = 0
        self.data.insert(locToInsert, self.data.pop(index))
        index += 1
        locToInsert += 1
        while self.data[index] != CENTER:
            self.data.insert(locToInsert, self.data.pop(index))
            index += 1
            locToInsert += 1


def getSumFit(genes):
    sum = 0
    for gene in genes:
        sum += gene.fit
    return sum


def getRandomGenes(size):
    genes = []
    for i in range(size):
        genes.append(Gene("Gene "+str(i)))
    return genes


# 計算適應度和
def getSumFit(genes):
    sumFit = 0
    for gene in genes:
        sumFit += gene.fit
    return sumFit


# 更新選擇機率
def updateChooseProb(genes):
    sumFit = getSumFit(genes)
    for gene in genes:
        gene.updateChooseProb(sumFit)

def updateUsedcar(genes):
    for gene in genes:
        gene.chooseCar()

# 計算累積概率
def getSumProb(genes):
    sum = 0
    for gene in genes:
        sum += gene.chooseProb
    return sum


# 選擇複製，選擇前 1/3
def choose(genes):
    num = int(geneNum/6) * 2
    key = lambda gene: gene.chooseProb
    genes.sort(reverse=True, key=key)
    return genes[0:num]


# 交叉一對
def crossPair(gene1, gene2, crossedGenes):
    gene1.moveRandSubPathLeft()
    gene2.moveRandSubPathLeft()
    newGene1 = []
    newGene2 = []
    centers = 0
    firstPos1 = 1
    for pos in gene1.data:
        firstPos1 += 1
        centers += (pos == CENTER)
        newGene1.append(pos)
        if centers >= 2:
            break
    centers = 0
    firstPos2 = 1
    for pos in gene2.data:
        firstPos2 += 1
        centers += (pos == CENTER)
        newGene2.append(pos)
        if centers >= 2:
            break
    
    for pos in gene2.data:
        if pos not in newGene1:
            newGene1.append(pos)
    for pos in gene1.data:
        if pos not in newGene2:
            newGene2.append(pos)

    newGene1.append(CENTER)
    newGene2.append(CENTER)
    
    key = lambda gene: gene.fit
    possible = []
 
    for i in range(7):
        newGene = newGene1.copy()
        newGene.insert(firstPos1, CENTER)
        possible.append(newGene)
        firstPos1 += 1
        


    P0index=[]
    cultime=0
    sumneed=0
    for i in range(len(possible)):
        x=0
        for j in possible[i]:
            if j==0:
                P0index.append(x)
                x+=1
                k=P0index[len(P0index)-2]+1
                lenk=len(possible[i])
                for k in range(lenk):
                    cultime += (t_ij[possible[i][k-1]][possible[i][k]])+20
                    sumneed+=t[possible[i][k]]
                    if possible[i][k]==0:
                        cultime=t_ij[0][possible[i][k]]
                        sumneed=t[possible[i][k]]
                        continue
                    elif  cultime>540 or sumneed>270:
                        possible[i].insert(k-1,CENTER)
                        cultime=t_ij[0][possible[i][k]]+20
                        sumneed=t[possible[i][k]]
                        lenk=len(possible[i])
                P0index=[]

    P0index=[]
    cultime=0
    sumneed=0
    for i in range(len(possible)):
        x=0
        for j in possible[i]:
            if j==0:
              P0index.append(x)
        k=P0index[len(P0index)-2]+1
        lenk=len(possible[i])
        for k in range(lenk):
            cultime += (t_ij[possible[i][k-1]][possible[i][k]])+20
            sumneed+=t[possible[i][k]]
            if possible[i][k]==0:
                cultime=t_ij[0][possible[i][k]]
                sumneed=t[possible[i][k]]
                continue
            elif  cultime>500 or sumneed>270:
                possible[i].insert(k-1,CENTER)
                cultime=t_ij[0][possible[i][k]]+20
                sumneed=t[possible[i][k]]
                              
        P0index=[]
    
    for i in range(len(possible)):
        possible[i] =  Gene(data=possible[i].copy())

    possible.sort(reverse=True, key=key)

    if len(possible)==0:
        print('bug1')
        possible.append(gene1)
    crossedGenes.append(possible[0])

    possible = []
    for i in range(7):
        newGene = newGene2.copy()
        newGene.insert(firstPos2, CENTER)
        possible.append(newGene)
        firstPos2 += 1
    
    P0index=[]
    cultime=0
    sumneed=0
    for i in range(len(possible)):
        x=0
        for j in possible[i]:
            if j==0:
              P0index.append(x)
            x+=1
        k=P0index[len(P0index)-2]+1
        lenk=len(possible[i])
        for k in range(lenk):
            cultime += (t_ij[possible[i][k-1]][possible[i][k]])+20
            sumneed+=t[possible[i][k]]
            if possible[i][k]==0:
                cultime=t_ij[0][possible[i][k]]
                sumneed=t[possible[i][k]]
                continue
            elif  cultime>540 or sumneed>270:
                possible[i].insert(k-1,CENTER)
                cultime=t_ij[0][possible[i][k]]+20
                sumneed=t[possible[i][k]]
                lenk=len(possible[i])
           
        P0index=[]
    P0index=[]
    cultime=0
    sumneed=0
    for i in range(len(possible)):
        x=0
        for j in possible[i]:
            if j==0:
              P0index.append(x)
            x+=1
        k=P0index[len(P0index)-2]+1
        lenk=len(possible[i])
        for k in range(lenk):
            cultime += (t_ij[possible[i][k-1]][possible[i][k]])+20
            sumneed+=t[possible[i][k]]
            if possible[i][k]==0:
                cultime=t_ij[0][possible[i][k]]
                sumneed=t[possible[i][k]]
                continue
            elif  cultime>500 or sumneed>270:
                possible[i].insert(k-1,CENTER)
                cultime=t_ij[0][possible[i][k]]+20
                sumneed=t[possible[i][k]]
                lenk=len(possible[i])                  

    for i in range(len(possible)):
        possible[i] =  Gene(data=possible[i].copy())
                 
    key = lambda gene: gene.fit
    possible.sort(reverse=True, key=key)  
     
    crossedGenes.append(possible[0])
    
   
# 交叉
def cross(genes):
    crossedGenes = []
    for i in range(0, len(genes), 2):
        crossPair(genes[i], genes[i+1], crossedGenes)  

    return crossedGenes



# 合併
def mergeGenes(genes, crossedGenes):
    key = lambda gene: gene.chooseProb
    genes.sort(reverse=True, key=key)
    pos = geneNum - 1  
    
    for gene in crossedGenes:
        genes[pos] = gene
        pos -= 1

    return  genes

# 變異一個
def varyOne(gene):
    varyNum = 10
    variedGenes = []
    for i in range(varyNum):
        p1, p2 = random.choices(list(range(1,len(gene.data)-2)), k=2)
        newGene = gene.data.copy()
        newGene[p1], newGene[p2] = newGene[p2], newGene[p1] 
        variedGenes.append(Gene(data=newGene.copy()))
    key = lambda gene: gene.fit
    variedGenes.sort(reverse=True, key=key)

    return variedGenes[0]


# 變異
def vary(genes):
    for index, gene in enumerate(genes):
        if index < 30:                  # 菁英主義，保留前三十
            continue
        if random.random() < VARY:
            genes[index] = varyOne(gene)

    return genes

def printGeneData(self):
    totalcost=timecost=fuelcost=0
        
    for i in range(len(self.carSubpath)):
        print('車輛編號',self.usedcar[i])
        print('車輛行走路徑',self.carSubpath[i])
        carFuel=a_k[self.usedcar[i]]
        for j in range(len(self.carSubpath[i])):
            subDist=self.subpathDist[self.carSubpath[i][j]]
            fuelcost+=subDist*carFuel*19.2
            timecost+=self.pathtimeadd20[self.carSubpath[i][j]]
        print('車輛路徑成本',fuelcost)
        print('車輛行駛時間',timecost)
        totalcost=fuelcost+timecost+966.667
        print('車輛總成本',totalcost)
        totalcost=fuelcost=timecost=0


if __name__ == "__main__" and not DEBUG:
    genes = getRandomGenes(geneNum) # 初始群種

    #迭代
    for i in tqdm(range(generationNum)):
        updateChooseProb(genes)
        sumProb = getSumProb(genes)
        chosenGenes = choose(deepcopy(genes))   # 選擇
        crossedGenes = cross(chosenGenes)   # 交叉
        genes = mergeGenes(genes, crossedGenes) # 複製交叉到子代
        genes = vary(genes)


    
    key = lambda gene: gene.fit
    genes.sort(reverse=True, key=key)

    print('\r\n')
    print('data:', genes[0].data)
    print('fit:', genes[0].genefit)
    
    # printGeneData(genes[0])
    

if DEBUG:
    print("START")
    gene = Gene()
    print(gene.data)
    gene.moveRandSubPathLeft()
    print(gene.data)


    print("FINISH")

