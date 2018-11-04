import numpy as np
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import cosine
# import pylab
import math
import scipy.cluster.hierarchy as sch

AUM=[[1,0,1,0],
     [0,1,1,0],
     [0,1,0,1],
     [0,0,1,1]]

# f = open ( 'aummatrix.txt' , 'r')
# AUM=[ map(int,line.split(',')) for line in f ]
#print (AUM)
AF=[45,5,75,3]
AFM=[[15,5,25,3],
    [20,0,25,0],
    [10,0,25,0]]
A=4
Q=4
S=3
sum1=0
AAM = np.zeros((A,A))
CA = np.zeros((A,A))

array = []
maxIndex = []
maxIndex2 = []
index = 0




#Calculate global bond energy
def cont(ai, ak, aj):
    if ak == aj:
        return 2 * bond(ai, ak) + 2 * bond(ak, aj)
    else:
        return 2 * bond(ai, ak) + 2 * bond(ak, aj) - 2 * bond(ai, aj)

#Calculate bond energy for 2 columns
def bond(ax, ay):
    result = 0
    if ax<0 or ay<0:
        return 0
    if ax == index:
        for i in range(0,A):
            result += array[i] * CA[i,ay]
        return result
    if ay == index:
        for i in range(0,A):
            result += array[i] * CA[i,ax]
        return result
    for i in range(0,A):
        result += CA[i,ax] * CA[i,ay]
    return result
    
def maxCont(res):
    maxValue = 0
    for i in range(len(res)):
        if res[i][3]>maxValue:
            maxValue = res[i][3]
            maxindex = res[i][2]
    return maxindex

def BEA(AA):
    #Add two columns to new matrix
    for i in range(0,A):
        global CA
        CA[i,0] = AA[i,0]
        CA[i,1] = AA[i,1]
    global index
    index = 2
    #Create array for storing the results
    s = [0,0,0,0]
    while (index<A):
        results = []
        global array
        array = AA[:,index]
        i = 0
        while (i<=index):
            result = cont(i-1,index,i)
            s[0]=i-1
            s[1]=index
            s[2]=i
            s[3]=result
            v=[]
            v=s[:]
            print(s)
            print(v)
            results.append(v)
            i+=1
        CA = np.insert(CA, maxCont(results), np.array((array)), 1)
        maxIndex.append(maxCont(results))
        maxIndex2.append(index)
        CA = np.delete(CA, (-1), 1)
        CA[:,[0, 1]] = CA[:,[1, 0]]
        index+=1
    return CA
    

print ("The no of attribute in the AUM matrix along Column and AUM matrix:")
print (A)
print (AUM)
print ("\nThe of queries are:")
print (Q)
print ("\nThe no of sites and access frequency matrix:")
print (S)
print (AF)
print(AFM)


#Calculating Attribute affinity matrix
for j in range(A):
    for i in range(A):
        for k in range(Q):
            if (i==j):
                if (AUM[k][j]==1):
                    for l in range(S):
                        sum1+=AFM[l][k]                    
            else:
                if ((AUM[k][j]==1) and (AUM[k][i]==1)):
                    for l in range(S):
                        sum1+=AFM[l][k]
        AAM[j][i]=sum1
        sum1=0


print ("\nThe Attribute affinity matrix is :\n")
print (AAM)

BEA(AAM)
#Rearrange the rows according to the order of columns
for i in range(len(maxIndex)):
    CA[[maxIndex[i], maxIndex2[i]],:] = CA[[maxIndex2[i], maxIndex[i]],:]
print ("\nThe clustered affinity matrix is:\n")
print(CA)





print ("\n\nThe partion algorithm stat from here:")
VA = np.zeros((A,Q))
SA = np.zeros((A,A))
print(VA)
print(SA)
# for i in range(A):
#     for j in range(Q):
#         if (AUM[i][j]==1):
#             VA[j][i]=AF[i]
#         else:
#             VA[j][i]=0

for i in range(A):
    for j in range(Q):
        sum=0
        for k in range(S):
            sum+=AFM[k][j]
        if (AUM[j][i]==1):
            VA[i][j]=sum
        else:
            VA[i][j]=0
print ("\n\nReference feature vectors table:")
print (VA)

print ("\n\nThe similarity matrix is given by:")
SA = 1-pairwise_distances(VA, metric="cosine")
print (SA)


def calculateInitialThresold(SA):
    DSA = pairwise_distances(VA, metric="cosine")
    print ("\n\n Dissimilarity :- \n")
    print (DSA)

    #DSA = 1 - pairwise_distances(VA, metric="cosine")
    w = 1
    Threshold = []
    print ("\n\n Threshold calculation :- ")
    for i in range(A):
        minimum_threshold = DSA[i][0]
        min_difference = 1000
        total_diff = 0
        for k in range(A):
            for j in range(A):
                total_diff += abs(DSA[i][j] - w*DSA[i][k])
            if (total_diff < min_difference):
                min_difference = total_diff
                minimum_threshold = DSA[i][k]
        Threshold.append(minimum_threshold)
    print ("\n\n Initial Threshold :- \n")
    print (Threshold)


F1=[]
F2=[]
print ("\n\nThe fragmentation on the basis of cosine similarity are:")
for i in range(A):
    for j in range(A):
        if(SA[i][j]>0 and (j not in F1) and (j not in F2)):
            F1.append(j)
        elif(SA[i][j]==0 and (j not in F1) and (j not in F2)):
            F2.append(j)
        else:
            continue

print ("\n\nFragment F1 and F2 are:")
print (F1)
print (F2)


# ========= Similarity based on clustred =====
print ("\n\nThe similarity matrix is given by after CAF:")
dissimalarity_CA = pairwise_distances(CA, metric="cosine")
print (dissimalarity_CA)
calculateInitialThresold(SA)


# print ("\nThe query executed in centralized environment\n\n")
# L=[]
# for i in range(Q):
#     for j in range(A):
#         if (AUM[i][j]==1):
#             L.append(j)
#     if (set(F1) & set(L)==set(L)):
#         print ("The Query {0} is executed by F1 only:".format(i))
#         print ("The time of query execution if one query executed in 1ms is:")
#         print (AF[i]*1)
#     elif (set(F2) & set(L)==set(L)):
#         print ("The Query {0} is executed by F2 only:".format(i))
#         print ("The time of query execution if one query executed in 1ms is:")
#         print (AF[i]*1)
#     else:
#         print ("The query {0} is executed with the use of both fragments:".format
#                (i))
#         print ("The time of query execution if one query executed in 1ms is:")
#         print (AF[i]*2)
#     L=[]


# print ("\nThe query executed in distributed environment\n\n")
# L=[]
# for i in range(Q):
#     for j in range(A):
#         if (AUM[i][j]==1):
#             L.append(j)
#     if (set(F1) & set(L)==set(L)):
#         print ("The Query {0} is executed by F1 only:".format(i))
#         print ("The time of query execution if one query executed in 1ms is:")
#         print (AF[i]*1)
#     elif (set(F2) & set(L)==set(L)):
#         print ("The Query {0} is executed by F2 only:".format(i))
#         print ("The time of query execution if one query executed in 1ms is:")
#         print (AF[i]*1)
#     else:
#         print ("The query {0} is executed with the use of both fragments:".format
#                (i))
#         print ("The time of query execution if one query executed in 1ms is:")
#         print (AF[i]*1)
#     L=[]

