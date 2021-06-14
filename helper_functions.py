import numpy as np
import random
import math
import copy
import itertools

from qiskit import *
from qiskit.visualization import plot_histogram
from qiskit.providers.aer.noise import pauli_error
from collections import defaultdict
from qiskit import Aer
import numpy as np
import matplotlib.pyplot as plt
import qiskit.quantum_info as qi

import scipy.linalg as la

def insert(current_circuit): #inserts random gate at random position
    types=["Rx","Ry","Rz","Rxx","Ryy","Rzz"]
    #types=["Rx","Ry","Rz"]
    gatetype=random.choice(types)
    theta= 2* math.pi * random.random()
    double_gate=True
    if gatetype in ["Rx","Ry","Rz"]:
        double_gate=False
        qubitnum,p=place_finder3(current_circuit)
        current_circuit=add_single_gate(qubitnum,p,gatetype,theta,current_circuit)
    else:
        [q1,p1],[q2,p2]=place_finder4(current_circuit)
        current_circuit=add_double_gate([q1,q2],p1, p2,gatetype,theta,current_circuit)
        
    return current_circuit
   
def delete(current_circuit): #deletes gate from random position
    if current_circuit==[[] for i in range(len(current_circuit))]:
        return current_circuit
    ind1,ind2=gate_picker(current_circuit)
    if len(current_circuit[ind1][ind2])==2: #means removing single gate
        del current_circuit[ind1][ind2]
        for i in current_circuit:
            for j in i:
                if len(j)==4 and j[2]==ind1 and j[3]>ind2:
                    j[3]-=1
    else: #double gate being removed
        q2=current_circuit[ind1][ind2][2]
        p2=current_circuit[ind1][ind2][3]
        del current_circuit[ind1][ind2]
        for i in current_circuit:
            for j in i:
                if len(j)==4 and j[2]==ind1 and j[3]>ind2:
                    j[3]-=1
                    
        del current_circuit[q2][p2]
        for i in current_circuit:
            for j in i:
                if len(j)==4 and j[2]==q2 and j[3]>p2:
                    j[3]-=1
    return current_circuit
                
def swap(current_circuit): #Combination of Delete and Insert at the same randomly chosen position
    cc=current_circuit
    if cc==[[] for i in range(len(cc))]:
        return cc
    
    ind1,ind2=gate_picker(cc)
    
    if len(cc[ind1][ind2])==2: #means single gate
        types=["Rx","Ry","Rz"]
        gatetype=random.choice(types)
        theta= 2* math.pi * random.random()
        cc[ind1][ind2][0]=gatetype
        cc[ind1][ind2][1]=theta
    else:
        types=["Rxx","Ryy","Rzz"]
        gatetype=random.choice(types)
        theta= 2* math.pi * random.random()
        q2=cc[ind1][ind2][2]
        p2=cc[ind1][ind2][3]
        cc[ind1][ind2][0]=gatetype
        cc[ind1][ind2][1]=theta
        cc[q2][p2][0]=gatetype
        cc[q2][p2][1]=theta
    return cc
    
def modify(current_circuit): #Modify parameter of randomly chosen gate
    cc=current_circuit
    if cc==[[] for i in range(len(cc))]:
        return cc
    
    ind1,ind2=gate_picker(cc)
    mu, sigma = 0, 0.1 # mean and standard deviation
    s = np.random.normal(mu, sigma, 1)[0]
    theta=cc[ind1][ind2][1]+s
    if len(cc[ind1][ind2])==2: #means single gate
        cc[ind1][ind2][1]=theta
    else:
        q2=cc[ind1][ind2][2]
        p2=cc[ind1][ind2][3]
        cc[ind1][ind2][1]=theta
        cc[q2][p2][1]=theta
        
    return cc
    
#helper function

def add_single_gate(qubitnum,ind,gatetype,gate_arg,current_circuit): #adds given single gate at given position
    for i in current_circuit:
        for j in i:
            if len(j)==4:
                if j[2]==qubitnum and j[3]>=ind:
                    j[3]+=1
    current_circuit[qubitnum].insert(ind,[gatetype,gate_arg])                
    return current_circuit 
            
def add_double_gate(qubitnums,ind1, ind2,gatetype,gate_arg,current_circuit): #adds given double gate at given position
    q1=qubitnums[0]
    q2=qubitnums[1]
    
    for i in current_circuit:
        for j in i:
            if len(j)==4:
                if j[2]==q1 and j[3]>=ind1:
                    j[3]+=1
    for i in current_circuit:
        for j in i:
            if len(j)==4:
                if j[2]==q2 and j[3]>=ind2:
                    j[3]+=1
    
    current_circuit[q1].insert(ind1,[gatetype,gate_arg,q2,ind2])
    current_circuit[q2].insert(ind2,[gatetype,gate_arg,q1,ind1])
    return current_circuit
    

def g(cc,q1,p1,q2,ans,done): #finding index of contact on left side
    #print("q1,p1,q2",q1,p1,q2)
    if (q1,p1) in done:
        return ans
    if p1==-1:
        return ans
    if len(cc[q1][p1])==2: #means its a single gate
        return g(cc,q1,p1-1,q2,ans,done) #move to next left gate
    else:
        if cc[q1][p1][2]==q2:
            ans.add(cc[q1][p1][3])
        #done.append(q1) #so that we don't enter this qubit multiple times
        done.append((q1,p1)) #so that we don't start exploring this cell while coming from some other path
        temp=set()
        temp=g(cc,cc[q1][p1][2],cc[q1][p1][3],q2,ans,done) #starting the whole thing for another qubit
        return temp.union(g(cc,q1,p1-1,q2,ans,done))
            
def h(cc,q1,p1,q2,ans,done): #finding index of contact on right side
    #print("q1,p1,q2",q1,p1,q2)
    if (q1,p1) in done:
        return ans
    if p1==len(cc[q1]):
        return ans
    if len(cc[q1][p1])==2:
        return h(cc,q1,p1+1,q2,ans,done)
    else:
        if cc[q1][p1][2]==q2:
            ans.add(cc[q1][p1][3])
        #done.append(q1) #so that we don't enter this qubit multiple times
        done.append((q1,p1)) #so that we don't start exploring this cell while coming from some other path
        temp=set()
        temp=h(cc,cc[q1][p1][2],cc[q1][p1][3],q2,ans,done) #starting the whole thing for another qubit
        return temp.union(h(cc,q1,p1+1,q2,ans,done))

    
def is_allowed(cc,start,end): #start=q1,p1 endq2,p2 #This function tells whether we can insert a double gate between 2 given positions
    #print("is_allowed,cc,start,end",cc,start,end)
    q1=start[0]
    p1=start[1]
    q2=end[0]
    p2=end[1]
    sett=set()
    done=[]
    a1=list(g(cc,q1,p1-1,q2,sett,done))
    sett=set()
    done=[]
    a2=list(h(cc,q1,p1,q2,sett,done))
    #print("a1",a1)
    for i in a1:
        if i>=p2:
            return False
    #print("a2",a2)
    for i in a2:
        if i<p2:
            return False
    return True


def place_finder1(current_circuit):  
    cc=current_circuit
    ct=0
    for i in cc:
        ct=ct+len(i)+1
    #print(ct)
    p=random.randrange(ct)
    #p=15 #REMOVE
    #print("p",p)
    j=0
    ans=[0,0]
    while True:
        if j==p:
            return ans
        if len(cc[ans[0]])==ans[1]:
            ans=[ans[0]+1,0]
        else:
            ans=[ans[0],ans[1]+1]
        #print("ans",ans)
        j=j+1
        
def place_finder2(cc):
    possible=[]
    for i in range(len(cc)):
        for j in range(i+1,len(cc)):
            for k in range(len(cc[i])+1):
                for l in range(len(cc[j])+1):
                    if is_allowed(cc,[i,k],[j,l]):
                        possible.append([[i,k],[j,l]])
    ans=random.choice(possible)
    return ans


def place_finder3(current_circuit): #place finder for single gate insertion
    cc=current_circuit
    qubitnum=random.randrange(len(cc))
    p=random.randrange(len(cc[qubitnum])+1)
    return [qubitnum,p]

def place_finder4(cc): #place finder for double gate insertion
    qubitnum=sorted(random.sample(set(range(len(cc))), 2))
    q1=qubitnum[0]
    q2=qubitnum[1]
    possible=[]
    for i in range(len(cc[q1])+1):
        for j in range(len(cc[q2])+1):
            if is_allowed(cc,[q1,i],[q2,j]):
                possible.append([[q1,i],[q2,j]])
    ans=random.choice(possible)
    return ans
    

def gate_picker(cc): #randomly picks a gate from the whole circuit
    al=[]
    for i in range(len(cc)):
        for j in range(len(cc[i])):
            al.append([i,j])
    for i in al:
        if len(cc[i[0]][i[1]])==4:
            al.remove([cc[i[0]][i[1]][2],cc[i[0]][i[1]][3]])
    ans=random.choice(al)
    return ans  

def gate_picker2(cc):
    al=[]
    qubitnum=random.randrange(len(cc))
    p=random.randrange(len(cc[qubitnum]))
    return [qubitnum,p]   



def add_gate(circuit,qb,gatetype,gate_arg): #add gate to actual qiskit circuit
    #print("add_gate function",gatetype,qb)
    if gatetype=="Rx":
        circuit.rx(gate_arg,qb)
    if gatetype=="Ry":
        circuit.ry(gate_arg,qb)
    if gatetype=="Rz":
        circuit.rz(gate_arg,qb)
    if gatetype=="Rxx":
        #print("qb",qb)
        circuit.rxx(gate_arg,qb[0],qb[1])
    if gatetype=="Ryy":
        #print("qb,gatearg",qb,gate_arg)
        circuit.ryy(gate_arg,qb[0],qb[1])
    if gatetype=="Rzz":
        #print("qb",qb)
        circuit.rzz(gate_arg,qb[0],qb[1])
    return circuit


def get_circuit(cc): #get actual qiskit circuit
    cq=QuantumCircuit(QuantumRegister(len(cc),'q'))
    pointer=[0]*len(cc)
    done=[False]*len(cc)
    #print(3)
    while False in done:
        for i in range(len(cc)):
            #print(i)
            #print(len(cc))
            if len(cc[i])==0:
                done[i]=True
            #print(4)
            for j in range(pointer[i],len(cc[i])):
                if j==len(cc[i])-1:
                    done[i]=True
                if len(cc[i][j])==4:
                    pointer[i]=j
                    break
                #print("add_gate",i,j)
                cq=add_gate(cq,i,cc[i][j][0],cc[i][j][1])
                pointer[i]+=1
        #print("pointers",pointer)
        pointer,cq=solve_dg(pointer,cc,cq)
        #print("pointer",pointer)
        #print("cc",cc)
    return cq
    

def solve_dg(pointer,cc,cq): #maintains the oorder in which double gates have to be added
    #print("solve_dg")
    l1=[]
    l2=[]
    for i in range(len(cc)):
        if pointer[i]==len(cc[i]):
            continue
        if ([i,pointer[i],cc[i][pointer[i]][2],cc[i][pointer[i]][3]] in l1) or ([cc[i][pointer[i]][2],cc[i][pointer[i]][3],i,pointer[i]]) in l1:
            l2.append([cc[i][pointer[i]][0],cc[i][pointer[i]][1],i,cc[i][pointer[i]][2]]) #[gatetype,arg,q1,q2]
        else:
            l1.append([i,pointer[i],cc[i][pointer[i]][2],cc[i][pointer[i]][3]])
    for i in l2:
        #print("l2",l2)
        pointer[i[2]]+=1
        pointer[i[3]]+=1
        cq=add_gate(cq,[i[2],i[3]],i[0],i[1])
    return [pointer,cq]


def get_min_circuit(H):
    pick=["i","i","i","i","i","d","s","m","m","m"]
    results=[]
    #pick=["i","i","d","s","m"]
    initial=[[] for i in range(10)]
    #initial=get_circuit(current_circuit)
    #initial=insert(initial)
    #sorted_circuits
    initial=[9999999999,initial]
    
    total_generations=300
    for i in range(total_generations):
        pick=["i","i","i","i","i","d","s","m","m","m"]
        offsprings=[initial]
        for j in range(4):
            cp=copy.deepcopy(initial[1])
            action=random.choice(pick)
            if action=="i":
                cp=insert(cp)
            elif action=="d":
                cp=delete(cp)
            elif action=="s":
                cp=swap(cp)
            else:
                cp=modify(cp)
            #print(1)
            qc=get_circuit(cp)
            starting_circuit=QuantumCircuit(QuantumRegister(10,'q'))
            for ii in range(10):
                #starting_circuit.h(ii)                
                pass
            qc2=starting_circuit+qc
            #print(2)
            #print("final_state",final_state(qc2))
            BR=bra_psi_U(final_state(qc2))
            KT=ket_psi_U(final_state(qc2))
            #print("BR",BR)
            #print("KT",KT)
            m1=np.matmul(BR,H)
            m2=np.matmul(m1,KT)
            offsprings.append([m2,cp])
            #print(m2)
        offsprings.sort()
        initial=offsprings[0]
        results.append(initial)
    return results
    
    
def final_state(circuit):
    svsim = Aer.get_backend('statevector_simulator')
    qobj = assemble(circuit)
    s = svsim.run(qobj).result().get_statevector() #q0 when empty circuit #psiU when circuit made
    return s
    
def bra_psi_U(s):
    temp=s.conjugate()
    return temp
    
def ket_psi_U(s):
    temp=[]
    for i in s:
        temp.append([i])
    return temp


def gamma(i, circuit):
    if i%2==0:
        for j in range(i//2):
            circuit.z(j)
        circuit.y(i//2)
    else:
        for j in range((i-1)//2):
            circuit.z(j)
        circuit.x((i-1)//2)
    return circuit

def gamma_quad_product(i1,i2,i3,i4):
    q=QuantumRegister(10,'q')
    circuit=QuantumCircuit(q)
    circuit=gamma(i1,circuit)
    circuit=gamma(i2,circuit)
    circuit=gamma(i3,circuit)
    circuit=gamma(i4,circuit)
    op=qi.Operator(circuit)
    array=op.data
    return array

def gamma_factor(i1,i2,i3,i4):
    ct=0
    for i in [i1,i2,i3,i4]:
        if i%2==0:
            ct=ct+1
    x = np.random.normal(loc=0, scale=((3/4000)**(0.5)))
    return 16*x*(-1)

def Hamiltonian_for_SYK():
    a20= [i for i in range(20)]
    a20=set(a20)
    a_set = a20

    data = list(set((itertools.combinations(a_set, 4))))
    data.sort()
    #subsets = set(data)
    for i in range(len(data)):
        data[i]=list(data[i])
        data[i].sort()
    #print(data)
    #print(type(data[0]))

    H=[[0 for i in range(1024)] for j in range(1024) ] #gives us Hamiltonian

    crr=1
    for i in data:
        #print(crr)
        #print(i[0],i[1],i[2],i[3])
        crr+=1
        arr=gamma_quad_product(i[0],i[1],i[2],i[3])
        fct=gamma_factor(i[0],i[1],i[2],i[3])
        for j in range(len(arr)):
            for k in range(len(arr[0])):
                arr[j][k]=fct*arr[j][k]
        H+=arr
    return H
    
    


