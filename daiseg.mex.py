import argparse
import argparse
import numpy as np
import HMMmex as hmm
import EMmex as EM
import sys
from functools import reduce
from operator import xor
from itertools import chain


def intersections(a,b):
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]

        if a_right < b_right:
            i += 1
        else:
            j += 1

        if a_right >= b_left and b_right >= a_left:
            end_pts = sorted([a_left, a_right, b_left, b_right])
            middle = [end_pts[1], end_pts[2]]
            ranges.append(middle)

    ri = 0
    while ri < len(ranges)-1:
        if ranges[ri][1] == ranges[ri+1][0]:
            ranges[ri:ri+2] = [[ranges[ri][0], ranges[ri+1][1]]]
        ri += 1
    return ranges

parser = argparse.ArgumentParser(description='DAIseg') 

parser.add_argument('--location', type=str, help='File with first-last positions on chr')
parser.add_argument('--gaps', type=str, help='File with gaps')

parser.add_argument('--obs_eu', type=str, help='File with observations with respect to Europeans')
parser.add_argument('--obs_na', type=str, help='File with observations with respect to Americans')
parser.add_argument('--obs_af', type=str, help='File with observations with respect to Africans')
parser.add_argument('--obs_archaic', type=str, help='File with observations with respect to Archaic reference genomes')

parser.add_argument('--EM', type=str, help='Whether or not to use EM algorithm')

parser.add_argument('--HMM_par', type= str, help='File with parameters')
parser.add_argument('--o_eu', type= str, help = 'Outfile for european-neanderthal ancestry' )
parser.add_argument('--o_na', type=str, help='Outfile for american-neanderthal ancestry')
parser.add_argument('--EM_est', type= str, help = 'Make estimation of the all parameters or only coalescent times' )



args = parser.parse_args()


N = 5 # number of hidden states

with open(args.location,'r') as f1:

    seq_start, seq_end = f1.readline().split(' ')
    seq_start = int(seq_start)
    seq_end = int(seq_end.replace('\n',''))    
    

with open(args.gaps,'r') as f:
    l=f.readline()
m=l.replace('\n','').replace('[','').replace(']','').split(',')

m=[[int(m[2*i]), int(m[2*i+1])] for i in range(int(len(m)/2))]

    


set2 = m
set1 = [[seq_start, seq_end]]

l = sorted((reduce(xor, map(set, chain(set1 , set2)))))
XOR=[l[i:i + 2] for i in range(0, len(l), 2)]


domain=intersections(XOR, set1)





    

f = open(args.HMM_par, 'r')
GEN_time = float(f.readline())
MU = float(f.readline())
RR = float(f.readline())
L = int(f.readline())

#seq_start, seq_end = f.readline().split(' ')
#seq_start = int(seq_start)
#seq_end = int(seq_end.replace('\n',''))



Lambda_0=np.zeros(11)
Lambda_0[1] = float(f.readline())/GEN_time*MU*L
Lambda_0[2] = float(f.readline())/GEN_time*MU*L
Lambda_0[0] = float(f.readline())/GEN_time*MU*L
Lambda_0[3] = float(f.readline())/GEN_time*MU* L
Lambda_0[4] = float(f.readline())/GEN_time*MU* L
Lambda_0[5] = float(f.readline())/GEN_time*MU* L
Lambda_0[6] = float(f.readline())/GEN_time*MU* L
Lambda_0[9] = float(f.readline())
Lambda_0[10] = Lambda_0[9]
Lambda_0[7] = float(f.readline())
Lambda_0[8] = float(f.readline())

f.close()

d=MU*L


seq0, seq1, seq2, seq3 = [], [], [], []
with open(args.obs_eu, "r") as f0, open(args.obs_na, "r") as f1, open(args.obs_af, "r") as f2,open(args.obs_archaic, "r") as f3:
    for line0, line1, line2, line3 in zip(f0, f1, f2, f3):
    
        row = line0.replace('\n','').split(' ')
        row = [int(i) for i in row]
        seq0.append(row)
        
        row = line1.replace('\n','').split(' ') 
        row = [int(i) for i in row]
        seq1.append(row)
        
        row = line2.replace('\n','').split(' ')
        row = [int(i) for i in row]
        seq2.append(row)

        row = line3.replace('\n','').split(' ')
        row = [int(i) for i in row]
        seq3.append(row)

seq0=np.array(seq0)
seq0 = np.transpose(seq0)        
seq1=np.array(seq1)
seq1 = np.transpose(seq1)
seq2=np.array(seq2)
seq2 = np.transpose(seq2)
seq3=np.array(seq3)
seq3 = np.transpose(seq3)



n0=seq0.max()
n1=seq1.max()
n2=seq2.max()
n3=seq3.max()


seq=[]
for i in range(len(seq1)):
    seq.append(np.column_stack((seq0[i],seq1[i], seq2[i], seq3[i])))   

SEQ=np.array(seq)
N_st=SEQ.max()+1

P=[0.4, 0.05, 0.4, 0.05, 0.1]
n_EM_steps = 10
epsilon = 1e-10
bnds=((0.5*Lambda_0[5], 2*Lambda_0[5]), (0.5*Lambda_0[6], 2*Lambda_0[6]))



def run_daiseg(lmbd_opt,seq, n_st, idx):
    d = MU * L       
    lmbd_opt=Lambda_0        
    a = hmm.initA(lmbd_opt[5]/d, lmbd_opt[6]/d, RR, L, lmbd_opt[7],  lmbd_opt[8],  lmbd_opt[9],  lmbd_opt[10])
    b=hmm.initB(MU, L, lmbd_opt[0:5], n_st) 
    tracts_HMM =  hmm.get_HMM_tracts(hmm.viterbi(seq [idx], P, a, b))
    for k in range(N):
        for j in range(len(tracts_HMM[k])):
            tracts_HMM[k][j][0]= L * tracts_HMM[k][j][0] +seq_start
            tracts_HMM[k][j][1]= L * (tracts_HMM[k][j][1]+1)-1 +seq_start

    return tracts_HMM
    
    
def EM_function(seq, lambda_0, n_st):

    Lambda_new = EM.EM_algorithm(P, seq, n_st, MU, RR, Lambda_0, epsilon, L, bnds)
    return Lambda_new
    
def EM2(seq, lambda_0,n_st):
    return EM.EM_common(P, seq, n_st, MU, RR, Lambda_0, epsilon, L, bnds)

tracts_HMM_result = []
if args.EM=='no':

    for idx in range(0, len(seq)):
        tracts_HMM_result.append(run_daiseg(Lambda_0, SEQ, N_st, idx))



        
if args.EM=='yes': 
               
    Lambda_opt = EM2(SEQ, Lambda_0, N_st) 
    for idx in range(0, len(seq)):       
        tracts_HMM_result.append(run_daiseg(Lambda_opt, SEQ, N_st, idx)) 
        

        
with open(args.o_eu, "w") as f:
   for idx in range(0, len(seq)):      

       f.write(str(intersections(tracts_HMM_result[idx][1], domain))+'\n') 
       
with open(args.o_na, "w") as f:
   for idx in range(0, len(seq)):
       

       f.write(str(intersections(tracts_HMM_result[idx][3], domain))+'\n') 
