import argparse
import numpy as np
import HMMmex as hmm
import EMmex as EM
import sys
import useful as usfl


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
parser.add_argument('--EM_steps', type=int, help='Number of steps in EM algorithm')
args = parser.parse_args()
N = 5 # number of hidden states

with open(args.location,'r') as f1:

    seq_start, seq_end = f1.readline().split(' ')
    seq_start = int(seq_start)
    seq_end = int(seq_end.replace('\n',''))    
    


f = open(args.HMM_par, 'r')
GEN_time = float(f.readline())
MU = float(f.readline())
RR = float(f.readline())
L = int(f.readline())

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





if args.gaps is not None:
    with open(args.gaps,'r') as f:
        l=f.readline()
    m=l.replace('\n','').replace('[','').replace(']','').split(',')
    m=[[int(m[2*i]), int(m[2*i+1])] for i in range(int(len(m)/2))]   
    domain=usfl.exclude_gaps([[seq_start, seq_end]], m)


    #list of gaps numbers consistent with windows and starting position
    gaps_numbers=[]
    seq_start_mas=[]
    seq_end_mas=[]

    len_mas=[]
    for i in range(len(domain)):
        if (domain[i][0] // L)*L + (domain[0][0]% L) >= domain[i][0]:
            seq_start_mas.append((domain[i][0] // L)*L + (domain[0][0]% L))
        else:
            seq_start_mas.append((domain[i][0] // L)*L + (domain[0][0]% L)+L)
        if (domain[i][1]//L)*L-1+(domain[0][0]% L) <= domain[i][1]:
            seq_end_mas.append((domain[i][1]//L)*L-1+(domain[0][0]% L))
        else:
            seq_end_mas.append((domain[i][1]//L)*L-1+(domain[0][0]% L)-L)        
        
        len_mas.append(int((domain[i][1]-domain[i][0])/1000))

        if i!=len(domain)-1:
            gaps_numbers.append([int((domain[i][1]-domain[0][0])/1000),int((domain[i+1][0]-domain[0][0])/1000)] )

    domain=[[seq_start_mas[i], seq_end_mas[i]] for i in range(len(domain))]


    SEQ_mas=[]
    for i in range(len(len_mas)):
        p1=int((seq_start_mas[i]-seq_start_mas[0])/1000)
        p2=int((seq_end_mas[i]-seq_start_mas[0])/1000)
        SEQ_mas.append(SEQ[:,p1:(p2+1)])


else:
    SEQ_mas=[SEQ]
    gaps_numbers=[[]]
    seq_start_mas=[seq_start]
    domain=[[seq_start, seq_end]]


P=[0.4, 0.05, 0.4, 0.05, 0.1]

epsilon = 1e-9
bnds=((0.2*Lambda_0[5], 100*Lambda_0[5]), (0.2*Lambda_0[6], 100*Lambda_0[6]))



def run_daiseg(lmbd_opt,seq, n_st, idx, start):
    d = MU * L       
    seq=np.array(seq)
    a = hmm.initA(lmbd_opt[5]/d, lmbd_opt[6]/d, RR, L, lmbd_opt[7],  lmbd_opt[8],  lmbd_opt[9],  lmbd_opt[10])
    b=hmm.initB(MU, L, lmbd_opt[0:5], n_st) 
    
    tracts_HMM =  hmm.get_HMM_tracts(hmm.viterbi(seq [idx], P, a, b))
    for k in range(N):
        for j in range(len(tracts_HMM[k])):
            tracts_HMM[k][j][0]= L * tracts_HMM[k][j][0] +start
            tracts_HMM[k][j][1]= L * (tracts_HMM[k][j][1]+1)-1 +start
            

    return tracts_HMM
    



def EM2_gaps(seq, lambda_0,n_st):
    if args.EM_est=='coal':
        return EM.EM_common_gaps(P, seq, n_st, MU, RR, Lambda_0, epsilon, L, bnds, args.EM_steps, gaps_numbers, True)
#    else:
#        return EM.EM_common_gaps(P, seq, n_st, MU, RR, Lambda_0, epsilon, L, bnds, args.EM_steps, gaps_numbers, False)
        
    
def run_daiseg_all(lmbd_0):
    tracts_HMM_mas=[]

    
    for idx in range(0, len(seq)):    
        tracts_HMM=[[],[],[],[],[]]
        for i in range(len(SEQ_mas)):
            tr=run_daiseg(lmbd_0, SEQ_mas[i], N_st, idx, seq_start_mas[i])
            for j in range(5):   
               for k in tr[j]:             
                   tracts_HMM[j].append( k )
 

        tracts_HMM_mas.append([tracts_HMM[j] for j in range(5)])
    return tracts_HMM_mas


if args.EM=='no': 
    Tracts_HMM_mas = run_daiseg_all(Lambda_0)

     
if args.EM=='yes': 
    tracts_HMM_mas=[]     

    Lambda_opt = EM2_gaps(SEQ, Lambda_0, N_st) 
    Tracts_HMM_mas = run_daiseg_all(Lambda_opt)
        

        
with open(args.o_eu, "w") as f:
   for idx in range(0, len(seq)): 
       f.write(str(Tracts_HMM_mas[idx][1])+'\n') 
       
with open(args.o_na, "w") as f:
   for idx in range(0, len(seq)):
       f.write(str(Tracts_HMM_mas[idx][3])+'\n') 
