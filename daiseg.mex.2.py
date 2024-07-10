import argparse
import argparse
import numpy as np
import sys
import useful2 as usfl
import HMMmex2 as HMM
import EMmex2 as EM

parser = argparse.ArgumentParser(description='DAIseg') 
parser.add_argument('--bed', type=str, help='Region bed file')
parser.add_argument('--EM', type=str, help='Whether or not to use EM algorithm')
parser.add_argument('--EM_steps', type=str, help='number of EMsteps')
parser.add_argument('--HMM_par', type= str, help='File with parameters')
parser.add_argument('--out_prefix', type= str, help = 'Prefix to output file' )

parser.add_argument('--EM_est', type= str, help = 'Make estimation of the all parameters or only coalescent times' )
parser.add_argument('--prepared_file', type=str, help='dlkfjgljk')
parser.add_argument('--arch_cover', type=str)
parser.add_argument('--obs_samples', type=str, help='File with samples names')
args = parser.parse_args()




dict_all = usfl.main_read2(args.prepared_file)
domain = usfl.read_bed(args.bed)
seq_start, seq_end = domain[0][0], domain[-1][1]
N = 5 # number of hidden states
GEN_time, MU, RR, L, Lambda_0 = usfl.read_par_HMM(args.HMM_par)
n_eu = len(dict_all[list(dict_all.keys())[0]]['Obs'])
epsilon=1e-10

P=[0.4,0.05, 0.4, 0.05, 0.1]
d=MU*L

#read the windows archaic covering 
with open(args.arch_cover,'r') as f:
    cover=f.readlines()
cover = [float(cover[i].replace('\n','')) for i in range(len(cover))]



gaps_numbers, seq_start_mas, seq_end_mas, len_mas=[], [], [], []
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
        
        
        
#creating observations sequence (in gaps is zero)
SEQ=[]
N_ST=[]
state_mas=[]
for ind in range(n_eu):
    o_eu=usfl.make_obs_ref(dict_all, domain, ind, L,  'EU')
    o_na=usfl.make_obs_ref(dict_all, domain, ind, L,  'NA'),
    o_af=usfl.make_obs_ref(dict_all, domain, ind, L,  'AF')
    o_nd=usfl.make_obs_ref(dict_all, domain, ind, L,  'Archaic')
    print(type(o_na))
    
    sq=np.vstack([o_eu, o_na, o_af, o_nd])

    state_mas.append([max(o_eu), max(o_na), max(o_af), max(o_nd)])
    
    
#    sq=np.vstack([usfl.make_obs_ref(dict_all, domain, ind, L,  'EU'), usfl.make_obs_ref(dict_all, domain, ind, L,  'NA'), usfl.make_obs_ref(dict_all, domain, ind, L,  'AF'), usfl.make_obs_ref(dict_all, domain, ind, L,  'Archaic')])
    
    sq=sq.transpose()
    
    n_st = sq.max()+1
    
    SEQ.append(sq)
    N_ST.append(n_st)
    
print(state_mas)
SEQ=np.array(SEQ)
N_st=SEQ.max()+1


#split observations by windows, removing gaps
SEQ_mas=[]
arch_cover=[]
for i in range(len(len_mas)):
    p1=int((seq_start_mas[i]-seq_start_mas[0])/1000)
    p2=int((seq_end_mas[i]-seq_start_mas[0])/1000)
    SEQ_mas.append(SEQ[:,p1:(p2+1)])
    arch_cover.append(cover[p1:(p2+1)])        
        
        
print('Observation sequences has been created')       
        
def run_daiseg(lmbd_opt,seq, n_st, idx, start, ar_cover, a, b_our_mas, b_Skov):



    tracts_HMM =  HMM.get_HMM_tracts(HMM.viterbi_modified(seq [idx], P, a, b_our_mas, b_Skov, ar_cover))

    for k in range(N):
       for j in range(len(tracts_HMM[k])):
           tracts_HMM[k][j][0]= L * tracts_HMM[k][j][0]+start
           tracts_HMM[k][j][1]= L * tracts_HMM[k][j][1]+start-1

    return tracts_HMM

def run_daiseg_all(lmbd_0):

    A = HMM.initA(lmbd_0[5]/d, lmbd_0[6]/d, RR, L, lmbd_0[7],  lmbd_0[8],  lmbd_0[9],  lmbd_0[10])    

    
    B_our_mas = np.array([HMM.initB_arch_cover( lmbd_0, N_st, 0.1+i*0.1) for i in range(10)])
    B_Skov = HMM.initBwN(lmbd_0[0:5], N_st)

    
    tracts_HMM_mas=[]

    
    for idx in range(0, n_eu):    
        tracts_HMM=[[] for i in range(N)]
        for i in range(len(SEQ_mas)):
            tr=run_daiseg(lmbd_0, SEQ_mas[i], N_st, idx, seq_start_mas[i], arch_cover[i], A, B_our_mas, B_Skov)
            for j in range(N):   
               for k in tr[j]:             
                   tracts_HMM[j].append( k )
 

        tracts_HMM_mas.append([tracts_HMM[j] for j in range(N)])
    return tracts_HMM_mas
    
def EM_gaps(seq, lambda_0, n_st, cover):
    return EM.EM_common_gaps(P, seq, n_st, MU, RR, lambda_0, epsilon, L, int(args.EM_steps), gaps_numbers, cover )   
    

if args.EM=='no': 
    Tracts_HMM_mas = run_daiseg_all(Lambda_0)  
    
    
if args.EM=='yes': 
    Lambda_opt = EM_gaps(SEQ[0:1], Lambda_0, N_st, cover)    
    Tracts_HMM_mas = run_daiseg_all(Lambda_opt)    
 
 
with open(args.obs_samples,'r') as f:
    names=f.readlines()
names=[str(names[i].replace('\n','')) for i in range(len(names))]           
        
#write into file arg.o results #Sample #Haplotype_number #Archaic tracts
with open(args.out_prefix+'.neand.eu.txt', "w") as f:
   for i in range(len(Tracts_HMM_mas)):
       if i % 2 ==0:
           f.write(names[int(i // 2)]+'\t0\t'+str(Tracts_HMM_mas[i][1])+'\n')
       else:
           f.write(names[int(i // 2)]+'\t1\t'+str(Tracts_HMM_mas[i][1])+'\n')      
 
        



with open(args.out_prefix+'.neand.na.txt', "w") as f:
   for i in range(len(Tracts_HMM_mas)):
       if i % 2 ==0:
           f.write(names[int(i // 2)]+'\t0\t'+str(Tracts_HMM_mas[i][3])+'\n')
       else:
           f.write(names[int(i // 2)]+'\t1\t'+str(Tracts_HMM_mas[i][3])+'\n')         
        
        
        
        
        
        
        
        
