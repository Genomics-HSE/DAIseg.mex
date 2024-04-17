import numpy as np



f_obs = "obs.chr22.ingroup.txt"

f_neand = "chr22.archaic.reference.txt"
f_yri = "chr22.yri.reference.txt"
f_na= "chr22.na.reference.txt"
f_eu = "chr22.eu.reference.txt"



def make_obs(lines, lines_ref, L, ind):

    start, end = lines_ref[0,0], lines_ref[-1,0]
    if float(end-start) % L != 0:
        T = int((end-start)/ L) + 1
    obs_ref = np.zeros(T, int)
    
    obs = lines[:,ind]
    
    for i in range(len(lines)):
        j=int((lines_ref[i][0]-start)/L)

        if obs[i]==0:
            if lines_ref[i][1]==-1:
                obs_ref[j]+=1
        if obs[i]==1:
            if lines_ref[i][2]==-1:
                obs_ref[j]+=1
    return obs_ref


with open(f_obs,'r') as f:
    
    
    lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].strip(' \n').split(' ')

        for j in range(len(lines[i])):
            lines[i][j]=int(lines[i][j])
lines = np.array(lines)  

    
    
   
with open(f_neand,'r') as f:
    lines_neand = f.readlines()
    
with open(f_yri,'r') as f:
    lines_yri = f.readlines()

with open(f_eu,'r') as f:
    lines_eu = f.readlines()
    
with open(f_na,'r') as f:
    lines_na = f.readlines()
        
    
for i in range(len(lines_yri)):
    lines_yri[i] = lines_yri[i].strip('\n').split('\t')
    for j in range(len(lines_yri[i])):
        lines_yri[i][j]=int(lines_yri[i][j])
        

        
        
        
for i in range(len(lines_eu)):
    lines_eu[i] = lines_eu[i].strip('\n').split('\t')
    for j in range(len(lines_eu[i])):
        lines_eu[i][j]=int(lines_eu[i][j])
        
for i in range(len(lines_na)):
    lines_na[i] = lines_na[i].strip('\n').split('\t')
    for j in range(len(lines_na[i])):
        lines_na[i][j]=int(lines_na[i][j])
        
for i in range(len(lines_neand)):
    lines_neand[i] = lines_neand[i].strip('\n').split('\t')
    for j in range(len(lines_neand[i])):
        lines_neand[i][j]=int(lines_neand[i][j])    
        
lines_yri = np.array(lines_yri)
lines_eu = np.array(lines_eu)
lines_na = np.array(lines_na)
lines_neand = np.array(lines_neand)





start, end = lines_yri[0,0], lines_yri[-1,0]

n_mex = len(lines[0])

SEQ=[]
N_ST=[]
L=1000

MU=1.25e-8
RR=1e-8
gen_time=29

T_NEAND_migration = 55000 #time of Neanderthal migration into Out_of_africa population
T_NEAND_AMH = 650000 # split time between AMH and Neanderthal
T_OOF_AF = 65700 # Out_of_Africa migration time
T_AS_EU_split = 41997 #Asian European split
T_AS_NA_split = 17000 # time of crossing Beringian (Nielsen et al., 2017; Mendes et al., 2020)
T_MEX_admix = 466 

with open('par.file.txt', "w") as file_pos:
    file_pos.write(str(gen_time)+'\n')
    file_pos.write(str(MU)+'\n')
    file_pos.write(str(RR)+'\n')
    file_pos.write(str(L)+'\n')
    file_pos.write(str(lines_yri[0,0])+' ' +str(lines_yri[-1,0])+'\n')
    file_pos.write(str( T_NEAND_AMH )+'\n')
    file_pos.write(str(T_OOF_AF)+'\n')
    file_pos.write(str( T_NEAND_migration )+'\n')
    file_pos.write(str( T_AS_EU_split )+'\n')
    file_pos.write(str( T_MEX_admix  )+'\n')
    file_pos.write(str( T_NEAND_migration )+'\n')
    file_pos.write(str( T_MEX_admix  )+'\n')
    file_pos.write(str(0.025)+'\n')
    file_pos.write(str(0.45)+'\n')
    file_pos.write(str(0.45)+'\n')
    file_pos.write(str(0.1))
        

for ind in range(n_mex):
    sq=np.vstack([make_obs(lines, lines_eu, L, ind ),make_obs(lines, lines_na, L, ind ), make_obs(lines, lines_yri, L, ind ),make_obs(lines, lines_neand, L, ind)])
    sq=sq.transpose()
    n_st = sq.max()+1
    SEQ.append(sq)
    N_ST.append(n_st)
SEQ=np.array(SEQ)





with open('obs.eu.txt', "w") as file1, open('obs.na.txt', "w") as file2, open('obs.yri.txt', "w") as file3,open('obs.neand.txt', "w") as file4:
    for j in range(len(SEQ[0])):
        s1, s2, s3, s4 = '','', '', ''
        for i in range(n_mex):
            s1 += str(SEQ[:,j,0][i])+' '
            s2 += str(SEQ[:,j,1][i])+' '
            s3 += str(SEQ[:,j,2][i])+' '
            s4 += str(SEQ[:,j,3][i])+' '                
        file1.write(s1[:-1]+'\n')
        file2.write(s2[:-1]+'\n')
        file3.write(s3[:-1]+'\n')
        file4.write(s4[:-1]+'\n')       
        
