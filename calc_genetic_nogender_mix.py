#%% 遺伝因子強度を計算する際のコントロールを変更
# ツイン間の性別を固定する
# 2021-03-15
# envggenetic.py を改変

cpuid = 0
#%% ライブラリ
import itertools
import pandas as pd
import pickle
import random
from scipy.stats import levene
import statistics
#%% Import Methylation data and monozygonal twin pairs
# 245 pairs ( 490 individuals)
# num. probes: (485,558) 

methyl_file = "data/methyl_T_intersection.pickle"

with open(methyl_file,'rb') as df_pickle:
    df_methyl = pickle.load(df_pickle)
    
with open("data/twin_mono_intersection.pickle",'rb') as df_pickle:
    df_twin = pickle.load(df_pickle)

del df_pickle


#%% dropna from methylation data
df_m = df_methyl.dropna(thresh=245,axis=1)

del methyl_file

#%% Generate Controls 
def all_comb_fromlist(person): #Get twin list and return non-twin pairs
    control_all = [(i,j) for i in person for j in person if i[:-1]  != j[:-1] ] #239120 組
    control_all = pd.DataFrame(control_all)
    control_all.columns = ['id0','id1']
    return(control_all)

def all_comb(twin_pair): #Get twin list and return all non-twin pairs
    person = list(twin_pair.id0)
    person.extend(list(twin_pair.id1))

    control_all = [(i,j) for i in person for j in person if i[:-1]  != j[:-1] ] #239120 組


    return(pd.DataFrame(control_all))


def shuffle_pair(twin_pair): #Get twin list and return non-twin list that is the same sizs of the twin list
    #In the list, the left and right side of the pairs are fixed
    id0 = list(twin_pair.id0)
    id1 = list(twin_pair.id1)
#    id0,id1 = zip(*twin_pair)
    
    
    x = list(id1[1:])
    x.append(id1[0])
    id1 = x
    
    new_pair = list(zip(id0,id1))
    new_pair = [ [i,j] for (i,j) in new_pair  ]
    
    num_pairs = len(twin_pair)   
    for i in range( num_pairs * 100):
        #右側だけ動かす。 id1のposition AとBを交換
        pos_A  = random.randrange(0,num_pairs)
        pos_B    = random.randrange(0,num_pairs)
        if pos_A == pos_B:
            next
            
        id0_A = new_pair[pos_A][0]
        id1_A = new_pair[pos_A][1]
    
        id0_B = new_pair[pos_B][0]
        id1_B = new_pair[pos_B][1]
        
        if id0_A[:-1] != id1_B[:-1] and id0_B[:-1] != id1_A[:-1] :
            new_pair[pos_A][1],new_pair[pos_B][1] = \
                new_pair[pos_B][1],new_pair[pos_A][1]
    
    new_pair = pd.DataFrame(new_pair)
    new_pair.columns = ['id0','id1']
    return(new_pair)

def shuffle(df): #Get twin data frame and shuffle the sleft and right side of the pair.
    pair = []
    for index,row in df.iterrows():
        if random.random() < 0.5:
            pair.append((row['id0'],row['id1']))
        else:
            pair.append((row['id1'],row['id0']))
            
    id0,id1 = zip(*pair)
    pair = pd.DataFrame({'id0': id0, 'id1': id1})
            
    return(pair)




"""
##実験群
#全双子の組み合わせ #0と1をシャッフルする。
twin_all = shuffle(df_twin)

#性別２種類        
twin_gender0 = shuffle(df_twin[df_twin.GENDER == 0])
twin_gender1 = shuffle(df_twin[df_twin.GENDER == 1])

#年齢２種類
twin_young = shuffle(df_twin[df_twin.AGE <= 52])
twin_old = shuffle(df_twin[df_twin.AGE >= 53])

twin_exp_pair = {"all": twin_all, "gender0": twin_gender0,
            "gender1": twin_gender1, "young": twin_young, "old": twin_old}
del twin_all,twin_gender0,twin_gender1,twin_young,twin_old

for key in twin_exp_pair.keys():
    twin_exp_pair[key] = pd.DataFrame(twin_exp_pair[key],columns=['id0','id1'])

for key in twin_exp_pair.keys():
    twin_exp_pair[key].id0 = twin_exp_pair[key].id0.astype(str)
    twin_exp_pair[key].id1 = twin_exp_pair[key].id1.astype(str)
    
with open("res/twin_exp_pair.pickle",'wb') as df_pickle:
    pickle.dump(twin_exp_pair,df_pickle)

#全組合せ- 双子組合せ
#人間全部
person = list(df_twin.ID0)
person.extend(list(df_twin.ID1))

control_all = [(i,j) for i in person for j in person if i /10  != j /10 ] #239120 組
"""

#save the list of pairs for the reproducibility
with open("res/twin_exp_pair.pickle",'rb') as df_pickle:
    twin_exp_pair = pickle.load(df_pickle)
    


#%% Functions for the statistical test
def set_info(info,df_data,prefix):
    
    info[prefix+'_num_pairs'] = len(df_data.x)
    
    data = df_data.x - df_data.y
    info[prefix+'_average'] = data.mean()
    info[prefix+'_max'] = data.max()
    info[prefix+'_min'] = data.min()
    info[prefix+'_std'] = data.std()  #unbiased standard deviation
    
    info[prefix+'_corr'] = df_data[['x','y']].corr().x.y

    return(info)
    
    

def set_methylation(twin_pair,df_data):
#    df_data.index = data.index.astype('int64')
    df_data.columns = ['x']
    m0  = pd.merge(twin_pair,df_data,left_on='id0',right_index=True,how='inner')
    df_data.columns = ['y']
    m01 = pd.merge(m0,df_data,left_on='id1',right_index=True,how='inner')
    
    m01.x.astype('float32') 
    m01.y.astype('float32') 

    return(m01)    

def set_levens(info,a_set,b_set,prefix):
    
    stat,p = levene(a_set.x- a_set.y, b_set.x - b_set.y,center="mean")
    
    
    
    info["leven_p_"+prefix] = p
    
    return(info)
    

def calc_statistics_nogendermix(methylation,exp_pair,probe_id): 
    # Calculation of genetic factor index
    # As the original data are type(str), we should conver to float.
    info = {}
    info['index'] = probe_id
        
    nul_flag = methylation.isnull()
#    if sum(nul_flag) != 0:
#        print("nul : ", sum(nul_flag) )
    info['num_valid'] = len(methylation) - sum(nul_flag)
    info['num_null'] = sum(nul_flag)
    data = methylation.dropna()
    data = data.astype('float32') #Use to calculate statistical values
    info['average'] = data.mean()
    info['max'] = data.max()
    info['min'] = data.min()
    info['std'] = data.std()     #Unbiased standard deviation

   
    df_data = pd.DataFrame(data)
    
    
    met_data_all = set_methylation(exp_pair['all'],df_data) 
    info = set_info(info,met_data_all,'all_twin')
    
    
    #Generate control pairs. The number of male-pairs and female-pairs in the control are the same as the twin-pairs.
    #We generate 11 controls.
    gender0_data = set_methylation(exp_pair['gender0'],df_data)
    gender1_data = set_methylation(exp_pair['gender1'],df_data)
    for i in range(11):
        lr_shuffle_gender0 = shuffle(gender0_data)
        lr_shuffle_gender1 = shuffle(gender1_data)
        shuffle_gender0 = shuffle_pair(lr_shuffle_gender0)
        shuffle_gender1 = shuffle_pair(lr_shuffle_gender1)
        
        control = pd.concat([shuffle_gender0,shuffle_gender1])
        control = set_methylation(control,df_data)
        
        info = set_info(info,control,'control_nogendermix' + str(i))
        
        info = set_levens(info,met_data_all,control,"all_twin_vs_control_nogendermix" + str(i))
        

    #make list to calculate medians.
    p_list = []
    for i in range(11):
        idx = 'leven_p_all_twin_vs_control_nogendermix' + str(i)
        p_list.append(info[idx])
    p_med = statistics.median(p_list)

    p_flag = False
    for i in range(11):
        idx = 'leven_p_all_twin_vs_control_nogendermix' + str(i)
        if info[idx] == p_med:
            info['leven_p_all_twin_vs_control_nogendermix_median'] = info[idx]
            info['control_nogendermix_std'] = info['control_nogendermix{}_std'.format(i)]
            info['control_nogendermix_corr'] = info['control_nogendermix{}_corr'.format(i)]
            p_flag=True

    if p_flag == False:
        print("Error to found median of p_value \nProbeid {}".format(probe_id))

    return(info)
    

    
#%% verification code
probe_id = "cg00004073" #nul 81個入り
methylation = df_methyl[probe_id]

info_all= calc_statistics_nogendermix(df_methyl[probe_id],twin_exp_pair,probe_id)


#%% Calc stats for each methyl probe
def set_init_dic_list(info_all):
    info_dic_list ={}
    
    for key in info_all.keys():
        info_dic_list[key] = [ info_all[key]]
    return(info_dic_list)
    
    
def append_dic_list(info_dic_list,info):
    for key in info.keys():
        info_dic_list[key].append( info[key])

    return(info_dic_list)


#%% Three Run in parallel

count = 0

start = [0,150000,300000]
end   = [150000,300000,500000]

start= start[cpuid]
end = end[cpuid]



if 'info_dic_list' in locals():
    del info_dic_list


for probe_id in df_methyl:
    if probe_id == "PAIRID":
        continue
    count += 1
    if count <= start:        
        continue
    if count > end:
        break
    

    print(count,probe_id)
    
    df_methyl[probe_id]
    
    info_all = calc_statistics_nogendermix(df_methyl[probe_id],twin_exp_pair,probe_id)

    if 'info_dic_list' not in locals():
        info_dic_list = set_init_dic_list(info_all)
    else:
        append_dic_list(info_dic_list,info_all)


#%% save the result of each run
df_res = pd.DataFrame(info_dic_list,index=info_dic_list['index'])

df_res.to_csv('res/genetic_noGenderMix_res'+ str(cpuid) + '.csv')
with open('res/genetic_noGenderMix_res'+str(cpuid) + '.pickle','wb') as df_pickle:
    pickle.dump(df_res,df_pickle)

#%% Merge the three resuts of each run. 
"""
df_res = []

for i in range(3):
    with open('res/genetic_noGenderMix_res'+str(i) + '.pickle','rb') as df_pickle:
        df = pickle.load(df_pickle)
        
        df_res.append(df)

df_concat = pd.concat([df_res[0],df_res[1],df_res[2]])

df_concat.to_csv('res/genetic_noGenderMix_res.csv')
with open('res/genetic_noGenderMix_res.pickle','wb') as df_pickle:
    pickle.dump(df_concat,df_pickle)
"""

#%% Merge all the result to one data frame.


# read the result of the three run
with open('res/genetic_noGenderMix_res.pickle','rb') as df_pickle:
    df_concat = pickle.load(df_pickle)

# Read the result of gender index and environmental index (this result is not published yet. )
with open('res/envGenetic_res.pickle','rb') as df_pickle:
    df_analysis = pickle.load(df_pickle)
df_analysis = pd.concat([df_analysis[0],df_analysis[1],df_analysis[2]])


col_ana = set(df_analysis.columns)
col_con = set(df_concat.columns)
need_col = list(col_con-col_ana)

df_analysis_res = pd.concat([df_analysis,df_concat[need_col]],join='inner',axis=1)


#Shave the probes that has less than 245 samples (less than half) or whose name does not start with 'ch'
df_nonch = df_analysis_res.query('not index.str.startswith("ch")',engine='python')
flag = df_nonch.num_valid >= 245  # 484279 probe
print(sum(flag))  #In fact, No probe was eliminated.
df_analysis_new = df_nonch[flag]
#number of probes whose name start with 'ch' was 3091,   481190 probes left. (the number of probes in the manuscript)

#Write the dataframe to a file
with open('res/df_analysis_res.pickle','wb') as df_pickle:
    pickle.dump(df_analysis_new,df_pickle)

#How to read the result file
with open('res/df_analysis_res.pickle','rb') as df_pickle:
    df_analysis_new = pickle.load(df_pickle) 






