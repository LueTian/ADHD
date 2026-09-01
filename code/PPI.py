############################################
### Core
############################################
'''
计算 PPI distance
'''
import pandas as pd
import numpy as np
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import os 
import random
import networkx as nx
#### Function 
# 定义函数 [计算距离]

def node_node_distance(G,node1,node2):
    '''
    :param G:
    :param node1: 
    :param node2:  
    :return: the length of the shortest path of node1 and node2
    '''
    outcome = nx.shortest_path_length(G,source = node1,target=node2)
    return outcome

class PPI_Compute():
    def __init__(self,G) -> None:
        self.G = G
        self.ref_dict = {}

    def PPI_distance(self,S):
        S1 = S[0]
        S2 = S[1]

        G = self.G
        ref_dict = self.ref_dict

        min_items_s2 = [] 
        # 对 S1 中每一个的 protein 到 S2 最近的距离，再对 S1 中每一个 点到集合的距离取均值
        for item_s in S1:
            min_items_itemt = []
            for item_t in S2:
                pair = str(sorted([item_s,item_t]))
                if pair in ref_dict.keys():
                    min_items_itemt.append(ref_dict[pair])
                else:
                    dis = node_node_distance(G,item_s,item_t)
                    min_items_itemt.append(dis)
                    ref_dict[pair] = dis
            min_items_s2.append(np.min(min_items_itemt))
        outcome = np.mean(min_items_s2)
        self.ref_dict = ref_dict
        return outcome
    
    def PPI_distance2(self,S):
        S1 = S[0]
        S2 = S[1]

        G = self.G
        ref_dict = self.ref_dict

        min_items_itemt = []
        for item_s in S1:
            for item_t in S2:
                pair = str(sorted([item_s,item_t]))
                if pair in ref_dict.keys():
                    min_items_itemt.append(ref_dict[pair])
                else:
                    dis = node_node_distance(G,item_s,item_t)
                    min_items_itemt.append(dis)
                    ref_dict[pair] = dis
        outcome = np.mean(min_items_itemt)
        self.ref_dict = ref_dict 
        return outcome 

CPI = pd.read_csv("./Result/Processed Data/Mapped_Data/CPI_LCC.csv")
DPI = pd.read_csv("./Result/Processed Data/Mapped_Data/DPI_LCC.csv")
DPI['EntrezID LCC'] = [[int(j) for j in eval(i)] for i in DPI['EntrezID LCC']]
CPI['EntrezID LCC'] = [[int(j) for j in eval(i)] for i in CPI['EntrezID LCC']] 
DPI_Dict = {i:j for i,j in zip(DPI['Disease Name'],DPI['EntrezID LCC'])}
CPI_Dict = {i:j for i,j in zip(CPI['Drug Name'],CPI['EntrezID LCC'])}
pairs = [[i,j] for i in DPI['Disease Name'] for j in CPI['Drug Name']]
Pairs = pd.DataFrame(pairs,columns=['Disease','Drug'])
Pairs['pairs'] = [str(i) for i in pairs]
Pairs['Disease Protein'] = [str(DPI_Dict[i]) for i in Pairs['Disease']]
Pairs['Drug Protein'] = [str(CPI_Dict[i]) for i in Pairs['Drug']]
Pairs['Disease Protein'] = [eval(i) for i in Pairs['Disease Protein']]
Pairs['Drug Protein'] = [eval(i) for i in Pairs['Drug Protein']]

PPI = pd.read_csv("./Data/PPI_LCC.csv")
G = nx.Graph()
G.add_edges_from(np.array(PPI[['proteinA_entrezid','proteinB_entrezid']]))

print("Positive Sample PPI distance computate")
Computer = PPI_Compute(G)
Pairs['PPI_Dis1'] = Pairs[['Disease Protein','Drug Protein']].apply(Computer.PPI_distance,axis =1)
Pairs['PPI_Dis2'] = Pairs[['Drug Protein','Disease Protein']].apply(Computer.PPI_distance,axis =1)
Pairs['PPI_Dis3'] = Pairs[['Disease Protein','Drug Protein']].apply(Computer.PPI_distance2,axis =1)
Pairs.to_csv("./PPI_Outcome/Pairs_True_v2.csv",index = False)
print("Negative Sample PPI distance computate")

file_path = "./Negative_Sample"
save_path = "./Negative_Outcome"
file_list = os.listdir(file_path)
for file_name in file_list:
    print(file_name)
    negative_pairs = pd.read_csv(os.path.join(file_path,file_name))
    negative_pairs.columns = ['edge', 'epoch', 'Drug Protein', 'Disease Protein']
    negative_pairs['Drug Protein'] = [eval(i) for i in negative_pairs['Drug Protein']]
    negative_pairs['Disease Protein'] = [eval(i) for i in negative_pairs['Disease Protein']]

    negative_pairs['PPI_Dis1'] = negative_pairs[['Disease Protein','Drug Protein']].apply(Computer.PPI_distance,axis =1)
    negative_pairs['PPI_Dis2'] = negative_pairs[['Drug Protein','Disease Protein']].apply(Computer.PPI_distance,axis =1)
    negative_pairs['PPI_Dis3'] = negative_pairs[['Disease Protein','Drug Protein']].apply(Computer.PPI_distance2,axis =1)
    negative_pairs.to_csv(os.path.join(save_path,file_name),index = False)
    
###############################################################################################################################
Pairs = pd.read_csv("./PPI_Outcome/Pairs_True_v2.csv") 
# 提取bg mean and bg std 
file_path = "./Negative_Outcome_v2"
file_names = os.listdir(file_path)
bg_dict = {}
for file_name in file_names:
    file_DF = pd.read_csv(os.path.join(file_path,file_name))
    edge_name = file_name.split('_')[0]
    PPI_Dis1_mean = file_DF['PPI_Dis1'].mean()
    PPI_Dis1_std  = file_DF['PPI_Dis1'].std()
    PPI_Dis2_mean = file_DF['PPI_Dis2'].mean()
    PPI_Dis2_std  = file_DF['PPI_Dis2'].std()
    PPI_Dis3_mean = file_DF['PPI_Dis3'].mean()
    PPI_Dis3_std  = file_DF['PPI_Dis3'].std()
    bg_dict[edge_name] = {'PPI_Dis1':{'mean':PPI_Dis1_mean,'std':PPI_Dis1_std,'value ls':file_DF['PPI_Dis1'].tolist()},
                          'PPI_Dis2':{'mean':PPI_Dis2_mean,'std':PPI_Dis2_std,'value ls':file_DF['PPI_Dis2'].tolist()},
                          'PPI_Dis3':{'mean':PPI_Dis3_mean,'std':PPI_Dis3_std,'value ls':file_DF['PPI_Dis3'].tolist()}}



def Z_score_compute(input,bg_dict,str):
    pair  = input[0]
    value = input[1]
    bg_mean = bg_dict[pair][str]['mean']
    bg_std  = bg_dict[pair][str]['std']
    z_score = (value-bg_mean)/bg_std
    return z_score

def P_value(input,bg_dict,str):
    pair = input[0]
    value = input[1]
    bg_ls = bg_dict[pair][str]['value ls']
    p_value = np.sum([1 if i<value else 0 for i in bg_ls])/len(bg_ls)
    return p_value



Pairs['PPI_Dis1_Z_score'] = Pairs[['pairs','PPI_Dis1']].apply(Z_score_compute,bg_dict=bg_dict,str ='PPI_Dis1',axis =1)
Pairs['PPI_Dis2_Z_score'] = Pairs[['pairs','PPI_Dis2']].apply(Z_score_compute,bg_dict=bg_dict,str ='PPI_Dis2',axis =1)
Pairs['PPI_Dis3_Z_score'] = Pairs[['pairs','PPI_Dis3']].apply(Z_score_compute,bg_dict=bg_dict,str ='PPI_Dis3',axis =1)
Pairs.to_csv("./PPI_Outcome/Pairs_True_bg_v2.csv",index = False)

# Pairs[Pairs['PPI_Dis1_Z_score']<-1.645].shape # 4
# Pairs[Pairs['PPI_Dis2_Z_score']<-1.645].shape # 21
# Pairs[Pairs['PPI_Dis3_Z_score']<-1.645].shape # 8

# Pairs['PPI_Dis1_p_value'] = [norm.cdf(i, loc=0, scale=1) for i in Pairs['PPI_Dis1_Z_score']]
# Pairs['PPI_Dis2_p_value'] = [norm.cdf(i, loc=0, scale=1) for i in Pairs['PPI_Dis2_Z_score']]
# Pairs['PPI_Dis3_p_value'] = [norm.cdf(i, loc=0, scale=1) for i in Pairs['PPI_Dis3_Z_score']]
# Pairs[(Pairs['PPI_Dis1_p_value']<0.05)&(Pairs['PPI_Dis2_p_value']<0.05)&(Pairs['PPI_Dis3_p_value']<0.05)]

Pairs = pd.read_csv("./PPI_Outcome/Pairs_True_bg_v2.csv")
Pairs = Pairs[Pairs['Disease']!= 'Any CVD'].reset_index(drop = True)
Pairs['p value'] = Pairs[['pairs','PPI_Dis2']].apply(P_value,bg_dict=bg_dict,str ='PPI_Dis2',axis =1)

# # adjust p value
# Pairs['PPI_Dis1_p_value_adjust'] = multipletests(Pairs['PPI_Dis1_p_value'].tolist(), method='bonferroni')[1]
# Pairs['PPI_Dis2_p_value_adjust'] = multipletests(Pairs['p value'].tolist(), method='bonferroni')[1]
# Pairs['PPI_Dis3_p_value_adjust'] = multipletests(Pairs['PPI_Dis3_p_value'].tolist(), method='bonferroni')[1]

Pairs['p value_adjust'] = multipletests(Pairs['p value'].tolist(), method='bonferroni')[1]
Pairs.to_csv("./PPI_Outcome/Pairs_True_bg_v3.csv",index = False)

Pairs[Pairs['PPI_Dis1_p_value_adjust'] < 0.05] # 1 
Pairs[Pairs['PPI_Dis2_p_value_adjust'] < 0.05] # 13
Pairs[Pairs['PPI_Dis3_p_value_adjust'] < 0.05] # 5
Pairs[Pairs['p value_adjust']<0.05]
# [Hypertension,Methylphenidate] & [Hypertension,Amphetamine]