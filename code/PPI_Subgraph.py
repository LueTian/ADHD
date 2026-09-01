###############################################
# Core 
###############################################
import pandas as pd 
import numpy as np 
import os 
###########################################################################################################################
# Function
def List2DF(input,color,type,size):
    outcome = pd.DataFrame(input,columns=['ID'])
    outcome['color'] = color
    outcome['type']  = type
    outcome['size']  = size 
    return outcome

def Overlap_Protein(input):
    S1 = eval(input[0])
    S2 = eval(input[1])
    outcome = sorted(list(set(S1).intersection(set(S2))))
    return str(outcome)

def New_Protein(input):
    S1 = eval(input[0])
    S2 = eval(input[1])
    outcome = sorted(list(set(S2)-set(S1)))
    return str(outcome)

###########################################################################################################################
# Data
Pairs = pd.read_csv("./PPI_Outcome/Pairs_True_v2.csv") 
Protein_Dict = pd.read_csv(".\GO_KEGG\Original\Protein_Map_Entrez_Symbol.csv")
Protein_Dict = {i:j for i,j in zip(Protein_Dict['Entrez ID'],Protein_Dict['Symbol'])}

Pairs['Overlap Protein'] = Pairs[['Disease Protein','Drug Protein']].apply(Overlap_Protein,axis = 1)
Pairs['Drug Protein New'] = Pairs[['Overlap Protein','Drug Protein']].apply(New_Protein,axis = 1)
Pairs['Disease Protein New'] = Pairs[['Overlap Protein','Disease Protein']].apply(New_Protein,axis = 1)

PPI_net = pd.read_csv("./Data/PPI_LCC.csv")
PPI_net.columns = ['source','target','database']


for idx in Pairs.index:
    diesase_pair = Pairs.loc[idx,'pairs']
    protein_disease = eval(Pairs.loc[idx,'Disease Protein New'])
    protein_drug = eval(Pairs.loc[idx,'Drug Protein New'])
    protein_disease_drug = eval(Pairs.loc[idx,'Overlap Protein'])
    protein_disease_old = eval(Pairs.loc[idx,'Disease Protein'])
    protein_drug_old = eval(Pairs.loc[idx,'Drug Protein'])
    protein_ls = set(protein_disease + protein_drug + protein_disease_drug)
    PPI_part = PPI_net[(PPI_net['source'].isin(protein_ls))&(PPI_net['target'].isin(protein_ls))].reset_index(drop = True)
    
    # 设置蛋白质节点属性
    node_disease = List2DF(protein_disease,'#007CD3','Disease Protein',10)
    node_drug = List2DF(protein_drug,'#C05683','Drug Protein',10)
    node_disease_drug = List2DF(protein_disease_drug,'#eba834','Disease Drug Overlapped Protein',10)

    # 设置药物，疾病节点属性
    disease = eval(diesase_pair)[0]
    drug = eval(diesase_pair)[1]

    node_add = [[disease,'#005983','Disease',50],[drug,'#C25E5E','Drug',50]]
    node_add = pd.DataFrame(node_add,columns=['ID','color','type','size'])
    node_add['extension'] = [disease,drug]
    node_add['polygon'] = 6

    # 增加节点的属性 Symbol
    node = pd.concat([node_disease,node_drug,node_disease_drug],axis = 0).reset_index(drop = True)
    node['extension'] = [Protein_Dict[i] for i in node['ID']]
    node['polygon'] = 1
    node = pd.concat([node_add,node],axis=0).reset_index(drop = True)

    # 边图表
    PPI_part2 = [[disease,i] for i in protein_disease_old]
    PPI_part3 = [[drug,i] for i in protein_drug_old]
    PPI_part2 = pd.DataFrame(PPI_part2,columns=['source','target'])
    PPI_part3 = pd.DataFrame(PPI_part3,columns=['source','target'])
    edge = pd.concat([PPI_part,PPI_part2,PPI_part3],axis =0).reset_index(drop = True)
    edge = edge[edge['source']!=edge['target']].reset_index(drop = True)
    path = f'.\\PPI_Subgraph\\{diesase_pair}_subgraph'
    if not os.path.exists(path):
        os.makedirs(path)
    edge.to_csv(os.path.join(path,'net_edge.csv'),index = False)
    node.to_csv(os.path.join(path,'net_node.csv'),index = False)