'''
生成 GO KEGG 分析数据 
'''
import pandas as pd 
import numpy as np 
import os 
import mygene 
#############################################################################################
# Function 
def entrez_id_to_gene_symbol(entrez_id):
    mg = mygene.MyGeneInfo()
    result = mg.getgene(entrez_id, fields='symbol')
    gene_symbol = result.get('symbol')
    return gene_symbol

def List2DF(list,col):
    DF = pd.DataFrame(list,columns=col)
    return DF
#############################################################################################
# Data
Pairs = pd.read_csv("./PPI_Outcome/Pairs_True_v2.csv") 
Protein_list = sorted(set([str(j) for i in Pairs['Disease Protein'] for j in eval(i)]+ [str(j) for i in Pairs['Drug Protein'] for j in eval(i)]))
len(Protein_list) # 1451
Protein_df = pd.DataFrame(Protein_list,columns=['Entrez ID'])
Protein_df['Symbol'] = [entrez_id_to_gene_symbol(i) for i in Protein_df['Entrez ID']]
Protein_df.to_csv(".\GO_KEGG\Original\Protein_Map_Entrez_Symbol.csv",index = False)

path1 = "./GO_KEGG/Original/Disease"
path2 = "./GO_KEGG/Original/Drug"

Pairs_Disease = Pairs[['Disease','Disease Protein']].drop_duplicates().reset_index(drop=True)
Pairs_Drug    = Pairs[['Drug','Drug Protein']].drop_duplicates().reset_index(drop = True)

for idx in Pairs_Disease.index:
    name_ = Pairs_Disease.loc[idx,'Disease']
    protein_ = [[name_,i] for i in eval(Pairs_Disease.loc[idx,'Disease Protein'])]
    protein_df = List2DF(protein_,['ID','EntrezID'])
    protein_df.to_csv(f"./GO_KEGG/Original/Disease/{name_}.csv",index = False)

for idx in Pairs_Drug.index:
    name_ = Pairs_Drug.loc[idx,'Drug']
    protein_ = [[name_,i] for i in eval(Pairs_Drug.loc[idx,'Drug Protein'])]
    protein_df = List2DF(protein_,['ID','EntrezID'])
    protein_df.to_csv(f"./GO_KEGG/Original/Drug/{name_}.csv",index = False)





