'''
分析 KEGG GO overlapped 部分 
'''
import pandas as pd
import numpy as np 
import os 

path1 = ".//GO_KEGG//Result//Disease"
path2 = ".//GO_KEGG//Result//Drug"

def Folder_Analysis(Path):
    '''
    收集文件夹内的数据 【edge,GO,KEGG】
    '''
    KEGG_path = Path + "\\KEGG"
    GO_path =  Path + "\\GO"
    #### Folder KEGG
    KEGG_List_DF = pd.DataFrame()
    Edge_KEGG_DF = []
    KEGG_file = os.listdir(KEGG_path)
    for KEGG_name in KEGG_file:
        KEGG_df = pd.read_csv(os.path.join(KEGG_path,KEGG_name))
        if KEGG_df.shape[0] == 0:
            continue
        KEGG_df = KEGG_df[KEGG_df['p.adjust']<0.05].reset_index(drop = True)
        KEGG_df['bg.gene'] = KEGG_df['BgRatio'].str.split('/').str.get(0).astype('int')
        KEGG_df = KEGG_df[(KEGG_df['bg.gene']>20)&(KEGG_df['bg.gene']<200)].reset_index(drop = True)
        Edge_KEGG_DF += [[KEGG_name,i] for i in KEGG_df['ID']]
        ### KEGG_detail_dataframe
        KEGG_List_DF = pd.concat([KEGG_List_DF,KEGG_df[['category','subcategory','ID','Description']]],axis = 0).drop_duplicates().reset_index(drop=True)

    #### Folder KEGG
    GO_List_DF = pd.DataFrame()
    Edge_GO_DF = []
    GO_file = os.listdir(GO_path)
    for GO_name in GO_file:
        GO_df = pd.read_csv(os.path.join(GO_path,GO_name))
        if GO_df.shape[0] == 0:
            continue
        GO_df = GO_df[GO_df['p.adjust']<0.05].reset_index(drop = True)
        GO_df['bg.gene'] = GO_df['BgRatio'].str.split('/').str.get(0).astype('int')
        # GO_df = GO_df[(GO_df['bg.gene']>20)&(GO_df['bg.gene']<200)].reset_index(drop = True)
        Edge_GO_DF += [[GO_name,i] for i in GO_df['ID']]
        ### GO_detail_dataframe
        GO_List_DF = pd.concat([GO_List_DF,GO_df[['ONTOLOGY','ID','Description']]],axis = 0).drop_duplicates().reset_index(drop=True)

    Edge_GO_DF = pd.DataFrame(Edge_GO_DF,columns= ['Edge','GO'])
    Edge_KEGG_DF = pd.DataFrame(Edge_KEGG_DF,columns= ['Edge','KEGG'])

    return Edge_GO_DF,Edge_KEGG_DF,GO_List_DF,KEGG_List_DF

Edge_GO_DF_Disease,Edge_KEGG_DF_Disease,GO_List_DF_Disease,KEGG_List_DF_Disease = Folder_Analysis(path1)
Edge_GO_DF_Drug,Edge_KEGG_DF_Drug,GO_List_DF_Drug,KEGG_List_DF_Drug = Folder_Analysis(path2)

GO_List = pd.concat([GO_List_DF_Disease,GO_List_DF_Drug],axis = 0).drop_duplicates().reset_index(drop = True)
KEGG_List = pd.concat([KEGG_List_DF_Disease,KEGG_List_DF_Drug],axis = 0).drop_duplicates().reset_index(drop =True)

Edge_GO_DF_Disease['Edge'] = Edge_GO_DF_Disease['Edge'].str.split('.').str.get(0).str.split('_').str.get(1)
Edge_KEGG_DF_Disease['Edge'] = Edge_KEGG_DF_Disease['Edge'].str.split('.').str.get(0).str.split('_').str.get(1)     

Edge_GO_DF_Drug['Edge'] = Edge_GO_DF_Drug['Edge'].str.split('.').str.get(0).str.split('_').str.get(1)
Edge_KEGG_DF_Drug['Edge'] = Edge_KEGG_DF_Drug['Edge'].str.split('.').str.get(0).str.split('_').str.get(1)     

Edge_GO_DF_Disease['Ref'] = Edge_GO_DF_Disease['Edge'].str.cat(Edge_GO_DF_Disease['GO'],sep = '_')
Edge_GO_DF_Drug['Ref'] = Edge_GO_DF_Drug['Edge'].str.cat(Edge_GO_DF_Drug['GO'],sep = '_')
Edge_KEGG_DF_Disease['Ref'] = Edge_KEGG_DF_Disease['Edge'].str.cat(Edge_KEGG_DF_Disease['KEGG'],sep = '_')
Edge_KEGG_DF_Drug['Ref'] = Edge_KEGG_DF_Drug['Edge'].str.cat(Edge_KEGG_DF_Drug['KEGG'],sep = '_')

def Overlap(DF1,DF2,str):
    outcome = []
    item1_ls = set(DF1['Edge'].tolist())
    item2_ls = set(DF2['Edge'].tolist())
    DF1_G = DF1.groupby('Edge')
    DF2_G = DF2.groupby('Edge')
    for item1 in item1_ls:
        for item2 in item2_ls:
            DF1_p = DF1_G.get_group(item1)
            DF2_p = DF2_G.get_group(item2)
            overlapped = set(DF1_p[str]).intersection(set(DF2_p[str]))
        outcome += [[item1,item2,i] for i in overlapped]
    outcome = pd.DataFrame(outcome,columns=['Disease','Drug',str])
    return outcome

KEGG_overlapped = Overlap(Edge_KEGG_DF_Disease,Edge_KEGG_DF_Drug,'KEGG')
GO_overlapped = Overlap(Edge_GO_DF_Disease,Edge_GO_DF_Drug,'GO')
KEGG_overlapped['Ref'] = KEGG_overlapped['Disease'].str.cat(KEGG_overlapped['Drug'],sep = '_')
GO_overlapped['Ref'] = GO_overlapped['Disease'].str.cat(GO_overlapped['Drug'],sep = '_')
KEGG_overlapped.columns = ['Disease', 'Drug', 'ID','Ref']
GO_overlapped.columns = ['Disease', 'Drug', 'ID', 'Ref']

KEGG_overlapped = pd.merge(KEGG_overlapped,KEGG_List,on = 'ID',how = 'left')
GO_overlapped = pd.merge(GO_overlapped,GO_List,on = 'ID',how = 'left')

KEGG_count = KEGG_overlapped.value_counts('Ref').reset_index(drop = False)
GO_count   = GO_overlapped.value_counts('Ref').reset_index(drop = False)

GO_count.to_csv(".//GO_KEGG_Analysis//GO_count.csv",index = False)
KEGG_count.to_csv(".//GO_KEGG_Analysis//KEGG_count.csv",index = False)
KEGG_overlapped.to_csv(".//GO_KEGG_Analysis//KEGG_overlapped.csv",index = False)
GO_overlapped.to_csv(".//GO_KEGG_Analysis//GO_overlapped.csv",index = False)