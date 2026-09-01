'''
1.  重复采样 生成负样本 
2.  基于 最短距离
3.  基于 Diffusion State Distance
'''
import pandas as pd 
import numpy as np 
import random
import networkx as nx
import os 
random.seed(0)

##########################
# Data Load
CPI = pd.read_csv("./Result/Processed Data/Mapped_Data/CPI_LCC.csv")
DPI = pd.read_csv("./Result/Processed Data/Mapped_Data/DPI_LCC.csv")
PPI_LCC = pd.read_csv(".\\Data\\PPI_LCC.csv")

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

def group_group_distance(G_S,G_T,Net,reference_dict): # 
    GG_dis = []
    NN_dis = []
    NN_dict = {}
    for node1 in set(G_T):
        for node2 in set(G_S):
            node_pair = sorted([node1,node2])
            node_pair_str = str(node_pair)
            if node_pair_str in reference_dict.keys():
                distance = reference_dict[node_pair_str]
            else:
                distance = node_node_distance(Net,node1,node2)
                reference_dict[node_pair_str] = distance
            NN_dis.append(distance)
            NN_dict[str(node_pair)] = distance
        NG_value = np.min(NN_dis)
        GG_dis.append(NG_value)

    for node1 in set(G_S):
        for node2 in set(G_T):
            node_pair = sorted([node1,node2])
            node_pair_str = str(node_pair)
            if node_pair_str in reference_dict.keys():
                distance = reference_dict[node_pair_str]
            else:
                distance = node_node_distance(Net,node1,node2)
                reference_dict[node_pair_str] = distance
            NN_dis.append(distance)
            NN_dict[str(node_pair)] = distance
        NG_value = np.min(NN_dis)
        GG_dis.append(NG_value)
    outcome = np.mean(GG_dis)
    return outcome,reference_dict


# Z-score 采样
def Network_construct(net_df):
    '''
    input PPI network data  return a network object used for PPI distance compute and a degree dataframe [node_ID,node_degree]
    '''
    G = nx.Graph()
    G.add_edges_from(np.array(net_df))
    c = max(nx.connected_components(G), key=len)
    G = G.subgraph(c)
    # following codes for sampling
    node_degree_df = pd.DataFrame([list(i) for i in G.degree()],columns = ['node','degree'])
    node_degree_df = node_degree_df.sort_values('degree').reset_index(drop = True)
    return G, node_degree_df

def Get_Degree_Details(node_degree_df):
    '''
    input node_degree_df return Counter Result [degree,num,nodes]
    return bins [according to nodes' degree, divide nodes into different sets under a set size restriction]
    '''
    df = node_degree_df
    degree_list = list(set(node_degree_df['degree'].tolist()))
    df_G = df.groupby('degree')
    degree_summary_df = []
    for i in degree_list:
        selected_df = df_G.get_group(i)
        selected_nodes = selected_df['node'].tolist()
        nodes_num = len(selected_nodes)
        degree_summary_df.append([i,nodes_num,str(selected_nodes)])
    degree_summary_df  = pd.DataFrame(degree_summary_df,columns= ['degree','num','nodes'])
    degree_summary_df = degree_summary_df.sort_values('degree').reset_index(drop = True)
    return degree_summary_df

def Get_Degree_Bins(degree_summary_df,min_bin_size):
    '''
    input degree_summary_df,min_bin_size     return: degree bins
    '''
    bins = []
    i = 0
    nodes_list = degree_summary_df['nodes'].tolist()
    degree_list = degree_summary_df['degree'].tolist()
    while i < len(degree_list):
        selected_nodes = eval(nodes_list[i]) # nodes list
        low = degree_list[i] # corresponding degree
        while len(selected_nodes) < min_bin_size:
            i += 1
            if i == len(degree_list):
                break
            add_selected_nodes = eval(nodes_list[i]) # nodes list
            selected_nodes += add_selected_nodes
        if i == len(degree_list):
            i -= 1
        high = degree_list[i]
        i += 1
        if len(selected_nodes) < min_bin_size:
            _low,_high,_node_list = bins[-1]
            bins[-1] = [_low,high,_node_list + selected_nodes]
        else:
            bins.append((low,high,selected_nodes))
    bins = bins
    return bins

def Random_Select_Node(G_node,node_degree_df,bins):
    '''
    for each ICD target protein, random selected porteins with similar degree in PPI network
    input: node_degree_df, G_node, bins  return: random select node
    '''
    select_nodes = []
    G_node_df = node_degree_df[node_degree_df['node'].isin(G_node)] # get node degree
    for l, h, binnode_list in bins:
        G_node_df_select = G_node_df[(G_node_df['degree'] >= l) & (G_node_df['degree'] <= h)].reset_index(drop = True)
        G_node_set = G_node_df_select['node'].tolist()
        count = len(G_node_set)
        left = list(set(binnode_list) - set(G_node_set))
        if len(left) < count:
            select_nodes += random.sample(binnode_list, count)
        else:
            select_nodes += random.sample(left, count)
    return set(select_nodes)

def Random_Select_Result(G_S,G_T,node_degree_df,bins):
    G_S = Random_Select_Node(G_S,node_degree_df,bins)
    G_T = Random_Select_Node(G_T,node_degree_df,bins)
    return G_S,G_T

min_bin_size = 100
repeat_num = 1000
G, node_degree_df = Network_construct(PPI_LCC[['proteinA_entrezid','proteinB_entrezid']])
degree_summary_df = Get_Degree_Details(node_degree_df)
bins = Get_Degree_Bins(degree_summary_df,min_bin_size)
#### repeat test for positive sample
random.seed(0)
DPI['EntrezID LCC'] = [[int(j) for j in eval(i)] for i in DPI['EntrezID LCC']]
CPI['EntrezID LCC'] = [[int(j) for j in eval(i)] for i in CPI['EntrezID LCC']] 

DPI_Dict = {i:j for i,j in zip(DPI['Disease Name'],DPI['EntrezID LCC'])}
CPI_Dict = {i:j for i,j in zip(CPI['Drug Name'],CPI['EntrezID LCC'])}
pairs = [[i,j] for i in DPI['Disease Name'] for j in CPI['Drug Name']]
len(pairs)

for pair in pairs:
    repeat_ls = []
    disease = pair[0]
    drug = pair[1]

    disease_protein = DPI_Dict[disease]
    drug_protein = CPI_Dict[drug]

    epoch_dis_ls = []
    for epoch in range(repeat_num):
        G_drug,G_disease= Random_Select_Result(drug_protein,disease_protein,node_degree_df,bins)
        repeat_ls.append([str(pair),epoch,str(sorted(list(G_drug))),str(sorted(list(G_disease)))])
    repeat_ls = pd.DataFrame(repeat_ls,columns=['edge','epoch','drug_protein','disease_protein'])
    repeat_ls.to_csv(f"./Result/Processed Data/PPI Distance Compute/Negative_Sample/{pair}_repeat.csv",index = False)
    print(pair)
