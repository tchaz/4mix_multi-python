#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

# CC BY-NC-SA original code by Eurogenes DESEUK1, modified by AJ14, python translation by tchaz


import pandas as pd
import numpy as np
from numpy import sqrt
from scipy.spatial import KDTree, cKDTree

working_dir = './'
   
# Get all admixture proportions of 4-way models. 
# Let's just calculate them out of 100 (ie percentages) as well for preformance
    
   
props = []
props_percent = []
        
for i in range(101):
    for j in range(101-i):
        for k in range(101-i-j):
            # props.append([i/100,j/100,k/100,1-((i+j+k)/100)])
            props_percent.append([i,j,k,100 - i -j -k])
            
props_percent = np.array(props_percent)
props = props_percent/100
   


def get_mix(working_dir = 'M:/4mix/', input_file = working_dir + 'K8avg.csv', target_file = working_dir + 'target.csv', pop1 = 'Bashkir', pop2 = 'PPNB', pop3 = 'Central_Greek', pop4 = 'Corded_Ware_LN', result_file_prefix = working_dir + 'results_'):

    # Load admixture and target values from files

    admixture_df = pd.read_csv(input_file)
    target_df = pd.read_csv(target_file)
 
    # Look for population names in input file

    pop_names = [pop1, pop2, pop3, pop4]
    pop_names.sort()
    pop_diff = [n for n in pop_names if n not in set(admixture_df.ix[:,0])]

    # Test if popdiffs is null character set; if not, error

    if len(pop_diff) == 0:
        # Get the relevant rows from admixture_df
        source_pops = admixture_df[admixture_df['Unnamed: 0'].isin(pop_names)]
    else:
        print("Error: The following population names are not in the input file:\n" + \
                         ' '.join(pop_diff) +"\n Correct the above population names and rerun.")

    # Perform matrix multiplication of props and sourcepops
       
    source_pops_matrix = [r[1:] for r in source_pops.values]
    target_matrix = [r[1:] for r in target_df.values]
    
    models = np.dot(props,source_pops_matrix)
    
    # for each row in the target_matrix find the nearest neighbour amongst the rows of the models matrix using cKDTree
    
    t = cKDTree(models)
    results = t.query(target_matrix, 1, 0, 2, 1)
    
    # make the results table
    
    result_table = [["Population", pop_names[0], pop_names[1], pop_names[2], pop_names[3], "D statistic"]]
    for i in range(len(results[0])):
        this_result = [target_df['Unnamed: 0'][i]] + [props_percent[results[1][i]][j] for j in range(4)] + ['%.6f' % results[0][i]]
        result_table.append(this_result)
    
    # write results_file
    
    headers = result_table.pop(0) 
    df = pd.DataFrame(result_table, columns=headers)
    result_file_suffix = "-".join(pop_names)
    df.to_csv(result_file_prefix + result_file_suffix + '.csv')
