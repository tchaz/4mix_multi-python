#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

# CC BY-NC-SA original code by Eurogenes DESEUK1, modified by AJ14, python translation, adaptation by tchaz


import pandas as pd
import numpy as np
import time
import itertools
from numpy import sqrt
from scipy.spatial import KDTree, cKDTree

import base64
import json

working_dir = 'E:/genetic anthropology/'

# working directory name should end with a /
if type(working_dir) != str:
    print('working_dir must be a string, please try again. Use:\n > working_dir = \'C:/my_directory/\'  \n replacing \'C:/my_directory/\' by the path you want')
else:
    if working_dir[-1] != '/':
        working_dir = working_dir +'/'
    else:
        pass



class NumpyAwareJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            if obj.ndim == 1:
                return obj.tolist()
            else:
                return [self.default(obj[i]) for i in range(obj.shape[0])]
        return json.JSONEncoder.default(self, obj)


def combs(c, r):
    """
    Return successive r-length combinations of elements in [0,c).
    Should produce the same output as array(list(combinations(a, r))), but 
    faster. Shouldn't it be faster with n rather than a? Doesn't seem to be.
    """
    a = np.arange(c)
    dt = np.dtype([('', a.dtype)]*r)
    b = np.fromiter(itertools.combinations(a, r), dt)
    return b.view(a.dtype).reshape(-1, r)    


def make_props_fast(n = 4, k=100):
    """
    Makes lists of all possible proportions of n many mixtures with granularity 100 / k - 
    So when n = 2 and k = 100 we list [0,100], [1,99], [2,98],...
    and when n=10 and k = 20 we list [0,...,100], [0,..,5,95], ...
    """
    # crashing memory for n = 6, k = 100 - n = 6, k = 90 takes 6Gb of mem and ~ 5 minutes to run
    # looks like a ~15 times speed up over original make_props for n = 5, k = 100.
    if n < 1 or n > k + 1:
        print('n must be between 1 and ' + str(k) + ' - please try again')
        return
    else:
        c = combs(k+n-1,n-1)
        m = np.zeros(c.shape)
        h = np.ones(c.shape[0])[:, np.newaxis] * k
        m[:,:] = np.arange(c.shape[1])
        r = c - m
        props = np.zeros((c.shape[0],n))
        props[:,0][:, np.newaxis] = r[:,0][:, np.newaxis]
        for i in range(n-2):
            props[:,1+i][:, np.newaxis] = r[:,i+1][:, np.newaxis] - r[:,i][:, np.newaxis]
        props[:,n-1][:, np.newaxis] = h - r[:,n-2][:, np.newaxis]
        return props
        
def get_mix(pop_names = ['Azeri_Dagestan', 'Baalberge_MN', 'Balkar', 'Bashkir', 'Basque_French', 'Bedouin','Bell_Beaker_LN'], input_file = working_dir + 'K8avg.csv', target_file = working_dir + 'target.csv', result_file_prefix = working_dir + 'results_'):
    """
    The function that searches for appropriate weightings of the sources populations 
    given in the pop_names list of names of populations in the input file
    to model the various populations in the target_file list
    """
   
    # Load admixture and target values from files

    admixture_df = pd.read_csv(input_file)
    target_df = pd.read_csv(target_file)
 
    # Look for population names in input file
    
    n = len(pop_names)
    
    if n <= 5:
        k = 100
    elif n == 6:
        k = 80
    elif n == 7:
        k = 50
    elif n == 8:
        k = 33
    elif 9 <= n <= 10:
        print('Warning: more than 8 source populations makes the estimates crude (as well as slow)')
        k = 20
    else:
        print('Warning: more than 10 source populations is problematic - please try with fewer')
        
    props_percent = make_props_fast(n, k)
    props = props_percent / k
    
    if n != len(props_percent[1]):
        print("Error: different number of population names from number of propostions")
    else:
        pass

    pop_names.sort()
    pop_diff = [na for na in pop_names if na not in set(admixture_df.ix[:,0])]

    # Test if popdiffs is null character set; if not, error

    if len(pop_diff) == 0:
        # Get the relevant rows from admixture_df
        source_pops = admixture_df[admixture_df['Unnamed: 0'].isin(pop_names)]
    else:
        print("Error: The following population names are not in the input file:\n" + \
                         ' '.join(pop_diff) +"\n Correct the above population names and rerun.")
        return

    # Perform matrix multiplication of props and sourcepops
       
    source_pops_matrix =  np.array([list(r[1:]) for r in source_pops.values])
    target_matrix = [list(r[1:]) for r in target_df.values]

    models =  np.einsum('ij,jk->ik', props, source_pops_matrix)
    
    
    # for each row in the target_matrix find the nearest neighbour amongst the rows of the models matrix using cKDTree
    
    t = cKDTree(models)
    results = t.query(target_matrix, 1, 0, 2, 1)
    
    # make the results table
    
    result_table = [["Population"] + pop_names + ["D statistic"]]
    for i in range(len(results[0])):
        this_result = [target_df['Unnamed: 0'][i]] + ['%.2f' % (props_percent[results[1][i]][j] * (100/k)) for j in range(n)] + ['%.6f' % results[0][i]]
        result_table.append(this_result)
    
    # write results_file
    
    headers = result_table.pop(0) 
    df = pd.DataFrame(result_table, columns=headers)
    result_file_suffix = "-".join(pop_names)
    df.to_csv(result_file_prefix + result_file_suffix + '.csv')
