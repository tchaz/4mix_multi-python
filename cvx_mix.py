import numpy as np
import cvxpy as cvx
import pandas as pd
from multiprocessing import Pool
import time

working_dir = 'E:/genetic anthropology/4mix/' # 'M:/4mix/'

def cvx_mix(data_file = working_dir + 'source.csv', target_file = working_dir + 'target.csv', result_file_prefix = working_dir + 'results_', dec_pl = 3, distances = True, norm = 2, solver='ECOS', verbose = False):
    
    p = norm 
    dec_places = '%.' + str(dec_pl) + 'f'
    results_table = []

    source_df = pd.read_csv(data_file)    
    target_df = pd.read_csv(target_file)
    
    # A will be the transpose of the matrix given by the source populations; pop_names is the names of the sources; 
    A = np.transpose(source_df.values[:,1:])
    pop_names = [r[0] for r in source_df.values]
    num_src_pops = A.shape[1]
    
    x = cvx.Variable(A.shape[1])
    
    # T is the matrix given by the target populations.
    T = target_df.values[:,1:] 
     
    # For the transpose of each row b of T we minimize || Ax - b ||_p  (using the L_p norm; this defaults to Euclidean: p = 2),
    # subject to x > 0 and || x ||_1 = 1
    
    for tr in range(T.shape[0]):
        
        b = T[tr,:]

        objective = cvx.Minimize(cvx.sum_entries(cvx.power((A*x - b),p)))
        constraints = [0 <= x, x <= 1, cvx.sum_entries(x) == 1]
    
        prob = cvx.Problem(objective, constraints)
    
        min_dist = np.power(prob.solve(solver=solver, verbose=verbose),1/p)
        
        y = 100 * x
        
        # our list comprehension should really be over x.value.shape[0], but this is always just equal to A.shape[1] = len(pop_names) = num_src_pops
        
        results_table.append([target_df.values[tr][0]] + [dec_places % (y.value[i][0,0]) for i in range(num_src_pops)] + ['%.6f' % min_dist])
               
    # write results_file
    result_file_suffix = time.strftime("%d-%m-%y_%H-%M-%S", time.gmtime()) 
    headers = ["Population"] + pop_names + ["Min dist^2"]
    df = pd.DataFrame(results_table, columns=headers)
    if distances == False:
        df = df.transpose()[:-1].transpose()
    df.to_csv(result_file_prefix + result_file_suffix + '.csv')
    
    return

















