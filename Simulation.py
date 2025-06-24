# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 13:13:48 2023

@author: Zhijie
"""

# Here, we assees whether the simulated species or the experimental species meet the first and second rules sensu Tilman.


import numpy as np
import pandas as pd
from itertools import combinations
import sympy as sp
import math
from tqdm import tqdm
import os





#%%
# function of calculating the impact vectors. Larger value indicate stronger impact on resource 1
def cal_imp(c_r1, s_r1, c_r2, s_r2, r1, r2 ): # r1 is the concentration of resource 1, and r2 the concentration of resource 2
    return (c_r1 * r1) / (s_r1 + r1) / (c_r2 * r2 / (s_r2 + r2))


# def cal_r_star(c, w, k, s, m):
#     return m * k * s/ (c * w - (k + c) * m)

# calculating the angle of the two impact vectors
def cal_angle(slope1, slope2):
    angle_radians = math.atan(abs((slope1 - slope2) / (1 + slope1 * slope2)))
    angle_degrees = math.degrees(angle_radians)
    return angle_degrees


#%%

# draw 12 random species
def f_create_paras(n_sp_final = 12):
    # -----------------prepare the parameters--------------------------#
    n_sp = n_sp_final *2 # simulate more species than wanted. Because in some cases, the R_star is negative, i.e., the algae can not grow when alone
    t = pd.DataFrame(np.random.uniform(0, 2, (n_sp, 10))) # random parameters for each species
    # Setting low mortality rate does not guarantee that birth > death rate. But at least in most of the cases.
    t[4] = t[5] = 0.1 # same mortality rate for all species
    t['r1_star'] = t[4]*t[6]*t[8]/(t[0] * t[2] - t[6] * t[4] - t[0]*t[4]) # m * k * s/ (c * w - (k + c) * m)
    t['r2_star'] = t[5]*t[7]*t[9]/(t[1] * t[3] - t[7] * t[5] - t[1]*t[5])
    t = t[(t['r1_star'] > 0) & (t['r2_star'] > 0)].head(n_sp_final) # only take n_sp_final species whose r_star is positive.
    t['sp'] = np.arange(1, len(t) + 1)
    return t


# assess the first and second rules
def f_simu(t, # the parameter table
           type = 'simulated', # whether simulted random species, or the experimental species (real species)
           export = 'summary'
           ):    
    # create an table to summarize the results
    t_meet_combinations = list(combinations(t['sp'], 2))
    t_meet = pd.DataFrame(t_meet_combinations, columns=['sp1', 'sp2'])
    t_meet['meet_req_ess'] = 0 # proportion of the species pairs meet the first rule when competing for essential resources
    t_meet['meet_imp_ess'] = 0 # second rule
    t_meet['angle_ess'] = 0    # the angle between two impact vectors
    t_meet['v1_ess'] = 0    # impact vector of species 1
    t_meet['v2_ess'] = 0    # impact vector of species 2
    t_meet['eq1_ess'] = 0   # equilibrium of resource 1
    t_meet['eq2_ess'] = 0   # equilbirum of resource 2
    
    t_meet['meet_req_sub'] = 0 # compting for substitutable resources
    t_meet['meet_imp_sub'] = 0
    t_meet['angle_sub'] = 0
    t_meet['v1_sub'] = 0
    t_meet['v2_sub'] = 0
    t_meet['eq1_sub'] = 0
    t_meet['eq2_sub'] = 0

    # check coexistence for each species pair
    for i, row in t_meet.iterrows():
        pars_sp1 = t[t['sp'] == row['sp1']].iloc[0][:12] # get the parameters of species 1 
        pars_sp2 = t[t['sp'] == row['sp2']].iloc[0][:12] # parameters of species 2
        
            
        
        c11, c12, w11, w12, m11, m12, k11, k12, s11, s12, r11_star, r12_star = pars_sp1 # first index is the species
        c21, c22, w21, w22, m21, m22, k21, k22, s21, s22, r21_star, r22_star = pars_sp2

           
        
        # ------------------- essential resources -------------------#
        # ------------It is actually easy because the rule one can be checked with r_star -------#
        # find the equalibrium
        
        ## In case the parameters are from the experiment
        if type == 'exp':
            pars_sp1 = t[t['sp'] == row['sp1']].iloc[0][:18]
            pars_sp2 = t[t['sp'] == row['sp2']].iloc[0][:18]
            c11, u11, m11, k11,  s11, r11_star,  c12, u12, m12, k12, s12, r12_star = pars_sp1[:12]
            c21, u21, m21, k21,  s21, r21_star,  c22, u22, m22, k22, s22, r22_star = pars_sp2[:12]

             
        eq1 = max(r11_star, r21_star) # resource 1
        eq2 = max(r12_star, r22_star) # resource 2
        
        # impact vector
        v1 = cal_imp(c11, s11, c12, s12, eq1, eq2)
        v2 = cal_imp(c21, s21, c22, s22, eq1, eq2)
        t_meet.loc[i, 'v1_ess'] = v1
        t_meet.loc[i, 'v2_ess'] = v2
        
        if r11_star < r21_star and r12_star > r22_star: # sp 1 more limited by resource 2 (this also guarantees the presence of solutions)
            t_meet.loc[i, 'meet_req_ess'] = 1
            t_meet.loc[i, 'eq1_ess']    = eq1
            t_meet.loc[i, 'eq2_ess']    = eq2
            if v1 < v2: # sp 1 has stronger impact on resource 2
                t_meet.loc[i, 'meet_imp_ess'] = 1
                t_meet.loc[i, 'angle_ess']    = cal_angle(v1, v2)
        
        if r11_star > r21_star and r12_star < r22_star:
            t_meet.loc[i, 'meet_req_ess'] = 1
            t_meet.loc[i, 'eq1_ess']    = eq1
            t_meet.loc[i, 'eq2_ess']    = eq2
            if v1 > v2:
                t_meet.loc[i, 'meet_imp_ess'] = 1
                t_meet.loc[i, 'angle_ess']    = cal_angle(v1, v2)
        
        
        # ----------------------substitutable resource---------------------------#. 
        # More difficult than the essental resource, because we need to solve the solutions for two nonlinear curves
        # however, for rule 1, it has similar result as the essential resource
        
        # find the equalibrium
        r1, r2 = sp.symbols('r1 r2')
        f1 = w11*c11*r1/(k11*s11 + k11* r1 + c11* r1) + w12*c12* r2/(k12*s12 + k12* r2 + c12*r2) - m11 #- m12
        f2 = w21*c21*r1/(k21*s21 + k21* r1 + c21* r1) + w22*c22* r2/(k22*s22 + k22* r2 + c22*r2) - m21 #- m22
        

        if type == 'exp':
            c11, u11, m11, k11,  s11, r11_star = pars_sp1[:6] 
            c12, u12, m12, k12, s12, r12_star  = pars_sp1[12:]
            c21, u21, m21, k21,  s21, r21_star = pars_sp2[:6] 
            c22, u22, m22, k22, s22, r22_star  = pars_sp2[12:]
            f1 = u11*r1/(k11 + r1) + u12 *r2/(k12 + r2) - m11 - m12 + math.log(2/3)/2
            f2 = u21*r1/(k21 + r1) + u22 *r2/(k22 + r2) - m21  -m22 + math.log(2/3)/2
            
        # whether there is a tradeoff in resource requirement, i.e. whether there is any positive solutions of the two functions
        solutions = sp.solve((f1, f2), (r1, r2))
        real_solutions = [arr for arr in solutions if all(sp.im(expr) == 0 for expr in arr)] # remove imaginary number
        positive_solutions = [arr for arr in real_solutions if min(arr) > 0] # only keep positive solutions

        v1 = cal_imp(c11, s11, c12, s12, eq1, eq2)
        v2 = cal_imp(c21, s21, c22, s22, eq1, eq2)
        t_meet.loc[i, 'v1_sub'] = v1
        t_meet.loc[i, 'v2_sub'] = v2


        if len(positive_solutions) > 1: # in very rare cases ( ~ 0.5%), there are two positive solutions due to the nonlinearity. Discard it.
            t_meet.loc[i, 'meet_req_ess'] = np.nan
            t_meet.loc[i, 'meet_imp_ess'] = np.nan
            t_meet.loc[i, 'angle_ess'] = np.nan            
            t_meet.loc[i, 'meet_req_sub'] = np.nan
            t_meet.loc[i, 'meet_imp_sub'] = np.nan
            t_meet.loc[i, 'angle_sub'] = np.nan
            continue
        
        if len(positive_solutions) == 1:
            t_meet.loc[i, 'meet_req_sub'] = 1
            # the equilibriums
            eq1 = positive_solutions[0][0] 
            eq2 = positive_solutions[0][1]
            t_meet.loc[i, 'eq1_sub']    = eq1
            t_meet.loc[i, 'eq2_sub']    = eq2
            # the slopes
            s11 = sp.diff(f1, r1).subs(r1, eq1) # resource 1 on sp1
            s12 = sp.diff(f1, r2).subs(r2, eq2)
            s21 = sp.diff(f2, r1).subs(r1, eq1)
            s22 = sp.diff(f2, r2).subs(r2, eq2)
            if s11/s12 > s21/s22 and v1 > v2: # sp1 more limited by resource 1, and sp1 has larger impact on resource 1
                t_meet.loc[i, 'meet_imp_sub'] = 1                
                t_meet.loc[i, 'angle_sub']    = cal_angle(v1, v2)
            if s11/s12 < s21/s22 and v1 < v2: # sp1 more limited by resource 1, and sp1 has larger impact on resource 1
                t_meet.loc[i, 'meet_imp_sub'] = 1
                t_meet.loc[i, 'angle_sub']    = cal_angle(v1, v2)
                
    # --------------calculate the percentages-------------------------#
    
    meet_req_ess  = np.nanmean(t_meet['meet_req_ess'])
    meet_req_sub  = np.nanmean(t_meet['meet_req_sub'])
    meet_imp_ess  = np.nansum(t_meet['meet_imp_ess']) / np.nansum(t_meet['meet_req_ess'])
    meet_imp_sub  = np.nansum(t_meet['meet_imp_sub']) / np.nansum(t_meet['meet_req_sub'])
    angle_ess     = np.nan if meet_imp_ess == 0 else np.nanmean(t_meet['angle_ess'][t_meet['meet_imp_ess'] == 1])
    angle_sub     = np.nan if meet_imp_sub == 0 else np.nanmean(t_meet['angle_sub'][t_meet['meet_imp_sub'] == 1])
    
    
    if export == 'summary': #only get the summary data
        return [meet_req_ess, meet_req_sub, 
                meet_imp_ess, meet_imp_sub,
                angle_ess ,   angle_sub
                ]
    if export == 'table': # get the data for each species pair
        return t_meet


#%%
# simulate 999 set of 12 species and calculate the proportions of meeting first and second rules
df = pd.DataFrame(columns= ['meet_req_ess', 'meet_req_sub', 
                            'meet_imp_ess', 'meet_imp_sub',
                            'angle_ess' ,   'angle_sub'
        ])

         

for i in tqdm(range(999), desc= 'Processing', unit= 'item'):
    t = f_create_paras()
    res = f_simu(t)
    df.loc[i,] = res
# df.to_csv('C:\\02_postdoc\\08_bottom_up\\main_nc\\df_simu.csv', index = False)



#%%

# same as above, but with the experimental data
# With the parameters from the experiment, calculate the proportions

t = pd.read_csv('C:\\02_postdoc\\08_bottom_up\\main_r2\\t_para.csv')  
res_obs = f_simu(t, type = 'exp', export = 'table')
res_obs = res_obs.drop(columns=['eq1_ess', 'eq2_ess', 'eq1_sub', 'eq2_sub'])
res_obs.to_csv('C:\\02_postdoc\\08_bottom_up\\main_r2\\df_obs.csv', index = False) # save to be anlyzed in R


#%%
# export the details to R, to test the case of multispecies
t_list = []
t_meet__list = []    
               
for i in tqdm(range(999), desc= 'Processing', unit= 'item'):
    t = f_create_paras()
    t_meet = f_simu(t, export='table')
    t_list.append(t.assign(i = i))
    t_meet__list.append(t_meet.assign(i = i))


pd.concat(t_list)
output_dir = 'C:\\02_postdoc\\08_bottom_up\\main_nc\\paras'
os.makedirs(output_dir, exist_ok=True)

pd.concat(t_list).to_csv(os.path.join(output_dir, 'para.csv'), index= False)
pd.concat(t_meet__list).to_csv(os.path.join(output_dir, 'meet.csv'), index= False)


