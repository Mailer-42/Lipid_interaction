#!/usr/bin/env python
# coding: utf-8

# In[11]:


#import lib
import sys
import pylab
import matplotlib.pyplot as plt
import numpy as np

import MDAnalysis as mda
from MDAnalysis.analysis.leaflet import LeafletFinder


# In[12]:


#insert file 
u = mda.Universe("md_final100ns_2.pdb", "md_final100ns_2.xtc")


# In[13]:


u


# In[14]:


len(u.trajectory)


# In[10]:


#selecting lower leaflet lipid around 6 Å protein 
def res(Input):
       #input == types of phospholipid
    if Input == 'POPC':
 
            res = []
            for ts in u.trajectory:
                atoms = u.select_atoms("prop z < 187.087 and resid 1154-1871 and resname POPC and around 6 protein")
                res += [len(atoms.residues)]
            res_traj = np.array(res)
            return res

    if Input == 'POPI':            
    
            res = []
            for ts in u.trajectory:
                atoms = u.select_atoms("prop z < 187.087 and resid 1154-1871 and resname POPI and around 6 protein")
                res += [len(atoms.residues)]
            res_traj = np.array(res)
            return res

    if Input == 'POPE':
    
            res = []
            for ts in u.trajectory:
                atoms = u.select_atoms("prop z < 187.087 and resid 1154-1871 and resname POPE and around 6 protein")
                res += [len(atoms.residues)]
            res_traj = np.array(res)
            return res
    
    if Input == 'POPS':
    
            res = []
            for ts in u.trajectory:
                atoms = u.select_atoms("prop z < 187.087 and resid 1154-1871 and resname POPS and around 6 protein")
                res += [len(atoms.residues)]
            res_traj = np.array(res)
            return res
 
    if Input == 'CHOL':
    
            res = []
            for ts in u.trajectory:
                atoms = u.select_atoms("prop z < 187.087 and resid 1154-1871 and resname CHOL and around 6 protein")
                res += [len(atoms.residues)]
            res_traj = np.array(res)
            return res
   
    if Input == 'DPSM':
    
            res = []
            for ts in u.trajectory:
                atoms = u.select_atoms("prop z < 187.087 and resid 1154-1871 and resname DPSM and around 6 protein")
                res += [len(atoms.residues)]
            res_traj = np.array(res)
            return res

    else:
        print("")
      


# In[15]:


#find all lipid
res_POPC = res('POPC')
res_POPC_c = sys.getsizeof(res_POPC)

res_POPI = res('POPI')
res_POPI_c = sys.getsizeof(res_POPI)

res_POPE = res('POPE')
res_POPE_c = sys.getsizeof(res_POPE)

res_POPS = res('POPS')
res_POPS_c = sys.getsizeof(res_POPS)

res_CHOL = res('CHOL')
res_CHOL_c = sys.getsizeof(res_CHOL)

res_DPSM = res('DPSM')
res_DPSM_c = sys.getsizeof(res_DPSM)

#find all lipid 
res_all = res_POPC + res_POPI + res_POPE + res_POPS + res_CHOL + res_DPSM 
res_all_c = sys.getsizeof(res_all)


# In[16]:


#calculation fraction (1 lipid type / all lipid)
def fraction(res,res_all):
    fractions = []
    for x in res:
        fractions.append(x / res_all_c)
        #fractions = [x / res_all for x in res]
    return fractions


# In[17]:


print(fraction(res('POPC'),res_all))


# In[18]:


# Enter data
f_POPC = fraction(res('POPC'),res_all)
f_POPI = fraction(res('POPI'),res_all)
f_POPE = fraction(res('POPE'),res_all)
f_POPS = fraction(res('POPS'),res_all)
f_CHOL = fraction(res('CHOL'),res_all)
f_DPSM = fraction(res('DPSM'),res_all)

# Calculate the average
mean_POPC = np.mean(f_POPC)
mean_POPI = np.mean(f_POPI)
mean_POPE = np.mean(f_POPE)
mean_POPS = np.mean(f_POPS)
mean_CHOL = np.mean(f_CHOL)
mean_DPSM = np.mean(f_DPSM)

# Calculate the standard deviation(SD, Error bars)
std_POPC = np.std(f_POPC)
std_POPI = np.std(f_POPI)
std_POPE = np.std(f_POPE)
std_POPS = np.std(f_POPS)
std_CHOL = np.std(f_CHOL)
std_DPSM = np.std(f_DPSM)


# In[19]:


#list for plot
Phospholipid = ['POPC', 'POPI', 'POPE', 'POPS', 'CHOL', 'DPSM']
Mean = [mean_POPC, mean_POPI, mean_POPE, mean_POPS, mean_CHOL, mean_DPSM]
Error = [std_POPC, std_POPI, std_POPE, std_POPS, std_CHOL, std_DPSM]
x_pos = np.arange(len(Phospholipid))


# In[24]:


# Build the plot
fig, ax = plt.subplots()
ax.bar(x_pos, Mean, yerr=Error, align='center', alpha=0.5, ecolor='black', color = 'c', capsize=10)
ax.set_ylabel('Interaction(%)')
ax.set_xticks(x_pos)
ax.set_xticklabels(Phospholipid)
ax.set_title('Time (Final 100ns)')
ax.yaxis.grid(True)


# In[ ]:




