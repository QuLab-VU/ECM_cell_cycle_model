# -*- coding: utf-8 -*-
"""
Created on Fri May  1 11:07:19 2020

@author: Grim.M
"""

#test random forest

path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"


import sys
sys.path.insert(0, path+"/cellular_cycle") #add a pathway to the model

import numpy as np

import matplotlib.pyplot as plt

import seaborn as sns
sns.set(style='white', context='poster')

import umap

import sobol_seq

from pysb.simulator import ScipyOdeSimulator

from peakfinder import peakfinder

from Austin_proliferation_model_concentration import model



 

"""
sobol_seq_set_accepted = np.concatenate([sobol_seq_set1_accepted, sobol_seq_set2_accepted])
sobol_seq_set_rejected = np.concatenate([sobol_seq_set1_rejected, sobol_seq_set2_rejected])

sobol_seq_set = np.concatenate([sobol_seq_set_accepted, sobol_seq_set_rejected])
target = [1 for i in range(len(sobol_seq_set_accepted))] + [0 for i in range(len(sobol_seq_set_rejected))]
"""


embedding = umap.UMAP().fit_transform(sobol_seq_set, y=target)

fig, ax = plt.subplots(1, figsize=(14, 10))
plt.scatter(*embedding.T, s=0.1, c=target, cmap='Spectral', alpha=1.0)
plt.setp(ax, xticks=[], yticks=[])
cbar = plt.colorbar(boundaries=np.arange(3)-0.5)
cbar.set_ticks(np.arange(10))
plt.title('sobolseq Embedded via UMAP using Labels');





mapper = umap.UMAP(n_neighbors=10).fit(sobol_seq_set, y=target)


test=sobol_seq.i4_sobol_generate(10, 8000, skip=202000)*2
test=np.random.rand(10000,10)*2
###############################################################################
pysb_sampled_parameter_names = ["V_s2skp","V_Fp","k_M","k_EMY","k_CDF","k_CDS","k_ED","K_S","K_MY","K_DF"]

t = np.linspace(0, 200, 120)
accepted=[]
rejected=[]

for parameter_vector in test :
    param_dict = {pname: pvalue for pname, pvalue in zip(pysb_sampled_parameter_names, parameter_vector)}
    for pname, pvalue in param_dict.items():
        #Change model parameter values to current location in parameter space
        model.parameters[pname].value = (pvalue)
    simres = ScipyOdeSimulator(model, t).run() 
    yout = simres.all
    a, maxima = peakfinder(yout['obsCdc20'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    b, maxima_b = peakfinder(yout['obsCycE'][40:120],[],0.2,includeEndpoints=1)   #find the peak of cdc20 - last protein to peak
    if len(a)>1 and len(b)>1:
        accepted+=[parameter_vector]
    else :
        rejected+=[parameter_vector]
    print(len(accepted)+len(rejected))

###############################################################################



mapper = umap.UMAP(n_neighbors=10).fit(sobol_seq_set, y=target)

test = np.concatenate([sobol_accepted_200000201000, sobol_rejected_200000201000, sobol_accepted_201000202000, sobol_rejected_201000202000])
test_labels = [1 for i in range(len(sobol_accepted_200000201000))] + [0 for i in range(len(sobol_rejected_200000201000))] + [1 for i in range(len(sobol_accepted_201000202000))] + [0 for i in range(len(sobol_rejected_201000202000))]


test_embedding = mapper.transform(test)




fig, ax = plt.subplots(1, figsize=(14, 10))
plt.scatter(*test_embedding.T, s=2, c=np.array(test_labels), cmap='Spectral', alpha=1.0)
plt.setp(ax, xticks=[], yticks=[])
cbar = plt.colorbar(boundaries=np.arange(3)-0.5)
cbar.set_ticks(np.arange(10))
plt.title('sobol seq Train Digits Embedded via UMAP');



from sklearn.ensemble import RandomForestRegressor



sobol_seq_set = np.concatenate([sobol_seq_set1_again_accepted, sobol_seq_set1_again_rejected, sobol_seq_set2_again_accepted, sobol_seq_set2_again_rejected, sobol_seq_set3_again_accepted, sobol_seq_set3_again_rejected, sobol_seq_set4_again_accepted, sobol_seq_set4_again_rejected])

target = [1 for i in range(len(sobol_seq_set1_again_accepted))] + [0 for i in range(len(sobol_seq_set1_again_rejected))] + [1 for i in range(len(sobol_seq_set2_again_accepted))] + [0 for i in range(len(sobol_seq_set2_again_rejected))] + [1 for i in range(len(sobol_seq_set3_again_accepted))] + [0 for i in range(len(sobol_seq_set3_again_rejected))] + [1 for i in range(len(sobol_seq_set4_again_accepted))] + [0 for i in range(len(sobol_seq_set4_again_rejected))]


test = np.concatenate([sobol_accepted_200000201000, sobol_rejected_200000201000, sobol_accepted_201000203000, sobol_rejected_201000203000, sobol_accepted_202000210000, sobol_rejected_202000210000])
test_labels = [1 for i in range(len(sobol_accepted_200000201000))] + [0 for i in range(len(sobol_rejected_200000201000))] + [1 for i in range(len(sobol_accepted_201000203000))] + [0 for i in range(len(sobol_rejected_201000203000))] + [1 for i in range(len(sobol_accepted_202000210000))] + [0 for i in range(len(sobol_rejected_202000210000))]

regressor = RandomForestRegressor(n_estimators = 60, random_state = 0)
regressor.fit(sobol_seq_set, target)

y_pred = regressor.predict(test)

error=0
false_positive=0
false_negative=0


for i in range(len(y_pred)) :
    if test_labels[i]==0 :
        if y_pred[i] >= 0.2 :
            error+=1
            false_positive+=1
        
    else:   
        if y_pred[i] < 0.2 :
            error+=1
            false_negative+=1
        
            
error_rate = error/len(y_pred)            
            
TPR = (sum(test_labels) - false_negative) / sum(test_labels) #true positive rate (true positive / all real positive) the proportion of the individuals with a known positive condition for which the test resultis positive.

TNR = ((len(y_pred)-sum(test_labels))- false_positive) / (len(y_pred)-sum(test_labels)) #true negative rate (true negative / all real negative)


#On veut avoir TPR = 1

################################################################################
#Hyperparameters
#nber of trees

{'bootstrap': True,
 'criterion': 'mse',
 'max_depth': None,
 'max_features': 'auto',
 'max_leaf_nodes': None,
 'min_impurity_decrease': 0.0,
 'min_impurity_split': None,
 'min_samples_leaf': 1,
 'min_samples_split': 2,
 'min_weight_fraction_leaf': 0.0,
 'n_estimators': 10,
 'n_jobs': 1,
 'oob_score': False,
 'random_state': 42,
 'verbose': 0,
 'warm_start': False}


TNR=[]
TPR=[]
error_rate=[]

for trees in np.linspace(32,80,25):
    print(trees)
    regressor = RandomForestRegressor(n_estimators = int(trees), random_state = 0)
    regressor.fit(sobol_seq_set, target)
    y_pred = regressor.predict(test)
    error=0
    false_positive=0
    false_negative=0
    
    
    for i in range(len(y_pred)) :
        if test_labels[i]==0 :
            if y_pred[i] >= 0.2 :
                error+=1
                false_positive+=1
            
        else:   
            if y_pred[i] < 0.2 :
                error+=1
                false_negative+=1
        
            
    error_rate += [error/len(y_pred)]           
                
    TPR += [(sum(test_labels) - false_negative) / sum(test_labels)] #true  positive rate (true positive / all real positive) the proportion of the individuals     with a known positive condition for which the test resultis positive.

    TNR += [((len(y_pred)-sum(test_labels))- false_positive) / (len(y_pred)-sum(test_labels))] #true negative rate (true negative / all real negative)


#thresholds


TNR=[]
TPR=[]
error_rate=[]

regressor = RandomForestRegressor(n_estimators = 60, random_state = 0)
regressor.fit(sobol_seq_set, target)

y_pred = regressor.predict(test)


for threshold in np.linspace(0.1,0.9,50):
    error=0
    false_positive=0
    false_negative=0
    for i in range(len(y_pred)) :
        if test_labels[i]==0 :
            if y_pred[i] >= threshold :
                error+=1
                false_positive+=1
        
        else:   
            if y_pred[i] < threshold :
                error+=1
                false_negative+=1
        
            
    error_rate += [error/len(y_pred)]            
            
    TPR += [(sum(test_labels) - false_negative) / sum(test_labels)] #true positive rate (true positive / all real positive) the proportion of the individuals with a known     positive condition for which the test resultis positive.

    TNR += [((len(y_pred)-sum(test_labels))- false_positive) / (len(y_pred)-sum(test_labels))] #true negative rate (true negative / all real negative)



indices = [i for i, x in enumerate(error_rate) if x == min(error_rate)]


fig = plt.figure(figsize=(12, 8))
line_weight = 3
alpha = .5

ax1 = fig.add_axes([0, 0, 1, 1])
ax2 = fig.add_axes()
ax2 = ax1.twinx()

lns1 = ax1.plot(np.linspace(0.1,0.9,50), error_rate, markersize=12, color='skyblue', linewidth=4, label="errors")
lns2 = ax2.plot(np.linspace(0.1,0.9,50), TPR, marker='', color='red', linewidth=2, label="TPR")
lns3 = ax2.plot(np.linspace(0.1,0.9,50), TNR, marker='', color='purple', linewidth=2, label="TNR")
ax1.axvspan(np.linspace(0.1,0.9,50)[indices[0]], np.linspace(0.1,0.9,50)[indices[-1]], alpha=0.5, color='grey')
leg = lns1 + lns2 + lns3
labs = [l.get_label() for l in leg]
ax1.legend(leg, labs, loc=0)

ax1.set_xticks(np.linspace(0.1,0.9,17))

fig.savefig(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\trees.pdf", dpi=fig.dpi, format="pdf", bbox_inches='tight')

fig.show()


np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\random_forest\TPR_randomforest_60trees_0.1-0.9threshold.npy",TPR)

np.save(r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code\fit_data\random_forest\error_rate_randomforest_60trees_0.1-0.9threshold.npy",error_rate)

##########################################################################
from sklearn.tree import export_graphviz
# Export as dot file

estimator = regressor.estimators_[5]

export_graphviz(estimator, out_file='tree.dot', 
                feature_names = iris.feature_names,
                class_names = iris.target_names,
                rounded = True, proportion = False, 
                precision = 2, filled = True)

# Convert to png using system command (requires Graphviz) 
from subprocess import call
call(['dot', '-Tpng', 'tree.dot', '-o', 'tree.png', '-Gdpi=600'])


import numpy
from sklearn.model_selection import train_test_split
from sklearn import metrics, datasets, ensemble

def print_decision_rules(rf):

    for tree_idx, est in enumerate(rf.estimators_):
        tree = est.tree_
        assert tree.value.shape[1] == 1 # no support for multi-output

        print('TREE: {}'.format(tree_idx))

        iterator = enumerate(zip(tree.children_left, tree.children_right, tree.feature, tree.threshold, tree.value))
        for node_idx, data in iterator:
            left, right, feature, th, value = data

            # left: index of left child (if any)
            # right: index of right child (if any)
            # feature: index of the feature to check
            # th: the threshold to compare against
            # value: values associated with classes            

            # for classifier, value is 0 except the index of the class to return
            class_idx = numpy.argmax(value[0])

            if left == -1 and right == -1:
                print('{} LEAF: return class={}'.format(node_idx, class_idx))
            else:
                print('{} NODE: if feature[{}] < {} then next={} else next={}'.format(node_idx, feature, th, left, right))    



data=[tree.children_left[8], tree.children_right[8], tree.feature[8], tree.threshold[8], tree.value[8]]


print_decision_rules(regressor)

import os
from sklearn.tree import export_graphviz

    estimator = regressor.estimators_[5]


export_graphviz(estimator, out_file='C:/users/Grim.M/Ecole/ens/Master/S2/Stage/report/figures/test_tree.dot', node_ids=True, filled=True, rounded=True)

#"dot -Tpdf test_tree.dot -o test_tree.pdf"copy in cmd.exe (trzanslate .dot to .pdf)

##################################################################################
#extracting subspaces from one tree
export_graphviz(estimator, out_file='C:/users/Grim.M/Ecole/ens/Master/S2/Stage/report/figures/test_tree.dot', node_ids=True, filled=True, rounded=True)


estimator = regressor.estimators_[5]
tree=estimator.tree_

iterator = enumerate(zip(tree.children_left, tree.children_right, tree.feature, tree.threshold, tree.value))

leave_0=[]
leave_1=[]
nodes=[]

for node_idx, data in iterator:
    left, right, feature, th, value = data
    if left == -1 and right == -1 and value[0]==1: #if it's a leaf corresponding to the good state
        leave_1+=[[data, node_idx]]
    elif left == -1 and right == -1 and value[0]==0:
        leave_0+=[[data, node_idx]]
    nodes+=[[data,node_idx]]
        

        
subspaces_0=[[[0,2] for i in range(10)] for j in range(len(leave_0))]    

for leaf in range(len(leave_0)):
    current_node = leave_0[leaf][-1]
    previous_node = len(nodes)-1 #on scanne tout
    while current_node != 0:

        left, right, feature, th, value = nodes[previous_node][0] #ggg
        
        if left == nodes[current_node][1]:
            print(feature)
            print("left")
            print(th)
            print(previous_node)
            current_node = previous_node
            if subspaces_0[leaf][feature][1]>th: #
                subspaces_0[leaf][feature][1]=th
            previous_node = previous_node - 1
            
        elif right == nodes[current_node][1]:
            print(feature)
            print("right")
            print(th)
            print(previous_node)
            current_node = previous_node
            if subspaces_0[leaf][feature][0]<th: #
                subspaces_0[leaf][feature][0]=th
            previous_node = previous_node - 1
        else : 
            previous_node = previous_node - 1




subspaces_1=[[[0,2] for i in range(10)] for j in range(len(leave_1))]    

for leaf in range(len(leave_1)):
    current_node = leave_1[leaf][-1]
    previous_node = len(nodes)-1 #on scanne tout
    while current_node != 0:

        left, right, feature, th, value = nodes[previous_node][0] #ggg
        
        if left == current_node:
            print(feature)
            print("left")
            print(th)
            current_node = previous_node
            if th<subspaces_1[leaf][feature][1]: #
                subspaces_1[leaf][feature][1]=th
            previous_node = previous_node - 1
            
        elif right == current_node:
            print(feature)
            print("right")
            print(th)
            current_node = previous_node
            if th>subspaces_1[leaf][feature][0]: #
                subspaces_1[leaf][feature][0]=th
            previous_node = previous_node - 1
        else : 
            previous_node = previous_node - 1

