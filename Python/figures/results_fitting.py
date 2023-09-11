# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 08:39:48 2021

@author: Grim.M
"""

import matplotlib.patheffects as mpe

######
#parameters for the simulation
init = 50
run_time = 102.3
observed_cells_time = 30
omega=200

#parameters of the frac prolif model to compare
parameter_names = ['a','d','q']

#######################################################################################
#Experimental
path=r"C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\Code"
all_events_exp={}

IMT_eoe_counts = pd.DataFrame(columns=['End_of_Expt', 'Death' ,'Birthtime','Lifespan'])

df = pd.read_csv(path + '/cell_fates/ds8_06_DMSO_fates.csv')
lines = {}
for col in df :
    lines[col]=df[col].dropna().tolist()

for (death, mitotic_time, end_time, lifetime) in zip(lines['Death'], lines['Mitotic Time'], lines['End Time'], lines['Lifetime']):
    tmp_counts = {'End_of_Expt' : end_time==102.3, 'Death' : death==1, 'Birthtime' : mitotic_time, 'Lifespan' : lifetime}
    IMT_eoe_counts = IMT_eoe_counts.append(tmp_counts, ignore_index=True)

all_events_exp['0.6 kPa DMSO'] = all_events(run_time, IMT_eoe_counts)


#calculating rate of fractional proliferation model

mitotic_lifespans_exp = all_events_exp['0.6 kPa DMSO']['Lifespan'][all_events_exp['0.6 kPa DMSO']['event']=='mitosis']
censored_lifespans_exp = all_events_exp['0.6 kPa DMSO']['Lifespan'][all_events_exp['0.6 kPa DMSO']['event']=='End_of_Expt']

number_of_deaths = len(all_events_exp['0.6 kPa DMSO']) - (all_events_exp['0.6 kPa DMSO']['event']=='apoptosis').value_counts()[0]
cell_time_total = sum(all_events_exp['0.6 kPa DMSO']['Lifespan'])


m = fp.mle(mitotic_lifespans_exp, censored_lifespans_exp)

parameters_exp={}
parameters_exp['0.6 kPa DMSO'] ={}

parameters_exp['0.6 kPa DMSO']['a'] = number_of_deaths/cell_time_total
parameters_exp['0.6 kPa DMSO']['mu'] = m['x'][0]
parameters_exp['0.6 kPa DMSO']['sigma'] = m['x'][1]
parameters_exp['0.6 kPa DMSO']['lambda'] = m['x'][2]
parameters_exp['0.6 kPa DMSO'].update(fp.q_rates(m))

parameters_exp['0.6 kPa DMSO']['q']=1e-5
#######################################


experimental_SMs={}
experimental_G1s={}
experimental_IMTs={}


df = pd.read_csv(path+'/IMT_times/DS8_SM_0_2_55.csv')

for col in df.columns :
    experimental_SMs[col]=df[col].dropna().tolist()


df = pd.read_csv(path+'/IMT_times/DS8_G1_0_2_55.csv')

for col in df.columns :
    experimental_G1s[col]=df[col].dropna().tolist()

for k in experimental_G1s.keys() :
    minimum=min(len(experimental_G1s[k]),len(experimental_SMs[k]))
    experimental_IMTs[k]=[experimental_G1s[k][i]+experimental_SMs[k][i] for i in range(minimum)]
    

#################################################################################
def RBS(x,x0,sigma) :
    observations , dimensions= np.shape(x0)
    s=0
    for i in range(observations):
        s+=math.exp(-np.linalg.norm(x-x0[i])**2/(2*sigma**2)) / (sigma*math.sqrt(2*math.pi))
    if observations != 0:
        return s/observations
    else:
        return s

################################################################################
mitotic_lifespans_exp.tolist()

sigma=5

x=np.linspace(0,120,360)
experimental_IMT_distribution=[RBS(i, np.transpose([mitotic_lifespans_exp.tolist()]), sigma) for i in x]



sigma=4
 
g1_distribution_exp = [RBS(i, np.transpose([experimental_G1s['55 kPa DMSO']]), sigma) for i in np.linspace(0,120,400)]
sm_distribution_exp = [RBS(i, np.transpose([experimental_SMs['55 kPa DMSO']]), sigma) for i in np.linspace(0,120,400)]

###################################################################################

g1_distribution_55 = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/g1_distribution_all_parameters_50_initial_cells_55.pkl")

sm_distribution_55 = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/sm_distribution_all_parameters_50_initial_cells_55.pkl")




IMT_distribution = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/IMT_distribution_all_parameters_50_initial_cells.pkl")


frac_prolif_parameters_006_55 = pd.read_pickle(path + r"/fit_data/sobol_seq/cost_function_model_4/cost_and_IMT_distribution/frac_prolif_parameters_all_parameters_006_55_kPa.pkl")
frac_prolif_parameters = frac_prolif_parameters_006_55



sum_column = frac_prolif_parameters["cost_006"] +frac_prolif_parameters["cost_frac_prolif_model_006"] + 0 * frac_prolif_parameters["cost_mmd_55"]/2

frac_prolif_parameters["cost_total"] = sum_column


sns.distplot(frac_prolif_parameters['cost_mmd_006'])
#sns.distplot(frac_prolif_parameters['cost_frac_prolif_model_006'])
sns.distplot(frac_prolif_parameters['cost_006'])

plt.scatter(frac_prolif_parameters['cost_mmd_006'],frac_prolif_parameters['cost_006'])

plt.scatter(frac_prolif_parameters['cost_mmd_55']/2,frac_prolif_parameters['cost_006'])


sns.distplot(frac_prolif_parameters['cost_total'])
sns.distplot(frac_prolif_parameters['cost_mmd_55']/2)
sns.distplot(frac_prolif_parameters['cost_006'])

sns.distplot(frac_prolif_parameters['cost_mmd_55']/3)
sns.distplot(frac_prolif_parameters['cost_006'])

np.mean(frac_prolif_parameters['cost_mmd_55']/3)
np.mean(frac_prolif_parameters['cost_006'])



##################################################################################


fn = Path(r'C:\Users\Grim.M\Ecole\ens\Master\S2\Stage\report\figures\sobol_seq\results_fitting\results_high_and_low_stiffness_cost_mmd_55_div_by_inf_best_60_costs.pdf').expanduser()

fig, axs = pl.subplots(ncols=2, nrows=2, sharex=False, sharey=False, figsize=(40,20), constrained_layout=True)

############### 1st subplot

minimum = frac_prolif_parameters.sort_values(by=['cost_total'])[0:50]['cost_total'].tolist()[2]


cost_limit = frac_prolif_parameters.sort_values(by=['cost_total'])[0:60]['cost_total'].tolist()[-1]



x=np.linspace(0,120,360)
colormap = cm.RdYlBu
normalize = mcolors.Normalize(vmin = minimum, vmax = cost_limit)
#or 


for idx in frac_prolif_parameters['index']:
    idx = int(idx)
    if frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_total']<cost_limit :
    #if frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_mmd']< inf :
        axs[0, 1].plot(x, IMT_distribution[idx].tolist()
,color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_total'])), lw=0.4)
#        axs[0, 1].grid(True)
        axs[0, 1].tick_params(axis="y", labelsize=24)
        axs[0, 1].tick_params(axis="x", labelsize=24)
#plotting original parameters and expt
axs[0, 1].plot(x, experimental_IMT_distribution, color='black', lw=6)
axs[0, 1].plot(x, IMT_distribution[1].tolist(), lw=5, color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==1.0].iloc[0]['cost_total'])), path_effects=[mpe.withStroke(linewidth=6, foreground='black')])


############### 2nd subplot


for idx in frac_prolif_parameters['index']:
    idx = int(idx)
    if frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_total']<cost_limit: 
        frac_prolif_parameters_temp = frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx]
        parameters_temp = frac_prolif_parameters_temp[['a_006','d_006','q_006']].iloc[0].as_matrix()
        
        axs[1, 1].plot(parameter_names, np.log10(parameters_temp), color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_total'])), lw=0.4)
#        axs[1, 1].grid(True)
        axs[1, 1].tick_params(axis="y", labelsize=24)
        axs[1, 1].tick_params(axis="x", labelsize=24)

#plotting original parameters and expt
axs[1, 1].plot(parameter_names,  np.log10([parameters_exp['0.6 kPa DMSO']['a'], parameters_exp['0.6 kPa DMSO']['d'], parameters_exp['0.6 kPa DMSO']['q']]), color='black', lw=6)
axs[1, 1].plot(parameter_names, np.log10(frac_prolif_parameters.loc[frac_prolif_parameters['index']==1.0][['a_006','d_006','q_006']].iloc[0].as_matrix()), lw=5, color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==1.0].iloc[0]['cost_total'])), path_effects=[mpe.withStroke(linewidth=6, foreground='black')])


############### 3rd subplot

x=np.linspace(0,120,400)

for idx in frac_prolif_parameters['index']:
    idx = int(idx)
    if frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_total']<cost_limit:
    #if frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_mmd']< inf :
        axs[0, 0].plot(x, g1_distribution_55[idx].tolist()
,color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_total'])), lw=0.4)
#        axs[0, 1].grid(True)
        axs[0, 0].tick_params(axis="y", labelsize=24)
        axs[0, 0].tick_params(axis="x", labelsize=24)
#plotting original parameters and expt
axs[0, 0].plot(x, g1_distribution_exp, color='black', lw=6)
#axs[0, 0].plot(x, g1_distribution_55[1].tolist(), lw=4,color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==1.0].iloc[0]['cost_total'])), path_effects = [mpe.Stroke(linewidth=8, foreground='black')])
#axs[0,0].plot(x, g1_distribution_55[1].tolist(), lw=5, color='red', path_effects=[mpe.withStroke(linewidth=8, foreground='black')], label="stroke no alpha")
axs[0,0].plot(x, g1_distribution_55[1].tolist(), lw=5, color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==1.0].iloc[0]['cost_total'])), path_effects=[mpe.withStroke(linewidth=6, foreground='black')])

############### 4th subplot


x=np.linspace(0,120,400)

for idx in frac_prolif_parameters['index']:
    idx = int(idx)
    if frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_total']<cost_limit:
    #if frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_mmd']< inf :
        axs[1, 0].plot(x, sm_distribution_55[idx].tolist()
,color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==idx].iloc[0]['cost_total'])), lw=0.4)
#        axs[0, 1].grid(True)
        axs[1, 0].tick_params(axis="y", labelsize=24)
        axs[1, 0].tick_params(axis="x", labelsize=24)
#plotting original parameters and expt
axs[1, 0].plot(x, sm_distribution_exp, color='black', lw=6)
axs[1, 0].plot(x, sm_distribution_55[1].tolist(), lw=4, color=colormap(normalize(frac_prolif_parameters.loc[frac_prolif_parameters['index']==1.0].iloc[0]['cost_total'])), path_effects=[mpe.withStroke(linewidth=6, foreground='black')])



scalarmappaple = cm.ScalarMappable(norm=normalize, cmap=colormap)
scalarmappaple.set_array([])

#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(scalarmappaple, ax=axs.ravel().tolist(), shrink=0.8)
cb.set_label(label='cost_total', fontsize=24)
cb.ax.tick_params(labelsize=12)

#pl.tight_layout()

#pl.tight_layout()

pl.savefig(fn, bbox_inches='tight')
