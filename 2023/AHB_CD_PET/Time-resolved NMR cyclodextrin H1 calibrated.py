'''
Time resolved NMR to investigate the effect of cyclodextrins on product inhibition in enzymatic PET-degradation using a Fusarium solani pisi cutinase
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



#calculate TSP
vol_D2O_TSP_uL = 5
total_vol_uL = 600
uL_to_L = 1e-6
D2O_density = 1.11 # 1 g/mL = 1 mg/uL
D2O_TSP_mg = vol_D2O_TSP_uL * D2O_density
TSP_wt_frac = 0.75/100
TSP_Mw = 172.27
TSP_mmol = TSP_wt_frac * D2O_TSP_mg / TSP_Mw
TSP_mM = TSP_mmol / (total_vol_uL * uL_to_L)
print('TSP:', '%.2f' % TSP_mM, 'mM')

#What is this?
#control = np.loadtxt('integrals_expno10.txt', usecols=1, delimiter=';')[:196]
time_step = 5 # min


fig, axs = plt.subplots(2,3, constrained_layout=True, figsize=(12,6))
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 12
size=2.4

c1 = 'cornflowerblue'
c2 = 'firebrick'
c3 = 'grey'


#Control
data = np.loadtxt('integrals_42.txt', delimiter=';', usecols=(1,2,3,4,5))
time = np.linspace(0,len(data)*time_step, num=len(data), dtype=int)
TSP = (data[:,4]/9)
BHET = 1e3 * TSP_mM * (data[:,0]/4)/TSP
MHET = 1e3 * TSP_mM * ((data[:,1]+data[:,2])/4)/TSP
TPA = 1e3 * TSP_mM * (data[:,3]/4)/TSP
axs[0,0].scatter(time,MHET, color=c1, s=size, label='MHET')
axs[0,0].scatter(time, TPA, s=size, color=c2, label='TPA')
axs[0,0].scatter(time,BHET,s=size, color=c3, alpha=0.4, label='BHET')
axs[0,0].set_title('Control')
axs[0,0].legend(loc='best', markerscale=3, labelspacing=0.1)
axs[0,0].set_ylabel(r'Concentration [µM]', fontsize='14')


axs[1,0].scatter(time,MHET/(MHET+TPA+BHET), color=c1, s=size)
axs[1,0].scatter(time, TPA/(MHET+TPA+BHET), s=size, color=c2)
axs[1,0].scatter(time,BHET/(MHET+TPA+BHET),s=size, color=c3, alpha=0.4)
axs[1,0].set_ylabel(r'Molar fraction', fontsize='14')
axs[1,0].set_ylim(0,1)

print(f'control: BHET: {max(BHET)}, MHET: {max(MHET)}, TPA: {max(TPA)}')
print(f'control, apparent rate (BHET): {max(BHET)/(max(time)/60)} µM/h')
print(f'control, apparent rate (MHET): {max(MHET)/(max(time)/60)} µM/h')
print(f'control, apparent rate (TPA): {max(TPA)/(max(time)/60)} µM/h')
'''

#Alpha-CD

data = np.loadtxt('integrals_45.txt', delimiter=';', usecols=(1,2,3,4,5))
time = np.linspace(0,len(data)*time_step, num=len(data), dtype=int)
TSP = (data[:,4]/9)
BHET = 1e3 * TSP_mM * (data[:,0]/4)/TSP
MHET = 1e3 * TSP_mM * ((data[:,1]+data[:,2])/4)/TSP
TPA = 1e3 * TSP_mM * (data[:,3]/4)/TSP
axs[0,1].scatter(time,MHET, color=c1, s=size, label='MHET')
axs[0,1].scatter(time, TPA, s=size, color=c2, label='TPA')
axs[0,1].scatter(time,BHET,s=size, color=c3, alpha=0.4, label='BHET')
axs[0,1].set_title('Alpha-CD 1.3 mM')

axs[1,1].scatter(time,MHET/(MHET+TPA+BHET), color=c1, s=size)
axs[1,1].scatter(time, TPA/(MHET+TPA+BHET), s=size, color=c2)
axs[1,1].scatter(time,BHET/(MHET+TPA+BHET),s=size, color=c3, alpha=0.4)

axs[1,1].set_ylim(0,1)
'''


#H1-cal alpha-CD
a_CD_mM = 1.547

data = np.loadtxt('integrals_45_alpha-CD-H1.txt', delimiter=';', usecols=(1,2,3,4,5))
time = np.linspace(0,len(data)*time_step, num=len(data), dtype=int)
H1 = (data[:,4]/6)
BHET = 1e3 * a_CD_mM * (data[:,0]/4)/H1
MHET = 1e3 * a_CD_mM * ((data[:,1]+data[:,2])/4)/H1
TPA = 1e3 * a_CD_mM * (data[:,3]/4)/H1
axs[0,1].scatter(time,MHET, color=c1, s=size, label='MHET')
axs[0,1].scatter(time, TPA, s=size, color=c2, label='TPA')
axs[0,1].scatter(time,BHET,s=size, color=c3, alpha=0.4, label='BHET')
axs[0,1].set_title('H1-cal. Alpha-CD 1.55 mM')


axs[1,1].scatter(time,MHET/(MHET+TPA+BHET), color=c1, s=size)
axs[1,1].scatter(time, TPA/(MHET+TPA+BHET), s=size, color=c2)
axs[1,1].scatter(time,BHET/(MHET+TPA+BHET),s=size, color=c3, alpha=0.4)

axs[1,1].set_ylim(0,1)

print(f'alpha: BHET: {max(BHET)}, MHET: {max(MHET)}, TPA: {max(TPA)}')
print(f'alpha, apparent rate (BHET): {max(BHET)/(max(time)/60)} µM/h')
print(f'alpha, apparent rate (MHET): {max(MHET)/(max(time)/60)} µM/h')
print(f'alpha, apparent rate (TPA): {max(TPA)/(max(time)/60)} µM/h')

'''

#Beta-CD

data = np.loadtxt('integrals_48.txt', delimiter=';', usecols=(1,2,3,4,5))
time = np.linspace(0,len(data)*time_step, num=len(data), dtype=int)
TSP = (data[:,4]/9)
BHET = 1e3 * TSP_mM * (data[:,0]/4)/TSP
MHET = 1e3 * TSP_mM * ((data[:,1]+data[:,2])/4)/TSP
TPA = 1e3 * TSP_mM * (data[:,3]/4)/TSP
axs[0,3].scatter(time,MHET, color=c1, s=size, label='MHET')
axs[0,3].scatter(time, TPA, s=size, color=c2, label='TPA')
axs[0,3].scatter(time,BHET,s=size, color=c3, alpha=0.4, label='BHET')
axs[0,3].set_title('Beta-CD 1.4 mM')




axs[1,3].scatter(time,MHET/(MHET+TPA+BHET), color=c1, s=size)
axs[1,3].scatter(time, TPA/(MHET+TPA+BHET), s=size, color=c2)
axs[1,3].scatter(time,BHET/(MHET+TPA+BHET),s=size, color=c3, alpha=0.4)

axs[1,3].set_ylim(0,1)
'''

#H1-cal beta_CD
b_CD_mM = 1.652

data = np.loadtxt('integrals_48_beta-CD-H1.txt', delimiter=';', usecols=(1,2,3,4,5))
time = np.linspace(0,len(data)*time_step, num=len(data), dtype=int)
H1 = (data[:,4]/7)
BHET = 1e3 * b_CD_mM * (data[:,0]/4)/H1
MHET = 1e3 * b_CD_mM * ((data[:,1]+data[:,2])/4)/H1
TPA = 1e3 * b_CD_mM * (data[:,3]/4)/H1
axs[0,2].scatter(time,MHET, color=c1, s=size, label='MHET')
axs[0,2].scatter(time, TPA, s=size, color=c2, label='TPA')
axs[0,2].scatter(time,BHET,s=size, color=c3, alpha=0.4, label='BHET')
axs[0,2].set_title('H1-cal. Beta-CD 1.65 mM')


axs[1,2].scatter(time,MHET/(MHET+TPA+BHET), color=c1, s=size)
axs[1,2].scatter(time, TPA/(MHET+TPA+BHET), s=size, color=c2)
axs[1,2].scatter(time,BHET/(MHET+TPA+BHET),s=size, color=c3, alpha=0.4)

axs[1,2].set_ylim(0,1)

print(f'beta: BHET: {max(BHET)}, MHET: {max(MHET)}, TPA: {max(TPA)}')
print(f'beta, apparent rate (BHET): {max(BHET)/(max(time)/60)} µM/h')
print(f'beta, apparent rate (MHET): {max(MHET)/(max(time)/60)} µM/h')
print(f'beta, apparent rate (TPA): {max(TPA)/(max(time)/60)} µM/h')



#fig.text(0.5, -0.03, r'Time \ min', fontsize='14', ha='center')
#plt.figtext(0.5, 0, 'Time [min]', fontsize='14')
#axs[1,1].text(0.5,-0.3, 'Time [min]', fontsize='14', ha='center')
axs[1,1].set_xlabel('Time [min]',fontsize='14', ha='center')


ymax=400
#Adjust y-axes for concentration plots
for i in range(0,3):
    axs[0,i].set_ylim(0,ymax)





#plt.savefig('integrals3.png', dpi=300)
plt.show()


