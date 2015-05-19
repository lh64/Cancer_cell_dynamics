from pysb import *
from pysb.integrate import odesolve
from pysb.bng import run_ssa
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import matplotlib.cm as cm
import matplotlib.colors as mplcolors
import numpy as np
import scipy.stats as ss
from sympy import sympify
import os

figs = os.path.join('.','FIGS')
if not os.path.exists(figs):
    os.makedirs(figs)

Model()

Monomer('Cell', ['d'], {'d' : ['0', '1']})
Monomer('Drug')

Parameter('kdiv0', 1.02)       # /hr
Parameter('kdth0', 1)          # /hr
Parameter('kon',   5e-2) #1e4)    # /uM-hr
Parameter('koff',  5e-4) #1e2)    # /hr
Parameter('kdiv1', 1.02)       # /hr
Parameter('kdth1', 1.025)      # /hr

Initial(Cell(d='0'), Parameter('Cell0_init', 1000)) # number
Initial(Cell(d='1'), Parameter('Cell1_init', 0)) # number
Initial(Drug(), Parameter('Drug0', 0)) # uM

Observable('Cell_total', Cell())
Observable('Cell_free', Cell(d='0'))
Observable('Cell_drug', Cell(d='1'))

def full_model():
    Rule('Division', Cell(d='0') >> Cell(d='0') + Cell(d='0'), kdiv0)
    Rule('Death', Cell(d='0') >> None, kdth0)

    # Rule('DrugActivation', Cell(d='0') + Drug() >> Cell(d='1') + Drug(), kf)
    Rule('DrugActivation', Cell(d='0') >> Cell(d='1'), Expression('kon_Drug', sympify("kon*Drug0")))
    Rule('DrugInactivation', Cell(d='1') >> Cell(d='0'), koff)
    
    Rule('DivisionInDrug', Cell(d='1') >> Cell(d='1') + Cell(d='1'), kdiv1)
    Rule('DeathInDrug', Cell(d='1') >> None, kdth1)

def reduced_model():
    Expression('kdiv_eff', sympify("(koff*kdiv0 + kon*Drug0*kdiv1) / (koff + kon*Drug0)"))
    Expression('kdth_eff', sympify("(koff*kdth0 + kon*Drug0*kdth1) / (koff + kon*Drug0)"))
    Rule('Division', Cell(d='0') >> Cell(d='0') + Cell(d='0'), kdiv_eff)
    Rule('Death', Cell(d='0') >> None, kdth_eff)

full_model()
# reduced_model()

def plot_timecourse(x, logD, clr, label):
    plt.figure('log(drug) = %g' % logD)
    plt.plot(tspan, np.log2(x['Cell_total']/x['Cell_total'][0]), lw=3, color=clr, label=label)
#     plt.plot(tspan[100:], np.log2(x['Cell_total'][100:]/x['Cell_total'][0]), 'r--', lw=3)
    plt.xlabel('time (h)')
    plt.ylabel('population doublings')
    plt.legend(loc=0)

for n,clr in zip([0,1,1], ['blue','green','darkorange']):
    kon.value  *= 10**n
    koff.value *= 10**n

    logD = np.arange(-4.0,2.1,0.1)
    DIP0 = kdiv0.value - kdth0.value
    DIPmax = kdiv1.value - kdth1.value
    colors = [cm.spectral(i) for i in np.linspace(0.2, 0.9, len(logD)+1)]
    sm = plt.cm.ScalarMappable(cmap=mplcolors.ListedColormap(colors), norm=plt.Normalize(vmin=int(logD[0]), vmax=int(logD[-1])))
    sm._A = []
    
    slopes_ode = []
    slopes_ssa = []
    tspan = np.linspace(0, 200, 201)
    
    ### Zero drug
    Drug0.value = 0
    x = odesolve(model, tspan, verbose=False)
    plt.figure('Total Cells (ODE) (kon,koff)=(%.2g,%.2g)' % (kon.value,koff.value))
    plt.plot(tspan, np.log2(x['Cell_total']/x['Cell_total'][0]), '--', lw=8, color='navy', label='Untreated', zorder=10)
    plt.legend(loc=0)
    
    for i,ld in enumerate(logD):
        
        print ld
        Drug0.value = 10.**ld
        
        # ODE simulations
        x = odesolve(model, tspan, verbose=False)
    #     x = run_ssa(model, tspan, verbose=False)
        #####
        if ld < -0.95 and ld > -1.05 or ld < 0.05 and ld > -0.05:
            plot_timecourse(x, round(ld*10.)/10., clr, '(kon,koff)=(%.2g,%.2g)' % (kon.value,koff.value))
            plt.savefig(os.path.join(figs,'timecourses_%.1g_uM.pdf' % (10**ld)))
        #####
        try:
            # throw out the first 100 hrs
            slope, intercept, r_value, p_value, std_err = ss.linregress(tspan[100:], np.log(x['Cell_total'][100:])) # Note: natural log
        except RuntimeWarning:
            pass
        else:
            slopes_ode.append(slope)
        
        plt.figure('Total Cells (ODE) (kon,koff)=(%.2g,%.2g)' % (kon.value,koff.value))
        plt.plot(tspan, np.log2(x['Cell_total']/x['Cell_total'][0]), lw=3, color=colors[i])
        
        fig = plt.figure('%% Drugged Cells (ODE) (kon,koff)=(%.2g,%.2g)' % (kon.value,koff.value))
        plt.plot(tspan[1:], x['Cell_drug'][1:]/x['Cell_total'][1:]*100., lw=3, color=colors[i])
        '''
        # SSA simulations (run each in triplicate)
        slopes_ssa.append(3*[None])
        for j in range(3):
            y = run_ssa(model, tspan, verbose=False)
            try:
                # throw out the first 10% of points (just to be safe)
                slope, intercept, r_value, p_value, std_err = ss.linregress(y['time'][100:], np.log(y['Cell_total'][100:])) # Note: natural log
            except RuntimeWarning:
                pass
            else:
                slopes_ssa[-1][j] = slope
                
        plt.figure('Total Cells (SSA)')
        plt.plot(y['time'], np.log2(y['Cell_total']/y['Cell_total'][0]), lw=3, color=colors[i])
        
        fig = plt.figure('% Drugged Cells (SSA)')
        plt.plot(y['time'][1:], y['Cell_drug'][1:]/y['Cell_total'][1:]*100., lw=3, color=colors[i])
        '''
    slopes_ode = np.array(slopes_ode)
    
    plt.figure('Total Cells (ODE) (kon,koff)=(%.2g,%.2g)' % (kon.value,koff.value))
    plt.xlabel('time (h)')
    plt.ylabel('population doublings')
    cb = plt.colorbar(sm, label='log10(drug)', orientation='horizontal', ticks=np.arange(-4,2.1))
    # cb.set_label("xxx", labelpad=50)
    
    plt.figure('%% Drugged Cells (ODE) (kon,koff)=(%.2g,%.2g)' % (kon.value,koff.value))
#     plt.xscale('log')
#     plt.yscale('log')
    plt.xlabel('time (h)')
    plt.ylabel('%% drugged cells')
    fig.get_axes()[0].get_xaxis().set_major_formatter(FormatStrFormatter('%g'))
    fig.get_axes()[0].get_yaxis().set_major_formatter(FormatStrFormatter('%g'))
    cb = plt.colorbar(sm, label='log10(drug)', orientation='horizontal', ticks=np.arange(-4,2.1))
    '''
    slopes_ssa = np.array(slopes_ssa)
    effects_ssa = (slopes_ssa-Emax)/(E0-Emax)
    mean_ssa = np.mean(effects_ssa, axis=1)
    sdev_ssa = np.std(effects_ssa, axis=1)
    
    plt.figure('Total Cells (SSA)')
    plt.xlabel('time (h)')
    plt.ylabel('population doublings')
    cb = plt.colorbar(sm, label='log10(drug)', orientation='horizontal')
    '''
    # plt.figure('% Drugged Cells (SSA)')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.xlabel('time')
    # plt.ylabel('% drugged cells')
    # fig.get_axes()[0].get_xaxis().set_major_formatter(FormatStrFormatter('%g'))
    # fig.get_axes()[0].get_yaxis().set_major_formatter(FormatStrFormatter('%g'))
    # cb = plt.colorbar(sm, label='log10(drug)', orientation='horizontal')
    
    plt.figure('Dose Response')
    # plt.plot(logD, (slopes_ode-Emax)/(E0-Emax), 'o', ms=10, mfc='none', mec='blue', mew=2, label='ODE')
    plt.plot(logD, slopes_ode, 'o', ms=10, mfc='none', mec=clr, mew=2, label='(kon,koff)=(%.2g,%.2g)' % (kon.value,koff.value))
    ####
    # plt.plot(logD, effects_ssa, 'x', ms=10, mfc='none', mec='k', mew=2, label='SSA')
    # plt.errorbar(logD, mean_ssa, sdev_ssa, marker='o', ls='none', ms=10, mfc='none', mec='g', mew=2, label='SSA (n=3)')
    ####

plt.figure('Dose Response')
ax = plt.subplot(111)
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
plt.plot(logD, (koff.value*DIP0 + kon.value*10**logD*DIPmax)/(koff.value + kon.value*10**logD), 'r', lw=3, label='analytical (PEA)')
plt.legend(loc=0, bbox_to_anchor=(0.6,1.))
plt.xlabel('log10(drug)')
plt.ylabel('DIP rate') #, fontsize=16)
plt.xlim(xmin=int(logD[0]), xmax=int(logD[-1]))

plt.savefig(os.path.join(figs,'two_state_dose_response_full_v_PEA.pdf'))

# np.savetxt('TEMP/logD.txt', logD)
# np.savetxt('TEMP/effect.txt', kr.value/(kr.value + kf.value*10**logD))

plt.show()
