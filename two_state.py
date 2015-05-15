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

Model()

Monomer('Cell', ['d'], {'d' : ['0', '1']})
Monomer('Drug')

Parameter('kdiv', 1.02)       # /hr
Parameter('kdth', 1)         # /hr
Parameter('kf', 1e4)         # /uM-hr
Parameter('kr', 100.0)       # /hr
Parameter('kdivDrug', 1.02)   # /hr
Parameter('kdthDrug', 1.025)   # /hr

Initial(Cell(d='0'), Parameter('Cell0', 1000)) # number
Initial(Drug(), Parameter('Drug0', 0)) # uM

Observable('Cell_total', Cell())
Observable('Cell_free', Cell(d='0'))
Observable('Cell_drug', Cell(d='1'))

def full_model():
    Rule('Division', Cell(d='0') >> Cell(d='0') + Cell(d='0'), kdiv)
    Rule('Death', Cell(d='0') >> None, kdth)

    # Rule('DrugActivation', Cell(d='0') + Drug() >> Cell(d='1') + Drug(), kf)
    Rule('DrugActivation', Cell(d='0') >> Cell(d='1'), Expression('kf_Drug', sympify("kf*Drug0")))
    Rule('DrugInactivation', Cell(d='1') >> Cell(d='0'), kr)
    
    Rule('DivisionInDrug', Cell(d='1') >> Cell(d='1') + Cell(d='1'), kdivDrug)
    Rule('DeathInDrug', Cell(d='1') >> None, kdthDrug)

def reduced_model():
    Expression('kdiv_eff', sympify("(kr*kdiv + kf*Drug0*kdivDrug) / (kr + kf*Drug0)"))
    Expression('kdth_eff', sympify("(kr*kdth + kf*Drug0*kdthDrug) / (kr + kf*Drug0)"))
    Rule('Division', Cell(d='0') >> Cell(d='0') + Cell(d='0'), kdiv_eff)
    Rule('Death', Cell(d='0') >> None, kdth_eff)

# full_model()
reduced_model()

logD = np.arange(-4.0,1.1,0.1)
E0 = kdiv.value - kdth.value
Emax = kdivDrug.value - kdthDrug.value
colors = [cm.spectral(i) for i in np.linspace(0.2, 0.9, len(logD)+1)]
sm = plt.cm.ScalarMappable(cmap=mplcolors.ListedColormap(colors), norm=plt.Normalize(vmin=int(logD[0]), vmax=int(logD[-1])))
sm._A = []

slopes_ode = []
slopes_ssa = []
tspan = np.linspace(0, 100, 1001)

### Zero drug
x = odesolve(model, tspan, verbose=False)
plt.figure('Total Cells (ODE)')
plt.plot(tspan, np.log2(x['Cell_total']/x['Cell_total'][0]), '--', lw=6, color='navy', label='Untreated')
plt.legend(loc=0)

for i,ld in enumerate(logD):
    
    print ld    
    Drug0.value = 10.**ld
    
    # ODE simulations
    x = odesolve(model, tspan, verbose=False)
#     x = run_ssa(model, tspan, verbose=False)
    try:
        # throw out the first 10% of points (just to be safe)
        slope, intercept, r_value, p_value, std_err = ss.linregress(tspan[100:], np.log(x['Cell_total'][100:])) # Note: natural log
    except RuntimeWarning:
        pass
    else:
        slopes_ode.append(slope)
    
    plt.figure('Total Cells (ODE)')
    plt.plot(tspan, np.log2(x['Cell_total']/x['Cell_total'][0]), lw=3, color=colors[i])
    
    fig = plt.figure('% Drugged Cells (ODE)')
    plt.plot(tspan[1:], x['Cell_drug'][1:]/x['Cell_total'][1:]*100., lw=3, color=colors[i])
    
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
    
slopes_ode = np.array(slopes_ode)

plt.figure('Total Cells (ODE)')
plt.xlabel('time')
plt.ylabel('population doublings')
cb = plt.colorbar(sm, label='log10(drug)', orientation='horizontal')
# cb.set_label("xxx", labelpad=50)

# plt.figure('% Drugged Cells (ODE)')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('time')
# plt.ylabel('% drugged cells')
# fig.get_axes()[0].get_xaxis().set_major_formatter(FormatStrFormatter('%g'))
# fig.get_axes()[0].get_yaxis().set_major_formatter(FormatStrFormatter('%g'))
# cb = plt.colorbar(sm, label='log10(drug)', orientation='horizontal')

slopes_ssa = np.array(slopes_ssa)
effects_ssa = (slopes_ssa-Emax)/(E0-Emax)
mean_ssa = np.mean(effects_ssa, axis=1)
sdev_ssa = np.std(effects_ssa, axis=1)

plt.figure('Total Cells (SSA)')
plt.xlabel('time')
plt.ylabel('population doublings')
cb = plt.colorbar(sm, label='log10(drug)', orientation='horizontal')

# plt.figure('% Drugged Cells (SSA)')
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel('time')
# plt.ylabel('% drugged cells')
# fig.get_axes()[0].get_xaxis().set_major_formatter(FormatStrFormatter('%g'))
# fig.get_axes()[0].get_yaxis().set_major_formatter(FormatStrFormatter('%g'))
# cb = plt.colorbar(sm, label='log10(drug)', orientation='horizontal')

plt.figure('Dose Response')
plt.plot(logD, (slopes_ode-Emax)/(E0-Emax), 'o', ms=10, mfc='none', mec='b', mew=2, label='ODE')
####
# plt.plot(logD, effects_ssa, 'x', ms=10, mfc='none', mec='k', mew=2, label='SSA')
plt.errorbar(logD, mean_ssa, sdev_ssa, marker='o', ls='none', ms=10, mfc='none', mec='g', mew=2, label='SSA (n=3)')
####
plt.plot(logD, kr.value/(kr.value + kf.value*10**logD), 'r', lw=3, label='analytical (PEA)')
plt.legend(loc=0)
plt.xlabel('log10(drug)')
plt.ylabel(r'$\frac{\mathrm{DIP-DIP_{max}}}{\mathrm{DIP_0-DIP_{max}}}$') #, fontsize=16)
plt.xlim(xmin=int(logD[0]), xmax=int(logD[-1]))

# np.savetxt('TEMP/logD.txt', logD)
# np.savetxt('TEMP/effect.txt', kr.value/(kr.value + kf.value*10**logD))

plt.show()
