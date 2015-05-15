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

# from matplotlib import rc
# 
# rc('text', usetex=True)

Model()

Monomer('Cell', ['d'], {'d' : ['0', '1', '2']})
Monomer('Drug')

Parameter('kdiv0', 1.02)       # /hr
Parameter('kdth0', 1)         # /hr
Parameter('kf1', 1e7) # 1e4)         # /uM-hr
Parameter('kr1', 1e3) #100.0)       # /hr
Parameter('kdiv1', 1.02)   # /hr
Parameter('kdth1', 1.01)   # /hr
Parameter('kf2', 1e4)         # /uM-hr
Parameter('kr2', 1e4) #100.0)       # /hr
Parameter('kdiv2', 1.02)   # /hr
Parameter('kdth2', 1.025)   # /hr

Initial(Cell(d='0'), Parameter('Cell0', 1000)) # number
Initial(Drug(), Parameter('Drug0', 0)) # uM

Observable('Cell_total', Cell())
Observable('Cell_free', Cell(d='0'))
Observable('Cell_inter', Cell(d='1'))
Observable('Cell_drugged', Cell(d='2'))

def full_model():
    Rule('Division0', Cell(d='0') >> Cell(d='0') + Cell(d='0'), kdiv0)
    Rule('Death0', Cell(d='0') >> None, kdth0)

    Rule('DrugActivation_0_1', Cell(d='0') >> Cell(d='1'), Expression('kf1_Drug', sympify("kf1*Drug0")))
    Rule('DrugInactivation_1_0', Cell(d='1') >> Cell(d='0'), kr1)
    
    Rule('Division1', Cell(d='1') >> Cell(d='1') + Cell(d='1'), kdiv1)
    Rule('Death1', Cell(d='1') >> None, kdth1)
    
    Rule('DrugActivation_1_2', Cell(d='1') >> Cell(d='2'), Expression('kf2_Drug', sympify("kf2*Drug0")))
    Rule('DrugInactivation_2_1', Cell(d='2') >> Cell(d='1'), kr2)
    
    Rule('Division2', Cell(d='2') >> Cell(d='2') + Cell(d='2'), kdiv2)
    Rule('Death2', Cell(d='2') >> None, kdth2)

# def reduced_model():
#     Expression('kdiv_eff', sympify("(kr*kdiv + kf*Drug0*kdivDrug) / (kr + kf*Drug0)"))
#     Expression('kdth_eff', sympify("(kr*kdth + kf*Drug0*kdthDrug) / (kr + kf*Drug0)"))
#     Rule('Division', Cell(d='0') >> Cell(d='0') + Cell(d='0'), kdiv_eff)
#     Rule('Death', Cell(d='0') >> None, kdth_eff)

full_model()
# reduced_model()

# logD = np.arange(-4.0,1.1,0.1)
logD = np.arange(-6.0,2.1,0.1)
min_logD = int(round(logD[0]))
max_logD = int(round(logD[-1]))
E0 = kdiv0.value - kdth0.value
Emax = kdiv2.value - kdth2.value
colors = [cm.spectral(i) for i in np.linspace(0.2, 0.9, len(logD)+1)]
sm = plt.cm.ScalarMappable(cmap=mplcolors.ListedColormap(colors), norm=plt.Normalize(vmin=min_logD, vmax=max_logD))
sm._A = []

slopes_ode = []
# slopes_ssa = []
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

slopes_ode = np.array(slopes_ode)

plt.figure('Total Cells (ODE)')
plt.xlabel('time')
plt.ylabel('population doublings')
cb = plt.colorbar(sm, label='log10(drug)', orientation='horizontal')
plt.savefig('trajectories.pdf')

plt.figure('Dose Response')
plt.plot(logD, (slopes_ode-Emax)/(E0-Emax), 'o', ms=10, mfc='none', mec='b', mew=2, label='ODE')
# plt.plot(logD, effects_ssa, 'x', ms=10, mfc='none', mec='k', mew=2, label='SSA')
# plt.errorbar(logD, mean_ssa, sdev_ssa, marker='o', ls='none', ms=10, mfc='none', mec='g', mew=2, label='SSA (n=3)')
#####
Kd1 = kr1.value / kf1.value
Kd2 = kr2.value / kf2.value
print "Kd1:", Kd1
print "Kd2:", Kd2
DIP_0 = kdiv0.value - kdth0.value
DIP_int = kdiv1.value - kdth1.value
DIP_max = kdiv2.value - kdth2.value
Gamma = (DIP_int - DIP_max) / (DIP_0 - DIP_max)
d = 10**logD
f = Kd1*Kd2 + Kd2*d + d*d
plt.plot(logD, (Kd1*Kd2 + Kd2*d*Gamma)/f, 'r', lw=3, label='analytical (PEA)')
#####
twostate_logD = np.loadtxt('TEMP/logD.txt')
twostate_effect = np.loadtxt('TEMP/effect.txt')
plt.plot(twostate_logD, twostate_effect, 'o', ms=8, mfc='none', mec='0.7', mew=2, label='two-state')
#####
# plt.plot(np.log10(np.sqrt(Kd1*Kd2)), 0.5, '*k', ms=30, color='gold')
# plt.annotate(r"$\sqrt{K_{d1} K_{d2}}$", (0.45,0.5), xycoords='axes fraction', fontsize=18)
#####
plt.legend(loc=0)
plt.xlabel('log10(drug)')
plt.ylabel(r'$\frac{\mathrm{DIP-DIP_{max}}}{\mathrm{DIP_0-DIP_{max}}}$')
print "max_logD:", max_logD
plt.xlim(xmin=min_logD, xmax=max_logD)
plt.savefig('dose_response_relative.pdf')
#####
plt.figure('DIP rate')
DIP = (Kd1*Kd2*DIP_0 + Kd2*d*DIP_int + d*d*DIP_max) / f
plt.plot(logD, DIP, 'r', lw=3, label='analytical (PEA)')
plt.legend(loc=0)
plt.xlabel('log10(drug)')
plt.ylabel('DIP') 
plt.xlim(xmin=min_logD, xmax=max_logD)
plt.savefig('dose_response_absolute.pdf')

plt.show()

