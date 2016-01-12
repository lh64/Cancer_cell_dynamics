from two_state import model
from pysb.integrate import odesolve
from pysb.bng import run_ssa
import numpy as np
import matplotlib.pyplot as plt

model.parameters['Drug0'].value = 0.63

# for p in model.parameters:
#     print p
# quit()

tspan = np.arange(120.1)

# ODE
plt.figure()
for drug in [0.6, 0.63, 0.7, 0.8, 0.9, 1.0]:
    model.parameters['Drug0'].value = drug
    x = odesolve(model, tspan, verbose=True)
    print x['Cell_total'][0]
    plt.plot(tspan, np.log2(x['Cell_total']/x['Cell_total'][0]), label=drug)
plt.legend(loc=0)

# SSA 
for cells in [100, 1000, 10000]:
    plt.figure()
    model.parameters['Cell0_init'].value = cells
    for drug in [0.6, 0.63, 0.7, 0.8, 0.9, 1.0]:
        model.parameters['Drug0'].value = drug
        x = run_ssa(model, t_end=tspan[-1], n_steps=len(tspan)-1, verbose=False)
        print x['Cell_total'][0]
        plt.plot(tspan, np.log2(x['Cell_total']/x['Cell_total'][0]), label=drug)
    plt.legend(loc=0)
    
plt.show()
