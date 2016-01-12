from pysb import *

Model()

Monomer('Cell', ['d'], {'d' : ['0', '1']})
Monomer('Drug')

Parameter('kdiv0', 1.06)      # /hr
Parameter('kdth0', 1.)        # /hr
Parameter('kon',   0.1)       # /uM-hr
Parameter('koff',  1e-3)      # /hr
Parameter('kdiv1', 0.97)      # /hr
Parameter('kdth1', 1.)        # /hr

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
    Rule('DrugActivation', Cell(d='0') >> Cell(d='1'), Expression('kon_Drug', kon*Drug0))
    Rule('DrugInactivation', Cell(d='1') >> Cell(d='0'), koff)
    
    Rule('DivisionInDrug', Cell(d='1') >> Cell(d='1') + Cell(d='1'), kdiv1)
    Rule('DeathInDrug', Cell(d='1') >> None, kdth1)

def reduced_model():
    Expression('kdiv_eff', (koff*kdiv0 + kon*Drug0*kdiv1) / (koff + kon*Drug0))
    Expression('kdth_eff', (koff*kdth0 + kon*Drug0*kdth1) / (koff + kon*Drug0))
    Rule('Division', Cell(d='0') >> Cell(d='0') + Cell(d='0'), kdiv_eff)
    Rule('Death', Cell(d='0') >> None, kdth_eff)

full_model()
# reduced_model()