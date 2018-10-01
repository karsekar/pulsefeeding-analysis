# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 05:24:39 2018

@author: karse
"""

import tellurium as te
import roadrunner

r = te.loada("""
model ftsZ_pulsing()
    // Reactions
    J0: -> FtsZ; alpha0 + f*alpha1
    J1: FtsZ -> ; Vmax*FtsZ/(Km + FtsZ)
    
    // Species initializations:
    FtsZ = 2000
    
    // Variable initialization:
    f = 0; Km = 600; Vmax = 10; alpha0 = 5.4; alpha1 = 12.9; 
end
""")

result = r.simulate(0, 120, 500)
r.plot(result)

sbml_file = open('ftsZ_pulsing.sbml','w')
print(r.getSBML(),file = sbml_file)
sbml_file.close()