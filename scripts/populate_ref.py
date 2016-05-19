#!/usr/bin/env python
"""populate_ref.py -YES -- Add the reference data for abinit/HGHk and wien2k/acc to FATMAN.

The reference data (V_0 [A^3/atom], B0 [GPa], B1) have been hardcoded here from the 
Lejaeghere, Science (2016) supporting information.
Execute the script on the FATMAN server to add these data to the TestResult table
and also add appropriate entries to the Method table.

The deltatests for all chemical elements have to be present in the Test table,
since the TestResult refers to them.

Only run once to prevent multiple entries in the DB.

Parameters:
    -h     Show this help.
    -YES   Call the program with this switch to add the entries. 
           Otherwise nothing happens to prevent accidental DB clusterfuck.
"""

from __future__ import print_function
from sys import argv

from fatman.models import *
from ase.units import Ry
from datetime import datetime


def abinit_ref():
    """Reference data for HGHk/abinit as given in the Lejaeghere, Science (2016) SI """
    mydata = { #V0        B0        B1
       "H":  (  17.422 ,  10.262  , 2.683  ), 
       "He": (  17.781 ,  0.865   , 6.201  ),
       "Li": (  20.210 ,  13.848  , 3.345  ),
       "Be": (  7.893  ,  122.959 , 3.332  ), 
       "B":  (  7.223  ,  234.668 , 3.435  ),
       "C":  (  11.633 ,  207.213 , 3.551  ),
       "N":  (  28.794 ,  53.568  , 3.666  ),
       "O":  (  18.681 ,  50.004  , 3.518  ),
       "F":  (  19.292 ,  33.769  , 3.730  ),
       "Ne": (  24.799 ,  1.730   , 9.142  ),
       "Na": (  37.586 ,  7.489   , 6.691  ),
       "Mg": (  23.201 ,  38.235  ,-12.582 ), 
       "Al": (  16.460 ,  77.517  , 4.861  ),
       "Si": (  20.355 ,  88.278  , 4.278  ),
       "P":  (  21.260 ,  68.722  , 4.313  ), 
       "S":  (  17.035 ,  84.313  , 4.049  ),
       "Cl": (  38.257 ,  19.302  , 4.372  ),
       "Ar": (  52.181 ,  0.755   , 7.308  ),
       "K":  (  73.602 ,  3.594   , 3.773  ),
       "Ca": (  42.276 ,  17.618  , 3.379  ), 
       "Sc": (  24.630 ,  54.548  , 3.392  ),
       "Ti": (  17.414 ,  111.829 , 3.603  ),
       "V":  (  13.473 ,  182.469 , 3.911  ),
       "Cr": (  12.175 ,  136.321 , 7.032  ),
       "Mn": (  12.006 ,  119.492 , 5.921  ),
       "Fe": (  11.464 ,  177.303 , 7.785  ),
       "Co": (  10.915 ,  210.008 , 5.081  ),
       "Ni": (  10.925 ,  193.949 , 4.730  ), 
       "Cu": (  11.955 ,  144.750 , 5.652  ),
       "Zn": (  15.234 ,  74.393  , 4.628  ),
       "Ga": (  20.439 ,  46.913  , 7.184  ),
       "Ge": (  24.071 ,  57.597  , 4.789  ),
       "As": (  22.565 ,  68.686  , 4.278  ),
       "Se": (  29.726 ,  47.199  , 4.445  ),
       "Br": (  39.336 ,  22.570  , 4.842  ),
       "Kr": (  65.867 ,  0.649   , 7.221  ), 
       "Rb": (  91.058 ,  2.794   , 3.799  ),
       "Sr": (  54.452 ,  11.304  , 5.145  ),
       "Y":  (  32.884 ,  41.228  , 3.132  ),
       "Zr": (  23.353 ,  93.903  , 3.270  ),
       "Nb": (  18.141 ,  170.147 , 3.696  ),
       "Mo": (  15.819 ,  259.942 , 4.350  ),
       "Tc": (  14.477 ,  299.499 , 4.529  ),
       "Ru": (  13.802 ,  312.577 , 4.875  ), 
       "Rh": (  14.090 ,  257.333 , 5.209  ),
       "Pd": (  15.369 ,  169.094 , 5.544  ),
       "Ag": (  18.022 ,  89.845  , 6.078  ),
       "Cd": (  22.869 ,  44.751  , 6.972  ),
       "In": (  27.664 ,  35.871  , 5.068  ),
       "Sn": (  37.094 ,  34.150  , 4.716  ),
       "Sb": (  31.977 ,  50.220  , 4.513  ),
       "Te": (  35.162 ,  44.717  , 4.702  ), 
       "I":  (  50.546 ,  18.567  , 5.055  ),
       "Xe": (  86.827 ,  0.542   , 7.182  ),
       "Cs": (  116.83 ,  1.962   , 3.608  ),
       "Ba": (  63.430 ,  8.747   , 2.063  ), 
       "Hf": (  22.470 ,  107.014 , 3.374  ),
       "Ta": (  18.151 ,  193.434 , 3.533  ),
       "W":  (  15.999 ,  303.367 , 4.101  ),
       "Re": (  14.818 ,  367.741 , 4.401  ), 
       "Os": (  14.137 ,  403.796 , 4.799  ),
       "Ir": (  14.342 ,  354.960 , 5.096  ),
       "Pt": (  15.603 ,  253.682 , 5.477  ),
       "Au": (  17.888 ,  143.464 , 6.011  ), 
       "Hg": (  29.405 ,  8.379   , 9.825  ),
       "Tl": (  31.484 ,  27.024  , 5.472  ),
       "Pb": (  32.047 ,  39.979  , 5.605  ),
       "Bi": (  36.970 ,  42.722  , 4.668  ),  
       "Po": (  37.627 ,  45.541  , 5.012  ),
       "Rn": (  93.350 ,  0.539   , 7.197  )                     
                      }
    
    basis  = BasissetFamily.get(BasissetFamily.name=='planewave')
    pseudo = PseudopotentialFamily.get(PseudopotentialFamily.name=='GTH-PBE')
        
    M, created = Method.create_or_get(basis_set        = basis, 
                         pseudopotential  = pseudo, 
                         code             = "ref-abinit", 
                         settings         = {   "cutoff_pw"  :  250.*Ry, 
                                                "fband"      : 1.5,
                                                "ecutsm"     : 0.0,
                                                "tolsym"     : 1e-12,
                                                "tsmear"     : 0.01,
                                                "xc"         : "PBE",
                                                "occopt"     : 3,
                                                "pawovlp"    : -1,
                                                "chksymbreak": 0,
                                                "prtwf"      : 0,
                                                "prtden"     : 0,
                                            }
                         )

    for e, data in mydata.items():
        v, b0, b1 =  data
        TestResult.create_or_get(ctime = datetime.now(),
                                 test = Test.get(Test.name == "deltatest_"+e),
                                 method = M,
                                 result_data = {"_status": "reference",
                                                "source" : "Science Manuscript, 2016",
                                                "V"      : v,
                                                "B0"     : b0,
                                                "B1"     : b1}
                                 )

def wien_ref():
    """Reference data for Wien2k/acc given in the Lejaeghere, Science (2016) SI """

    mydata = {
    "H":  ( 17.388   ,  10.284     ,  2.711  ), 
    "He": ( 17.771   ,  0.847      ,  7.708  ),
    "Li": ( 20.219   ,  13.839     ,  3.336  ),
    "Be": ( 7.910    ,  122.903    ,  3.036  ), 
    "B":  ( 7.240    ,  237.290    ,  3.468  ),
    "C":  ( 11.637   ,  208.991    ,  3.579  ),
    "N":  ( 28.885   ,  54.220     ,  3.724  ),
    "O":  ( 18.559   ,  51.378     ,  3.895  ),
    "F":  ( 19.167   ,  34.325     ,  3.935  ),
    "Ne": ( 24.249   ,  1.406      ,  14.437 ),
    "Na": ( 37.469   ,  7.472      ,  3.771  ),
    "Mg": ( 22.936   ,  35.933     ,  4.067  ), 
    "Al": ( 16.480   ,  78.077     ,  4.570  ),
    "Si": ( 20.453   ,  88.545     ,  4.306  ),
    "P":  ( 21.471   ,  68.208     ,  4.348  ), 
    "S":  ( 17.184   ,  83.407     ,  4.261  ),
    "Cl": ( 38.889   ,  19.081     ,  4.343  ),
    "Ar": ( 52.385   ,  0.743      ,  7.256  ),
    "K":  ( 73.679   ,  3.574      ,  4.593  ),
    "Ca": ( 42.199   ,  17.114     ,  3.312  ), 
    "Sc": ( 24.620   ,  54.393     ,  3.424  ),
    "Ti": ( 17.390   ,  112.213    ,  3.583  ),
    "V":  ( 13.452   ,  181.674    ,  3.745  ),
    "Cr": ( 11.773   ,  183.899    ,  7.158  ),
    "Mn": ( 11.447   ,  118.632    , -0.206  ),
    "Fe": ( 11.344   ,  197.652    ,  5.801  ),
    "Co": ( 10.864   ,  216.489    ,  4.363  ),
    "Ni": ( 10.891   ,  199.970    ,  5.006  ), 
    "Cu": ( 11.957   ,  141.121    ,  4.845  ),
    "Zn": ( 15.195   ,  74.572     ,  5.271  ),
    "Ga": ( 20.307   ,  49.223     ,  5.384  ),
    "Ge": ( 23.915   ,  59.128     ,  4.988  ),
    "As": ( 22.589   ,  68.285     ,  4.225  ),
    "Se": ( 29.744   ,  47.070     ,  4.441  ),
    "Br": ( 39.447   ,  22.415     ,  4.870  ),
    "Kr": ( 65.658   ,  0.671      ,  9.857  ), 
    "Rb": ( 90.809   ,  2.787      ,  5.798  ),
    "Sr": ( 54.527   ,  11.256     ,  3.490  ),
    "Y":  ( 32.844   ,  41.593     ,  3.016  ),
    "Zr": ( 23.385   ,  93.684     ,  3.207  ),
    "Nb": ( 18.137   ,  171.270    ,  3.548  ),
    "Mo": ( 15.786   ,  258.928    ,  4.332  ),
    "Tc": ( 14.437   ,  299.149    ,  4.459  ),
    "Ru": ( 13.762   ,  312.502    ,  4.953  ), 
    "Rh": ( 14.040   ,  257.824    ,  5.321  ),
    "Pd": ( 15.310   ,  168.629    ,  5.562  ),
    "Ag": ( 17.847   ,  90.148     ,  5.420  ),
    "Cd": ( 22.835   ,  44.082     ,  6.969  ),
    "In": ( 27.471   ,  34.937     ,  4.781  ),
    "Sn": ( 36.817   ,  36.030     ,  4.637  ),
    "Sb": ( 31.730   ,  50.367     ,  4.516  ),
    "Te": ( 34.977   ,  44.787     ,  4.691  ), 
    "I":  ( 50.233   ,  18.654     ,  5.046  ),
    "Xe": ( 86.681   ,  0.548      ,  6.344  ),
    "Cs": ( 117.080  ,  1.982      ,  2.141  ),
    "Ba": ( 63.140   ,  8.677      ,  3.771  ), 
   # "Lu": ( 29.054   ,  46.384     ,  2.943  ),
    "Hf": ( 22.532   ,  107.004    ,  3.498  ),
    "Ta": ( 18.286   ,  195.147    ,  3.714  ),
    "W":  ( 16.139   ,  301.622    ,  4.279  ),
    "Re": ( 14.958   ,  362.850    ,  4.517  ), 
    "Os": ( 14.280   ,  397.259    ,  4.844  ),
    "Ir": ( 14.500   ,  347.680    ,  5.179  ),
    "Pt": ( 15.642   ,  248.711    ,  5.465  ),
    "Au": ( 17.975   ,  139.109    ,  5.757  ), 
    "Hg": ( 29.612   ,  8.055      ,  8.899  ),
    "Tl": ( 31.390   ,  26.865     ,  5.489  ),
    "Pb": ( 32.003   ,  39.544     ,  4.533  ),
    "Bi": ( 36.905   ,  42.630     ,  4.705  ),  
    "Po": ( 37.587   ,  45.458     ,  4.926  ),
    "Rn": ( 92.685   ,  0.564      ,  8.618  )                     
                }
    
    basis  = BasissetFamily.get(BasissetFamily.name=='planewave')
    pseudo = PseudopotentialFamily.get(PseudopotentialFamily.name=='all-electron')
        
    M, created = Method.create_or_get(basis_set        = basis, 
                         pseudopotential  = pseudo, 
                         code             = "ref-wien2k", 
                         settings         = {  
                                            }
                         )

    for e, data in mydata.items():
        v, b0, b1 =  data
        TestResult.create_or_get(ctime = datetime.now(),
                                 test = Test.get(Test.name == "deltatest_"+e),
                                 method = M,
                                 result_data = {"_status": "reference",
                                                "source" : "Science Manuscript, 2016",
                                                "V"      : v,
                                                "B0"     : b0,
                                                "B1"     : b1}
                                 )

if __name__ == "__main__":
    if len(argv) == 2 and argv[1] == "-YES":
        abinit_ref()
        wien_ref()
    else:
        print (__doc__)
