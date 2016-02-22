#!/usr/bin/python
import numpy as np
import os, json

from ase.calculators.abinit import Abinit
from ase.calculators.cp2k import CP2K
from ase import Atoms
from ase.test.tasks import dcdft
from ase.units import Ry

class abinitHandler():
    def __init__(self):
        _deltaRunner.workdir_prefix = "/data/ralph/deltatests/abinit"

    def getDeltaCalc(self, ecut, kpts=[16,16,16],pps='hgh.k', basis="none"):
        """
        Set the required environment variables to run abinit (unless already set).
        Return a calculator object with all required parameters set to high accuracy
        """
        if pps=="hgh.k.nlcc2015":
            os.environ['ABINIT_PP_PATH'] = '/users/ralph/work/abinit/pseudos/HGH-NLCC2015'
            pps="hgh.k"
        elif pps=="hgh.k.nlcc":
            os.environ['ABINIT_PP_PATH'] = '/users/ralph/work/abinit/pseudos/HGH-NLCC'
            pps="hgh.k"
        else:
            os.environ['ABINIT_PP_PATH'] = '/users/ralph/work/abinit/pseudos/HGH'

        if not 'ABINIT_ASE_ABINIT_COMMAND' in os.environ.keys():
            os.environ['ASE_ABINIT_COMMAND'] = 'module load gcc-suite >/dev/null; nice -n 2 mpirun -np 4 /users/ralph/work/abinit/abinit-7.10.5/tmp9/abinit < PREFIX.files > PREFIX.log'


        kptdensity = 16.0  #in the original, MP kpts are generated somehow automatically, we just use what they have determined as converged
        # this is converged
        width = 0.01
        ecutsm = 0.0
        fband = 1.5
        tolsym = 1.e-12


        calc = Abinit(label='deltatest',
                      pps=pps, # uses highest valence hgh.k pps
                      xc='PBE',
                      kpts=kpts,
                      ecut=ecut,
                      occopt=3,
                      tsmear=width,
                      ecutsm=ecutsm,
                      toldfe=1.0e-6,
                      nstep=900,
                      pawovlp=-1, # bypass overlap check
                      fband=fband,
                      # http://forum.abinit.org/viewtopic.php?f=8&t=35
                      chksymbreak=0,
                      tolsym=tolsym,
                      prtwf=0,
                      prtden=0,
                  )
        return calc


class cp2kHandler():
    def __init__(self,structure, settings = {}, workdir_prefix = "/data/ralph/deltatests/cp2k"):
        self.workdir_prefix = workdir_prefix
        self.settings = settings
        self.structure = structure
        self.input_template = """
            @SET DO_GAPW {do_gapw:s}
            &GLOBAL
            &END GLOBAL
            &FORCE_EVAL
               &DFT
                   &QS
                    @IF ${{DO_GAPW}}==TRUE
                      METHOD       GAPW
                     !EPS_GVG      1.0E-8
                     !EPS_PGF_ORB  1.0E-8
                     !QUADRATURE   GC_LOG
                     !EPSFIT       1.E-4
                     !EPSISO       1.0E-12
                     !EPSRHO0      1.E-8
                     !LMAXN0       2
                     !LMAXN1       6
                     !ALPHA0_H     10
                    @ENDIF
             
                      EXTRAPOLATION USE_GUESS
                      EPS_DEFAULT 1.0E-10
                   &END QS
                   &KPOINTS
                        FULL_GRID .TRUE.
                        SYMMETRY .FALSE.
                        SCHEME MACDONALD {kpts1:d} {kpts2:d} {kpts3:d} 0.5 0.5 0.5
                        PARALLEL_GROUP_SIZE 2
                   &END KPOINTS
                   &MGRID
                        REL_CUTOFF [Ry] 100
                   &END MGRID
                   &SCF
                      SCF_GUESS  ATOMIC
                      ADDED_MOS  50
                      CHOLESKY OFF
                      &SMEAR  ON
                          METHOD FERMI_DIRAC
                          ELECTRONIC_TEMPERATURE [K] 30
                      &END SMEAR
                      &DIAGONALIZATION
                           ALGORITHM STANDARD
                      &END DIAGONALIZATION
                      &MIXING
                           METHOD BROYDEN_MIXING
                           ALPHA   0.1
                           BETA    0.5
                           NBROYDEN  8
                      &END
                   &END SCF
               &END DFT
            &END FORCE_EVAL
        """

    def getCalculator(self):
        """
        Set the required environment variables to run CP2K (unless already set).
        Return a calculator object with all required parameters set to high accuracy
        """
        pps = self.settings["pseudopotential"]
        basis = self.settings["basis_set"]
        ecut  = self.settings["settings"]["cutoff_rho"]
        kpts  = self.settings["kpoints"]

        magmoms = self.structure.get_initial_magnetic_moments()

        #an ugly hack because the cp2k ase calculator can't set these parameters, we have to do a search/replace in an input template
        formatparams = {"kpts1":kpts[0],"kpts2": kpts[1],"kpts3": kpts[2], "do_gapw": "FALSE"}

        if pps=="GTH-PBE":
            pp_file = "/users/ralph/cp2k/cp2k-code/data/GTH_POTENTIALS"
        elif pps=="GTH-NLCC-PBE":
            pp_file = "/users/ralph/cp2k/cp2k-code/data/NLCC_POTENTIALS"
        elif pps=="GTH-NLCC2015-PBE":
            pp_file = "/users/ralph/cp2k/cp2k-code/data/NLCC2015_POTENTIALS"
        elif pps=="ALL":
            #deal with GAPW calculations
            pp_file = "/users/ralph/cp2k/cp2k-code/data/POTENTIAL"
            formatparams["do_gapw"] = "TRUE"

        inp_template = self.input_template.format(**formatparams)

        calc = CP2K (
                label          = "test",
                command        = "module load gcc-suite; mpirun -np 4 /users/ralph/cp2k/cp2k-code/exe/Linux-x86-64-gfortran-local/cp2k_shell.popt",
                basis_set_file = "/users/ralph/cp2k/cp2k-code/data/BASIS_SETS_inc_MOLOPT",
                potential_file = pp_file,
                basis_set      = basis,
                pseudo_potential = pps,
                xc             ="PBE",
                cutoff         = ecut,
                max_scf        = 200,
                inp            = inp_template,
                uks            = magmoms.any(),
                debug=True,
                )

        return calc


    def runOne(self):
        struct = self.structure
        identifier = struct.info["key_value_pairs"]["identifier"]

        workdir = os.path.join(self.workdir_prefix, "method_{:04d}/{:s}".format(self.settings["id"], identifier))
        print workdir
#       try:
#           os.makedirs(workdir)
#       except OSError:
#           pass

#       os.chdir(workdir)

        deltacalc = self.getCalculator()


        struct.set_calculator(deltacalc)

        print deltacalc._generate_input()

#       try:
#           e = struct.get_potential_energy()
#       except RuntimeError:
#           return 0,0.

        return 0

def HandlerFactory(structure, methodsettings = {"code":"cp2k"}):
    if methodsettings["code"] =="abinit":
        return abinitHandler(structure, settings = methodsettings)

    elif methodsettings["code"] =="cp2k":
        return cp2kHandler(structure, settings = methodsettings)

    else:
        raise RuntimeError ("Asked for unknown code")

if __name__=="__main__":
    print "DONT RUN ME!"

