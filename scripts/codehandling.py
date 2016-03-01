#!/usr/bin/env python
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
            &GLOBAL
            &END GLOBAL
            &FORCE_EVAL
               &DFT
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
        ecut  = self.settings["settings"]["cutoff_rho"]
        kind_settings = self.settings["kind_settings"]

        if "kpoints" in self.settings.keys():
            kpoints  = list(self.settings["kpoints"])
            if "kpoint_shift" in self.settings["settings"].keys():     #apply k-point origin shift
                kpoints.extend(self.settings["settings"]["kpoint_shift"])
        else:
            kpoints  = None

        if "qs_settings" in self.settings["settings"].keys():
            qs_settings = self.settings["settings"]["qs_settings"]
        else:
            qs_settings = {}

        if "rel_settings" in self.settings["settings"].keys():
            rel_settings = self.settings["settings"]["rel_settings"]
        else:
            rel_settings = {}

        magmoms = self.structure.get_initial_magnetic_moments()

        #an ugly hack because the cp2k ase calculator can't set these parameters, we have to do a search/replace in an input template
        #formatparams = {"kpts1":kpts[0],"kpts2": kpts[1],"kpts3": kpts[2], "do_gapw": "FALSE"}

        inp_template = self.input_template

        calc = CP2K (
                label          = "test",
                command        = "module load gcc-suite/5.3.0; mpirun -np 4 /users/ralph/cp2k/cp2k-code/exe/Linux-x86-64-gfortran-local/cp2k_shell.popt",
                kind_settings  = kind_settings,
                basis_set_file = False,
                potential_file = False,
                basis_set      = False,
                pseudo_potential = False,
                kpoints        = kpoints,
                xc             ="PBE",
                cutoff         = ecut,
                max_scf        = 200,
                inp            = inp_template,
                uks            = magmoms.any(),
                qs_settings    = qs_settings,
                rel_settings   = rel_settings,
                )

        return calc


    def runOne(self):
        struct = self.structure
        identifier = struct.info["key_value_pairs"]["identifier"]

        workdir = os.path.join(self.workdir_prefix, "method_{:04d}/{:s}".format(self.settings["id"], identifier))

        try:
            os.makedirs(workdir)
        except OSError:
            pass

        os.chdir(workdir)

        deltacalc = self.getCalculator()

        struct.set_calculator(deltacalc)

        e = struct.get_potential_energy()
        
        return e, os.path.join(workdir, "test.out")

def HandlerFactory(structure, methodsettings = {"code":"cp2k"}):
    if methodsettings["code"] =="abinit":
        return abinitHandler(structure, settings = methodsettings)

    elif methodsettings["code"] =="cp2k":
        return cp2kHandler(structure, settings = methodsettings)

    else:
        raise RuntimeError ("Asked for unknown code")

if __name__=="__main__":
    print("DONT RUN ME!")

