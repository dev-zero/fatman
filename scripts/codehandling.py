#!/usr/bin/env python
import numpy as np
import os, json

from ase.calculators.abinit import Abinit
from ase.calculators.cp2k import CP2K
from ase import Atoms
from ase.test.tasks import dcdft
from ase.units import Ry

import espresso

class HandlerParent():
    """
    TODO: docstring
    """
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




class abinitHandler(HandlerParent):
    #THIS IS CURRENTLY UNTESTED CODE, MAKE SURE IT WORKS!
    """TODO: docstring here"""
    def __init__(self,structure, settings = {}, workdir_prefix = "/data/ralph/deltatests/abinit"):
        self.workdir_prefix = workdir_prefix
        self.settings = settings
        self.structure = structure

    def getCalculator(self):
        """
        TODO: docstring
        """

        if "tolsym" in self.settings["settings"].keys():
            tolsym = self.settings["settings"]["tolsym"]
        else:
            tolsym = 1.0e-12

        if "fband" in self.settings["settings"].keys():
            fband = self.settings["settings"]["fband"]
        else:
            fband = 1.5

        if "xc" in self.settings["settings"].keys():
            xc = self.settings["settings"]["xc"]
        else:
            xc = 'PBE'

        if "ecut" in self.settings["settings"].keys():
            ecut = self.settings["settings"]["cutoff"]
        else:
            ecut = 100.

        if "ecutsm" in self.settings["settings"].keys():
            ecutsm = self.settings["settings"]["ecutsm"]
        else:
            ecutsm = 0.0

        if "tsmear" in self.settings["settings"].keys():
            tsmear = self.settings["settings"]["tsmear"]
        else:
            tsmear = 0.01

        if "kpoints" in self.settings.keys():
            kpts = list(self.settings["kpoints"])
        else:
            kpts = [12,12,12]

        pseudopotential = self.settings["pseudopotential"]

        if pseudopotential=="hgh.k.nlcc2015":
            os.environ['ABINIT_PP_PATH'] = '/users/ralph/work/abinit/pseudos/HGH-NLCC2015'
            pseudopotential="hgh.k"
        elif pseudopotential=="hgh.k.nlcc":
            os.environ['ABINIT_PP_PATH'] = '/users/ralph/work/abinit/pseudos/HGH-NLCC'
            pseudopotential="hgh.k"
        else:
            os.environ['ABINIT_PP_PATH'] = '/users/ralph/work/abinit/pseudos/HGH'

        if not 'ABINIT_ASE_ABINIT_COMMAND' in os.environ.keys():
            os.environ['ASE_ABINIT_COMMAND'] = 'module load gcc-suite >/dev/null; nice -n 2 mpirun -np 4 /users/ralph/work/abinit/abinit-7.10.5/tmp9/abinit < PREFIX.files > test.out'


        calc = Abinit(label='test',
                      pps=pseudopotential, # uses highest valence hgh.k pps
                      xc='PBE',
                      kpts=kpts,
                      ecut=ecut,
                      occopt=3,
                      tsmear=tsmear,
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

    def createOne(self):
        print "NOT IMPLEMENTED"
        quit()


class cp2kHandler(HandlerParent):
    """TODO: docstring here"""
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
        TODO: UPDATE Docstring
        """
        ecut  = self.settings["settings"]["cutoff_rho"]
        kind_settings = self.settings["kind_settings"]

        if 'pseudopotential' in self.settings.keys() and self.settings["pseudopotential"]=="ALL":
            pseudo = "ALL"
        else:
            pseudo = False

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

        if "multiplicity" in self.settings.keys():
            multiplicity = self.settings["multiplicity"]
        else:
            multiplicity = 1

        if "charge" in self.settings.keys():
            charge = self.settings["charge"]
        else:
            charge = 0

        magmoms = self.structure.get_initial_magnetic_moments()

        inp_template = self.input_template

        calc = CP2K (
                label          = "test",
                command        = "module load gcc-suite/5.3.0; mpirun -np 4 /users/ralph/cp2k/cp2k-code/exe/Linux-x86-64-gfortran-local/cp2k_shell.popt",
                kind_settings  = kind_settings,
                basis_set_file = False,
                potential_file = False,
                basis_set      = False,
                pseudo_potential = pseudo,
                kpoints        = kpoints,
                xc             ="PBE",
                cutoff         = ecut,
                max_scf        = 200,
                inp            = inp_template,
                uks            = magmoms.any() or (multiplicity>1),
                qs_settings    = qs_settings,
                rel_settings   = rel_settings,
                charge         = charge,
                multiplicity   = multiplicity
                )

        return calc


    def createOne(self):
        print "NOT IMPLEMENTED"
        quit()



class espressoHandler():
    def __init__(self,structure, settings = {}, workdir_prefix = "/data/ralph/deltatests/espresso"):
        self.workdir_prefix = workdir_prefix
        self.settings = settings
        self.structure = structure

    def getCalculator(self, input_only = None, outdir = "./"):
        """
        Set the required environment variables to run QuantumEspresso (unless already set).
        Return a calculator object with all required parameters set to high accuracy
        """
        
        os.popen("module load quantum-espresso 2> /dev/null")
        os.popen("module load gcc-suite/5.3.0 2> /dev/null")

        ecut     = self.settings["settings"]["cutoff_pw"]
        xc       = self.settings["settings"]["xc"]
        sigma    = self.settings["settings"]["sigma"]
        smearing = str(self.settings["settings"]["smearing"])
        kpts     = list(self.settings["kpoints"])
        pseudopotential = self.settings["pseudopotential"]

        calc = espresso.espresso(
                        mode="scf",
                        pw=ecut,                    #density cutoff automatically 4*ecut ?
                        kpts=kpts,
                        xc='PBE',
                        outdir=outdir,
                        smearing=smearing,
                        sigma=sigma,
                        psppath = "/users/ralph/work/espresso/" + pseudopotential ,
                        output = {'removewf': True, 'removesave': True},
                        parflags='-npool 4',
                        onlycreatepwinp = input_only
                        )
            #            onlycreatepwinp=True,
            #            outdir="./pwinp")
        return calc

    def runOne(self):
        struct = self.structure
        identifier = struct.info["key_value_pairs"]["identifier"]

        workdir = os.path.join(self.workdir_prefix, "method_{:04d}/{:s}".format(self.settings["id"], identifier))
        deltacalc = self.getCalculator()

        try:
            os.makedirs(workdir)
        except OSError:
            pass

        deltacalc.outdir = workdir

        if np.sum(abs(struct.get_initial_magnetic_moments()))>0:
            deltacalc.spinpol = True
        os.chdir(workdir)

        struct.set_calculator(deltacalc)

        e = struct.get_potential_energy()
        return e, os.path.join(workdir, "log")


    def createOne(self):
        struct = self.structure
        identifier = struct.info["key_value_pairs"]["identifier"]

        workdir = os.path.join("/tmp/espresso-remote-task")
        deltacalc = self.getCalculator(input_only=True)

        try:
            os.makedirs(workdir)
        except OSError:
            pass

        deltacalc.outdir = workdir
        deltacalc.pwinp = os.path.join(workdir,"pw.inp")

        struct.set_calculator(deltacalc)
        deltacalc.initialize(struct) 

        return workdir



def HandlerFactory(structure, methodsettings = {"code":"cp2k"}):
    if methodsettings["code"] =="abinit":
        return abinitHandler(structure, settings = methodsettings)

    elif methodsettings["code"] =="cp2k":
        return cp2kHandler(structure, settings = methodsettings)

    elif methodsettings["code"] =="espresso":
        return espressoHandler(structure, settings = methodsettings)

    else:
        raise RuntimeError ("Asked for unknown code")

if __name__=="__main__":
    print("DONT RUN ME!")
