#!/usr/bin/env python
"""codehandling.py
Provides functionality to run computational codes with input data defined in FATMAN.

Get a universally usable codehandler by calling HandlerFactory(), supplying the
structure for which the calculation is to be run and the method settings. The returned
class then provides the methods runOne() and createOne(), which either execute the code
and return the energy and name of the output file, or create an input file and return
its location.

Public Methods:
    - HandlerFactory(structure, settings)
      Returns a codehandler object that provides method to run calculations or create inputs.

Classes:
    - HandlerParent
      Base class for all the codehandler objects. Provides a generic runOne() method to be inherited

    - AbinitHandler
      Handler for the abinit code package, uses ASE to interface with the code.

    - Cp2kHandler
      Handler for the cp2k package, uses ASE to interface with the code.

    - EspressoHandler
      Handler for Quantum Espresso, uses the ase-espresso module (https://github.com/vossjo/ase-espresso/).
"""

import numpy as np
import os
import json

from ase.calculators.abinit import Abinit
from ase.calculators.cp2k import CP2K
from ase import Atoms
from ase.test.tasks import dcdft
from ase.units import Ry

import espresso


class HandlerParent():
    """Base class for the codehandler.

    Currently only provides the runOne() method which all derived classes inherit.
    """
    def __init__(self):
        self.output_filename = "test.out"

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

        return e, os.path.join(workdir, self.output_filename)


# THIS IS CURRENTLY UNTESTED CODE, MAKE SURE IT WORKS!
class AbinitHandler(HandlerParent):
    """Class to run the Abinit code with FATMAN data through ASE.

    Instantiated by the HandlerFactory for a particular ASE chemical structure
    and a dictionary of settings supplied from FATMAN. The workdir_prefix defines a
    writable 'scratch' directory, in which subdirectories will be created for the
    particular calculation.

    Private methods:
        - getCalculator()
          return an ASE Calculator object with the correct settings.

    Public methods:
        - runOne()
          inherited; create a calculator, run the job and return total energy and output path
        - createOne()
          create a calculator and make it write an input file. NOT IMPLEMENTED YET.
    """

    def __init__(self, structure, settings={}, workdir_prefix="/data/ralph/deltatests/abinit"):
        self.workdir_prefix = workdir_prefix
        self.settings = settings
        self.structure = structure
        self.output_filename = "PREFIX.out"

    def getCalculator(self):
        """Mangle the supplied settings and return an Abinit ASE Calculator object."""

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
            kpts = [12, 12, 12]

        pseudopotential = self.settings["pseudopotential"]

        if pseudopotential == "hgh.k.nlcc2015":
            os.environ['ABINIT_PP_PATH'] = '/users/ralph/work/abinit/pseudos/HGH-NLCC2015'
            pseudopotential = "hgh.k"
        elif pseudopotential == "hgh.k.nlcc":
            os.environ['ABINIT_PP_PATH'] = '/users/ralph/work/abinit/pseudos/HGH-NLCC'
            pseudopotential = "hgh.k"
        else:
            os.environ['ABINIT_PP_PATH'] = '/users/ralph/work/abinit/pseudos/HGH'

        if 'ABINIT_ASE_ABINIT_COMMAND' not in os.environ.keys():
            os.environ['ASE_ABINIT_COMMAND'] = 'module load gcc-suite >/dev/null; nice -n 2 mpirun -np 4 /users/ralph/work/abinit/abinit-7.10.5/tmp9/abinit < PREFIX.files > PREFIX.out'

        calc = Abinit(label='test',
                      pps=pseudopotential,  # uses highest valence hgh.k pps
                      xc='PBE',
                      kpts=kpts,
                      ecut=ecut,
                      occopt=3,
                      tsmear=tsmear,
                      ecutsm=ecutsm,
                      toldfe=1.0e-6,
                      nstep=900,
                      pawovlp=-1,  # bypass overlap check
                      fband=fband,
                      # http://forum.abinit.org/viewtopic.php?f=8&t=35
                      chksymbreak=0,
                      tolsym=tolsym,
                      prtwf=0,
                      prtden=0,
                      )
        return calc

    def createOne(self):
        raise NotImplementedError("This method is not implemented")


class Cp2kHandler(HandlerParent):
    """Class to run CP2K with FATMAN data through ASE.

    Instantiated by the HandlerFactory for a particular ASE chemical structure
    and a dictionary of settings supplied from FATMAN. The workdir_prefix defines a
    writable 'scratch' directory, in which subdirectories will be created for the
    particular calculation.

    The CP2K ASE interface is somewhat bare, so we have to provide an input template
    that will be extended with specific settings and the system definition. Stored
    as a string attribute of the class.

    Private methods:
        - getCalculator()
          return an ASE Calculator object with the correct settings.

    Public methods:
        - runOne()
          inherited; create a calculator, run the job and return total energy and output path
        - createOne()
          create a calculator and make it write an input file. NOT IMPLEMENTED YET.
    """

    def __init__(self, structure, settings={}, workdir_prefix="/data/ralph/deltatests/cp2k"):
        self.workdir_prefix = workdir_prefix
        self.settings = settings
        self.structure = structure
        self.output_filename = "test.out"
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
        """Mangle the supplied settings and return an Abinit ASE Calculator object."""

        ecut = self.settings["settings"]["cutoff_rho"]
        kind_settings = self.settings["kind_settings"]

        if 'pseudopotential' in self.settings.keys() and self.settings["pseudopotential"] == "ALL":
            pseudo = "ALL"
        else:
            pseudo = False

        if "kpoints" in self.settings.keys():
            kpoints = list(self.settings["kpoints"])
            if "kpoint_shift" in self.settings["settings"].keys():     # apply k-point origin shift
                kpoints.extend(self.settings["settings"]["kpoint_shift"])
        else:
            kpoints = None

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

        calc = CP2K(
                label = "test",
                command = "module load gcc-suite/5.3.0; mpirun -np 4 /users/ralph/cp2k/cp2k-code/exe/Linux-x86-64-gfortran-local/cp2k_shell.popt",
                kind_settings = kind_settings,
                basis_set_file = False,
                potential_file = False,
                basis_set = False,
                pseudo_potential = pseudo,
                kpoints = kpoints,
                xc ="PBE",
                cutoff = ecut,
                max_scf = 200,
                inp = inp_template,
                uks = magmoms.any() or (multiplicity > 1),
                qs_settings = qs_settings,
                rel_settings = rel_settings,
                charge = charge,
                multiplicity = multiplicity
                )

        return calc

    def createOne(self):
        """Attach the cp2k calculator to the structure, write an input file and return its path.
        """
        raise NotImplementedError("This method is not implemented")


class EspressoHandler():
    """Class to run Quantum Espresso with FATMAN data through ASE and the ASE-Espresso extension.

    Instantiated by the HandlerFactory for a particular ASE chemical structure
    and a dictionary of settings supplied from FATMAN. The workdir_prefix defines a
    writable 'scratch' directory, in which subdirectories will be created for the
    particular calculation.

    The ASE Espresso interface is a demanding princess. Needs to be installed separately, on top
    of ASE, and configured through a system-specific espsite.py where settings are hardcoded.
    In particular a scratch directory has to be defined, which espresso will fill with a lot of garbage.

    ========
    CAUTION:
    ========
      - Each job leaves data in its scratch directory, which can quickly grow to 10s-100s of GiB.
      - The pseudopotential has to be written to disk for Q.E. to find/use it. Avoid collisions!

    Private methods:
        - getCalculator()
          return an ASE Calculator object with the correct settings. Optional parameter to only create .inp.

    Public methods:
        - runOne()
          inherited; create a calculator, run the job and return total energy and output path
        - createOne()
          create a calculator and make it write an input file.
    """

    def __init__(self, structure, settings={}, workdir_prefix="/data/ralph/deltatests/espresso"):
        self.workdir_prefix = workdir_prefix
        self.settings = settings
        self.structure = structure
        self.output_filename = "log"

    def getCalculator(self, input_only=None):
        """Mangle the supplied settings and return an Abinit ASE Calculator object.

        The pseudopotential path is semi-hardcoded. Because crap.

        Arguments:
            - input_only: set to True if you want to only write a Q.E. input file, not calculate.
        """

        os.popen("module load quantum-espresso 2> /dev/null")
        os.popen("module load gcc-suite/5.3.0 2> /dev/null")

        ecut = self.settings["settings"]["cutoff_pw"]
        xc = self.settings["settings"]["xc"]
        sigma = self.settings["settings"]["sigma"]
        smearing = str(self.settings["settings"]["smearing"])
        kpts = list(self.settings["kpoints"])
        pseudopotential = self.settings["pseudopotential"]

        calc = espresso.espresso(
                        mode="scf",
                        pw=ecut,                    # density cutoff automatically 4*ecut ?
                        kpts=kpts,
                        xc='PBE',
                        outdir='./',
                        smearing=smearing,
                        sigma=sigma,
                        psppath="/users/ralph/work/espresso/" + pseudopotential,
                        output={'removewf': True, 'removesave': True},
                        parflags='-npool 4',
                        onlycreatepwinp=input_only
                        )
        return calc

    def runOne(self):
        """Attach the ASE Calculator (espresso) to the structure, run the energy calc, and return its results.

        Create working directory for the calculation based on the workdir_prefix, get the espresso calculator,
        set it to work in that directory and run the energy calculation.
        If the structure specifies magnetic moments, the option `spinpol` has to be set for the calculator.

        Returns:
            - e: the total energy from the DFT calculation in eV.
            - the path of the output file (`log`) which should then eventually be uploaded to the server.
        """
        struct = self.structure
        identifier = struct.info["key_value_pairs"]["identifier"]

        workdir = os.path.join(self.workdir_prefix, "method_{:04d}/{:s}".format(self.settings["id"], identifier))
        deltacalc = self.getCalculator()

        try:
            os.makedirs(workdir)
        except OSError:
            pass

        deltacalc.outdir = workdir

        if np.sum(abs(struct.get_initial_magnetic_moments())) > 0:
            deltacalc.spinpol = True
        os.chdir(workdir)

        struct.set_calculator(deltacalc)

        e = struct.get_potential_energy()
        return e, os.path.join(workdir, "log")

    def createOne(self):
        """Attach the ASE Calculator (espresso) to the structure, write a Q.E. input file and return its path.

        Create output directory for the input file base on workdir_prefix, get the espresso calculator set to
        only produce an input file, write that input file and return the pathname where it was created.
        If the structure specifies magnetic moments, the option `spinpol` has to be set for the calculator.

        Returns:
            - the path of the input file (`pw.inp`) for Q. Espresso.
        """

        struct = self.structure
        identifier = struct.info["key_value_pairs"]["identifier"]

        workdir = os.path.join("/tmp/espresso-remote-task")
        deltacalc = self.getCalculator(input_only=True)

        try:
            os.makedirs(workdir)
        except OSError:
            pass

        deltacalc.outdir = workdir
        deltacalc.pwinp = os.path.join(workdir, "pw.inp")

        if np.sum(abs(struct.get_initial_magnetic_moments())) > 0:
            deltacalc.spinpol = True

        struct.set_calculator(deltacalc)
        deltacalc.initialize(struct)

        return deltacalc.pwinp


def HandlerFactory(structure, methodsettings={"code": "cp2k"}):
    """Depending on the selected code, return the appropriate handler object"""

    if methodsettings["code"] == "abinit":
        return AbinitHandler(structure, settings=methodsettings)

    elif methodsettings["code"] == "cp2k":
        return Cp2kHandler(structure, settings=methodsettings)

    elif methodsettings["code"] == "espresso":
        return EspressoHandler(structure, settings=methodsettings)

    else:
        raise RuntimeError("Asked for unknown code")

if __name__ == "__main__":
    print("DONT RUN ME, I'M ONLY A LIBRARY!")
