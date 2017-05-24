
import numpy as np

from ase.eos import EquationOfState
from ase.units import kJ  # pylint: disable=locally-disabled, no-name-in-module

from . import energies_for_test

ATOMIC_ELEMENTS = {
    "Cu": { "sym": "Cu", "name": "copper", "num": 29 },
    "Pt": { "sym": "Pt", "name": "platinum", "num": 78 },
    "Hf": { "sym": "Hf", "name": "hafnium", "num": 72 },
    "Cl": { "sym": "Cl", "name": "chlorine", "num": 17 },
    "Ir": { "sym": "Ir", "name": "iridium", "num": 77 },
    "Hg": { "sym": "Hg", "name": "mercury", "num": 80 },
    "Ag": { "sym": "Ag", "name": "silver", "num": 47 },
    "Ge": { "sym": "Ge", "name": "germanium", "num": 32 },
    "Li": { "sym": "Li", "name": "lithium", "num": 3 },
    "Mn": { "sym": "Mn", "name": "manganese", "num": 25 },
    "K":  { "sym": "K",  "name": "potassium", "num": 19 },
    "Be": { "sym": "Be", "name": "beryllium", "num": 4 },
    "Rn": { "sym": "Rn", "name": "radon", "num": 86 },
    "Md": { "sym": "Md", "name": "mendelevium", "num": 101 },
    "Pm": { "sym": "Pm", "name": "promethium", "num": 61 },
    "Dy": { "sym": "Dy", "name": "dysprosium", "num": 66 },
    "Np": { "sym": "Np", "name": "neptunium", "num": 93 },
    "W":  { "sym": "W",  "name": "tungsten", "num": 74 },
    "I":  { "sym": "I",  "name": "iodine", "num": 53 },
    "Tm": { "sym": "Tm", "name": "thulium", "num": 69 },
    "Nb": { "sym": "Nb", "name": "niobium", "num": 41 },
    "N":  { "sym": "N",  "name": "nitrogen", "num": 7 },
    "La": { "sym": "La", "name": "Lanthanum", "num": 57 },
    "At": { "sym": "At", "name": "astatine", "num": 85 },
    "Og": { "sym": "Og", "name": "oganesson", "num": 118 },
    "Nd": { "sym": "Nd", "name": "neodymium", "num": 60 },
    "Ar": { "sym": "Ar", "name": "argon", "num": 18 },
    "Cs": { "sym": "Cs", "name": "caesium", "num": 55 },
    "Rb": { "sym": "Rb", "name": "rubidium", "num": 37 },
    "Hs": { "sym": "Hs", "name": "hassium", "num": 108 },
    "P":  { "sym": "P",  "name": "phosphorus", "num": 15 },
    "O":  { "sym": "O",  "name": "oxygen", "num": 8 },
    "Rg": { "sym": "Rg", "name": "roentgenium", "num": 111 },
    "Mg": { "sym": "Mg", "name": "magnesium", "num": 12 },
    "Sb": { "sym": "Sb", "name": "antimony", "num": 51 },
    "Th": { "sym": "Th", "name": "thorium", "num": 90 },
    "Cd": { "sym": "Cd", "name": "cadmium", "num": 48 },
    "F":  { "sym": "F",  "name": "fluorine", "num": 9 },
    "Yb": { "sym": "Yb", "name": "ytterbium", "num": 70 },
    "Lu": { "sym": "Lu", "name": "lutetium", "num": 71 },
    "Pd": { "sym": "Pd", "name": "palladium", "num": 46 },
    "Zr": { "sym": "Zr", "name": "zirconium", "num": 40 },
    "Ca": { "sym": "Ca", "name": "calcium", "num": 20 },
    "Er": { "sym": "Er", "name": "erbium", "num": 68 },
    "Es": { "sym": "Es", "name": "einsteinium", "num": 99 },
    "Rf": { "sym": "Rf", "name": "rutherforium", "num": 104 },
    "Ds": { "sym": "Ds", "name": "darmstadtium", "num": 110 },
    "Ho": { "sym": "Ho", "name": "holmium", "num": 67 },
    "Xe": { "sym": "Xe", "name": "xenon", "num": 54 },
    "Mt": { "sym": "Mt", "name": "meitnerium", "num": 109 },
    "Kr": { "sym": "Kr", "name": "krypton", "num": 36 },
    "Se": { "sym": "Se", "name": "selenium", "num": 34 },
    "Am": { "sym": "Am", "name": "americium", "num": 95 },
    "Gd": { "sym": "Gd", "name": "gadolinium", "num": 64 },
    "Fl": { "sym": "Fl", "name": "flerovium", "num": 114 },
    "Ga": { "sym": "Ga", "name": "gallium", "num": 31 },
    "Ac": { "sym": "Ac", "name": "actinium", "num": 89 },
    "Lr": { "sym": "Lr", "name": "lawrencium", "num": 103 },
    "Pa": { "sym": "Pa", "name": "protactinium", "num": 91 },
    "Nh": { "sym": "Nh", "name": "nihonium", "num": 113 },
    "Cn": { "sym": "Cn", "name": "copernicium", "num": 112 },
    "Cf": { "sym": "Cf", "name": "californium", "num": 98 },
    "C":  { "sym": "C",  "name": "carbon", "num": 6 },
    "V":  { "sym": "V",  "name": "vandium", "num": 23 },
    "Na": { "sym": "Na", "name": "sodium", "num": 11 },
    "He": { "sym": "He", "name": "helium", "num": 2 },
    "Pr": { "sym": "Pr", "name": "praseodymium", "num": 59 },
    "Bh": { "sym": "Bh", "name": "bohrium", "num": 107 },
    "Tl": { "sym": "Tl", "name": "thalium", "num": 81 },
    "Mo": { "sym": "Mo", "name": "molybdenum", "num": 42 },
    "Bk": { "sym": "Bk", "name": "berkelium", "num": 97 },
    "Co": { "sym": "Co", "name": "cobalt", "num": 27 },
    "Sm": { "sym": "Sm", "name": "samarium", "num": 62 },
    "Pu": { "sym": "Pu", "name": "plutonium", "num": 94 },
    "Ru": { "sym": "Ru", "name": "ruthenium", "num": 44 },
    "Sn": { "sym": "Sn", "name": "tin", "num": 50 },
    "Sr": { "sym": "Sr", "name": "strontium", "num": 38 },
    "Po": { "sym": "Po", "name": "polonium", "num": 84 },
    "Rh": { "sym": "Rh", "name": "rhodium", "num": 45 },
    "No": { "sym": "No", "name": "nobelium", "num": 102 },
    "Ne": { "sym": "Ne", "name": "neon", "num": 10 },
    "S":  { "sym": "S",  "name": "sulfur", "num": 16 },
    "Br": { "sym": "Br", "name": "bromine", "num": 35 },
    "Tb": { "sym": "Tb", "name": "terbium", "num": 65 },
    "Os": { "sym": "Os", "name": "osmium", "num": 76 },
    "Sc": { "sym": "Sc", "name": "scandium", "num": 21 },
    "Al": { "sym": "Al", "name": "aluminium", "num": 13 },
    "Si": { "sym": "Si", "name": "silicon", "num": 14 },
    "As": { "sym": "As", "name": "arsenic", "num": 33 },
    "Bi": { "sym": "Bi", "name": "bismuth", "num": 83 },
    "Pb": { "sym": "Pb", "name": "lead", "num": 82 },
    "Eu": { "sym": "Eu", "name": "europium", "num": 63 },
    "Lv": { "sym": "Lv", "name": "livermorium", "num": 116 },
    "Y":  { "sym": "Y",  "name": "yttrium", "num": 39 },
    "Fr": { "sym": "Fr", "name": "francium", "num": 87 },
    "Ni": { "sym": "Ni", "name": "nickel", "num": 28 },
    "Fe": { "sym": "Fe", "name": "iron", "num": 26 },
    "Tc": { "sym": "Tc", "name": "technetium", "num": 43 },
    "Sg": { "sym": "Sg", "name": "seaborgium", "num": 106 },
    "Cr": { "sym": "Cr", "name": "chromium", "num": 24 },
    "Cm": { "sym": "Cm", "name": "curium", "num": 96 },
    "Ce": { "sym": "Ce", "name": "cerium", "num": 58 },
    "Fm": { "sym": "Fm", "name": "fermium", "num": 100 },
    "Zn": { "sym": "Zn", "name": "zinc", "num": 30 },
    "U":  { "sym": "U",  "name": "uranium", "num": 92 },
    "H":  { "sym": "H",  "name": "hydrogen", "num": 1 },
    "In": { "sym": "In", "name": "indium", "num": 49 },
    "Re": { "sym": "Re", "name": "rhenium", "num": 75 },
    "Mc": { "sym": "Mc", "name": "moscovium", "num": 115 },
    "Ra": { "sym": "Ra", "name": "radium", "num": 88 },
    "Ba": { "sym": "Ba", "name": "barium", "num": 56 },
    "Ta": { "sym": "Ta", "name": "tantalum", "num": 73 },
    "B":  { "sym": "B",  "name": "boron", "num": 5 },
    "Db": { "sym": "Db", "name": "dubnium", "num": 105 },
    "Au": { "sym": "Au", "name": "gold", "num": 79 },
    "Ts": { "sym": "Ts", "name": "tennessine", "num": 117 },
    "Ti": { "sym": "Ti", "name": "titanium", "num": 22 },
    "Te": { "sym": "Te", "name": "tellurium", "num": 52 },
  }

def eos(V0,B0, B1, E0=0.):
    B0 = B0 * 1e9 / 1.602176565e-19 / 1e30
    rng = np.linspace(0.93*V0, 1.07*V0, 40)
    E = [ E0 + 9./16. * V0 * B0 *  ( ((V0/v)**(2./3.) -1)**3 * B1 +
                                     ((V0/v)**(2./3.) -1)**2 * (6-4*(V0/v)**(2./3.)) ) for v in rng]
    return rng, np.array(E)


class FullEquationOfState(EquationOfState):

    """Fit equation of state for bulk systems.

    Salvaged from https://gitlab.com/ase/ase/blob/3.8.1/ase/test/tasks/dcdft.py
    Based on eosfit.py from http://molmod.ugent.be/DeltaCodesDFT

    """

    def __init__(self, volumes, energies, eos='birch'):
        assert eos == 'birch', eos + ' eos not available.'
        self.v = np.array(volumes)
        self.e = np.array(energies)
        self.eos_string = 'birch'

        self.v0 = None

    def fit(self):
        """Calculate volume (v0), energy (e0), bulk modulus (B0), and
        bulk modulus pressure derivative (B1).

        Returns v0, e0, B0, B1, fit residuals.

        Notice that the ASE units for the bulk modulus is
        eV/Angstrom^3 - to get the value in GPa, do this::

          v0, e0, B0, B1, R = eos.fit()
          print B0 / kJ * 1.0e24, 'GPa'

        """

        fitdata = np.polyfit(self.v**(-2./3.), self.e, 3, full=True)
        ssr = fitdata[1]
        sst = np.sum((self.e - np.average(self.e))**2.)
        residuals0 = ssr/sst
        deriv0 = np.poly1d(fitdata[0])
        deriv1 = np.polyder(deriv0, 1)
        deriv2 = np.polyder(deriv1, 1)
        deriv3 = np.polyder(deriv2, 1)

        self.v0 = None
        for x in np.roots(deriv1):
            if x > 0 and deriv2(x) > 0:
                self.v0 = x**(-3./2.)
                break

        if self.v0 is None:
            raise ValueError('No minimum!')

        derivV2 = 4./9. * x**5. * deriv2(x)
        derivV3 = (-20./9. * x**(13./2.) * deriv2(x) -
                   8./27. * x**(15./2.) * deriv3(x))
        bulk_modulus0 = derivV2 / x**(3./2.)
        bulk_deriv0 = -1 - x**(-3./2.) * derivV3 / derivV2

        self.e0 = deriv0(x)
        self.B0 = bulk_modulus0
        self.B1 = bulk_deriv0

        return self.v0, self.e0, self.B0, self.B1, residuals0


def deltatest_ev_curve(volumes, energies):
    eos = FullEquationOfState(volumes, energies)

    try:
        v, e, B0, B1, R = eos.fit()
    except ValueError:
        print("failure")
        return "fail", "fail", "fail", "fail", "fail"
    else:
        return v, e, B0/kJ * 1.0e24, B1, R[0]


def calcDelta(data_f, data_w):
    """
    Calculate the Delta using the data in data_f, data_w on
    element in eloverlap
    This is taken from Delta_v3-0.zip from the Ghent website [23.11.2015]
    """
    v0w = data_f[0]
    b0w = data_f[1] * 10.**9. / 1.602176565e-19 / 10.**30.
    b1w = data_f[2]

    v0f = data_w[0]
    b0f = data_w[1] * 10.**9. / 1.602176565e-19 / 10.**30.
    b1f = data_w[2]

    vref = 30.
    bref = 100. * 10.**9. / 1.602176565e-19 / 10.**30.

    Vi = 0.94 * (v0w + v0f) / 2.        #is this really correct? we are using smaller strain ranges
    Vf = 1.06 * (v0w + v0f) / 2.        #probably okay, since we are supplying the fitted V, B0 and B1 here,
                                        #and the integration range should not make a difference

    a3f = 9. * v0f**3. * b0f / 16. * (b1f - 4.)
    a2f = 9. * v0f**(7./3.) * b0f / 16. * (14. - 3. * b1f)
    a1f = 9. * v0f**(5./3.) * b0f / 16. * (3. * b1f - 16.)
    a0f = 9. * v0f * b0f / 16. * (6. - b1f)

    a3w = 9. * v0w**3. * b0w / 16. * (b1w - 4.)
    a2w = 9. * v0w**(7./3.) * b0w / 16. * (14. - 3. * b1w)
    a1w = 9. * v0w**(5./3.) * b0w / 16. * (3. * b1w - 16.)
    a0w = 9. * v0w * b0w / 16. * (6. - b1w)

    x = [0, 0, 0, 0, 0, 0, 0]

    x[0] = (a0f - a0w)**2
    x[1] = 6. * (a1f - a1w) * (a0f - a0w)
    x[2] = -3. * (2. * (a2f - a2w) * (a0f - a0w) + (a1f - a1w)**2.)
    x[3] = -2. * (a3f - a3w) * (a0f - a0w) - 2. * (a2f - a2w) * (a1f - a1w)
    x[4] = -3./5. * (2. * (a3f - a3w) * (a1f - a1w) + (a2f - a2w)**2.)
    x[5] = -6./7. * (a3f - a3w) * (a2f - a2w)
    x[6] = -1./3. * (a3f - a3w)**2.

    y = [0, 0, 0, 0, 0, 0, 0]

    y[0] = (a0f + a0w)**2 / 4.
    y[1] = 3. * (a1f + a1w) * (a0f + a0w) / 2.
    y[2] = -3. * (2. * (a2f + a2w) * (a0f + a0w) + (a1f + a1w)**2.) / 4.
    y[3] = -(a3f + a3w) * (a0f + a0w) / 2. - (a2f + a2w) * (a1f + a1w) / 2.
    y[4] = -3./20. * (2. * (a3f + a3w) * (a1f + a1w) + (a2f + a2w)**2.)
    y[5] = -3./14. * (a3f + a3w) * (a2f + a2w)
    y[6] = -1./12. * (a3f + a3w)**2.

    Fi = 0.
    Ff = 0.

    Gi = 0.
    Gf = 0.

    for n in range(7):
        Fi = Fi + x[n] * Vi**(-(2.*n-3.)/3.)
        Ff = Ff + x[n] * Vf**(-(2.*n-3.)/3.)

        Gi = Gi + y[n] * Vi**(-(2.*n-3.)/3.)
        Gf = Gf + y[n] * Vf**(-(2.*n-3.)/3.)

    Delta = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi))
   #Deltarel = 100. * np.sqrt((Ff - Fi) / (Gf - Gi))
   #if useasymm:
   #    Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
   #             / v0w / b0w * vref * bref
   #else:
   #    Delta1 = 1000. * np.sqrt((Ff - Fi) / (Vf - Vi)) \
   #             / (v0w + v0f) / (b0w + b0f) * 4. * vref * bref

    return Delta

def test():
    #testing the delta computation
    #          X:  v       b0         b1      REF: v         b0       b1        REQ. VALUE
    testdata = [[17.388,  10.284,  2.711,         17.388,  10.284,  2.711,      0.000 ],   #hydrogen wien2k vs. wien2k
                [17.388,  10.284,  2.711,         17.422,  10.262,  2.683,      0.07 ],    #hydrogen abinit vs. wien2k
                [36.817,  36.030,  4.637,         37.094,  34.150,  4.716,      2.10 ]]    #Sn abinit vs. wien2k

    epsilon = 0.01   #tolerance of 0.01 meV

    for line in testdata:
        data_f = line[:3]
        data_w = line[3:6]
        ref = line[6]

        delta = calcDelta(data_f, data_w)
        assert abs(delta-ref)<epsilon

    print("Accuracy of Delta-factor calculation: apparently no issue")


if __name__=="__main__":
    #if called execute "unit" tests
    test()
