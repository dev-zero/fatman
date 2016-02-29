#!/usr/bin/env python

def Atoms2Json(structure, additional_information={}):
    """Serialize an ASE Structure definition to JSON and return it as a string"""

    import json,os
    from ase.db.row import AtomsRow
    from ase.db.core import now
    from ase.io.jsonio import MyEncoder as AseJsonEncoder

    row = AtomsRow(structure)       #this is what ASE would store in its DB.
    row.ctime = mtime = now()       #the Row object has an attribute ctime, but not mtime -- we have to wiggle it into the dict later.
    row.user = os.getenv("USER")

    dct = row.__dict__.copy()
    del dct["_keys"], dct["_data"],dct["_constraints"]   #containing useless default entries that shouldn't be stored
    dct["mtime"] = mtime
    dct["key_value_pairs"] = additional_information 

    return json.dumps(dct, sort_keys=True, cls=AseJsonEncoder)


def Json2Atoms(jsonstring):
    """Read a JSON string and return an Atoms object"""

    from ase.io.jsonio import decode
    from ase.db.row import AtomsRow

    dct = decode(jsonstring)
    row = AtomsRow(dct)

    return row.toatoms(attach_calculator=False, add_additional_information=True)


def calcDelta(data_f, data_w):
    """
    Calculate the Delta using the data in data_f, data_w on
    element in eloverlap
    This is taken from Delta_v3-0.zip from the Ghent website [23.11.2015]
    """
    import numpy as np

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






def test1():
    from ase import Atoms
    teststructure = Atoms('N2',[(0,0,0),(0,0,1.1)])

    teststring = Atoms2Json(teststructure, additional_information={"abc":"def"}) #convert to the string

    newstructure = Json2Atoms(teststring)  #and convert back

    assert (newstructure.get_atomic_numbers() == teststructure.get_atomic_numbers()).all()
    assert (newstructure.get_positions()      == teststructure.get_positions()).all()
    assert (newstructure.get_cell()           == teststructure.get_cell()).all()
    assert (newstructure.get_pbc()            == teststructure.get_pbc()).all()
    assert  newstructure.info["key_value_pairs"]["abc"] == "def"  

    print("JSON Conversion in both directions: apparently no issue")


def test2():
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
    test1()
    test2()
