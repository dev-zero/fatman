#!/usr/bin/env python3

from fatman.models import Calculation, TestResult2

for calc in Calculation.query:
    if len(calc.testresults) > 1:
        print("Calculation {c.id} ({c.collection.name}) linked to:".format(c=calc))
        for tres in calc.testresults:
            print("    TestResult {r.id} in {r.collections}".format(r=tres))

