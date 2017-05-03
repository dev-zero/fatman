#!/usr/bin/env python3

from fatman.models import db, TestResult2

count = (TestResult2.query
         .filter(TestResult2.collections == None)
         .delete(synchronize_session=False)
         )

db.session.commit()
print("Dropped %i test results" % count)
