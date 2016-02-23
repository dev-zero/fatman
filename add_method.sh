#!/bin/bash
psql fatman <<EOF

INSERT INTO method (code,pseudopotential,basis_set,settings) VALUES  ('cp2k', 'GTH-PBE', 'DZVP-MOLOPT-GTH', '{"cutoff_rho":5000, "GAPW": "false", "relativistic": "false"}');

EOF
