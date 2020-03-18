#!/bin/bash

perl BadiRate.pl -treefile genome.root.nwk -sizefile Orthogroups.tbl -anc -bmodel FR -rmodel BDI -ep CSP > badirate_CSP.out
