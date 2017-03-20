#!/bin/bash

BISON_DIR='/Users/calsrw/bison'
METHOD='oprof'
N=10

# ./run_tests --oprof --re ifba_he_production.base_ifba_only_smeared --heavy 
TEST_DIR='tests/ifba_he_production'
BISON_FLAGS='-i ifba_only_template.i Postprocessors/he_prod/zrb2_load=1.181e-4 Postprocessors/he_prod/ifba_len=1.0e-2 Postprocessors/he_prod/b10_enrich=0.5 Postprocessors/he_prod/zrb2_rel_dens=0.7 Postprocessors/he_prod/model=burnup Postprocessors/he_prod/u235_enrich=0.045 Postprocessors/he_prod/burnup=average_burnup'

for i in $(seq $N); do
    (cd "$BISON_DIR/$TEST_DIR" && time ${BISON_DIR}/bison-${METHOD} ${BISON_FLAGS}) 2>&1 | tail -n3 | awk 'NR==0 {print "";} {printf " %s", $2}'
done

