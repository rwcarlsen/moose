#!/bin/bash

APP='phase_field'
APP_DIR='/Users/calsrw/moose/modules/phase_field'
TEST_DIR='examples/grain_growth'
APP_FLAGS='-i grain_growth_2D_graintracker.i Executioner/end_time=110'

METHOD='opt'
N=10

# ./run_tests --oprof --re ifba_he_production.base_ifba_only_smeared --heavy 

for i in $(seq $N); do
    (cd "$APP_DIR/$TEST_DIR" && time ${APP_DIR}/${APP}-${METHOD} ${APP_FLAGS}) 2>&1 | tail -n3 | awk 'NR==1 {print "";} {printf " %s", $2}'
done

#sed -e 's/[ms]/ /g' | awk '{print $1*60+$2" "$3*60+$4" "$5*60+$6}'

