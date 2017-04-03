#!/bin/bash

APP='combined'
APP_DIR='/Users/calsrw/moose/modules/combined'
TEST_DIR='tests/contact_verification/hertz_cyl/half_symm_q4' # relative to APP_DIR
APP_FLAGS='-i hertz_cyl_half_1deg_template1_sm.i Contact/interface/model=frictionless Contact/interface/formulation=kinematic Outputs/file_base=hertz_cyl_half_1deg_frictionless_kin_out Outputs/chkfile/file_base=hertz_cyl_half_1deg_frictionless_kin_check Outputs/chkfile2/file_base=hertz_cyl_half_1deg_frictionless_kin_check2 Outputs/chkfile/start_time=3.49'

METHOD='oprof'
N=10

# ./run_tests --oprof --re ifba_he_production.base_ifba_only_smeared --heavy 

for i in $(seq $N); do
    (cd "$APP_DIR/$TEST_DIR" && time ${APP_DIR}/${APP}-${METHOD} ${APP_FLAGS}) 2>&1 | tail -n3 | awk 'NR==1 {print "";} {printf " %s", $2}'
done

#sed -e 's/[ms]/ /g' | awk '{print $1*60+$2" "$3*60+$4" "$5*60+$6}'

