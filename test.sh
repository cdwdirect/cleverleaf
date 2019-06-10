export CALI_LOG_VERBOSITY=0
export CALI_AGGREGATE_KEY="APOLLO_time_flush"
export CALI_SOS_TRIGGER_ATTR="APOLLO_time_flush"
export CALI_SERVICES_ENABLE="sos,timestamp"
export CALI_TIMER_SNAPSHOT_DURATION="false"
export CALI_SOS_ITER_PER_PUBLISH="1"

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/g/g17/wood67/src/apollo/install/lib \
    ~/src/cleverleaf/apollo-test/RelWithDebInfo/install/cleverleaf/bin/cleverleaf /g/g17/wood67/src/apollo/jobs/cleaf_test.in
