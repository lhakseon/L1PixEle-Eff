def makeBatchConfigFile( job_dir):

    config='#!/bin/sh'
    config+='\n'
    config+='#$ -S /bin/bash \n'
    config+='cd /share/apps/root_v5-34-32/root/ \n'
    config+='. bin/thisroot.sh \n'
    config+='cd ' + job_dir +  ' \n'
    config+='rm -rf test_C.so \n'
    config+='rm -rf test_C.d \n'
    config+='rm -rf result.log \n'
    config+='root -l -b < x_test.C \n'
    config+='rm -rf test_C.so \n'
    config+='rm -rf test_C.d \n'
    return config
