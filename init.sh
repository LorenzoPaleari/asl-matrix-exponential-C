#!/bin/bash

# Sets up environment variables - Done in .bashrc
# export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/2023.1.0/lib/intel64:$LD_LIBRARY_PATH

# Sets up git hooks
chmod +x git_hooks/pre-commit
chmod +x test_runner.sh
chmod +x test_runner_flop.sh
chmod +x test_runner_data_collection.sh
# git config core.hooksPath git_hooks

# Sets up advisor environment
output=$(env | grep ADVISOR)
if [ -z "$output" ]; then
    . /opt/intel/oneapi/setvars.sh
fi