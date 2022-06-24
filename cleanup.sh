#!/bin/bash

find /tmp/uncurl -maxdepth 1 -type f -mtime +30 -exec rm {} +

# delete directories that do not have an sc_analysis.json file that are older than 30 days
# i.e. the analysis was not run or was already finished
find /tmp/uncurl -maxdepth 1 -mtime +30 -type d '!' -exec test -e "{}/sc_analysis.json" ';' -print | xargs rm -rf

# delete all test upload files
ls -d /tmp/uncurl/*-yjz-test-upload | xargs rm -rf
