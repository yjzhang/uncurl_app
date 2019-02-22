#!/bin/bash

killall -9 redis-server
nohup redis-server > nohup_redis.out &

# kill all gunicorn processes at port 8888
pids=`ps ax | grep gunicorn | grep "8888" | awk '{split($0,a," "); print a[1]}'`
for pid in $pids; do
    kill -9 $pid
    echo "killed gunicorn process $pid"
done
# 20000 second timeout
nohup gunicorn --workers 4 --threads 1 --bind 0.0.0.0:8888 --log-level debug -t 20000 wsgi:app >> nohup_uncurl.out &
