#!/bin/bash

killall -9 redis-server
nohup redis-server > nohup_redis.out &

killall -9 gunicorn
gunicorn --workers 4 --threads 1 --bind 0.0.0.0:8888 --log-level debug -t 20000 wsgi:app
