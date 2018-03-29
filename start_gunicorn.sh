#!/bin/bash

killall -9 gunicorn
nohup gunicorn --workers 1 --threads 1 --bind 0.0.0.0:8888 --log-level debug -t 100 wsgi:app &
