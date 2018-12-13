#!/bin/sh

wget https://github.com/mozilla/geckodriver/releases/download/v0.20.1/geckodriver-v0.20.1-linux64.tar.gz
tar -xzf geckodriver-v0.20.1-linux64.tar.gz

export PATH=$PATH:`pwd`
echo $PATH
