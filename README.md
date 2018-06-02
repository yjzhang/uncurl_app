Flask web app for UNCURL
=======================

## Installation

Requirements: see `requirements.txt`. To install requirements, run `pip install -r requirements.txt`

To deploy locally: run `python run.py`

To deploy on a server (requires redis): run `sh start-gunicorn.sh`

Main code is located in `app/`

## Building with Docker

(still working on this...)

`docker build . -t uncurl-app`

`docker run uncurl-app -p 8888:<port>`

This exposes the given port.

## Testing

`python test_app.py`

`python test_frontend.py` (requires selenium and the Firefox webdriver from https://github.com/mozilla/geckodriver/releases)

## Usage

See `user_guide.md`


## Included resources

loading gif from http://www.ajaxload.info/

Includes plotly, jquery, poppler
