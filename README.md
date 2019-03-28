Flask web app for UNCURL
=======================

Try it out: the current url is at http://uncurl-app.yjzhang.com:8888/

## Installation

Requirements: see `requirements.txt`. To install requirements, run `pip install -r requirements.txt`

To deploy locally: run `python run.py`

To deploy on a server (requires redis): run `sh start-gunicorn.sh`

Main code is located in `app/`

## Building/Deploying with Docker

To build a new image from the repository's root directory:

`docker build . -t uncurl-app`

To run the server:

`docker run uncurl-app -p 8888:<port>`

This exposes the given port, and then the uncurl-app website can be visited in the browser at http://your-ip-address:port. To stop the server, run `sudo killall gunicorn` in another terminal.

Alternatively, we have built reasonably up-to-date images at ayuezhang27/uncurl-app. To run the server using these images (does not require cloning this repository):

    docker pull ayuezhang27/uncurl-app
    docker run -p 8888:8888 -p 6379:6379  ayuezhang27/uncurl-app

When deploying on AWS, make sure to allow HTTP requests to and from the selected port in the security group.


## Testing

`python test_app.py`

`python test_frontend.py` (requires selenium and the Firefox webdriver from https://github.com/mozilla/geckodriver/releases)

## Usage

See `user_guide.md`, or http://uncurl-app.yjzhang.com:8888/help


## Included resources

loading gif from http://www.ajaxload.info/

Javascript libraries: plotly, jquery, poppler
