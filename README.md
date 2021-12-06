UNCURL-App
=======================

Try it out: the current url is at <https://uncurl.cs.washington.edu>

## Installation

Requirements: see `requirements.txt`. To install requirements, run `pip install -r requirements.txt`

To install: run `pip install -e .` - this is necessary to use the `uncurl_app_split_seq` script (see `user_guide.md`).

To deploy locally: run `python run.py`

To deploy on a server (requires redis): run `sh start-gunicorn.sh`

Main code is located in `uncurl_app/`

By default, data is stored at `/tmp/uncurl/`.

## Building/Deploying with Docker

To build a new image from the repository's root directory:

`docker build . -t uncurl-app`

To run the server:

`docker run -p 6379:6379 -p <port>:8888 uncurl-app`

This exposes the given port, and then the uncurl-app website can be visited in the browser at http://your-ip-address:port. To stop the server, run `sudo killall gunicorn` in another terminal.

Alternatively, we have built reasonably up-to-date images at [ayuezhang27/uncurl-app](https://hub.docker.com/repository/docker/ayuezhang27/uncurl-app). To run the server using these images (does not require cloning this repository):

    docker pull ayuezhang27/uncurl-app
    docker run -p 8888:8888 -p 6379:6379  ayuezhang27/uncurl-app

When deploying on AWS, make sure to allow HTTP requests to and from the selected port in the security group.

Pushing to docker-hub:

    docker build . -t uncurl-app
    docker tag uncurl-app ayuezhang27/uncurl-app
    docker push ayuezhang27/uncurl-app


## Testing

`python test/test_app.py`

`python test/test_frontend.py` (requires selenium and the Firefox webdriver from <https://github.com/mozilla/geckodriver/releases>)


## Usage

See <http://uncurl.cs.washington.edu/help>

In order to change the file upload size limits when running locally, set the `MAX_CONTENT_LENGTH` environment variable to the desired value (in bytes). 


## Command-line usage

Most of the features of uncurl-app are provided through the [uncurl_analysis](https://github.com/yjzhang/uncurl_analysis) package. The [SCAnalysis](https://github.com/yjzhang/uncurl_analysis/blob/master/uncurl_analysis/sc_analysis.py#L33) class represents the data analysis.

Cell type databases include [cellmarker_python](https://github.com/yjzhang/cellmarker_python) and [cellmesh](https://github.com/yjzhang/cellmesh).


## Included resources

loading gif from http://www.ajaxload.info/

Javascript libraries: [plotly](https://plot.ly/), [jquery](https://jquery.com), [popper.js](https://popper.js.org/)
