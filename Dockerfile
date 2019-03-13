FROM python:3.6
MAINTAINER Yue Zhang <yjzhang@cs.washington.edu>

RUN apt-get update && apt-get install -y redis-server && apt-get install -y libhdf5-dev
RUN pip install cython
RUN pip install numpy

WORKDIR /app

ADD requirements.txt /app

RUN pip install -r requirements.txt

ADD . /app

EXPOSE 8888

ENV NAME World

CMD ["sh", "start_gunicorn_docker.sh"]
