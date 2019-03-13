FROM python:3.6
MAINTAINER Yue Zhang <yjzhang@cs.washington.edu>

RUN apt-get update && apt-get install -y redis-server
RUN pip install cython
RUN pip install numpy

WORKDIR /app

ADD . /app

RUN pip install -r requirements.txt

EXPOSE 8888

ENV NAME World

CMD ["sh", "start_gunicorn_docker.sh"]
