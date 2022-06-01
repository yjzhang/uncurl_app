FROM python:3.8
MAINTAINER Yue Zhang <yjzhang@cs.washington.edu>

RUN apt-get update && apt-get install -y redis-server && apt-get install -y libhdf5-dev && apt-get install -y cmake && apt-get install -y nginx
RUN pip install cython==0.29.27
RUN pip install numpy==1.19.4

WORKDIR /app

ADD requirements.txt /app

RUN pip install -r requirements.txt

ADD . /app

EXPOSE 8888

ENV NAME World

ENV MAX_CONTENT_LENGTH 9999999999
ENV SHOW_ALL_RESULTS true

CMD ["sh", "start_gunicorn_docker.sh"]
