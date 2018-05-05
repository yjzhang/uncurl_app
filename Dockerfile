FROM python:2.7

WORKDIR /app

ADD . /app

RUN apt-get update && apt-get install -y redis-server

RUN pip install cython
RUN pip install numpy

RUN pip install -r requirements.txt

EXPOSE 8888

ENV NAME World

CMD ["sh", "start_gunicorn_docker.sh"]
