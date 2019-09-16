#FROM alpine:latest
### I tried python:3.7-alpine, but no success during installation
FROM ubuntu:latest
MAINTAINER Jin Cui <cuijinjincui4@gmail.com>


RUN apt-get update -y
RUN apt-get install -y python3.6 python3-pip



RUN mkdir /app
WORKDIR /app
ADD requirements.txt /app/

RUN pip3 install --upgrade pip
RUN pip3 install -r requirements.txt
ADD . /app/


ENTRYPOINT ["python3"]
CMD ["dashapp.py"]

## docker build -t dashcircos-docker .
## docker run -it --rm -p 8000:8050 dashcircos-docker
## http://localhost:8000