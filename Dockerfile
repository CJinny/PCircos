FROM ubuntu:latest
LABEL org.opencontainers.image.authors="antony.lebechec@chru-strasbourg.fr"


RUN apt-get update -y
RUN apt-get install -y python3.6 python3-pip

RUN mkdir /app
WORKDIR /app
ADD requirements.txt /app/

RUN pip3 install --upgrade pip
ADD . /app/
RUN python3 -m pip install -e .

ENTRYPOINT ["vcf2circos"]
