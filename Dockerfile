
FROM ubuntu:16.04


COPY . /sortmerna

RUN apt-get update && apt-get install -y make g++ zlib1g-dev python2.7

RUN cd /sortmerna && ./configure && make && make install
