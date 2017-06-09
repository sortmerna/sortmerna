
FROM ubuntu:16.04


COPY . /sortmerna

RUN apt-get update && apt-get install -y make g++ zlib1g-dev python2.7 python-pip python-numpy python-scipy python-tk

RUN cd /sortmerna && ./configure && make && make install

RUN pip install --upgrade pip && pip install scikit-bio==0.2.3
