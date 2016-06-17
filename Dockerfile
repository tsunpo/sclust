FROM ubuntu

# File Author / Maintainer
MAINTAINER Tsun-Po Yang <tyang2@uni-koeln.de>

RUN apt-get update && apt-get install -y gfortran build-essential \
make gcc build-essential

WORKDIR /opt

RUN git clone https://github.com/tsunpo/sclust-smc-het.git && cd sclust-smc-het && git checkout smchet
