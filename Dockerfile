FROM ubuntu

# File Author / Maintainer
MAINTAINER Tsun-Po Yang <tyang2@uni-koeln.de>

RUN apt-get update \
make gcc build-essential

WORKDIR /opt

RUN git clone https://github.com/tsunpo/sclust-smc-het.git && cd sclust-smc-het && git checkout smchet
