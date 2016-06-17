FROM ubuntu

# File Author / Maintainer
MAINTAINER Tsun-Po Yang <tyang2@uni-koeln.de>

RUN apt-get update && apt-get install -y r-base-core

WORKDIR /opt

RUN git clone https://github.com/tsunpo/sclust-smc-het.git
RUN cd sclust-smc-het && R BATCH -f install_quadprog.R