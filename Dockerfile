FROM ubuntu

# File Author / Maintainer
MAINTAINER Tsun-Po Yang <tyang2@uni-koeln.de>

RUN sudo apt-get install r-base-core
RUN sudo R BATCH -f /opt/galaxy/tools/sclust-smc-het/inst_qp.R

WORKDIR /opt
