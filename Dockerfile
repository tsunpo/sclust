FROM ubuntu

# File Author / Maintainer
MAINTAINER Tsun-Po Yang <tyang2@uni-koeln.de>

RUN apt-get install r-base
RUN R BATCH -f /opt/galaxy/tools/sclust-smc-het/inst_qp.R

WORKDIR /opt
