FROM ubuntu

# File Author / Maintainer
MAINTAINER Tsun-Po Yang <tyang2@uni-koeln.de>

WORKDIR /opt

RUN apt-get install -y r-base-core

ADD /opt/galaxy/tools/sclust-smc-het/inst_qp.R /opt/inst_qp.R
RUN R BATCH -f /opt/inst_qp.R
