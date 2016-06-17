FROM ubuntu

# File Author / Maintainer
MAINTAINER Tsun-Po Yang <tyang2@uni-koeln.de>

RUN apt-get update && apt-get install -y r-base-core
RUN R BATCH -f /opt/galaxy/tools/sclust-smc-het/inst_qp.R

WORKDIR /opt
