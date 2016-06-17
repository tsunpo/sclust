FROM ubuntu

# File Author / Maintainer
MAINTAINER Tsun-Po Yang <tyang2@uni-koeln.de>

RUN apt-get update

WORKDIR /opt

ADD /opt/galaxy/tools/sclust-smc-het/Sclust /opt/Sclust
