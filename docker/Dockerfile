FROM ubuntu:20.04
LABEL org.opencontainers.image.authors="antony.lebechec@chru-strasbourg.fr"


# ARGS

ARG CLOUDGENE_RELEASE="2.5.3"

##########
# Python #
##########
RUN apt-get update -y
RUN apt-get install -y openjdk-11-jre
RUN apt-get install -y python3.6 python3-pip curl docker 
#docker-ce docker-ce-cli


#############
# CloudGene #
#############

ENV TOOLS=/tools

# TOOL INFO
ENV TOOL_NAME="cloudgene"
ENV TOOL_VERSION=$CLOUDGENE_RELEASE
ENV TOOL_TARBALL="v$TOOL_VERSION.tar.gz"
ENV TOOL_SOURCE_EXTERNAL="https://github.com/genepi/cloudgene/archive/refs/tags/$TOOL_TARBALL"
ENV PATH=$TOOLS/$TOOL_NAME/current/bin:$PATH
# TOOL PARAMETERS


RUN echo "BUILD"

# TOOL INSTALLATION
RUN if [ -z $TOOL_VERSION ]; then TOOL_VERSION_FOLDER=latest; else TOOL_VERSION_FOLDER=$TOOL_VERSION; fi && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin && \
    (cd $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin && \
    curl -s install.cloudgene.io | bash -s $TOOL_VERSION && \
	ln -s $TOOL_VERSION_FOLDER $TOOLS/$TOOL_NAME/current && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/config/apps && \
    ln -s config/apps $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/apps && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/config/pages && \
    mv $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/sample/pages $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/sample/pages.orignal && \
    ln -s ../config/pages $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/sample/pages && \
    cp -pR $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/sample/pages.orignal/*.stache config/pages/ && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/data/tmp && \
    ln -s data/tmp $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/tmp && \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/data/workspace && \
    ln -s data/workspace $TOOLS/$TOOL_NAME/$TOOL_VERSION_FOLDER/bin/workspace)


# Service
ENV CLOUDGENE_SERVICE_NAME="Cloudgene vcf2circos"
ENV CLOUDGENE_HELP_PAGE=https://github.com/bioinfo-chru-strasbourg/vcf2circos
ENV START_CLOUDGENE=true
ENV START_HADOOP=true

# ENDPOINT 
# EXPOSE 8082
# ENTRYPOINT [ "tools/cloudgene/current/bin/launch.sh" ]



##############
# vcf2circos #
##############

RUN mkdir /app
WORKDIR /app

RUN pip3 install --upgrade pip
ADD ./README.md /app/
ADD ./LICENSE /app/
ADD ./setup* /app/
ADD ./vcf2circos /app/vcf2circos
ADD ./demo_data /app/demo_data
ADD ./docs /app/docs
RUN python3 -m pip install -e .


####################
# Cloudgene Config #
####################

# Scripts
# ADD ./docker/cloudgene/conf/* $TOOLS/$TOOL_NAME/current/bin/config/
# ADD ./docker/cloudgene/scripts/launch.sh $TOOLS/$TOOL_NAME/current/bin/

ADD ./docker/cloudgene/conf/log4j2.xml $TOOLS/$TOOL_NAME/current/bin/config/
ADD ./docker/cloudgene/conf/settings.yaml $TOOLS/$TOOL_NAME/current/bin/config/
ADD ./docker/cloudgene/conf/apps $TOOLS/$TOOL_NAME/current/bin/config/apps
ADD ./docker/cloudgene/conf/pages $TOOLS/$TOOL_NAME/current/bin/config/pages
ADD ./docker/cloudgene/scripts/launch.sh $TOOLS/$TOOL_NAME/current/bin/


##############
# ENTRYPOINT #
##############

ENTRYPOINT ["vcf2circos"]
