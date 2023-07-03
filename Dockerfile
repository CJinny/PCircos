FROM python:3.9.16

LABEL maintainer="Jean-Baptiste Lamouche"
LABEL name="jb.lamouche@unistra.fr"

#Install vcf2circos
RUN pip install git+https://github.com/JbaptisteLam/vcf2circos@manuscript

#Download config
RUN wget https://www.lbgi.fr/~lamouche/vcf2circos/config_vcf2circos_29032023.tar.gz

#Untar config
RUN tar -xzf config_vcf2circos_29032023.tar.gz

#Set configuration path
RUN sed -i 's,\"Static\": \"/enadisk/maison/lamouche/dev_vcf2circos/Static\"\,,\"Static\": \"/Static\"\,,' /Static/options.json

ENTRYPOINT ["vcf2circos"]

