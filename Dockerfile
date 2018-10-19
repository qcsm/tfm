FROM rocker/shiny:latest
MAINTAINER Carlos Quijano

## Update and install git:
RUN apt-get update && apt-get install -y git

## Install Secure Sockets Layer (SSL) from OpenSSL
RUN apt-get install -y libssl-dev &&  \
    apt-get clean && \ 
    rm -rf /var/lib/apt/lists/ && \ 
    rm -rf /tmp/downloaded_packages/ /tmp/*.rds


## Install R Packages for Shiny app:
RUN sudo R -e "install.packages(c('caret','tidyverse','RColorBrewer','scales'), repos='http://cran.us.r-project.org')"

## Install R Packages for R ggplot2 theming:
RUN sudo R -e "install.packages(c('extrafont'))"

## Install R Packages for WGCNA App
RUN Rscript -e "install.packages('ape')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('shiny')"
RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('WGCNA')"

## Copy Shiny application into Docker container
COPY app.R /srv/shiny-server/app.R

## Copy Shiny Configuration into Docker container
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /srv/shiny-server/shiny-server.sh

## Copy data into Docker container
COPY BLCA_tumor.Rda /srv/shiny-server/BLCA_tumor.Rda
COPY BRCA_tumor.Rda /srv/shiny-server/BRCA_tumor.Rda
COPY COAD_tumor.Rda /srv/shiny-server/COAD_tumor.Rda
COPY HNSC_tumor.Rda /srv/shiny-server/HNSC_tumor.Rda
COPY KICH_tumor.Rda /srv/shiny-server/KICH_tumor.Rda
COPY KIRC_tumor.Rda /srv/shiny-server/KIRC_tumor.Rda
COPY KIRP_tumor.Rda /srv/shiny-server/KIRP_tumor.Rda
COPY LIHC_tumor.Rda /srv/shiny-server/LIHC_tumor.Rda
COPY LUAD_tumor.Rda /srv/shiny-server/LUAD_tumor.Rda
COPY LUSC_tumor.Rda /srv/shiny-server/LUSC_tumor.Rda
COPY PRAD_tumor.Rda /srv/shiny-server/PRAD_tumor.Rda
COPY THCA_tumor.Rda /srv/shiny-server/THCA_tumor.Rda

RUN chown -R shiny.shiny /srv/shiny-server

EXPOSE 3838
EXPOSE 80

CMD ["/usr/bin/shiny-server.sh"]
