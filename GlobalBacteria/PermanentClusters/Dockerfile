# Install R version 3.5
#FROM r-base:3.4.4


# get shiny server plus tidyverse packages image
FROM rocker/shiny-verse:latest

# Install Ubuntu packages

RUN apt-get update && apt-get install -y \
  apt-utils debconf debconf-utils \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-openssl-dev \
    libssl-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \ 
    libjpeg-dev \
    libcairo2-dev \
    libxt-dev \
    libssh2-1-dev \
  ncbi-blast+ \
  libxml2-dev \
  libmariadbclient-dev \
  libbz2-dev \
  libpcre3-dev \
  liblzma-dev \
  zlib1g-dev \  
  openjdk-8-jdk-headless

# install R packages required 
# (change it dependeing on the packages you need)
RUN R -e "install.packages('shiny', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinydashboard', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyjs', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('data.table', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('leaflet', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('DT', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ECharts2Shiny', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('RMySQL', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dplyr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('leaflet.extras', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('plyr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('digest', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stringr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('markdown', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('readxl', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('stringi', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('lazyeval', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('sp', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyBS', repos='http://cran.rstudio.com/')"
RUN R CMD javareconf
RUN R -e "install.packages('mailR', repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_github('RedOakStrategic/geoshaper')"

RUN mkdir --parents /home/fungal/databases

# Copy configuration files into the Docker image
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf

# copy the app to the image
COPY /app /srv/shiny-server/

# select port
EXPOSE 3838

# Copy further configuration files into the Docker image
COPY shiny-server.sh /usr/bin/shiny-server.sh

RUN adduser --system --disabled-password --group shiny

# allow permission
RUN sudo chown -R shiny:shiny /srv/shiny-server

#RUN apt-get install -y mariadb-client

# run app
RUN ["chmod", "+x", "/usr/bin/shiny-server.sh"]
CMD "/usr/bin/shiny-server.sh"
