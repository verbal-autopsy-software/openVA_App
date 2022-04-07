# Install rocker
FROM rocker/shiny:4.1.3

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libncurses5-dev \
    libgdbm-dev \
    libnss3-dev \
    libssl-dev \
    libreadline-dev \
    libffi-dev \
    wget \
    #python3-pip \
    #python3-dev \
    liblzma-dev \
    libbz2-dev \
    libicu-dev \
    openjdk-8-jdk \
    git

# Install R packages that are required
WORKDIR /opt
RUN wget https://www.python.org/ftp/python/3.7.5/Python-3.7.5.tgz
RUN tar -xf Python-3.7.5.tgz
WORKDIR Python-3.7.5
RUN ./configure --enable-optimizations --enable-shared
RUN sudo make altinstall
RUN mkdir /home/shiny/GitHub
RUN git -C /home/shiny/GitHub/ clone https://github.com/verbal-autopsy-software/pyCrossVA -b openVA_App 
RUN sudo ldconfig /usr/local/lib
RUN pip3.7 install pip setuptools wheel --upgrade
RUN pip3.7 install -r /home/shiny/GitHub/pyCrossVA/requirements.txt
WORKDIR /home/shiny/GitHub/pyCrossVA
RUN python3.7 setup.py install
RUN R CMD javareconf
RUN R -e "install.packages(c('glue', 'shinyjs', 'openVA'), repos='http://cran.rstudio.com/')"
#RUN git -C /home/shiny/GitHub/ clone https://github.com/verbal-autopsy-software/InterVA5
#RUN R CMD INSTALL /home/shiny/GitHub/InterVA5/InterVA5_1.1.1.tar.gz
#RUN git -C /home/shiny/GitHub/ clone https://github.com/verbal-autopsy-software/InSilicoVA
#RUN Rscript -e "install.packages('/home/shiny/GitHub/InSilicoVA/InSilicoVA', repose = NULL, type = 'source')"
#RUN R CMD INSTALL /home/shiny/GitHub/InSilicoVA/InSilicoVA
RUN git -C /home/shiny/GitHub/ clone https://github.com/verbal-autopsy-software/openVA_App
RUN R CMD INSTALL /home/shiny/GitHub/openVA_App/pkg

# Set up SmartVA-Analyze
RUN wget -P /opt/SmartVA https://github.com/ihmeuw/SmartVA-Analyze/releases/download/v2.1.0/smartva
RUN chmod 755 /opt/SmartVA/smartva
ENV PATH="/opt/SmartVA:${PATH}"

# Set up app
RUN mkdir /srv/shiny-server/openVA_App
RUN echo "appDir <- system.file('app', package = 'openVAapp'); options(shiny.maxRequestSize = 100*1024^2); library(openVAapp); shinyAppDir(appDir)" > /srv/shiny-server/openVA_App/app.R
RUN chown -R shiny:shiny /srv/shiny-server
RUN chown -R shiny:shiny /usr/local/lib/R/site-library/openVAapp

# Start container with:
# docker run --rm --user shiny -p 3838:3838 -v /srv/shinylog/:/var/log/shiny-server/ openVA_App

# Open app at:
# localhost:3838/openVA_App
