# Install rocker
FROM rocker/shiny:3.6.0

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
RUN cd /opt; wget https://www.python.org/ftp/python/3.7.5/Python-3.7.5.tgz
RUN cd /opt; tar -xf Python-3.7.5.tgz; cd Python-3.7.5; ./configure --enable-optimizations --enable-shared; sudo make altinstall
RUN mkdir /home/shiny/GitHub
RUN git -C /home/shiny/GitHub/ clone https://github.com/verbal-autopsy-software/pyCrossVA
RUN sudo ldconfig /usr/local/lib
RUN pip3.7 install pip setuptools wheel --upgrade
RUN pip3.7 install -r /home/shiny/GitHub/pyCrossVA/requirements.txt
RUN cd /home/shiny/GitHub/pyCrossVA; python3.7 setup.py install
RUN R CMD javareconf
RUN R -e "install.packages(c('glue', 'shinyjs', 'openVA', 'CrossVA'), repos='http://cran.rstudio.com/')"
RUN git -C /home/shiny/GitHub/ clone https://github.com/verbal-autopsy-software/shinyVA
RUN R CMD INSTALL /home/shiny/GitHub/shinyVA/pkg

## NEED TO CHANGE WORKING DIRECTORY
RUN mkdir /srv/shiny-server/shinyVA
RUN echo "appDir <- system.file('app', package = 'shinyVA'); library(shinyVA); shinyAppDir(appDir)" > /srv/shiny-server/shinyVA/app.R
RUN chown -R shiny:shiny /srv/shiny-server
RUN chown -R shiny:shiny /usr/local/lib/R/site-library/shinyVA

# Start container with:
# docker run --rm --user shiny -p 3838:3838 -v /srv/shinylog/:/var/log/shiny-server/ shiny_app

# Open app at:
# localhost:3838/shinyVA

# Open app at:
# localhost:3838/shinyVA
