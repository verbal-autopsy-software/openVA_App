# Install rocker
FROM rocker/shiny:3.6.0

# Install Ubuntu packages
RUN apt-get update && apt-get install -y \
    python3-pip \
    python3-dev \
    liblzma-dev \
    libbz2-dev \
    libicu-dev \
    openjdk-8-jdk \
    git

# Install R packages that are required
RUN R CMD javareconf
RUN R -e "install.packages(c('glue', 'shinyjs', 'openVA', 'CrossVA'), repos='http://cran.rstudio.com/')"
RUN mkdir /home/shiny/GitHub
RUN git -C /home/shiny/GitHub/ clone https://github.com/verbal-autopsy-software/shinyVA
RUN cd /home/shiny/GitHub/shinyVA/pkg
RUN R CMD INSTALL /home/shiny/GitHub/shinyVA/pkg

## NEED TO CHANGE WORKING DIRECTORY
RUN mkdir /srv/shiny-server/shinyVA
RUN echo "appDir <- system.file('app', package = 'shinyVA'); library(shinyVA); shinyAppDir(appDir)" > /srv/shiny-server/shinyVA/app.R
RUN chown -R shiny:shiny /srv/shiny-server
RUN chmod 777 /srv/shiny-server/shinyVA


# Start container with
# docker run --rm --user shiny -p 3838:3838 -v /srv/shinylog/:/var/log/shiny-server/ shiny_app
# and open app with at localhost:3838/shinyVA
