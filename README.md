# OpenVA_RShiny

Analyze Verbal Autopsy (VA) data with the InSilicoVA & InterVA algorithms using `shiny-openVA.R`, which  builds a web app with the 
**R** package [Shiny](https://cran.r-project.org/web/packages/shiny/index.html).  The app can process VA data from the 2016 WHO VA
instrument (version 1.4.1).  Data from an Open Data Kit (ODK) Briefcase export are also accepted formats for the app.

## Installation

The dependencies for the OpenVA_RShiny app include Java, and R, as well as several **R** packages.  

To run the OpenVA_RShiny app: (1) download the files `shiny-openVA.R`; (2) install Java Development Kit; 
(3) install [**R**](https://cran.r-project.org/); and (4) start **R** and install the necessary packages with the following command:

```r
install.packages(c("openVA", "shinyjs", "CrossVA"), dependencies = TRUE)
```

*Notes* 

- The app depends on [**CrossVA** ](https://cran.r-project.org/package=CrossVA) **version 0.9.3** as well as version 1.2.5
of the **R** package [InSilicoVA](https://github.com/verbal-autopsy-software/InSilicoVA) package. To update a previously
installed versions of these packages, run the following command at the **R** prompt (which may ask you to select a CRAN 
mirror and to confirm the update of each package for which you have administrative privileges):

  ```r
  ## update all packages and have R ask about updates for each package:
  update.packages()
  ## update all packages without confirmation
  ## update.packages(ask = FALSE)
  ```

## Usage

After installing the dependencies, open the file `shiny-openVA.R` in R, select all of the code and execute it.  This repo 
has two example ODK Briefcase exported data sets (_odkBriefcaseExport.csv_) and (_who151_va_output.csv_), which can be used 
by the Shiny app.  When running InterVA5 or InSilico2016, users can specify the version of the 2016 WHO VA Questionnaire
used to collected the data.  The file (_odkBriefcaseExport.csv_) is an example of version 1.4.1 of the 2016
WHO VA Questionnaire, while the file (_who151_va_output.csv_) is an example of version 1.5.1.  The shiny app will appear 
in your default web browser.  See the [vignette](https://github.com/verbal-autopsy-software/shinyVA/blob/master/shiny-openVA-vignette.pdf) for more details.  _Note: The vignette needs to be updated, and may not reflect the current state of the shiny app._


## Troubleshooting

InSilicoVA depends on the **R** package rJava.  It is common to run into problems with loading rJava (and thus InSilicoVA).  

- For linux and Macs, it may help to open a terminal run the command ```R CMD javareconf```, and then try loading InSilicoVA. 

- For Macs R.app will look for jdk-9.  If you are using a different version of Java (e.g., jdk-8 or jdk10) then you need to tell R.app where to look for your Java installation.  In the terminal, run the following command with **your computer-specific paths**
  ```bash
  install_name_tool -change \
  /Library/Java/JavaVirtualMachines/jdk-9.jdk/Contents/Home/lib/server/libjvm.dylib \
  /Library/Java/JavaVirtualMachines/jdk-10.0.2.jdk/Contents/Home/lib/server/libjvm.dylib \ 
  /Library/Frameworks/R.framework/Resources/library/rJava/libs/rJava.so
  ```
  **you will probably need to modify the third line above to reference the path to jdk-10 (or jdk-8) on your machine.** See [this discussion](https://github.com/s-u/rJava/issues/151) for more details.
  
- For Macs, you may want to try using [RStudio](https://www.rstudio.com/) instead of R.app, since RStudio does a better job finding the appropriate version of Java.
  
- On Windows, you may have
  luck with the following commands at the **R** prompt:
  ```r
  options(java.home = "C:\\Path\\to\\Java\\jdk")
  library(InSilicoVA)
  ## Option 2:
  ## Sys.setenv("JAVA_HOME" = "C:\\Path\\to\\Java\\jdk")
  ```
  
  Another alternative is to create a new Environment Variable with variable name "JAVA_HOME" and variable
  value "C:\Program Files\Java\jdk-8u181\jre" (or wherever Java JDK is installed on your computer), and then restart your computer.  For
  more information on setting Environment Variables in windows see: 
  [https://www.java.com/en/download/help/path.xml](https://www.java.com/en/download/help/path.xml).

## To Do

- update vignette
- add support for Tariff and WHO VA instrument

## Licence
GNU General Public License v3.0
