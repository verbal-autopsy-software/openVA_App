# openVA app

Analyze Verbal Autopsy (VA) data with the InSilicoVA, InterVA5, and Tariff2 algorithms using the openVA app, available 
from this repository in the form of an R package (in the pkg folder).  The app can process VA data from the 2016 WHO VA 
instrument (versions 1.4.1 and 1.5.1) and the PHMRC Shortened Questionnaire.  Results are rendered in the app and can be 
saved in CSV and PDF formats.  The openVA app and all the dependencies are also available as a Docker image at
[https://hub.docker.com/r/openvateam/openva_app](https://hub.docker.com/r/openvateam/openva_app).

## Installing and Running openVA app

The recommended way to install and use the openVA app is by installing Docker Desktop and then downloading and running
the GUI located in [releases](https://github.com/verbal-autopsy-software/openVA_App/releases).

## Non-Docker Usage

The dependencies for the openVA app include R, Java, Python (version >=3.6 *and* version 2.7 if you wish to run the 
Tariff2 algorithm), and [pyCrossVA](https://github.com/verbal-autopsy-software/pyCrossVA).  If you wish to run the Tariff2 algorithm, [SmartVA-Analyze](https://github.com/ihmeuw/SmartVA-Analyze/releases) (command line version) must also be installed []()

After satisfying the dependencies, install and run the openVA app with the following commands in R:

```r
install.packages(c("openVA", "shinyjs", "devtools"), dependencies = TRUE)
devtools::install_github('verbal-autopsy-software/openVA_App', subdir = 'pkg', INSTALL_opts=c('--no-multiarch'))
library(openVAapp)
launchApp()
```

R may ask if you would like to update some of the packages openVA app depends on; you can type in the number for
the option you choose (e.g., 1: All, 2: CRAN packages only; 3: None). 


### Troubleshooting installation of R package

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

- Another error message on Windows complains about the inability to install the 32bit version of InSilico (i386),
  even though you are using the 64-bit version.  A potential solution is to use the option that tells R not to
  install the 32bit version:
  ```r
  install.packages("openVA", INSTALL_opts = "--no-multiarch")
  ```

## Video Tutorials
Installation Guide (Windows): https://youtu.be/C2EPOpTzJTk

Analysis Guide: https://youtu.be/K1wkSbTwxkg

## Licence
GNU General Public License v3.0
