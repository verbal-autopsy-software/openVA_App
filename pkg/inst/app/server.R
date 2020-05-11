#' @import shiny
server <- function(input, output, session) {

  ## Read in data
  rvFile <- reactiveValues(clear = 0)
  getData <- reactive({
    vaData <- input$readIn
    if (is.null(vaData)) {
      return (NULL)
    }
    read.csv(input$readIn$datapath, stringsAsFactors = FALSE)
  })
  output$csvCheck <- renderText({
    if (is.null(input$readIn$datapath)) {
      return ("")
    }
    #vaData <- input$readIn
    validate(
      need(tools::file_ext(input$readIn$datapath) %in% c(
        "text/csv", "text/comma-separated-values,text/plain", "csv"), 
        "Wrong file format, please use CSV file")
      )
  })
  output$fileUploaded <- reactive({
    return(!is.null(getData()))
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)
  selectedAlgorithm <- reactive({
    switch(input$algorithm,
           "InSilicoVA" = 1, "InterVA5" = 2, "Tariff2" = 3)
  })

  ## Run model
  rv <- reactiveValues()
  rv$male     <- TRUE
  rv$female   <- TRUE
  rv$neonate  <- TRUE
  rv$child    <- TRUE
  rv$adult    <- TRUE
  rv$mNeonate  <- TRUE
  rv$mChild    <- TRUE
  rv$mAdult    <- TRUE
  rv$fNeonate  <- TRUE
  rv$fChild    <- TRUE
  rv$fAdult    <- TRUE


  observeEvent(input$processMe, {

    if (names(dev.cur()) != "null device") dev.off()
    fileNames <- dir(system.file('app', package = 'shinyVA'))
    filesToRemove <- grepl("plot|*warnings*|VA_results", fileNames)
    if (sum(filesToRemove) > 0) {
      file.remove(paste(system.file('app', package = 'shinyVA'),
                        fileNames[filesToRemove],
                        sep = "/"))
    }

    rv$fitAll     <- NULL
    rv$fitMale    <- NULL
    rv$fitFemale  <- NULL
    rv$fitNeonate <- NULL
    rv$fitChild   <- NULL
    rv$fitAdult   <- NULL
    rv$fitMNeonate  <- NULL
    rv$fitMChild    <- NULL
    rv$fitMAdult    <- NULL
    rv$fitFNeonate  <- NULL
    rv$fitFChild    <- NULL
    rv$fitFAdult    <- NULL

    shinyjs::disable("processMe")
    shinyjs::disable("algorithm")
    shinyjs::disable("downloadAgeDist")
    shinyjs::disable("downloadMetadata")
    shinyjs::disable("downloadAgeDist")
    shinyjs::disable("downloadMetadata")
    shinyjs::disable("downloadCOD1")
    shinyjs::disable("downloadCOD2")
    shinyjs::disable("downloadCOD3")
    shinyjs::disable("downloadCOD4")
    shinyjs::disable("downloadCOD5")
    shinyjs::disable("downloadCOD6")
    shinyjs::disable("downloadCOD7")
    shinyjs::disable("downloadCOD8")
    shinyjs::disable("downloadCOD9")
    shinyjs::disable("downloadCOD10")
    shinyjs::disable("downloadCOD11")
    shinyjs::disable("downloadCOD12")
    shinyjs::disable("downloadData1")
    shinyjs::disable("downloadData2")
    shinyjs::disable("downloadData3")
    shinyjs::disable("downloadData4")
    shinyjs::disable("downloadData5")
    shinyjs::disable("downloadData6")
    shinyjs::disable("downloadData7")
    shinyjs::disable("downloadData8")
    shinyjs::disable("downloadData9")
    shinyjs::disable("downloadData10")
    shinyjs::disable("downloadData11")
    shinyjs::disable("downloadData12")
    shinyjs::disable("downloadPlot1")
    shinyjs::disable("downloadPlot2")
    shinyjs::disable("downloadPlot3")
    shinyjs::disable("downloadPlot4")
    shinyjs::disable("downloadPlot5")
    shinyjs::disable("downloadPlot6")
    shinyjs::disable("downloadPlot7")
    shinyjs::disable("downloadPlot8")
    shinyjs::disable("downloadPlot9")
    shinyjs::disable("downloadPlot10")
    shinyjs::disable("downloadPlot11")
    shinyjs::disable("downloadPlot12")
    shinyjs::disable("downloadWarnings1")
    shinyjs::disable("downloadWarnings2")
    shinyjs::disable("downloadWarnings3")
    shinyjs::disable("downloadWarnings4")
    shinyjs::disable("downloadWarnings5")
    shinyjs::disable("downloadWarnings6")
    shinyjs::disable("downloadWarnings7")
    shinyjs::disable("downloadWarnings8")
    shinyjs::disable("downloadWarnings9")
    shinyjs::disable("downloadWarnings10")
    shinyjs::disable("downloadWarnings11")
    shinyjs::disable("downloadWarnings12")
    shinyjs::disable("downloadEverything")

    progress <- shiny::Progress$new()
    msg <- paste("Starting openVA...(this may take a while)")
    if (input$algorithm == "Tariff2") {
        msg <- paste("Starting SmartVA-Analyze...(this may take a while)")
    }
    progress$set(message = msg, value = 1/6)

    # set up objects needed to loop through all the calls to openVA
    if (input$algorithm == "InSilicoVA") {
      burn <- round(input$simLength / 2)
      modelArgs <- list(model = "InSilicoVA",
                        data.type = "WHO2016",
                        Nsim = isolate(input$simLength),
                        burnin = burn,
                        warning.write = TRUE) # include jump.scale and extend
    } else {
      modelArgs <- list(model = "InterVA",
                        version = "5.0",
                        data.type = "WHO2016",
                        HIV=input$HIV,
                        Malaria=input$Malaria,
                        directory=getwd(),
                        filename="VA_result",
                        returnCheckedData = TRUE)
    }
    # read in data
    pyCallStdout <- ""
    badData <- 0
    badConversion <- FALSE
    if (input$odkBC & input$algorithm != "Tariff2") {
      ## records <- CrossVA::odk2openVA(getData())
      ## records$ID <- getData()$meta.instanceID
      write.csv(getData(), file = 'tmpOut.csv', row.names = FALSE)
      pyAlg <- ifelse(input$algorithm == "InSilicoVA", "InsillicoVA", "InterVA5")
      pyCall <- paste0("/usr/local/bin/pycrossva-transform AUTODETECT ",
                       "InterVA5 tmpOut.csv --dst pyOut.csv")
      pyCallStdout <- system(pyCall, intern = TRUE)
      records <- read.csv("pyOut.csv", stringsAsFactors = FALSE)
      nCells <- nrow(records[,-1]) * ncol(records[,-1])
      if (sum(records[,-1] == ".") == nCells) badConversion <- TRUE
    } else {
      records <- getData()
      data(RandomVA5)
      if (ncol(RandomVA5) != ncol(records)) badData <- badData + 1
      vaColNames <- names(RandomVA5)
      if (badData == 0) badData <- badData + sum(vaColNames != names(records))
    }
    names(records) <- tolower(names(records))
    whoData <- "i004a" %in% names(records)
    all <- rep(TRUE, nrow(records))
    male <- rep(FALSE, length(records$i019a))
    male[records$i019a == "y"] <- TRUE
    female <- rep(FALSE, length(records$i019b))
    female[records$i019b == "y"] <- TRUE
    neonate <- rep(FALSE, length(records$i022g))
    neonate[records$i022g == "y"] <- TRUE
    child  <- rep(FALSE, length(records$i022f))
    child[records$i022f == "y" | records$i022e == "y" | records$i022d == "y"] <- TRUE
    adult <- rep(FALSE, length(records$i022a))
    adult[records$i022a == "y" | records$i022b == "y" | records$i022c == "y"] <- TRUE
    mNeonate <- rep(FALSE, length(records$i022g))
    mNeonate[male & neonate] <- TRUE
    mChild <- rep(FALSE, length(records$i022f))
    mChild[male & child] <- TRUE
    mAdult <- rep(FALSE, length(records$i022a))
    mAdult[male & adult] <- TRUE
    fNeonate <- rep(FALSE, length(records$i022g))
    fNeonate[female & neonate] <- TRUE
    fChild <- rep(FALSE, length(records$i022f))
    fChild[female & child] <- TRUE
    fAdult <- rep(FALSE, length(records$i022a))
    fAdult[female & adult] <- TRUE
    namesAll <- c("all", "male", "female", "neonate", "child", "adult",
                   "mNeonate", "mChild", "mAdult", "fNeonate", "fChild", "fAdult")
    namesRuns <- c("all", "male", "female", "neonate", "child", "adult",
                   "mNeonate", "mChild", "mAdult", "fNeonate", "fChild", "fAdult")
    namesNumericCodes <- 1:12
    names(namesNumericCodes) <- namesRuns
    includeRuns <- c(input$byAll, rep(input$bySex, 2), rep(input$byAge, 3),
                     rep(input$byAgeSex, 6))
    if (whoData) {
      nonZero <- c(sum(adult) > 0, sum(male) > 0, sum(female) > 0,
                   sum(neonate) > 0, sum(child) > 0, sum(adult) > 0,
                   sum(mNeonate) > 0, sum(mChild) > 0, sum(mAdult) > 0,
                   sum(fNeonate) > 0, sum(fChild) > 0, sum(fAdult) > 0)
      nGE100 <- c(sum(adult) >= 100, sum(male) >= 100, sum(female) >= 100,
                  sum(neonate) >= 100, sum(child) >= 100, sum(adult) >= 100,
                  sum(mNeonate) >= 100, sum(mChild) >= 100, sum(mAdult) >= 100,
                  sum(fNeonate) >= 100, sum(fChild) >= 100, sum(fAdult) >= 100)
      nRuns <- sum(includeRuns & nonZero)
      if (input$algorithm == "InSilicoVA") {
        namesRuns <- namesRuns[includeRuns & nonZero & nGE100]
      }
      if (input$algorithm == "InterVA") {
        namesRuns <- namesRuns[includeRuns & nonZero]
      }
    } else {
      nRuns <- sum(includeRuns)
      namesRuns <- namesRuns[includeRuns]
    }
    tmpDirResults <- paste0(getwd(), "/__tmp__", namesRuns)
    lapply(tmpDirResults, function (i) {
      if (dir.exists(i)) unlink(i, recursive = TRUE, force = TRUE)
    })
    lapply(tmpDirResults, function (i) {dir.create(i)})

    if (input$odkBC & input$algorithm != "Tariff2" & grepl("PHMRC", pyCallStdout[1])) {
      progress$close()
      rvFile$clear <- 1
      observe(rvFile$clear, shinyjs::reset("readIn"))
      msg <- "openVA App does not run InSilicoVA or InterVA5 with PHMRC
      data due to large number of missing items.  We may enable this in future versions."
      showNotification(msg, duration = NULL, type = "error", closeButton = FALSE,
                       action = a(href = "javascript:history.go(0)", "reset?"))
    } else if (badConversion & input$algorithm != "Tariff2") {
      progress$close()
      rvFile$clear <- 1
      observe(rvFile$clear, shinyjs::reset("readIn"))
      msg <- "Problem with data conversion (all cells are missing).  Are data really from ODK?"
      showNotification(msg, duration = NULL, type = "error", closeButton = FALSE,
                       action = a(href = "javascript:history.go(0);", "reset?"))
    } else if (badData > 0 & input$algorithm != "Tariff2") {
      progress$close()
      rvFile$clear <- 1
      observe(rvFile$clear, shinyjs::reset("readIn"))
      msg <- "Column names or number of columns do not match what openVA is expecting.  Unable to process data."
      showNotification(msg, duration = NULL, type = "error", closeButton = FALSE,
                       action = a(href = "javascript:history.go(0);", "reset?"))
    } else {

      # InSilico & InterVA5
      if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
        ageGroup <- rep(NA, length(records$i022a))
        ageGroup[neonate] <- "neonate"
        ageGroup[child]  <- "child"
        ageGroup[adult]  <- "ages >14"

        counts <- c(length(male[male]), length(female[female]),
                    length(neonate[neonate]), length(child[child]),
                    length(adult[adult]),
                    length(records$ID[records$i022a=="." & records$i022b=="." & records$i022c=="." &
                                        records$i022d=="." & records$i022e=="." & records$i022f=="." &
                                        records$i022g=="."]),
                    nrow(records))
        if(file.exists("plotAgeDist.pdf")) file.remove("plotAgeDist.pdf")
        pdf("plotAgeDist.pdf")
        barplot(table(ageGroup), horiz = TRUE, main="Age Distribution", xlab="Counts")
        dev.off()
        output$downloadAgeDist <- downloadHandler(
          filename = "plotAgeDist.pdf",
          content = function(file) {
            file.copy("plotAgeDist.pdf", file)
          }
        )
        # run codeVA()
        if (input$algorithm == "InSilicoVA") {
          nCores <- min(parallel::detectCores() - 1, length(namesRuns))
          if (nCores == 0) nCores <- 1
          cl <- parallel::makeCluster(nCores)
          parallel::clusterEvalQ(cl, {
            library(shiny)
            library(openVA)
          })
          seeds <- sample(1:2^15, length(namesRuns))
          vaOut <- list(all = NULL, male = NULL, female = NULL,
                        neonate = NULL, child = NULL, adult = NULL,
                        mNeonate = NULL, mChild = NULL, mAdult = NULL,
                        fNeonate = NULL, fChild = NULL, fAdult = NULL)
          parallel::clusterExport(cl,
                                  varlist = c("namesRuns", "tmpDirResults", "modelArgs", "seeds",
                                              "records", "all", "male", "female",
                                              "neonate", "child", "adult",
                                              "mNeonate", "mChild", "mAdult",
                                              "fNeonate", "fChild", "fAdult"),
                                  envir = environment())
          vaOut <- parallel::parLapply(cl, 1:length(namesRuns), function (i) {
            set.seed(seeds[i])
            modelArgs$data <- records[get(namesRuns[i]), ]
            modelArgs$inputID <- modelArgs$data[,1]
            modelArgs$directory <- tmpDirResults[i]
            okRun <- try(
              do.call(openVA::codeVA, modelArgs)
            )
            okRun$ID_orig <- records[get(namesRuns[i]), 1]
            okRun
          })
          parallel::stopCluster(cl)
          progress$set(message = "done with analyses", value = 5/6)
          names(vaOut) <- namesRuns
        } else {
          modelArgs$data <- records
          modelArgs$inputID <- modelArgs$data[,1]
          okRun <- do.call(openVA::codeVA, modelArgs)
          okRun$ID_orig <- records[, 1]
          vaOut <- sepVAResults(okRun)
          sepVALog(okRun, namesRuns, "errorlogV5.txt")
          progress$set(message = "done with analyses", value = 5/6)
        }
        lapply(1:length(namesRuns), function (i) {

          tmpNameRun <- namesRuns[i]
          groupName <- gsub("^(.)", "\\U\\1", tmpNameRun, perl = TRUE)
          rvName <- paste0("fit", groupName)
          rv[[rvName]] <- vaOut[[tmpNameRun]]

          titleDescriptiveStats <- paste0("titleDescriptiveStats", groupName)
          output[[titleDescriptiveStats]] <- renderText({
            "Counts of Deaths by Sex & Age"
          })
          descriptiveStatsName <- paste0("descriptiveStats", groupName)
          output[[descriptiveStatsName]] <- renderTable({
            if (!is.null(counts)) {
              matrix(counts, nrow=1, ncol=7,
                     dimnames = list(c("# of Deaths"),
                                     c("Male", "Female",
                                       "Neonate", "Child", "Ages >14",
                                       "Age is Missing", "Total")))
            }
          })

          # produce outputs
          if (!is.null(rv[[rvName]])) {
            rvNameIndivCOD <- paste0("indivCOD", groupName)            
            rv[[rvNameIndivCOD]] <- indivCOD(rv[[rvName]], top = 3)
            
            if (input$algorithm == "InSilicoVA" ) {
              orderedCSMF <- summary(rv[[rvName]])$csmf.ordered[, 1]
            } else {
              orderedCSMF <- summary(rv[[rvName]])$csmf.ordered[, 2]
            }
            newTop <- min(input$topDeaths, sum(orderedCSMF > 0))

            # CSMF Summary
            titleSummary <- paste0("titleSummary", groupName)
            output[[titleSummary]] <- renderText({
              paste("Summary of Results using", groupName, "Records")
            })
            summaryGrp <- paste0("summary", groupName)
            output[[summaryGrp]] <- renderPrint({
              if (!is.null(rv[[rvName]])) {
                print(summary(rv[[rvName]], top = newTop))
              }
            })
            # CSMF Plot
            plotName <- paste0("plot-", tmpNameRun, "-", input$algorithm, "-", Sys.Date(), ".pdf")
            if (file.exists(plotName)) file.remove(plotName)
            plotVA(rv[[rvName]], top = newTop)
            ggsave(plotName, device="pdf")
            downloadPlot <- paste0("downloadPlot", namesNumericCodes[tmpNameRun])
            output[[downloadPlot]] <- downloadHandler(
              filename = plotName,
              content = function(file) {
                if(!is.null(rv[[rvName]])){
                  file.copy(plotName, file)
                }
              }
            )
            titlePlot <- paste0("titlePlot", groupName)
            output[[titlePlot]] <- renderText({
              paste("CSMF Plot for", groupName, "Records")
            })
            plotGrp <- paste0("plot", groupName)
            output[[plotGrp]] <- renderPlot({
              if (!is.null(rv[[rvName]])) {
                plotVA(rv[[rvName]], top = newTop)
              }
            })
            # Download individual cause assignments
            downloadCOD <- paste0("downloadCOD", namesNumericCodes[tmpNameRun])
            output[[downloadCOD]] <- downloadHandler(
              filename = paste0("individual-causes-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv"),
              content = function(file) {
                if(!is.null(rv[[rvName]])){
                  write.csv(rv[[rvNameIndivCOD]], file = file, row.names = FALSE)
                }
              }
            )
            downloadData <- paste0("downloadData", namesNumericCodes[tmpNameRun])
            rvNameCSMFSummary <- paste0("csmfSummary", groupName)
            rv[[rvNameCSMFSummary]] <- csmfSummaryCSV(summary(rv[[rvName]], top = newTop))
            output[[downloadData]] <- downloadHandler(
              filename = paste0("results-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv"),
              content = function(file) {
                if(!is.null(rv[[rvName]])){
                  write.table(rv[[rvNameCSMFSummary]], file = file, row.names = FALSE, col.names = FALSE)
                }
              }
            )
            downloadWarnings <- paste0("downloadWarnings", namesNumericCodes[tmpNameRun])
            warningFileName <- paste0(input$algorithm, "-warnings-", namesRuns[i], ".txt")
            if(file.exists(warningFileName)) file.remove(warningFileName)
            file.create(warningFileName)
            cat("Warnings and Errors from", input$algorithm, "\t", namesRuns[i], "\t", date(),
                "\n", file = warningFileName)
            ovaLogFileName <- ifelse(input$algorithm == "InSilicoVA",
                                     paste0(tmpDirResults[i], "/errorlog_insilico.txt"),
                                     paste0(tmpDirResults[i], "/errorlogV5.txt"))
            file.append(warningFileName, ovaLogFileName)
            #unlink(ovaLogFileName, recursive = TRUE, force = TRUE)
            output[[downloadWarnings]] <- downloadHandler(
              filename = paste0(input$algorithm, "-warnings-", namesRuns[i], ".txt"),
              content = function(file) {
                if(!is.null(rv[[rvName]])){
                  file.copy(warningFileName, file)
                }
              }
            )
          }
            if (is.null(rv[[rvName]])) {
              titleSummary <- paste0("titleSummary", groupName)
              output[[titleSummary]] <- renderText({
                paste("No Summary for", groupName, "(not enough deaths for analysis)")
              })
            }
          if(is.null(rv[[rvName]])) rv[[tmpNameRun]] <- NULL
        })
        lapply(1:length(namesAll), function (i) {
          tmpNameRun <- namesAll[i]
          groupName <- gsub("^(.)", "\\U\\1", tmpNameRun, perl = TRUE)
          rvName <- paste0("fit", groupName)
          if (is.null(rv[[rvName]])) {
            titleSummary <- paste0("titleSummary", groupName)
            output[[titleSummary]] <- renderText({
              paste("No Summary for", groupName, "(results not requested or not enough deaths for analysis)")
            })
          }
        })
        output$downloadEverything <- downloadHandler(
          filename = "results.zip",
          content = function (file) {
            files <- NULL
            for (i in 1:length(namesRuns)) {
              tmpNameRun <- namesRuns[i]
              groupName <- gsub("^(.)", "\\U\\1", tmpNameRun, perl = TRUE)
              rvName <- paste0("fit", groupName)
              if (input$algorithm == "InSilicoVA" ) {
                orderedCSMF <- summary(rv[[rvName]])$csmf.ordered[, 1]
              } else {
                orderedCSMF <- summary(rv[[rvName]])$csmf.ordered[, 2]
              }
              newTop <- min(input$topDeaths, sum(orderedCSMF > 0))
              rvNameIndivCOD <- paste0("indivCOD", groupName)
              rvNameCSMFSummary <- paste0("csmfSummary", groupName)
              if(!is.null(rv[[rvName]])){
                # individual COD
                tmpFileName <- paste0("individual-causes-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv")
                write.csv(rv[[rvNameIndivCOD]], file = tmpFileName, row.names = FALSE)
                files <- c(files, tmpFileName)
                # summary
                tmpFileName <- paste0("results-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv")
                write.table(rv[[rvNameCSMFSummary]], file = tmpFileName, row.names = FALSE, col.names = FALSE)
                files <- c(files, tmpFileName)
                # plot
                tmpFileName <- paste0("plot-", tmpNameRun, "-", input$algorithm, "-", Sys.Date(), ".pdf")
                if (file.exists(tmpFileName)) file.remove(tmpFileName)
                plotVA(rv[[rvName]], top = newTop)
                ggsave(tmpFileName, device="pdf")
                files <- c(files, tmpFileName)
                # warnings
                tmpFileName <- paste0(input$algorithm, "-warnings-", namesRuns[i], ".txt")
                if(file.exists(tmpFileName)) file.remove(tmpFileName)
                file.create(tmpFileName)
                cat("Warnings and Errors from", input$algorithm, "\t", namesRuns[i], "\t", date(),
                    "\n", file = tmpFileName)
                ovaLogFileName <- ifelse(input$algorithm == "InSilicoVA",
                                         paste0(tmpDirResults[i], "/errorlog_insilico.txt"),
                                         paste0(tmpDirResults[i], "/errorlogV5.txt"))
                file.append(tmpFileName, ovaLogFileName)
                files <- c(files, tmpFileName)
              }
            }
            zip(file, files)
          }
        )
      }
      # Tariff2
      if (input$algorithm == "Tariff2") {

        file.remove(grep("plot-.*-Tariff2", dir(), value = TRUE))
        if (dir.exists("svaOut")) unlink("svaOut", recursive = TRUE, force = TRUE)
        dir.create("svaOut")
        if (file.exists("tmpOut.csv")) file.remove("tmpOut.csv")
        tmpOut <- getData()
        names(tmpOut) <- gsub("\\.", ":", names(tmpOut))
        write.csv(tmpOut, file = "tmpOut.csv", na = "", row.names = FALSE)
        svaCall <- paste("smartva", "--country", input$svaCountry,
                         "--hiv", ifelse(input$svaHIV, "True", "False"),
                         "--malaria", ifelse(input$svaMalaria, "True", "False"),
                         "--hce", ifelse(input$svaHCE, "True", "False"),
                         "--freetext", ifelse(input$svaFreeText, "True", "False"),
                         "--figures False",
                         "tmpOut.csv", "svaOut")
        system(svaCall)

        # render demographic table
        indCOD <- read.csv("svaOut/1-individual-cause-of-death/individual-cause-of-death.csv",
                           stringsAsFactors = FALSE)
        all <- rep(TRUE, nrow(indCOD))
        adult <- rep(FALSE, nrow(indCOD))
        adult[indCOD$age >= 12] <- TRUE
        child <- rep(FALSE, nrow(indCOD))
        child[indCOD$age < 12 & indCOD$age >= 29/365.25] <- TRUE
        neonate <- rep(FALSE, nrow(indCOD))
        neonate[indCOD$age < 29/365.25] <- TRUE
        ageGroup <- rep(NA, nrow(indCOD))
        ageGroup[neonate] <- "neonate"
        ageGroup[child]  <- "child"
        ageGroup[adult]  <- "ages >11"
        female <- rep(FALSE, nrow(indCOD))
        female[indCOD$sex == 2] <- TRUE
        male <- rep(FALSE, nrow(indCOD))
        male[indCOD$sex == 1] <- TRUE
        mNeonate <- rep(FALSE, nrow(indCOD))
        mNeonate[male & neonate] <- TRUE
        mChild <- rep(FALSE, nrow(indCOD))
        mChild[male & child] <- TRUE
        mAdult <- rep(FALSE, nrow(indCOD))
        mAdult[male & adult] <- TRUE
        fNeonate <- rep(FALSE, nrow(indCOD))
        fNeonate[female & neonate] <- TRUE
        fChild <- rep(FALSE, nrow(indCOD))
        fChild[female & child] <- TRUE
        fAdult <- rep(FALSE, nrow(indCOD))
        fAdult[female & adult] <- TRUE

        counts <- c(length(male[male]), length(female[female]),
                    length(neonate[neonate]), length(child[child]),
                    length(adult[adult]),
                    nrow(records) - nrow(indCOD),
                    nrow(records))
        if(file.exists("plotAgeDist.pdf")) file.remove("plotAgeDist.pdf")
        pdf("plotAgeDist.pdf")
        barplot(table(ageGroup), horiz = TRUE, main="Age Distribution", xlab="Counts")
        dev.off()
        output$downloadAgeDist <- downloadHandler(
          filename = "plotAgeDist.pdf",
          content = function(file) {
            file.copy("plotAgeDist.pdf", file)
          }
        )
        # render results
        svaCSMF <- read.csv("svaOut/2-csmf/csmf.csv", stringsAsFactors = FALSE)
        if (input$byAll & nrow(indCOD) > 0) {
          rv$fitAll <- svaCSMF[order(svaCSMF[, "all"], decreasing = TRUE),
                               c("cause34", "all")]
          names(rv$fitAll) <- c("cause34", "csmf")
          rownames(rv$fitAll) <- NULL
        }
        if (nrow(indCOD) == 0) rv$fitAll <- NULL
        if (input$bySex & length(male[male]) > 0) {
          rv$fitMale <- svaCSMF[order(svaCSMF[, "male"], decreasing = TRUE),
                                c("cause34", "male")]
          names(rv$fitMale) <- c("cause34", "csmf")
          rownames(rv$fitMale) <- NULL
        }
        if (length(male[male]) == 0) rv$fitMale <- NULL
        if (input$bySex & length(female[female]) > 0) {
          rv$fitFemale <- svaCSMF[order(svaCSMF[, "female"], decreasing = TRUE),
                                  c("cause34", "female")]
          names(rv$fitFemale) <- c("cause34", "csmf")
          rownames(rv$fitFemale) <- NULL
        }
        if (length(female[female]) == 0) rv$fitFemale <- NULL
        if (input$byAge & length(neonate[neonate]) > 0) {
          svaCSMFNeo <- read.csv("svaOut/2-csmf/neonate-csmf.csv", stringsAsFactors = FALSE)
          rv$fitNeonate <- svaCSMFNeo[order(svaCSMFNeo[, "all"], decreasing = TRUE),
                                      c("cause34", "all")]
          names(rv$fitNeonate) <- c("cause34", "csmf")
          rownames(rv$fitNeonate) <- NULL
        }
        if (length(neonate[neonate]) == 0) rv$fitNeonate <- NULL
        if (input$byAge & length(child[child]) > 0) {
          svaCSMFChild <- read.csv("svaOut/2-csmf/child-csmf.csv", stringsAsFactors = FALSE)
          rv$fitChild <- svaCSMFChild[order(svaCSMFChild[, "all"], decreasing = TRUE),
                                      c("cause34", "all")]
          names(rv$fitChild) <- c("cause34", "csmf")
          rownames(rv$fitChild) <- NULL
        }
        if (length(child[child]) == 0) rv$fitChild <- NULL
        if (input$byAge & length(adult[adult]) > 0) {
          svaCSMFAdult <- read.csv("svaOut/2-csmf/adult-csmf.csv", stringsAsFactors = FALSE)
          rv$fitAdult <- svaCSMFAdult[order(svaCSMFAdult[, "all"], decreasing = TRUE),
                                      c("cause34", "all")]
          names(rv$fitAdult) <- c("cause34", "csmf")
          rownames(rv$fitAdult) <- NULL
        }
        if (length(adult[adult]) == 0) rv$fitAdult <- NULL
        if (input$byAgeSex & length(neonate[neonate]) > 0) {
          svaCSMFNeonate <- read.csv("svaOut/2-csmf/neonate-csmf.csv", stringsAsFactors = FALSE)
          rv$fitMNeonate <- svaCSMFNeonate[order(svaCSMFNeonate[, "male"], decreasing = TRUE),
                                           c("cause34", "male")]
          names(rv$fitMNeonate) <- c("cause34", "csmf")
          rownames(rv$fitMNeonate) <- NULL
          rv$fitFNeonate <- svaCSMFNeonate[order(svaCSMFNeonate[, "female"], decreasing = TRUE),
                                            c("cause34", "female")]
          names(rv$fitFNeonate) <- c("cause34", "csmf")
          rownames(rv$fitFNeonate) <- NULL
        }
        if (length(neonate[neonate & male]) == 0) rv$fitMNeonate <- NULL
        if (length(neonate[neonate & female]) == 0) rv$fitFNeonate <- NULL
        if (input$byAgeSex & length(child[child]) > 0) {
          svaCSMFChild <- read.csv("svaOut/2-csmf/child-csmf.csv", stringsAsFactors = FALSE)
          rv$fitMChild <- svaCSMFChild[order(svaCSMFChild[, "male"], decreasing = TRUE),
                                       c("cause34", "male")]
          names(rv$fitMChild) <- c("cause34", "csmf")
          rownames(rv$fitMChild) <- NULL
          rv$fitFChild <- svaCSMFChild[order(svaCSMFChild[, "female"], decreasing = TRUE),
                                       c("cause34", "female")]
          names(rv$fitFChild) <- c("cause34", "csmf")
          rownames(rv$fitFChild) <- NULL
        }
        if (length(child[child & male]) == 0) rv$fitMChild <- NULL
        if (length(child[child & female]) == 0) rv$fitFChild <- NULL
        if (input$byAgeSex & length(adult[adult]) > 0) {
          svaCSMFAdult <- read.csv("svaOut/2-csmf/adult-csmf.csv", stringsAsFactors = FALSE)
          rv$fitMAdult <- svaCSMFAdult[order(svaCSMFAdult[, "male"], decreasing = TRUE),
                                       c("cause34", "male")]
          names(rv$fitMAdult) <- c("cause34", "csmf")
          rownames(rv$fitMAdult) <- NULL
          rv$fitFAdult <- svaCSMFAdult[order(svaCSMFAdult[, "female"], decreasing = TRUE),
                                       c("cause34", "female")]
          names(rv$fitFAdult) <- c("cause34", "csmf")
          rownames(rv$fitFAdult) <- NULL
        }
        if (length(adult[adult & male]) == 0) rv$fitMAdult <- NULL
        if (length(adult[adult & female]) == 0) rv$fitFAdult <- NULL

        lapply(1:length(namesRuns), function (i) {
          tmpNameRun <- namesRuns[i]
          groupName <- gsub("^(.)", "\\U\\1", tmpNameRun, perl = TRUE)
          rvName <- paste0("fit", groupName)

          titleDescriptiveStats <- paste0("titleDescriptiveStats", groupName)
          output[[titleDescriptiveStats]] <- renderText({
            "Counts of Deaths by Sex & Age"
          })
          descriptiveStatsName <- paste0("descriptiveStats", groupName)
          output[[descriptiveStatsName]] <- renderTable({
            if (!is.null(counts)) {
              matrix(counts, nrow=1, ncol=7,
                     dimnames = list(c("# of Deaths"),
                                     c("Male", "Female",
                                       "Neonate", "Child", "Ages >11",
                                       "Age is Missing", "Total")))
            }
          })
          if (!is.null(rv[[rvName]])) {
            newTop <- min(input$topDeaths, sum(rv[[rvName]]$csmf > 0))
            # CSMF Summary
            titleSummary <- paste0("titleSummary", groupName)
            titleSummaryText <- paste("Summary of Results using", groupName, "Records")
            output[[titleSummary]] <- renderText({
              paste("Tariff2: CSMF calculated using", groupName, "Records")
            })
            summaryGrp <- paste0("summary", groupName)
            output[[summaryGrp]] <- renderPrint({
              print(rv[[rvName]][1:newTop,], right = FALSE, row.names = FALSE)
            })
            # CSMF Plot
            plotName <- paste0("plot-", tmpNameRun, "-", input$algorithm, "-", Sys.Date(), ".pdf")
            if (file.exists(plotName)) file.remove(plotName)
            #pdf(plotName)
            marOld <- par()$mar
            par(mar = c(5, 10, 4, 2))
            barplot.df <- data.frame(Probability = rev(rv[[rvName]]$csmf[1:newTop]),
                                     Causes = rev(rv[[rvName]]$cause34[1:newTop]))
            g <- ggplot2::ggplot(barplot.df,
                                 ggplot2::aes(x = stats::reorder(Causes, seq(1:length(Causes))),
                                              y = Probability,
                                              fill = stats::reorder(Causes, seq(1:length(Causes))))) +
              ggplot2::geom_bar(stat="identity") + ggplot2::xlab("") + ggplot2::ylab("")
            g <- g + ggplot2::coord_flip() +
              ggplot2::scale_fill_grey(start = 0.8, end = 0.2) +
              ggplot2::theme(legend.position = "none")
            ggsave(plotName, device="pdf")
            par(mar = marOld)
            downloadPlot <- paste0("downloadPlot", namesNumericCodes[tmpNameRun])
            output[[downloadPlot]] <- downloadHandler(
              filename = plotName,
              content = function(file) {
                if(!is.null(rv[[rvName]])){
                  file.copy(plotName, file)
                }
              }
            )
            titlePlot <- paste0("titlePlot", groupName)
            output[[titlePlot]] <- renderText({
              paste("CSMF Plot for", groupName, "Records")
            })
            plotGrp <- paste0("plot", groupName)
            if (!is.null(rv[[rvName]])) {
              output[[plotGrp]] <- renderPlot({
                marOld <- par()$mar
                par(mar = c(5, 10, 4, 2))
                barplot.df <- data.frame(Probability = rev(rv[[rvName]]$csmf[1:newTop]),
                                         Causes = rev(rv[[rvName]]$cause34[1:newTop]))
                g <- ggplot2::ggplot(barplot.df,
                                     ggplot2::aes(x = stats::reorder(Causes, seq(1:length(Causes))),
                                                  y = Probability,
                                                  fill = stats::reorder(Causes, seq(1:length(Causes))))) +
                  ggplot2::geom_bar(stat="identity") + ggplot2::xlab("") + ggplot2::ylab("")
                g <- g + ggplot2::coord_flip() +
                  ggplot2::scale_fill_grey(start = 0.8, end = 0.2) +
                  ggplot2::theme(legend.position = "none")
                print(g)
                par(mar = marOld)
              })
            }
            # Download individual cause assignments
            downloadCOD <- paste0("downloadCOD", namesNumericCodes[tmpNameRun])
            output[[downloadCOD]] <- downloadHandler(
              filename = paste0("individual-causes-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv"),
              content = function(file) {
                write.csv(indCOD[get(namesRuns[i]),], file = file, row.names = FALSE)
              }
            )
            downloadData <- paste0("downloadData", namesNumericCodes[tmpNameRun])
            output[[downloadData]] <- downloadHandler(
              filename = paste0("results-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv"),
              content = function(file) {
                write.csv(print(rv[[rvName]]), file = file)
              }
            )
          }
          titleSummary <- paste0("titleSummary", groupName)
          emptySummary <- paste0("emptySummary", groupName)
          if (is.null(rv[[rvName]])) {
            output[[titleSummary]] <- renderText({
              paste("No Summary for", groupName, "(not enough deaths for analysis)")
            })
          }
          #if(is.null(rv[[rvName]])) rv[[tmpNameRun]] <- NULL
        })
        output$downloadEverything <- downloadHandler(
          filename = "results.zip",
          content = function (file) {
            files <- NULL
            for (i in 1:length(namesRuns)) {
              tmpNameRun <- namesRuns[i]
              groupName <- gsub("^(.)", "\\U\\1", tmpNameRun, perl = TRUE)
              rvName <- paste0("fit", groupName)
              newTop <- min(input$topDeaths, sum(rv[[rvName]]$csmf > 0))
              rvNameIndivCOD <- paste0("indivCOD", groupName)
              rvNameCSMFSummary <- paste0("csmfSummary", groupName)
              if(!is.null(rv[[rvName]])){
                # individual COD
                tmpFileName <- paste0("individual-causes-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv")
                write.csv(indCOD[get(namesRuns[i]),], file = tmpFileName, row.names = FALSE)
                files <- c(files, tmpFileName)
                # summary
                tmpFileName <- paste0("results-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv")
                write.csv(print(rv[[rvName]]), file = tmpFileName)
                files <- c(files, tmpFileName)
                # plot
                tmpFileName <- paste0("plot-", tmpNameRun, "-", input$algorithm, "-", Sys.Date(), ".pdf")
                if (file.exists(tmpFileName)) file.remove(tmpFileName)
                marOld <- par()$mar
                par(mar = c(5, 10, 4, 2))
                barplot.df <- data.frame(Probability = rev(rv[[rvName]]$csmf[1:newTop]),
                                         Causes = rev(rv[[rvName]]$cause34[1:newTop]))
                g <- ggplot2::ggplot(barplot.df,
                                     ggplot2::aes(x = stats::reorder(Causes, seq(1:length(Causes))),
                                                  y = Probability,
                                                  fill = stats::reorder(Causes, seq(1:length(Causes))))) +
                  ggplot2::geom_bar(stat="identity") + ggplot2::xlab("") + ggplot2::ylab("")
                g <- g + ggplot2::coord_flip() +
                  ggplot2::scale_fill_grey(start = 0.8, end = 0.2) +
                  ggplot2::theme(legend.position = "none")
                ggsave(tmpFileName, device="pdf")
                par(mar = marOld)
                files <- c(files, tmpFileName)
              }
            }
            zip(file, files)
          }
        )
        if (dir.exists("fontconfig")) unlink("fontconfig", recursive = TRUE, force = TRUE)
      }
      progress$close()
      # download for metadata
      metaData <- matrix("", nrow = 5, ncol = 2)
      metaData[1,] <- c("Metadata for openVA_App", as.character(Sys.time()))
      if (input$odkBC & input$algorithm != "Tariff2") {
        metaData[2,] <- c("pyCrossVA info", pyCallStdout[1])
        if (input$algorithm == "InSilicoVA") pkgVersion <- packageVersion("InSilicoVA")
        if (input$algorithm == "InterVA5") pkgVersion <- packageVersion("InterVA5")
        metaData[3,] <- c("Algorithm Package Version", as.character(pkgVersion))
      }
      if (input$algorithm == "Tariff2") {
        algCall <- svaCall
      } else {
        algCall <- modelArgs
      }

      algCall <- unlist(algCall)
      metaData[4,] <- c("openVA Arguments",
                        paste(names(algCall), algCall, sep = ": ", collapse = "; "))
      if (input$algorithm == "InSilicoVA") {
        metaData[5,] <- c("Random seeds",
                          paste(namesRuns, seeds, sep = ": ", collapse = "; "))
      }
      colnames(metaData) <- c("", "")
      output$downloadMetadata <- downloadHandler(
        filename = paste0("metadata", "-", input$algorithm, "-", Sys.Date(), ".csv"),
        content = function(file) {
          write.csv(metaData, file = file, row.names = FALSE)
        }
      )
      # download warning messages
      output$downloadWarnings <- downloadHandler(
        filename = warningFileName,
        content = function(file) {
          file.copy(filename, file)
        }
      )

      shinyjs::enable("processMe")
      shinyjs::enable("algorithm")
      shinyjs::enable("downloadAgeDist")
      shinyjs::enable("downloadMetadata")
      shinyjs::enable("downloadEverything")
      # HERE -- but this in an apply statement
      if (input$byAll) {
        shinyjs::enable("downloadPlot1")
        shinyjs::enable("downloadCOD1")
        shinyjs::enable("downloadData1")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings1")
        }
      }
      enableDownloads2 <- input$bySex & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
        (input$algorithm == "InSilicoVA" & length(male[male]) >= 100))
      #if (input$bySex & length(male[male]) > 0) {
      if (enableDownloads2) {
        shinyjs::enable("downloadPlot2")
        shinyjs::enable("downloadCOD2")
        shinyjs::enable("downloadData2")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings2")
        }
      }
      enableDownloads3 <- input$bySex & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
           (input$algorithm == "InSilicoVA" & length(female[female]) >= 100))
      #if (input$bySex & length(female[female]) > 0) {
      if (enableDownloads3) {
        shinyjs::enable("downloadPlot3")
        shinyjs::enable("downloadCOD3")
        shinyjs::enable("downloadData3")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings3")
        }
      }
      enableDownloads4 <- input$byAge & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
           (input$algorithm == "InSilicoVA" & length(neonate[neonate]) >= 100))
      #if (input$byAge & length(neonate[neonate]) > 0) {
      if (enableDownloads4) {
        shinyjs::enable("downloadPlot4")
        shinyjs::enable("downloadCOD4")
        shinyjs::enable("downloadData4")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings4")
        }
      }
      enableDownloads5 <- input$byAge & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
           (input$algorithm == "InSilicoVA" & length(child[child]) >= 100))
      #if (input$byAge & length(child[child]) > 0) {
      if (enableDownloads5) {
        shinyjs::enable("downloadPlot5")
        shinyjs::enable("downloadCOD5")
        shinyjs::enable("downloadData5")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings5")
        }
      }
      enableDownloads6 <- input$byAge & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
           (input$algorithm == "InSilicoVA" & length(adult[adult]) >= 100))
      #if (input$byAge & length(adult[adult]) > 0) {
      if (enableDownloads6) {
        shinyjs::enable("downloadPlot6")
        shinyjs::enable("downloadCOD6")
        shinyjs::enable("downloadData6")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings6")
        }
      }
      enableDownloads7 <- input$byAgeSex & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
           (input$algorithm == "InSilicoVA" & length(mNeonate[mNeonate]) >= 100))
      #if (input$byAgeSex & length(mNeonate[mNeonate]) > 0) {
      if (enableDownloads7) {
        shinyjs::enable("downloadPlot7")
        shinyjs::enable("downloadCOD7")
        shinyjs::enable("downloadData7")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings7")
        }
      }
      enableDownloads8 <- input$byAgeSex & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
           (input$algorithm == "InSilicoVA" & length(mChild[mChild]) >= 100))
      #if (input$byAgeSex & length(mChild[mChild]) > 0) {
      if (enableDownloads8) {
        shinyjs::enable("downloadPlot8")
        shinyjs::enable("downloadCOD8")
        shinyjs::enable("downloadData8")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings8")
        }
      }
      enableDownloads9 <- input$byAgeSex & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
           (input$algorithm == "InSilicoVA" & length(mAdult[mAdult]) >= 100))
      #if (input$byAgeSex & length(mAdult[mAdult]) > 0) {
      if (enableDownloads9) {
        shinyjs::enable("downloadPlot9")
        shinyjs::enable("downloadCOD9")
        shinyjs::enable("downloadData9")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings9")
        }
      }
      enableDownloads10 <- input$byAgeSex & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
           (input$algorithm == "InSilicoVA" & length(fNeonate[fNeonate]) >= 100))
      #if (input$byAgeSex & length(fNeonate[fNeonate]) > 0) {
      if (enableDownloads10) {
        shinyjs::enable("downloadPlot10")
        shinyjs::enable("downloadCOD10")
        shinyjs::enable("downloadData10")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings10")
        }
      }
      enableDownloads11 <- input$byAgeSex & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
           (input$algorithm == "InSilicoVA" & length(fChild[fChild]) >= 100))
      #if (input$byAgeSex & length(fChild[fChild]) > 0) {
      if (enableDownloads11) {
        shinyjs::enable("downloadPlot11")
        shinyjs::enable("downloadCOD11")
        shinyjs::enable("downloadData11")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings11")
        }
      }
      enableDownloads12 <- input$byAgeSex & 
        (input$algorithm == "Tariff2" | input$algorithm == "InterVA5" |
           (input$algorithm == "InSilicoVA" & length(fAdult[fAdult]) >= 100))
      #if (input$byAgeSex & length(fAdult[fAdult]) > 0) {
      if (enableDownloads12) {
        shinyjs::enable("downloadPlot12")
        shinyjs::enable("downloadCOD12")
        shinyjs::enable("downloadData12")
        if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
          shinyjs::enable("downloadWarnings12")
        }
      }
    }
  })

  # disable download button on page load
  shinyjs::disable("downloadAgeDist")
  shinyjs::disable("downloadMetadata")
  shinyjs::disable("downloadCOD1")
  shinyjs::disable("downloadCOD2")
  shinyjs::disable("downloadCOD3")
  shinyjs::disable("downloadCOD4")
  shinyjs::disable("downloadCOD5")
  shinyjs::disable("downloadCOD6")
  shinyjs::disable("downloadCOD7")
  shinyjs::disable("downloadCOD8")
  shinyjs::disable("downloadCOD9")
  shinyjs::disable("downloadCOD10")
  shinyjs::disable("downloadCOD11")
  shinyjs::disable("downloadCOD12")
  shinyjs::disable("downloadData1")
  shinyjs::disable("downloadData2")
  shinyjs::disable("downloadData3")
  shinyjs::disable("downloadData4")
  shinyjs::disable("downloadData5")
  shinyjs::disable("downloadData6")
  shinyjs::disable("downloadData7")
  shinyjs::disable("downloadData8")
  shinyjs::disable("downloadData9")
  shinyjs::disable("downloadData10")
  shinyjs::disable("downloadData11")
  shinyjs::disable("downloadData12")
  shinyjs::disable("downloadPlot1")
  shinyjs::disable("downloadPlot2")
  shinyjs::disable("downloadPlot3")
  shinyjs::disable("downloadPlot4")
  shinyjs::disable("downloadPlot5")
  shinyjs::disable("downloadPlot6")
  shinyjs::disable("downloadPlot7")
  shinyjs::disable("downloadPlot8")
  shinyjs::disable("downloadPlot9")
  shinyjs::disable("downloadPlot10")
  shinyjs::disable("downloadPlot11")
  shinyjs::disable("downloadPlot12")
  shinyjs::disable("downloadWarnings1")
  shinyjs::disable("downloadWarnings2")
  shinyjs::disable("downloadWarnings3")
  shinyjs::disable("downloadWarnings4")
  shinyjs::disable("downloadWarnings5")
  shinyjs::disable("downloadWarnings6")
  shinyjs::disable("downloadWarnings7")
  shinyjs::disable("downloadWarnings8")
  shinyjs::disable("downloadWarnings9")
  shinyjs::disable("downloadWarnings10")
  shinyjs::disable("downloadWarnings11")
  shinyjs::disable("downloadWarnings12")
  shinyjs::disable("downloadEverything")
}
