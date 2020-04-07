#' @import shiny
server <- function(input, output, session) {

  ## Read in data
  getData <- reactive({
    vaData <- input$readIn
    if(is.null(vaData)){
      return(NULL)
    }
    read.csv(vaData$datapath, stringsAsFactors = FALSE)
  })
  output$fileUploaded <- reactive({
    return(!is.null(getData()))
  })
  outputOptions(output, "fileUploaded", suspendWhenHidden = FALSE)

  selectedAlgorithm <- reactive({
    # validate(
    #   need(compareVersion("0.9.3", packageDescription("CrossVA")$Version) <= 0,
    #        "Please update the CrossVA package: update.packages()")
    # )
    switch(input$algorithm,
           "InSilicoVA" = 1, "InterVA5" = 2)
  })

  ## Run model
  rv <- reactiveValues()
  rv$male     <- TRUE
  rv$female   <- TRUE
  rv$neonate  <- TRUE
  rv$child    <- TRUE
  rv$adult    <- TRUE

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

    shinyjs::disable("processMe")
    shinyjs::disable("algorithm")
    shinyjs::disable("downloadAgeDist")
    shinyjs::disable("downloadCOD1")
    shinyjs::disable("downloadCOD2")
    shinyjs::disable("downloadCOD3")
    shinyjs::disable("downloadCOD4")
    shinyjs::disable("downloadCOD5")
    shinyjs::disable("downloadCOD6")
    shinyjs::disable("downloadData1")
    shinyjs::disable("downloadData2")
    shinyjs::disable("downloadData3")
    shinyjs::disable("downloadData4")
    shinyjs::disable("downloadData5")
    shinyjs::disable("downloadData6")
    shinyjs::disable("downloadPlot1")
    shinyjs::disable("downloadPlot2")
    shinyjs::disable("downloadPlot3")
    shinyjs::disable("downloadPlot4")
    shinyjs::disable("downloadPlot5")
    shinyjs::disable("downloadPlot6")
    shinyjs::disable("downloadWarnings1")
    shinyjs::disable("downloadWarnings2")
    shinyjs::disable("downloadWarnings3")
    shinyjs::disable("downloadWarnings4")
    shinyjs::disable("downloadWarnings5")
    shinyjs::disable("downloadWarnings6")

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
                        #directory=getwd(),
                        filename="VA_result")
    }
    namesRuns <- c("all", "male", "female", "neonate", "child", "adult")
    namesNumericCodes <- 1:6
    names(namesNumericCodes) <- namesRuns
    includeRuns <- c(input$byAll, rep(input$bySex, 2), rep(input$byAge, 3))
    nRuns <- sum(includeRuns)
    namesRuns <- namesRuns[includeRuns]
    tmpDirResults <- paste0(getwd(), "/__tmp__", namesRuns)
    lapply(tmpDirResults, function (i) {
      if (dir.exists(i)) unlink(i, recursive = TRUE, force = TRUE)
      })
    lapply(tmpDirResults, function (i) {dir.create(i)})
    # read in data
    if(input$odkBC & input$algorithm != 'Tariff2'){
      ## records <- CrossVA::odk2openVA(getData())
      ## records$ID <- getData()$meta.instanceID
      write.csv(getData(), file = 'tmpOut.csv', row.names = FALSE)
      pyAlg <- ifelse(input$algorithm == "InSilicoVA", "InsillicoVA", "InterVA5")
      pyCall <- paste0('/usr/local/bin/pycrossva-transform AUTODETECT ',
                       'InterVA5', ' tmpOut.csv --dst pyOut.csv')
      system(pyCall)
      records <- read.csv('pyOut.csv', stringsAsFactors = FALSE)
    } else{
      records <- getData()
    }

    # object needed to render table of demographic variables
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
    child[records$i022f == "y"] <- TRUE
    child[records$i022e == "y"] <- TRUE
    child[records$i022d == "Y"] <- TRUE
    adult <- rep(FALSE, length(records$i022a))
    adult[records$i022a == "y"] <- TRUE
    adult[records$i022b == "y"] <- TRUE
    adult[records$i022c == "y"] <- TRUE

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

    # InSilico & InterVA5
    if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
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
      nCores <- min(parallel::detectCores() - 1, length(namesRuns))
      if (nCores == 0) nCores <- 1
      cl <- parallel::makeCluster(nCores)
      parallel::clusterEvalQ(cl, {
        library(shiny)
        library(openVA)
      })
      vaOut <- list(all = NULL, male = NULL, female = NULL,
                    neonate = NULL, child = NULL, adult = NULL)
      parallel::clusterExport(cl,
                              varlist = c("namesRuns", "tmpDirResults", "modelArgs",
                                          "records", "all", "male", "female",
                                          "neonate", "child", "adult"),
                              envir = environment())
      vaOut <- parallel::parLapply(cl, 1:length(namesRuns), function (i) {
        modelArgs$data <- records[get(namesRuns[i]), ]
        modelArgs$directory <- tmpDirResults[i]
        okRun <- try(
          do.call(openVA::codeVA, modelArgs)
        )
        okRun
      })
      parallel::stopCluster(cl)
      progress$set(message = "done with analyses", value = 5/6)
      names(vaOut) <- namesRuns

      lapply(1:length(namesRuns), function (i) {

        tmpNameRun <- namesRuns[i]
        groupName <- gsub('^(.)', '\\U\\1', tmpNameRun, perl = TRUE)
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

          rv$indivCOD <- indivCOD(rv[[rvName]], top = 3)

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
          if (input$algorithm == "InSilicoVA") {
            plotVA(rv[[rvName]], top = newTop)
            ggsave(plotName, device="pdf")
          } else {
            pdf(plotName)
            plotVA(rv[[rvName]], top = newTop)
            dev.off()
          }
          downloadPlot <- paste0('downloadPlot', namesNumericCodes[tmpNameRun])
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
          downloadCOD <- paste0('downloadCOD', namesNumericCodes[tmpNameRun])
          output[[downloadCOD]] <- downloadHandler(
            filename = paste0("individual-causes-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv"),
            content = function(file) {
              if(!is.null(rv[[rvName]])){
                write.csv(rv$indivCOD, file = file, row.names = FALSE)
              }
            }
          )
          downloadData <- paste0('downloadData', namesNumericCodes[tmpNameRun])
          output[[downloadData]] <- downloadHandler(
            filename = paste0("results-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv"),
            content = function(file) {
              if(!is.null(rv[[rvName]])){
                write.csv(print(summary(rv[[rvName]], top = newTop)), file = file)
              }
            }
          )
          downloadWarnings <- paste0('downloadWarnings', namesNumericCodes[tmpNameRun])
          warningFileName <- paste(input$algorithm, "warnings", namesRuns[i], ".txt", sep = "-")
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
            filename = paste(input$algorithm, "warnings", namesRuns[i], ".txt", sep = "-"),
            content = function(file) {
              if(!is.null(rv[[rvName]])){
                file.copy(warningFileName, file)
              }
            }
          )
        }
        titleSummary <- paste0("titleSummary", groupName)
        output[[titleSummary]] <- renderText({
          if (is.null(rv[[rvName]])) {
            paste("No Summary for", groupName, "(not enough deaths for analysis)")
          }
        })
        if(is.null(rv[[rvName]])) rv[[tmpNameRun]] <- NULL
      })
    }
    # Tariff2
    if (input$algorithm == "Tariff2"){
      file.remove(grep('plot-.*-Tariff2', dir(), value = TRUE))
      if (dir.exists('svaOut')) unlink('svaOut', recursive = TRUE, force = TRUE)
      dir.create('svaOut')
      if (file.exists('tmpOut.csv')) file.remove('tmpOut.csv')
      tmpOut <- getData()
      names(tmpOut) <- gsub("\\.", ":", names(tmpOut))
      write.csv(tmpOut, file = 'tmpOut.csv', na = "", row.names = FALSE)
      svaCall <- paste('smartva', '--country', input$svaCountry,
                       '--hiv', ifelse(input$svaHIV, 'True', 'False'),
                       '--malaria', ifelse(input$svaMalaria, 'True', 'False'),
                       '--hce', ifelse(input$svaHCE, 'True', 'False'),
                       '--freetext', ifelse(input$svaFreeText, 'True', 'False'),
                       '--figures False',
                       'tmpOut.csv', 'svaOut')
      system(svaCall)

      # render demographic table
      indCOD <- read.csv('svaOut/1-individual-cause-of-death/individual-cause-of-death.csv',
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
      svaCSMF <- read.csv('svaOut/2-csmf/csmf.csv', stringsAsFactors = FALSE)
      if (input$byAll & nrow(indCOD) > 0) {
        rv$fitAll <- svaCSMF[order(svaCSMF[, 'all'], decreasing = TRUE),
                             c('cause34', 'all')]
        names(rv$fitAll) <- c('cause34', 'csmf')
        rownames(rv$fitAll) <- NULL
      }
      if (input$bySex & length(male[male]) > 0) {
        rv$fitMale <- svaCSMF[order(svaCSMF[, 'male'], decreasing = TRUE),
                             c('cause34', 'male')]
        names(rv$fitMale) <- c('cause34', 'csmf')
        rownames(rv$fitMale) <- NULL
      }
      if (input$bySex & length(female[female]) > 0) {
        rv$fitFemale <- svaCSMF[order(svaCSMF[, 'female'], decreasing = TRUE),
                                c('cause34', 'female')]
        names(rv$fitFemale) <- c('cause34', 'csmf')
        rownames(rv$fitFemale) <- NULL
      }
      if (input$byAge & length(neonate[neonate]) > 0) {
        svaCSMFNeo <- read.csv('svaOut/2-csmf/neonate-csmf.csv', stringsAsFactors = FALSE)
        rv$fitNeonate <- svaCSMFNeo[order(svaCSMFNeo[, 'all'], decreasing = TRUE),
                                    c('cause34', 'all')]
        names(rv$fitNeonate) <- c('cause34', 'csmf')
        rownames(rv$fitNeonate) <- NULL
      }
      if (input$byAge & length(child[child]) > 0) {
        svaCSMFChild <- read.csv('svaOut/2-csmf/child-csmf.csv', stringsAsFactors = FALSE)
        rv$fitChild <- svaCSMFChild[order(svaCSMFChild[, 'all'], decreasing = TRUE),
                                  c('cause34', 'all')]
        names(rv$fitChild) <- c('cause34', 'csmf')
        rownames(rv$fitChild) <- NULL
      }
      if (input$byAge & length(adult[adult]) > 0) {
        svaCSMFAdult <- read.csv('svaOut/2-csmf/adult-csmf.csv', stringsAsFactors = FALSE)
        rv$fitAdult <- svaCSMFAdult[order(svaCSMFAdult[, 'all'], decreasing = TRUE),
                                    c('cause34', 'all')]
        names(rv$fitAdult) <- c('cause34', 'csmf')
        rownames(rv$fitAdult) <- NULL
      }

      lapply(1:length(namesRuns), function (i) {

        tmpNameRun <- namesRuns[i]
        groupName <- gsub('^(.)', '\\U\\1', tmpNameRun, perl = TRUE)
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
          pdf(plotName)
          barplot(height = rev(rv[[rvName]]$csmf[1:newTop]), horiz = TRUE,
                  names = gsub(' |\\/', '\n', rev(rv[[rvName]]$cause34[1:newTop])),
                  col = grey.colors(length(rv[[rvName]]$csmf[1:newTop])),
                  las = 1, cex.names = .75, tcl = -0.2,
                  mar = c(5, 25, 4, 2), mgp = c(3, .25, 0))
          dev.off()
          downloadPlot <- paste0('downloadPlot', namesNumericCodes[tmpNameRun])
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
            barplot(height = rev(rv[[rvName]]$csmf[1:newTop]), horiz = TRUE,
                    names = gsub(' |\\/', '\n', rev(rv[[rvName]]$cause34[1:newTop])),
                    col = grey.colors(length(rv[[rvName]]$csmf[1:newTop])),
                    las = 1, cex.names = .75, tcl = -0.2,
                    mar = c(5, 25, 4, 2), mgp = c(3, .25, 0))
          })
          # Download individual cause assignments
          downloadCOD <- paste0('downloadCOD', namesNumericCodes[tmpNameRun])
          output[[downloadCOD]] <- downloadHandler(
            filename = paste0("individual-causes-", namesRuns[i], "-", input$algorithm, "-", Sys.Date(), ".csv"),
            content = function(file) {
              write.csv(indCOD[get(namesRuns[i]),], file = file, row.names = FALSE)
            }
          )
          downloadData <- paste0('downloadData', namesNumericCodes[tmpNameRun])
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
        if(is.null(rv[[rvName]])) rv[[tmpNameRun]] <- NULL
      })
      if (dir.exists('fontconfig')) unlink('fontconfig', recursive = TRUE, force = TRUE)
    }
    progress$close()

    shinyjs::enable("processMe")
    shinyjs::enable("algorithm")
    shinyjs::enable("downloadAgeDist")
    if (input$byAll) {
      shinyjs::enable("downloadPlot1")
      shinyjs::enable("downloadCOD1")
      shinyjs::enable("downloadData1")
      if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
        shinyjs::enable("downloadWarnings1")
      }
    }
    if (input$bySex & length(male[male]) > 0) {
      shinyjs::enable("downloadPlot2")
      shinyjs::enable("downloadCOD2")
      shinyjs::enable("downloadData2")
      if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
        shinyjs::enable("downloadWarnings2")
      }
    }
    if (input$bySex & length(female[female]) > 0) {
      shinyjs::enable("downloadPlot3")
      shinyjs::enable("downloadCOD3")
      shinyjs::enable("downloadData3")
      if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
        shinyjs::enable("downloadWarnings3")
      }
    }
    if (input$byAge & length(neonate[neonate]) > 0) {
      shinyjs::enable("downloadPlot4")
      shinyjs::enable("downloadCOD4")
      shinyjs::enable("downloadData4")
      if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
        shinyjs::enable("downloadWarnings4")
      }
    }
    if (input$byAge & length(child[child]) > 0) {
      shinyjs::enable("downloadPlot5")
      shinyjs::enable("downloadCOD5")
      shinyjs::enable("downloadData5")
      if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
        shinyjs::enable("downloadWarnings5")
      }
    }
    if (input$byAge & length(adult[adult]) > 0) {
      shinyjs::enable("downloadPlot6")
      shinyjs::enable("downloadCOD6")
      shinyjs::enable("downloadData6")
      if (input$algorithm == "InSilicoVA" | input$algorithm == "InterVA5") {
        shinyjs::enable("downloadWarnings6")
      }
    }
  })

  output$downloadWarnings <- downloadHandler(
    filename = warningFileName,
    content = function(file) {
        file.copy(filename, file)
    }
  )
  # disable download button on page load
  shinyjs::disable("downloadAgeDist")
  shinyjs::disable("downloadCOD1")
  shinyjs::disable("downloadCOD2")
  shinyjs::disable("downloadCOD3")
  shinyjs::disable("downloadCOD4")
  shinyjs::disable("downloadCOD5")
  shinyjs::disable("downloadCOD6")
  shinyjs::disable("downloadData1")
  shinyjs::disable("downloadData2")
  shinyjs::disable("downloadData3")
  shinyjs::disable("downloadData4")
  shinyjs::disable("downloadData5")
  shinyjs::disable("downloadData6")
  shinyjs::disable("downloadPlot1")
  shinyjs::disable("downloadPlot2")
  shinyjs::disable("downloadPlot3")
  shinyjs::disable("downloadPlot4")
  shinyjs::disable("downloadPlot5")
  shinyjs::disable("downloadPlot6")
  shinyjs::disable("downloadWarnings1")
  shinyjs::disable("downloadWarnings2")
  shinyjs::disable("downloadWarnings3")
  shinyjs::disable("downloadWarnings4")
  shinyjs::disable("downloadWarnings5")
  shinyjs::disable("downloadWarnings6")
}
