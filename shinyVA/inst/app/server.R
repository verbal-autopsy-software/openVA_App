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
    validate(
      need(compareVersion("0.9.3", packageDescription("CrossVA")$Version) <= 0,
           "Please update the CrossVA package: update.packages()")
    )
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
    shinyjs::disable("downloadWarnings")

    withProgress(value=0, {

      setProgress(message = paste("Starting analysis of data (this may take a while)"))

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
                          filename="VA_result")
      }
      namesRuns <- c("all", "male", "female", "neonate", "child", "adult")
      namesNumericCodes <- 1:6
      names(namesNumericCodes) <- namesRuns
      includeRuns <- c(input$byAll, rep(input$bySex, 2), rep(input$byAge, 3))
      nRuns <- sum(includeRuns)
      namesRuns <- namesRuns[includeRuns]

      # read in data
      if(input$odkBC){
        records <- CrossVA::odk2openVA(getData())
        records$ID <- getData()$meta.instanceID
      } else{
        records <- getData()
      }

      # object needed to render table of demographic variables
      # HERE -- create helper function that does table for age and sex
      names(records) <- tolower(names(records))
      all <- rep(TRUE, nrow(records))
      male <- rep(FALSE, length(records$i019a))
      male[records$i019a=="y"] <- TRUE
      female <- rep(FALSE, length(records$i019b))
      female[records$i019b=="y"] <- TRUE
      neonate <- rep(FALSE, length(records$i022g))
      neonate[records$i022g=="y"] <- TRUE
      child  <- rep(FALSE, length(records$i022f))
      child[records$i022f=="y"] <- TRUE
      child[records$i022e=="y"] <- TRUE
      child[records$i022d=="Y"] <- TRUE
      adult <- rep(FALSE, length(records$i022a))
      adult[records$i022a=="y"] <- TRUE
      adult[records$i022b=="y"] <- TRUE
      adult[records$i022c=="y"] <- TRUE

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

      output$descriptiveStats <- renderTable({
        if(!is.null(counts)){
          matrix(counts, nrow=1, ncol=7,
                 dimnames = list(c("# of Deaths"),
                                 c("Male", "Female", "Neonate", "Child", "Ages >14",
                                   "Age is Missing", "Total")))
        }
      })
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
      warningFileName <- paste(input$algorithm, "warnings.txt", sep = "-")
      if(file.exists(warningFileName)) file.remove(warningFileName)
      file.create(warningFileName)
      cat("Warnings and Errors from", input$algorithm, "\t", date(),
          "\n", file = warningFileName)
      ovaLogFileName <- ifelse(input$algorithm == "InSilicoVA",
                               "errorlog_insilico.txt",
                               "errorlogV5.txt")
      for (i in namesRuns) {
        # HERE -- create a cleanup function
        # HERE -- only run data check once (if possible)
        modelArgs$data <- records[get(i), ]
        rvName <- paste0("fit", gsub('^(.)', '\\U\\1', i, perl = TRUE))
        okRun <- try(
          rv[[rvName]] <- do.call(openVA::codeVA, modelArgs)
          )
        incProgress(1/length(namesRuns),
                    detail = paste("Completed analysis using", namesRuns[i], "deaths")
                    )
        # produce outputs
        if (!is.null(rv[[rvName]])) {
          rv$indivCOD <- indivCOD(rv[[rvName]], top = 3)
          cat("Analysis with ", i, "\t", date(), "\n", file = warningFileName)
          file.append(warningFileName, ovaLogFileName)
          file.remove(ovaLogFileName)
          plotName <- paste0("plot-", i, "-", input$algorithm, "-", Sys.Date(), ".pdf")
          if (file.exists(plotName)) file.remove(plotName)
          plotVA(rv[[rvName]], top = input$topDeaths); ggsave(plotName, device="pdf")

          downloadPlot <- paste0('downloadPlot', namesNumericCodes[i])
          output[[downloadPlot]] <- downloadHandler(
            filename = plotName,
            content = function(file) {
              if(!is.null(rv[[rvName]])){
                file.copy(plotName, file)
              }
            }
          )
          downloadCOD <- paste0('downloadCOD', namesNumericCodes[i])
          output[[downloadCOD]] <- downloadHandler(
            filename = paste0("individual-causes-", i, "-", input$algorithm, "-", Sys.Date(), ".csv"),
            content = function(file) {
              if(!is.null(rv[[rvName]])){
                write.csv(rv$indivCOD, file = file, row.names = FALSE)
              }
            }
          )
          downloadData <- paste0('downloadData', namesNumericCodes[i])
          output[[downloadData]] <- downloadHandler(
            filename = paste0("results-", i, "-", input$algorithm, "-", Sys.Date(), ".csv"),
            content = function(file) {
              if(!is.null(rv[[rvName]])){
                write.csv(print(summary(rv[[rvName]], top = input$topDeaths)), file = file)
              }
            }
          )
        }
        if(is.null(rv[[rvName]])) rv[[i]] <- NULL
      }
    })
    shinyjs::enable("processMe")
    shinyjs::enable("algorithm")
    shinyjs::enable("downloadAgeDist")
    if(input$byAll){
      shinyjs::enable("downloadPlot1")
      shinyjs::enable("downloadCOD1")
      shinyjs::enable("downloadData1")
    }
    if(input$bySex & length(male[male])>0){
      shinyjs::enable("downloadPlot2")
      shinyjs::enable("downloadCOD2")
      shinyjs::enable("downloadData2")
    }
    if(input$bySex & length(female[female])>0){
      shinyjs::enable("downloadPlot3")
      shinyjs::enable("downloadCOD3")
      shinyjs::enable("downloadData3")
    }
    if(input$byAge & length(neonate[neonate])>0){
      shinyjs::enable("downloadPlot4")
      shinyjs::enable("downloadCOD4")
      shinyjs::enable("downloadData4")
    }
    if(input$byAge & length(child[child])>0){
      shinyjs::enable("downloadPlot5")
      shinyjs::enable("downloadCOD5")
      shinyjs::enable("downloadData5")
    }
    if(input$byAge & length(adult[adult])>0){
      shinyjs::enable("downloadPlot6")
      shinyjs::enable("downloadCOD6")
      shinyjs::enable("downloadData6")
    }
  })

  ## Print warning messages -- HERE

  ## Output table with descriptive statistics
  output$titleDescriptiveStats <- renderText({
    ## if(!is.null(rv$counts)){
    "Counts of Deaths by Sex & Age"
    ## }
  })

  ## output$descriptiveStats <- renderTable({
  ##     if(!is.null(rv$counts)){
  ##         ## matrix(rv$counts, nrow=1, ncol=11, dimnames = list(c("# of Deaths"),
  ##         ##                                                    c("Male", "Female", "Neonate", "Neonate",
  ##         ##                                                      "Age 1-4", "Age 5-14", "Age 14-49",
  ##         ##                                                      "Age 50-64", "Age 65+", "Age is Missing", "Total")))
  ##         matrix(rv$counts, nrow=1, ncol=7, dimnames = list(c("# of Deaths"),
  ##                                                            c("Male", "Female", "Neonate (0-1)", "Child (1-4)", "Ages 5+",
  ##                                                              "Age is Missing", "Total")))
  ##     }
  ## })

  output$downloadWarnings <- downloadHandler(
    filename = "warnings-openVA.txt",
    content = function(file) {
        file.copy(warningFileName, file)
    }
  )

  ## All
  #### Summarize and print output (all)
  output$titleSummaryAll <- renderText({
    ## if(!is.null(rv$fitAll)){
    "Summary of Results using All Records"
    ## }
  })
  output$summaryAll <- renderPrint({
    # validate(
    #   need(!is.null(rv$fitALL), "No Results")
    # )
    algorithm <- selectedAlgorithm()
    if(!is.null(rv$fitAll)){
      if(algorithm==1 & is.null(rv$fitAll$HIV)){
        ## printIndiv(rv$fitAll, rv$agg.csmf)
        print(summary(rv$fitAll, top=input$topDeaths))
      }
      if(algorithm==2 & !is.null(rv$fitAll$HIV)){
        print(summary(rv$fitAll, top=input$topDeaths))
      }
      if(algorithm==3 & is.null(rv$fitAll$HIV)){
        print(summary(rv$fitAll, top=input$topDeaths))
      }
      if(algorithm==4 & !is.null(rv$fitAll$HIV)){
        print(summary(rv$fitAll, top=input$topDeaths))
      }
    }
  })
  #### Create plot (all)
  output$titlePlotAll <- renderText({
    ## if(!is.null(rv$fitAll)){
    "CSMF Plot for Total Population"
    ## }
  })
  output$plotAll <- renderPlot({
    #      validate(
    #        need(!is.null(rv$fitALL), "No Results")
    #      )
    if(!is.null(rv$fitAll)){
      if(input$algorithm=="InterVA4" & !is.null(rv$fitAll$HIV)){
        ## CSMF(rv$fitAll, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF2(rv$fitAll, top=input$topDeaths)
      }
      if(input$algorithm=="InterVA5" & !is.null(rv$fitAll$HIV)){
        ## CSMF(rv$fitAll, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF5(rv$fitAll, top.plot=input$topDeaths, InterVA.rule = TRUE)
      }
      if(input$algorithm=="InSilicoVA_2012" & is.null(rv$fitAll$HIV)){
        plotVA(rv$fitAll, top=input$topDeaths)
      }
      if(input$algorithm=="InSilicoVA_2016" & is.null(rv$fitAll$HIV)){
        plotVA(rv$fitAll, top=input$topDeaths)
      }
    }
  })
  ## Male
  #### Summarize and print output (male)
  output$titleSummaryMale <- renderText({
    ## if(!is.null(rv$fitMale)){
    "Summary of Results for Males"
    ## }
  })
  output$emptySummaryMale <- renderText({
    if(is.null(rv$male)){
      "No Summary for Males (not enough deaths for analysis)"
    }
  })
  output$summaryMale <- renderPrint({

    algorithm <- selectedAlgorithm()
    if(!is.null(rv$fitMale)){
      if(algorithm==1 & is.null(rv$fitMale$HIV)){
        ## sumMale <- printIndiv(rv$fitMale, rv$agg.csmfMale)
        print(summary(rv$fitMale, top=input$topDeaths))
      }
      if(algorithm==2 & !is.null(rv$fitMale$HIV)){
        print(summary(rv$fitMale, top=input$topDeaths))
      }
      if(algorithm==3 & is.null(rv$fitMale$HIV)){
        ## sumMale <- printIndiv(rv$fitMale, rv$agg.csmfMale)
        print(summary(rv$fitMale, top=input$topDeaths))
      }
      if(algorithm==4 & !is.null(rv$fitMale$HIV)){
        print(summary(rv$fitMale, top=input$topDeaths))
      }
    }
  })
  #### Create plot (male)
  output$titlePlotMale <- renderText({
    ## if(!is.null(rv$fitMale)){
    "CSMF Plot for Males"
    ## }
  })
  output$emptyPlotMale <- renderText({
    if(is.null(rv$male)){
      "No Plot for Males (not enough deaths for analysis)"
    }
  })
  output$plotMale <- renderPlot({
    if(!is.null(rv$fitMale)){
      if(input$algorithm=="InterVA4" & !is.null(rv$fitMale$HIV)){
        ## CSMF(rv$fitMale, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF2(rv$fitMale, top=input$topDeaths)
      }
      if(input$algorithm=="InterVA5" & !is.null(rv$fitMale$HIV)){
        ## CSMF(rv$fitMale, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF5(rv$fitMale, top.plot=input$topDeaths, InterVA.rule = TRUE)
      }
      if(input$algorithm=="InSilicoVA_2012" & is.null(rv$fitMale$HIV)){
        ## indivplot(rv$agg.csmfMale, top=20, title="Aggregated COD distribution")
        plotVA(rv$fitMale, top=input$topDeaths)
      }
      if(input$algorithm=="InSilicoVA_2016" & is.null(rv$fitMale$HIV)){
        ## indivplot(rv$agg.csmfMale, top=20, title="Aggregated COD distribution")
        plotVA(rv$fitMale, top=input$topDeaths)
      }
    }
  })

  ## Female
  #### Summarize and print output (female)
  output$titleSummaryFemale <- renderText({
    ## if(!is.null(rv$fitFemale)){
    "Summary of Results for Females"
    ## }
  })
  output$emptySummaryFemale <- renderText({
    if(is.null(rv$female)){
      "No Summary for Females (not enough deaths for analysis)"
    }
  })
  output$summaryFemale <- renderPrint({
    algorithm <- selectedAlgorithm()
    if(!is.null(rv$fitFemale)){
      if(algorithm==1 & is.null(rv$fitFemale$HIV)){
        ## sumFemale <- printIndiv(rv$fitFemale, rv$agg.csmfFemale)
        print(summary(rv$fitFemale, top=input$topDeaths))
      }
      if(algorithm==2 & !is.null(rv$fitFemale$HIV)){
        print(summary(rv$fitFemale, top=input$topDeaths))
      }
      if(algorithm==3 & is.null(rv$fitFemale$HIV)){
        ## sumFemale <- printIndiv(rv$fitFemale, rv$agg.csmfFemale)
        print(summary(rv$fitFemale, top=input$topDeaths))
      }
      if(algorithm==4 & !is.null(rv$fitFemale$HIV)){
        print(summary(rv$fitFemale, top=input$topDeaths))
      }
    }
  })
  #### Create plot (female)
  output$titlePlotFemale <- renderText({
    ## if(!is.null(rv$fitFemale)){
    "CSMF Plot for Females"
    ## }
  })
  output$emptyPlotFemale <- renderText({
    if(is.null(rv$female)){
      "No Plot for Females (not enough deaths for analysis)"
    }
  })
  output$plotFemale <- renderPlot({
    if(!is.null(rv$fitFemale)){
      if(input$algorithm=="InterVA4" & !is.null(rv$fitFemale$HIV)){
        ## CSMF(rv$fitFemale, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF2(rv$fitFemale, top=input$topDeaths)
      }
      if(input$algorithm=="InterVA5" & !is.null(rv$fitFemale$HIV)){
        ## CSMF(rv$fitFemale, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF5(rv$fitFemale, top.plot=input$topDeaths, InterVA.rule = TRUE)
      }
      if(input$algorithm=="InSilicoVA_2012" & is.null(rv$fitFemale$HIV)){
        ## indivplot(rv$agg.csmfFemale, top=20, title="Aggregated COD distribution")
        plotVA(rv$fitFemale, top=input$topDeaths)
      }
      if(input$algorithm=="InSilicoVA_2016" & is.null(rv$fitFemale$HIV)){
        ## indivplot(rv$agg.csmfFemale, top=20, title="Aggregated COD distribution")
        plotVA(rv$fitFemale, top=input$topDeaths)
      }
    }
  })

  ## Neonate
  #### Summarize and print output
  output$titleSummaryNeonate <- renderText({
    ## if(!is.null(rv$fitNeonate)){
    "Summary of Results for Neonates"
    ## }
  })
  output$emptySummaryNeonate <- renderText({
    if(is.null(rv$neonate)){
      "No Summary for Neonates (not enough deaths for analysis)"
    }
  })
  output$summaryNeonate <- renderPrint({
    algorithm <- selectedAlgorithm()
    if(!is.null(rv$fitNeonate)){
      if(algorithm==1 & is.null(rv$fitNeonate$HIV)){
        ## sumNeonate <- printIndiv(rv$fitNeonate, rv$agg.csmfNeonate)
        print(summary(rv$fitNeonate, top=input$topDeaths))
      }
      if(algorithm==2 & !is.null(rv$fitNeonate$HIV)){
        print(summary(rv$fitNeonate, top=input$topDeaths))
      }
      if(algorithm==3 & is.null(rv$fitNeonate$HIV)){
        ## sumNeonate <- printIndiv(rv$fitNeonate, rv$agg.csmfNeonate)
        print(summary(rv$fitNeonate, top=input$topDeaths))
      }
      if(algorithm==4 & !is.null(rv$fitNeonate$HIV)){
        print(summary(rv$fitNeonate, top=input$topDeaths))
      }
    }
  })
  #### Create plot (neonates)
  output$titlePlotNeonate <- renderText({
    ## if(!is.null(rv$fitNeonate)){
    "CSMF Plot for Neonates"
    ## }
  })
  output$emptyPlotNeonate <- renderText({
    if(is.null(rv$neonate)){
      "No Plot for Neonates (not enough deaths for analysis)"
    }
  })
  output$plotNeonate <- renderPlot({
    if(!is.null(rv$fitNeonate)){
      if(input$algorithm=="InterVA4" & !is.null(rv$fitNeonate$HIV)){
        ## CSMF(rv$fitNeonate, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF2(rv$fitNeonate, top=input$topDeaths)
      }
      if(input$algorithm=="InterVA5" & !is.null(rv$fitNeonate$HIV)){
        ## CSMF(rv$fitNeonate, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF5(rv$fitNeonate, top.plot=input$topDeaths, InterVA.rule = TRUE)
      }
      if(input$algorithm=="InSilicoVA_2012" & is.null(rv$fitNeonate$HIV)){
        ## indivplot(rv$agg.csmfNeonate, top=20, title="Aggregated COD distribution")
        plotVA(rv$fitNeonate, top=input$topDeaths)
      }
      if(input$algorithm=="InSilicoVA_2016" & is.null(rv$fitNeonate$HIV)){
        ## indivplot(rv$agg.csmfNeonate, top=20, title="Aggregated COD distribution")
        plotVA(rv$fitNeonate, top=input$topDeaths)
      }
    }
  })

  ## Child
  #### Summarize and print output
  output$titleSummaryChild <- renderText({
    ## if(!is.null(rv$fitChild)){
    "Summary of Results for Children"
    ## }
  })
  output$emptySummaryChild <- renderText({
    if(is.null(rv$child)){
      "No Summary for Children (not enough deaths for analysis)"
    }
  })
  output$summaryChild <- renderPrint({
    algorithm <- selectedAlgorithm()
    if(!is.null(rv$fitChild)){
      if(algorithm==1 & is.null(rv$fitChild$HIV)){
        ## sumChild <- printIndiv(rv$fitChild, rv$agg.csmfChild)
        print(summary(rv$fitChild, top=input$topDeaths))
      }
      if(algorithm==2 & !is.null(rv$fitChild$HIV)){
        print(summary(rv$fitChild, top=input$topDeaths))
      }
      if(algorithm==3 & is.null(rv$fitChild$HIV)){
        ## sumChild <- printIndiv(rv$fitChild, rv$agg.csmfChild)
        print(summary(rv$fitChild, top=input$topDeaths))
      }
      if(algorithm==4 & !is.null(rv$fitChild$HIV)){
        print(summary(rv$fitChild, top=input$topDeaths))
      }
    }
  })
  #### Create plot (children)
  output$titlePlotChild <- renderText({
    ## if(!is.null(rv$fitChild)){
    "CSMF Plot for Children"
    ## }
  })
  output$emptyPlotChild <- renderText({
    if(is.null(rv$child)){
      "No Plot for Children (not enough deaths for analysis)"
    }
  })
  output$plotChild <- renderPlot({
    if(!is.null(rv$fitChild)){
      if(input$algorithm=="InterVA4" & !is.null(rv$fitNeonate$HIV)){
        ## CSMF(rv$fitChild, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF2(rv$fitChild, top=input$topDeaths)
      }
      if(input$algorithm=="InterVA5" & !is.null(rv$fitNeonate$HIV)){
        ## CSMF(rv$fitChild, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF5(rv$fitChild, top.plot=input$topDeaths, InterVA.rule = TRUE)
      }
      if(input$algorithm=="InSilicoVA_2012" & is.null(rv$fitNeonate$HIV)){
        ## indivplot(rv$agg.csmfChild, top=20, title="Aggregated COD distribution")
        plotVA(rv$fitChild, top=input$topDeaths)
      }
      if(input$algorithm=="InSilicoVA_2016" & is.null(rv$fitNeonate$HIV)){
        ## indivplot(rv$agg.csmfChild, top=20, title="Aggregated COD distribution")
        plotVA(rv$fitChild, top=input$topDeaths)
      }
    }
  })

  ## Adult
  #### Summarize and print output
  output$titleSummaryAdult <- renderText({
    ## if(!is.null(rv$fitAdult)){
    "Summary of Results for Adults"
    ## }
  })
  output$emptySummaryAdult <- renderText({
    if(is.null(rv$adult)){
      "No Summary for Adults (not enough deaths for analysis)"
    }
  })
  output$summaryAdult <- renderPrint({
    algorithm <- selectedAlgorithm()
    if(!is.null(rv$fitAdult)){
      if(algorithm==1 & is.null(rv$fitAdult$HIV)){
        ## sumAdult <- printIndiv(rv$fitAdult, rv$agg.csmfAdult)
        print(summary(rv$fitAdult, top=input$topDeaths))
      }
      if(algorithm==2 & !is.null(rv$fitAdult$HIV)){
        print(summary(rv$fitAdult, top=input$topDeaths))
      }
      if(algorithm==3 & is.null(rv$fitAdult$HIV)){
        ## sumAdult <- printIndiv(rv$fitAdult, rv$agg.csmfAdult)
        print(summary(rv$fitAdult, top=input$topDeaths))
      }
      if(algorithm==4 & !is.null(rv$fitAdult$HIV)){
        print(summary(rv$fitAdult, top=input$topDeaths))
      }
    }
  })
  #### Create plot
  output$titlePlotAdult <- renderText({
    ## if(!is.null(rv$fitAdult)){
    "CSMF Plot for Adults"
    ## }
  })
  output$emptyPlotAdult <- renderText({
    if(is.null(rv$adult)){
      "No Plot for Adults (not enough deaths for analysis)"
    }
  })
  output$plotAdult <- renderPlot({
    if(!is.null(rv$fitAdult)){
      if(input$algorithm=="InterVA4" & !is.null(rv$fitAdult$HIV)){
        ## CSMF(rv$fitAdult, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF2(rv$fitAdult, top=input$topDeaths)
      }
      if(input$algorithm=="InterVA5" & !is.null(rv$fitAdult$HIV)){
        ## CSMF(rv$fitAdult, top.plot=input$topDeaths, InterVA.rule=TRUE, min.prob=.001)
        CSMF5(rv$fitAdult, top.plot=input$topDeaths, InterVA.rule = TRUE)
      }
      if(input$algorithm=="InSilicoVA_2012" & is.null(rv$fitAdult$HIV)){
        ## indivplot(rv$agg.csmfAdult, top=20, title="Aggregated COD distribution")
        plotVA(rv$fitAdult, top=input$topDeaths)
      }
      if(input$algorithm=="InSilicoVA_2016" & is.null(rv$fitAdult$HIV)){
        ## indivplot(rv$agg.csmfAdult, top=20, title="Aggregated COD distribution")
        plotVA(rv$fitAdult, top=input$topDeaths)
      }
    }
  })

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
  shinyjs::disable("downloadWarnings")
}
