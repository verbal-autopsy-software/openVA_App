#' @import shiny
server <- function(input, output, session) {

  options(width = 100)
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
    shinyjs::disable("downloadWarnings")

    withProgress(value = 0, {  ## HERE move this to inside parLapply

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
        # records <- CrossVA::odk2openVA(getData())
        # records$ID <- getData()$meta.instanceID
        write.csv(getData(), file = 'tmpOut.csv', row.names = FALSE)
        pyAlg <- ifelse(input$algorithm == "InSilicoVA", "InsillicoVA", "InterVA5")
        pyCall <- paste0('/usr/local/bin/pycrossva-transform AUTODETECT ',
                         pyAlg, 'tmpOut.csv --dst pyOut.csv')
        system(pyCall)
        records <- read.csv('pyOut.csv', stringsAsFactors = FALSE)
      } else{
        records <- getData()
      }

      # object needed to render table of demographic variables
      # HERE -- create helper function that does table for age and sex
      names(records) <- tolower(names(records))
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
      
      # nCores <- min(parallel::detectCores() - 1, 6)
      # if (nCores == 0) nCores <- 1
      # cl <- parallel::makeCluster(nCores)
      #for (i in 1:length(namesRuns)) {
      lapply(1:length(namesRuns), function (i) {
      # parallel::clusterEvalQ(cl, {
      #   library(shiny)
      #   library(openVA)
      # })
      # parallel::clusterExport(cl, 
      #                         varlist = c("rv", "namesRuns", "records"), 
      #                         envir = environment())
      # parallel::parLapply(cl, 1:length(namesRuns), function (i) {
        # HERE -- create a cleanup function
        # HERE -- only run data check once (if possible)
        tmpNameRun <- namesRuns[i]
        groupName <- gsub('^(.)', '\\U\\1', tmpNameRun, perl = TRUE)
        
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
        modelArgs$data <- records[get(namesRuns[i]), ]
        rvName <- paste0("fit", groupName)
        incProgress(0.5/length(namesRuns),
                    detail = paste("Analyzing", groupName, "deaths")
                    )
        okRun <- try(
          rv[[rvName]] <- do.call(openVA::codeVA, modelArgs)
          )
        incProgress(0.5/length(namesRuns),
                    detail = paste("Completed analysis using", groupName, "deaths")
                    )
        # produce outputs
        if (!is.null(rv[[rvName]])) {

          rv$indivCOD <- indivCOD(rv[[rvName]], top = 3)
          cat("Analysis with ", groupName, "\t", date(), "\n", file = warningFileName)
          file.append(warningFileName, ovaLogFileName)
          file.remove(ovaLogFileName)
          
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
          emptySummary <- paste0("emptySummary", groupName)
          output[[emptySummary]] <- renderText({
            if (is.null(rv[[rvName]])) {
              paste("No Summary for", groupName, "(not enough deaths for analysis)")
            }
          })
          summaryGrp <- paste0("summary", groupName)
          output[[summaryGrp]] <- renderPrint({
            if (!is.null(rv[[rvName]])) {
              print(summary(rv[[rvName]], top = newTop))
              #print(summary(rv[[rvName]], top = input$topDeaths))
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
        }
        if(is.null(rv[[rvName]])) rv[[tmpNameRun]] <- NULL
      })
      # parallel::stopCluster(cl)
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
  shinyjs::disable("downloadWarnings")
}
