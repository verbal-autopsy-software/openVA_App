library(shinytest)
app <- ShinyDriver$new('..')
app$snapshotInit('testInterVA5', screenshot = FALSE)

app$setInputs(
  algorithm = 'InterVA5',
  HIV = 'v',
  Malaria = 'v')

app$uploadFile(readIn = 'small.csv')

app$setInputs(processMe = 'click', wait_ = FALSE, values_ = FALSE)

app$waitForValue("downloadCOD1", iotype = 'output', ignore = list(NULL))

vals <- app$getAllValues()
str(vals)
