
modelTempsSubsurfaceDates = function(distances,
                                     directory,
                                     fileNamePrefix,
                                     fileNameSuffix,
                                     mdyHMS) {
  file.names = paste0(directory, fileNamePrefix, distances, fileNameSuffix)
  theData = lapply(
    file.names,
    read.table,
    header = F,
    stringsAsFactors = F,
    skip = 3
  )
  theData = lapply(theData,
                   `names<-`,
                   value = c("Time", "C", "X", "Y", "Z", "Node")
  )
  times = theData[[1]][["Time"]]
  times = times + lubridate::mdy_hms(mdyHMS)
  temperatureData = lapply(
    theData,
    `[[`, i = "C"
  )
  names(temperatureData) = paste0("temp", distances)
  tempTS = zoo::zoo(
    as.data.frame(temperatureData),
    order.by = times
  )
}
