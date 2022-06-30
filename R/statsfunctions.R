#' Get some stats after reading in data
#'
#' @param IntLimObject output of ReadData()
#' @return data.frame with some # of samples, features, etc.
#' @export
ShowStats <- function(IntLimObject) {
  if(length(IntLimObject@analyteType1)>0 && length(IntLimObject@analyteType2)>0){
    mystats <- data.frame(Num_Analytes_Type1 = nrow(IntLimObject@analyteType1),
                          Num_Analytes_Type2 = nrow(IntLimObject@analyteType2),
                          Num_Samples = nrow(IntLimObject@sampleMetaData))
  } else if(length(IntLimObject@analyteType1)>0){
    mystats <- data.frame(Num_Analytes_Type1 = nrow(IntLimObject@analyteType1),
                          Num_Samples = nrow(IntLimObject@sampleMetaData))      
  } else if(length(IntLimObject@analyteType2)>0){
    mystats <- data.frame(Num_Analytes_Type2 = nrow(IntLimObject@analyteType2),
                          Num_Samples = nrow(IntLimObject@sampleMetaData))
  }
  print(mystats)
}
