methods::setGeneric(
	"getStats", 
	function(IntLimObject, ...)
    	base::standardGeneric("getStats")
)


methods::setMethod(
       f = "getStats",
       signature = c("MultiDataSet"),
       definition = function(IntLimObject, ...) {
		incommon<-MultiDataSet::commonSamples(IntLimObject)
		mystats <- NULL
		mytypes <- names(Biobase::assayData(IntLimObject))
    if(any(mytypes=="expression") && any(mytypes=="metabolite")){
      mystats <- data.frame(Num_Genes = nrow(Biobase::fData(IntLimObject[["expression"]])),
                            Num_Metabolites = nrow(Biobase::fData(IntLimObject[["metabolite"]])),
                            Num_Samples_withGeneExpression = ncol(Biobase::assayDataElement(IntLimObject[["expression"]], 
                                                                                            'exprs')),
                            Num_Samples_withMetabolomics = ncol(Biobase::assayDataElement(IntLimObject[["metabolite"]], 
                                                                                          'metabData')),
                            Num_Samples_inCommon = ncol(Biobase::assayDataElement(incommon[["expression"]], 'exprs'))
                            
		  )
    } else if(any(mytypes=="metabolite")){
      mystats <- data.frame(Num_Metabolites = nrow(Biobase::fData(IntLimObject[["metabolite"]])),
                            Num_Samples_withMetabolomics = ncol(Biobase::assayDataElement(
                              IntLimObject[["metabolite"]], 'metabData'))      
      )
    } else if(any(mytypes=="expression")){
      mystats <- data.frame(Num_Genes = nrow(Biobase::fData(IntLimObject[["expression"]])),
                            Num_Samples_withGeneExpression = ncol(Biobase::assayDataElement(IntLimObject[["expression"]], 
                                                                                            'exprs'))
      )
    }
    return(mystats)
})
