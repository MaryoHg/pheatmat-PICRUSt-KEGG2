## THIS FUNCTION TRANSPOSE DATA FROM WIDE TO "TRANSPOSE WIDE FORMAT".
## HOW TO USE: transpose.table (data.frame)
## The format use has input is as follow:
##      Taxon   Sample_A        Sample_B        Sample_C        Sample_n
##      taxoni  count           count           count           count_n
##	...	count		count		count		count_n
##	taxon_j	count		count		count		count_n


transpose.any.table <- function (input.tab) {
        
        ## new data as input
        input.table <- input.tab 
        
        ## take the Taxon names from first colum
        taxon_names <- input.table[,1]
        
        ## Transposing the matrix (only counts):
        the.counts <- as.data.frame (t(input.table[,-1]))
        
        ## Rename the header column with the appropiate Taxon variable,
        colnames (the.counts) <- taxon_names
        
        ## print new dataframe:
        print (the.counts)
        
        #Save the new table: This data is ordered by rownames   
        #print (the.counts[order(rownames(the.counts)),])
}
