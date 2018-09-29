library(RTCGA)
dir.create('data')

# check available files
# change LUAD for corresponding tumor label (from labels used in TCGA Data Portal repository)
datasets <- checkTCGA(what = 'DataSets', cancerType = 'LUAD', date = NULL)

# downloading rnaseq data
downloadTCGA( cancerTypes = 'LUAD', 
              dataSet = 'rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level',
              destDir = 'data' )

library(dplyr)
# shortening paths and directories
list.files( 'data/') %>% 
    file.path( 'data', .) %>%
    file.rename( to = substr(.,start=1,stop=50))

# reading data
list.files( 'data/') %>% 
    file.path( 'data', .) -> folder

folder %>%
    list.files %>%
    file.path( folder, .) %>%
    grep( pattern = 'illuminahiseq', x = ., value = TRUE) -> pathRNA
readTCGA( path = pathRNA, dataType = 'rnaseq' ) -> expression_matrix

# expression_matrix is the data to use in our Shiny app
# It has to be persisted in RDA format (it is done to speed up the app load)
# output file: <TUMOR_LABEL>_tumor.Rda
