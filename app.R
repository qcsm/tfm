library(shiny)
library(WGCNA)
library(ape)

# First we load alll the tumor expression data
# ETL process got the data from TCGA
# Data was preprocessed and stored in a persistent R data structure
load("BRCA_tumor.Rda")
load("BLCA_tumor.Rda")
load("COAD_tumor.Rda")
load("PRAD_tumor.Rda")
load("LUAD_tumor.Rda")
load("LUSC_tumor.Rda")
load("HNSC_tumor.Rda")
load("KICH_tumor.Rda")
load("KIRC_tumor.Rda")
load("KIRP_tumor.Rda")
load("LIHC_tumor.Rda")
load("THCA_tumor.Rda")

# This code block can be used to log the tumor data load
# but it is pretty fast to load from the Rda datastructure
# and so it was commented out
# file.names <- list.files(pattern = "*_tumor.Rda")
# for( i in 1:length(file.names)){
#   print(paste0("Loading data file: ", file.names[i]))
#   load(file.names[i])
# }

# This are important WGCNA settings
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

# Functon to search for saved WGCNA networks in the script directory
# Matrix and Network data is stored in an RDS R datastructure
# Matrix and Network data is stored with the parameter information in the filename
# filename: [matrix|bwnet]_[ngenes]_[counts_filter]_[soft_threshold].rds
# bwnet_LUAD_448_80_4.rds
# matrix_LUAD_448_80.rds
get_saved_bwnets <- function(){
  list.files(pattern = "^bwnet_.*rds$")
}

# Function to get the soft threshold from file name
getSft <- function(x){
  f <- strsplit(x, "[_\\.]")
  f[[1]][5]
}

# Function to get the tumor label (shortname)
getTumorLabel <- function(x){
  f <- strsplit(x, "_")
  f[[1]][2]
} 

getTumorName <- function(x){
  f <- strsplit(x, "_")
  f[[1]][2]
}

# Functoin to get a Matrix from a persistent RDS file
getWGCNAMatrix <- function(x){
  f <- strsplit(x, "_")
  y <- paste0(paste("matrix", f[[1]][2], f[[1]][3], f[[1]][4], sep = "_"), ".rds")
  matrix <- readRDS(y)
  matrix
}

# Functoin to get a Network from a persistent RDS file
getWGCNANetwork <- function(x){
  net <- readRDS(x)
  net
}

# Function to add links to NCBI for accessing gene information
addNCBILinks <- function(x){
  lapply(x, function(y) {
    gene_name_id <- strsplit(y, "\\.")
    name <- gene_name_id[[1]][1]
    id <- gene_name_id[[1]][2]
    a(name, href=paste("https://www.ncbi.nlm.nih.gov/gene/",id, sep=""), target="_blank")
  })
}

# Shiny user interface
ui <- fluidPage(
  
  tags$head(
    tags$style(HTML("hr {border-top: 1px solid #000000;}"))
  ),
  
  # Sidebar layout with a input and output definitions
  sidebarLayout(
    
    # Inputs
    sidebarPanel(
      
      h4("Create Weighted Gene Co-expression Networks"),
      h5("Select your parameters for filtering and soft-threshold. You can also perform a soft-thresold analysis if needed"),
      hr(),
      selectInput(inputId = "tumorName", 
                  label = "Select Tumor:",
                  choices = c("Breast Invasive Carcinoma" = "BRCA",
                              "Bladder Urothelial Carcinoma" = "BLCA",
                              "Colon Adenocarcinoma" = "COAD",
                              "Prostate Adenocarcinoma" = "PRAD", 
                              "Head and Neck Squamous Cell Carcinoma" = "HNSC",
                              "Lung Adenocarcinoma" = "LUAD", 
                              "Lung Squamous Cell Carcinoma" = "LUSC",
                              "Kidney Chromophobe" = "KICH",
                              "Kidney Renal Clear Cell Carcinoma" = "KIRC",
                              "Kidney Renal Papillary Cell Carcinoma" = "KIRP",
                              "Liver Hepatocellular Carcinoma" = "LIHC",
                              "Thyroid Carcinoma" = "THCA"
                  )
      ),
      sliderInput("ngenes", label = h5("Reduce to the most variable genes"), min = 100, 
                  max = 100, value = 1000, step = 1),
      sliderInput("pctfilter", label = h5("Filter genes with with counts in less than x% of cases"), min = 10, 
                  max = 100, value = 80, step = 1),
      actionButton('str','Calculate soft-threshold'),
      sliderInput("sft", label = h5("Choose a soft-threshold"), min = 1, 
                  max = 20, value = 1, step = 1),
      actionButton('wgcna','Calculate WGCNA network'),
      hr(),
      h4("Consenus Analysis"),
      h5("You need to create two or more Weighted Gene Co-expression Networs to perform a consensus analysis"),
      hr(),
      selectInput(inputId = "consenseTumors", 
                  label = "Select WGCNA networks for consensus analysis:",
                  multiple = TRUE,
                  choices = get_saved_bwnets()
      ),
      actionButton('cons','Perform consensus analysis'),
      hr(),
      h4("Export genes"),
      h5("After consensus analysis, shared and specific genes between any 2 modules can be exported for further analysis"),
      hr(),
      selectInput(inputId = "moduleA", 
                  label = "Select Tumor Module:",
                  choices = "None"),
      selectInput(inputId = "moduleCons", 
                  label = "Select Consensus Module:",
                  choices = "None"),
      actionButton('mods','Get genes!')
    ),
    # Output
    mainPanel(
      uiOutput(outputId = "tumorSelected"),
      plotOutput("str"),
      plotOutput("str2"),
      plotOutput("phylo"),
      htmlOutput("title1"),
      htmlOutput("mod1genes"),
      htmlOutput("title2"),
      htmlOutput("mod2genes"),
      htmlOutput("title3"),
      htmlOutput("mod3genes"),
      plotOutput("consensusAnalysis")
    )
  )
)

# Shiny server side
server = function(input, output, session){
  outVar = reactive({
    nrow(get(paste(input$tumorName,"_tumor",sep = "")))
  })

  observe({
    updateSliderInput(session, "ngenes",
                      max = outVar(),
                      value = 1000
    )})
  observe({
    if(input$str > 0) {
      isolate({
        print(input$tumorName)
        tumorData <- get(paste(input$tumorName,"_tumor",sep = ""))
        print(dim(tumorData))
        tumorData = tumorData[apply(tumorData,1,function(x) sum(x==0))<ncol(tumorData)*input$pctfilter/100,]
        WGCNA_matrix = t(tumorData[order(apply(tumorData,1,mad), decreasing = T)[1:input$ngenes],])
        print(dim(WGCNA_matrix))
        
        saveRDS(WGCNA_matrix, file = paste0(paste("matrix",input$tumorName,input$ngenes,input$pctfilter,sep="_"), ".rds"))
        
        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)
        cex1 = 0.9;
        
        output$str <- renderPlot({
          # Scale-free topology fit index as a function of the soft-thresholding power
          plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
               xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
               main = paste("Scale independence"));
          text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
               labels=powers,cex=cex1,col="red");
          # this line corresponds to using an R^2 cut-off of h
          abline(h=0.90,col="red")
        })
        output$str2 <- renderPlot({
          # Mean connectivity as a function of the soft-thresholding power
          plot(sft$fitIndices[,1], sft$fitIndices[,5],
               xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
               main = paste("Mean connectivity"))
          text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
        })
      })
    }
  })
  observe({
    if(input$wgcna > 0) {
      isolate({
        print(input$tumorName)
        tumorData <- get(paste(input$tumorName,"_tumor",sep = ""))
        print(dim(tumorData))
        tumorData = tumorData[apply(tumorData,1,function(x) sum(x==0))<ncol(tumorData)*input$pctfilter/100,]
        WGCNA_matrix = t(tumorData[order(apply(tumorData,1,mad), decreasing = T)[1:input$ngenes],])
        print(dim(WGCNA_matrix))
        
        saveRDS(WGCNA_matrix, file = paste0(paste("matrix",input$tumorName,input$ngenes,input$pctfilter,sep="_"), ".rds"))
        
        print(paste0("Matrix dimension: ",dim(WGCNA_matrix)))
        
        bwnet = blockwiseModules(WGCNA_matrix, maxBlockSize = 5000,
                                 power = input$sft, TOMType = "unsigned", minModuleSize = 30,
                                 reassignThreshold = 0, mergeCutHeight = 0.25,
                                 numericLabels = TRUE,
                                 verbose = 3)
        
        saveRDS(bwnet, file = paste0(paste("bwnet",input$tumorName,input$ngenes,input$pctfilter,input$sft,sep="_"), ".rds"))
        # saveRDS(WGCNA_matrix, file = paste0(paste("matrix",input$tumorName,input$ngenes,input$pctfilter,sep="_"), ".rds"))
        updateSelectInput(session,
                    inputId = "consenseTumors", 
                    choices = get_saved_bwnets()
        )
        
        MEs <- bwnet$MEs
        MEDiss = 1-cor(MEs)
        METree = hclust(as.dist(MEDiss), method = 'average')
        
        output$phylo <- renderPlot({
            plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
            tiplabels(frame = 'circle',col='black', text=unique(bwnet$colors), bg = levels(as.factor(labels2colors(bwnet$colors))), cex = 1.7 - input$ngenes / nrow(get(paste(input$tumorName,"_tumor",sep = ""))))
        })
      })
    }
  })
  observe({
    if(input$cons > 0) {
      isolate({
        selections <- strsplit(input$consenseTumors, split = " ")
        
        tumor_A <- getTumorName( selections[[1]] )
        tumor_B <- getTumorName( selections[[2]] )
        tumor_A_matrix <- getWGCNAMatrix(selections[[1]])
        tumor_A_network <- getWGCNANetwork(selections[[1]])
        tumor_A_genes <- colnames(tumor_A_matrix)
        tumor_B_name <- getTumorName(selections[[2]])
        tumor_B_data <- get( paste0(tumor_B_name,'_tumor') )
        tumor_B_matrix <- t(tumor_B_data)
        # Genes in both matrices must be the same
        # To do this, one can take advantange of precalculated tumor_A matrix and network
        # Then take the whole tumor_B expression dataset and produce a corresponding matrix
        tumor_A_genes <- tumor_A_genes[!is.na(tumor_A_genes)]
        tumor_B_matrix <- tumor_B_matrix[,tumor_A_genes]
        # The filtering process may have done some NAs for genes with low counts
        # We need to add the same NAs in tumor_B matrix in order to avoid differences in dimension
        # This differences may affect later during the Fisher test, so we solve any difference beforehand
        gap_length <- ncol(tumor_A_matrix) - ncol(tumor_B_matrix)
        tumor_B_matrix <- cbind(tumor_B_matrix, matrix(rep(0,nrow(tumor_B_matrix)*gap_length), nrow=nrow(tumor_B_matrix)))
        write.csv(x = colnames(tumor_B_matrix), file = "tumor_B_genes.2")
        print(length(colnames(tumor_A_matrix)))
        print(length(colnames(tumor_B_matrix)))
        
        nSets = 2;
        
        setLabels = c(getTumorLabel(selections[[1]]), 
                      getTumorLabel(selections[[2]]))
        shortLabels = c(getTumorLabel(selections[[1]]), 
                        getTumorLabel(selections[[2]]))
        
        multiExpr = vector(mode = "list", length = nSets)
        multiExpr[[1]] = list(data = as.data.frame(tumor_A_matrix))
        multiExpr[[2]] = list(data = as.data.frame(tumor_B_matrix))
        
        exprSize = checkSets(multiExpr)
        
        gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
        print(gsg$allOK)
        
        sampleTrees = list()
        for (set in 1:nSets)
        {
          sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
        }
        
        powers <- c(getSft(selections[[1]]), getSft(selections[[2]]))
        print(powers)
        
        net = blockwiseConsensusModules(
          multiExpr, power = min(powers), minModuleSize = 30, deepSplit = 2,
          pamRespectsDendro = FALSE,
          #mergeCutHeight = 0.25, 
          numericLabels = TRUE,
          minKMEtoStay = 0,
          # You may want to save the TOM matrices
          # They can be easily converted to vertex data
          #saveTOMs = TRUE, 
          verbose = 5)
        
        saveRDS(net, file = "current_consensus_network.rds")
        
        consMEs = net$multiMEs
        moduleLabels = net$colors
        # Convert the numeric labels to color labels
        moduleColors = labels2colors(moduleLabels)
        consTree = net$dendrograms[[1]]
        
        tumorAMEs <- tumor_A_network$MEs
        tumorAcolors = tumor_A_network$colors
        tumorAMEs = orderMEs(tumorAMEs, greyName = "ME0")
        
        tumorAModuleLabels = substring(names(tumorAMEs), 3)
        consModuleLabels = substring(names(consMEs[[1]]$data), 3)
        # Convert the numeric module labels to color labels
        tumorAModules = labels2colors(as.numeric(tumorAModuleLabels))
        consModules = labels2colors(as.numeric(consModuleLabels))
        # Numbers of female and consensus modules
        nTumorAMods = length(tumorAModules)
        nConsMods = length(consModules)
        # Initialize tables of p-values and of the corresponding counts
        pTable = matrix(0, nrow = nTumorAMods, ncol = nConsMods)
        CountTbl = matrix(0, nrow = nTumorAMods, ncol = nConsMods)
        # Execute all pairwaise comparisons
        for (fmod in 1:nTumorAMods)
          for (cmod in 1:nConsMods)
          {
            tumor_aMembers = (labels2colors(tumorAcolors) == tumorAModules[fmod])
            consMembers = (moduleColors == consModules[cmod])
            print(paste(length(tumor_aMembers), length(consMembers), sep = ","))
            pTable[fmod, cmod] = -log10(fisher.test(tumor_aMembers, consMembers, alternative = "greater")$p.value)
            CountTbl[fmod, cmod] = sum(labels2colors(tumorAcolors) == tumorAModules[fmod] & moduleColors ==
                                         consModules[cmod])
          }
        
        # Truncate p values smaller than 10^{-50} to 10^{-50}
        pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
        pTable[pTable>50 ] = 50 ;
        # Marginal counts (really module sizes)
        tumor_aModTotals = apply(CountTbl, 1, sum)
        consModTotals = apply(CountTbl, 2, sum)
        
        updateSelectInput(session,
                          inputId = "moduleA", 
                          choices = tumorAModules
        )
        updateSelectInput(session,
                          inputId = "moduleCons", 
                          choices = consModules
        )
        
        # Actual plotting
        output$consensusAnalysis <- renderPlot({
            # sizeGrWindow(10,7 )
            par(mfrow=c(1,1))
            par(cex = 1.0)
            par(mar=c(8, 10.4, 2.7, 1)+0.8)
            # Use function labeledHeatmap to produce the color-coded table with all the trimmings
            labeledHeatmap(Matrix = pTable,
                           xLabels = paste(" ", consModules),
                           yLabels = paste(" ", tumorAModules),
                           colorLabels = TRUE,
                           xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
                           ySymbols = paste(getTumorLabel(selections[[1]]), tumorAModules, ": ", tumor_aModTotals, sep=""),
                           textMatrix = CountTbl,
                           colors = greenWhiteRed(100)[50:100],
                           main = paste("Correspondence of ",tumor_A," set-specific and ",tumor_A,"-",tumor_B," consensus modules", sep=""),
                           cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)
        }, height = length(tumorAModules)*80)
      })
    }
  })
  observe({
    if(input$mods > 0) {
      
      isolate({
        if( input$moduleA != "None" && input$moduleCons != "None" ){
          selections <- strsplit(input$consenseTumors, split = " ")
          
          tumor_A <- getTumorName(selections[[1]])
          tumor_A_network <- getWGCNANetwork(selections[[1]])
          tumor_A_matrix <- getWGCNAMatrix(selections[[1]])
          GeneNames <- colnames(tumor_A_matrix)
          net = readRDS("current_consensus_network.rds")
          
          consMEs = net$multiMEs
          moduleLabels = net$colors
          # Convert the numeric labels to color labels
          moduleColors = labels2colors(moduleLabels)
          tumorAMEs <- tumor_A_network$MEs
          tumorAcolors = tumor_A_network$colors
          tumorAMEs = orderMEs(tumorAMEs, greyName = "ME0")
          
          tumorAModuleLabels = substring(names(tumorAMEs), 3)
          consModuleLabels = substring(names(consMEs[[1]]$data), 3)
          # Convert the numeric module labels to color labels
          tumorAModules = labels2colors(as.numeric(tumorAModuleLabels))
          consModules = labels2colors(as.numeric(consModuleLabels))
          
          tumor_aMembers = (labels2colors(tumorAcolors) == input$moduleA)
          consMembers = (moduleColors == input$moduleCons)
          
          output$title1 <- renderText({
             HTML(paste("<b>Tumor:", tumor_A, "Module:", input$moduleA, "</b> - Gene list below:", sep = " "))
          })
  
          output$mod1genes <- renderUI({
            addNCBILinks( GeneNames[tumor_aMembers] )
          })
          
          output$title2 <- renderText({
            HTML(paste("<b>Consensus module:", input$moduleCons, "</b> - Gene list below:", sep = " "))
          })
          
          output$mod2genes <- renderUI({
            addNCBILinks( GeneNames[consMembers] )
          })
        }
      })
    }
  })
}

# Create a Shiny app object
shinyApp(ui = ui, server = server)
