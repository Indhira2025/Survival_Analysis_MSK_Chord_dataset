#install.packages(c("shiny","shinydashboard","data.table","maftools","pheatmap","circlize", "survival","survminer","DT" ))



#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("maftools")


library(shiny)
library(shinydashboard)
library(data.table)
library(ggplot2)
library(maftools)
library(pheatmap)
library(circlize)
library(survival)
library(survminer)
library(DT)
library(dplyr)

getwd()


clin_patient <- fread("Data/data_clinical_patient.txt", skip = 4)
clin_sample  <- fread("Data/data_clinical_sample.txt", skip = 4)
muta <- fread("Data/data_mutations.txt")
cna  <- fread("Data/data_cna.txt")
sv   <- fread("Data/data_sv.txt")


# Create survival binary variable
clin_patient$OS_STATUS_BIN <- 
  ifelse(grepl("DECEASED", clin_patient$OS_STATUS), 1, 0)
clin_patient <- as.data.table(clin_patient)  # make sure it's a data.table

##UI

ui <- dashboardPage(
  dashboardHeader(title = "MSK-IMPACT Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      selectInput("cancer_type",
                  "Select Cancer Type:",
                  choices = sort(unique(clin_sample$CANCER_TYPE)),
                  selected = unique(clin_sample$CANCER_TYPE)[1]),
      
      menuItem("Overview", tabName = "overview"),
      menuItem("Oncoplot", tabName = "oncoplot"),
      menuItem("CNA", tabName = "cna"),
      menuItem("Fusions", tabName = "fusion"),
      menuItem("TMB", tabName = "tmb"),
      menuItem("Survival", tabName = "survival")
    )
  )
  ,
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "overview",
              fluidRow(
                valueBoxOutput("totalPatients"),
                valueBoxOutput("totalSamples"),
                valueBoxOutput("totalMutations")
              ),
              fluidRow(
                box(width = 12,
                    DTOutput("cancerTable"))
              )
      ),
      tabItem(tabName = "oncoplot",
              box(width = 12,
                  plotOutput("oncoplotPlot", height = 700))
      ),
      tabItem(tabName = "cna",
              box(width = 12,
                  plotOutput("cnaPlot", height = 800))
      ),
      tabItem(tabName = "fusion",
              box(width = 12,
                  plotOutput("fusionPlot", height = 700))
      ),
      tabItem(tabName = "tmb",
              box(width = 12,
                  plotOutput("tmbPlot", height = 500))
      ),
      tabItem(tabName = "survival",
              box(width = 12,
                  plotOutput("survivalPlot", height = 600))
      )
    )
  )
)



 server <- function(input, output) {
    
    # Reactive filtered samples
    filtered_samples <- reactive({
      clin_sample[CANCER_TYPE == input$cancer_type]
    })
    
    # Total Patients
    output$totalPatients <- renderValueBox({
      valueBox(
        length(unique(filtered_samples()$PATIENT_ID)),
        "Patients",
        icon = icon("users"),
        color = "blue"
      )
    })
    
    # Total Samples
    output$totalSamples <- renderValueBox({
      valueBox(
        nrow(filtered_samples()),
        "Samples",
        icon = icon("vial"),
        color = "green"
      )
    })
    
    # Total Mutations
    output$totalMutations <- renderValueBox({
      
      n_mut <- muta[
        Tumor_Sample_Barcode %in% filtered_samples()$SAMPLE_ID,
        .N
      ]
      
      valueBox(
        n_mut,
        "Mutations",
        icon = icon("dna"),
        color = "red"
      )
    })
    
    # Data Table
    output$cancerTable <- renderDT({
      datatable(filtered_samples(),
                options = list(pageLength = 10,
                               scrollX = TRUE))
    })
    output$oncoplotPlot <- renderPlot({
      
      samples <- filtered_samples()$SAMPLE_ID
      
      muta_sub <- muta[
        Tumor_Sample_Barcode %in% samples
      ]
      
      if(nrow(muta_sub) == 0) return(NULL)
      
      maf_object <- read.maf(muta_sub)
      
      oncoplot(maf = maf_object, top = 20)
      
    })
    
    output$tmbPlot <- renderPlot({
      
      data <- clin_sample[
        CANCER_TYPE == input$cancer_type
      ]
      
      ggplot(data, aes(y = TMB_NONSYNONYMOUS)) +
        geom_boxplot() +
        theme_bw()
      
    })
  
    
    output$cnaPlot <- renderPlot({
      samples <- filtered_samples()$SAMPLE_ID
      samples <- samples[samples %in% colnames(cna)]  # keep only matching
      
      if(length(samples) == 0) return(NULL)
      
      cna_sub <- cna[, c("Hugo_Symbol", samples), with = FALSE]
      
      mat <- as.matrix(cna_sub[, -1, with = FALSE])
      rownames(mat) <- cna_sub$Hugo_Symbol
      
      var_genes <- head(order(apply(mat, 1, var), decreasing = TRUE), 30)
      
      pheatmap(mat[var_genes, ],
               show_colnames = FALSE,
               main = paste("Top 30 variable genes for", input$cancer_type))
    })
    
    
    output$fusionPlot <- renderPlot({
      
      sv_sub <- sv[Sample_Id %in% filtered_samples()$SAMPLE_ID]
      
      if(nrow(sv_sub) == 0) return(NULL)
      
      top_pairs <- sv_sub[, .N, by = .(Site1_Hugo_Symbol, Site2_Hugo_Symbol)]
      top_pairs <- top_pairs[order(-N)]
      top_pairs <- top_pairs[1:min(20, nrow(top_pairs))]  # max 20
      
      top_pairs <- top_pairs[!is.na(Site1_Hugo_Symbol) & !is.na(Site2_Hugo_Symbol)]
      if(nrow(top_pairs) == 0) return(NULL)
      
      chordDiagram(top_pairs, grid.col = NULL, transparency = 0.5)
      
    })
    
    
    output$survivalPlot <- renderPlot({
      
      # Get patient IDs for the selected cancer type
      patients <- filtered_samples()$PATIENT_ID
      
      # Filter clinical patient table for these patients
      surv_data <- clin_patient[PATIENT_ID %in% patients]  # THIS MUST BE INSIDE renderPlot
      print(head(surv_data))  # debug
      # If no patients, return NULL
      if(nrow(surv_data) == 0) return(NULL)
      
      # Fit Kaplan-Meier curve
      fit <- survfit(Surv(OS_MONTHS, OS_STATUS_BIN) ~ 1, data = surv_data)
      
      # Plot
      ggsurvplot(fit, risk.table = TRUE)$plot
      
    })
    
    
    
    
  
    }
  
 

shinyApp(ui, server)


