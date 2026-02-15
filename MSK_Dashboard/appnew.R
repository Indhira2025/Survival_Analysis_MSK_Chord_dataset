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

# --- Load Data ---
clin_patient <- fread("Data/data_clinical_patient.txt", skip = 4)
clin_sample  <- fread("Data/data_clinical_sample.txt", skip = 4)
muta <- fread("Data/data_mutations.txt")
cna  <- fread("Data/data_cna.txt")
sv   <- fread("Data/data_sv.txt")

# Ensure survival variable
clin_patient$OS_STATUS_BIN <- ifelse(grepl("DECEASED", clin_patient$OS_STATUS), 1, 0)
clin_patient <- as.data.table(clin_patient)
clin_patient <- as.data.table(clin_patient)

# Merge by PATIENT_ID
merged_data <- merge(clin_sample, clin_patient, by = "PATIENT_ID", all.x = TRUE)


# --- UI ---
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
  ),
  
  dashboardBody(
    tabItems(
      tabItem(tabName = "overview",
              fluidRow(
                valueBoxOutput("totalPatients"),
                valueBoxOutput("totalSamples"),
                valueBoxOutput("totalMutations")
              ),
              fluidRow(
                box(width = 12, DTOutput("cancerTable"))
              )
      ),
      tabItem(tabName = "oncoplot",
              box(width = 12, plotOutput("oncoplotPlot", height = 700))
      ),
      tabItem(tabName = "cna",
              box(width = 12, plotOutput("cnaPlot", height = 800))
      ),
      tabItem(tabName = "fusion",
              box(width = 12, plotOutput("fusionPlot", height = 700))
      ),
      tabItem(tabName = "tmb",
              box(width = 12, plotOutput("tmbPlot", height = 500))
      ),
      tabItem(tabName = "survival",
              box(width = 12, plotOutput("survivalPlot", height = 600))
      )
    )
  )
)

# --- Server ---
server <- function(input, output) {
  
  # --- Reactive filtered samples for selected cancer type ---
  filtered_samples <- reactive({
    merged_data[CANCER_TYPE == input$cancer_type]
  })
  
  # --- Overview Value Boxes ---
  output$totalPatients <- renderValueBox({
    n_pat <- length(unique(filtered_samples()$PATIENT_ID))
    valueBox(n_pat, "Patients", icon = icon("users"), color = "blue")
  })
  
  output$totalSamples <- renderValueBox({
    n_samp <- nrow(filtered_samples())
    valueBox(n_samp, "Samples", icon = icon("vial"), color = "green")
  })
  
  output$totalMutations <- renderValueBox({
    samp_ids <- filtered_samples()$SAMPLE_ID
    n_mut <- muta[Tumor_Sample_Barcode %in% samp_ids, .N]
    valueBox(n_mut, "Mutations", icon = icon("dna"), color = "red")
  })
  
  output$cancerTable <- renderDT({
    datatable(filtered_samples(), options = list(pageLength = 10, scrollX = TRUE))
  })
  
  # --- Oncoplot ---
  output$oncoplotPlot <- renderPlot({
    samp_ids <- filtered_samples()$SAMPLE_ID
    muta_sub <- muta[Tumor_Sample_Barcode %in% samp_ids]
    if(nrow(muta_sub) == 0) return(NULL)
    maf_obj <- read.maf(muta_sub)
    oncoplot(maf_obj, top = 20)
  })
  
  # --- CNA Heatmap ---
  output$cnaPlot <- renderPlot({
    samp_ids <- filtered_samples()$SAMPLE_ID
    samp_ids <- samp_ids[samp_ids %in% colnames(cna)]
    if(length(samp_ids) == 0) return(NULL)
    
    cna_sub <- cna[, c("Hugo_Symbol", samp_ids), with = FALSE]
    mat <- as.matrix(cna_sub[, -1, with = FALSE])
    rownames(mat) <- cna_sub$Hugo_Symbol
    var_genes <- head(order(apply(mat,1,var), decreasing = TRUE), 30)
    
    pheatmap(mat[var_genes, ], show_colnames = FALSE,
             main = paste("Top 30 variable genes for", input$cancer_type))
  })
  
  # --- Fusion Chord Diagram ---
  output$fusionPlot <- renderPlot({
    sv_sub <- sv[Sample_Id %in% filtered_samples()$SAMPLE_ID]
    if(nrow(sv_sub) == 0) return(NULL)
    
    top_pairs <- sv_sub[, .N, by = .(Site1_Hugo_Symbol, Site2_Hugo_Symbol)]
    top_pairs <- top_pairs[order(-N)][1:min(20, .N)]
    top_pairs <- top_pairs[!is.na(Site1_Hugo_Symbol) & !is.na(Site2_Hugo_Symbol)]
    if(nrow(top_pairs) == 0) return(NULL)
    
    chordDiagram(top_pairs, grid.col = NULL, transparency = 0.5)
  })
  
  # --- Tumor Mutation Burden ---
  output$tmbPlot <- renderPlot({
    data <- filtered_samples()
    if(nrow(data) == 0) return(NULL)
    
    ggplot(data, aes(y = TMB_NONSYNONYMOUS)) +
      geom_boxplot(fill = "steelblue") +
      theme_bw() +
      ylab("TMB (Non-synonymous)") +
      ggtitle(paste("TMB for", input$cancer_type))
  })
  
  # --- Survival Plot ---
  output$survivalPlot <- renderPlot({
    datasurv <- merged_data[merged_data$CANCER_TYPE == input$cancer_type, ]  
    
    if(nrow(datasurv) == 0) return(NULL)
    
    fit <- do.call(survfit, list(
      formula = Surv(as.numeric(OS_MONTHS), OS_STATUS_BIN) ~ 1, 
      data = as.data.frame(datasurv)
    ))
    
    ggsurvplot(fit,
               risk.table = TRUE,
               conf.int = TRUE,
               xlab = "Months",
               ylab = "Survival Probability",
               title = paste("Overall Survival for", input$cancer_type))$plot
  })
}



# --- Run App ---
shinyApp(ui, server)

