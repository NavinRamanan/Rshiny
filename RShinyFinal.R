## Author: Navin Ramanan
## BF591
## Final Project



library(shiny)
library(ggplot2)
library(colourpicker) 
library(tidyverse)
library(rlang)
library(dplyr)
library(DT) 
library(psych) 
library(genefilter)
library('RColorBrewer')
library(matrixStats)
library(pheatmap)
library('biomaRt')
library('fgsea')
library(abind)

###############################

#library(shinythemes)
options(shiny.maxRequestSize=30*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(

  titlePanel(h5("BF591 Final Project", align = "center")),
  titlePanel(h4("Navin Ramanan", align = "center")),
  titlePanel(h3("Nested Tab Input Control Example", align = "center")),
  
  # Creating Input Tabs
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(width = 0
    ),
    mainPanel(width = 12,
         tabsetPanel(type = "tabs",
              ##############  Tab1 : Samples  #############
              tabPanel(title = "Samples", fluid = TRUE,
                  sidebarLayout(
                      sidebarPanel(width=3,
                                  fileInput("Samples_file", "Sample File", accept = ".csv", placeholder = "sample_info.csv")
                                  ), #sidebarPanel
                      mainPanel(width = 9,
                          tabsetPanel(type = "tabs",
                              tabPanel("Summary", verbatimTextOutput("Samples_summary")),
                              #tabPanel("Table", tableOutput("Samples_table")),
                              tabPanel("Table", DT::dataTableOutput("Samples_table")),
                              tabPanel("Plots", 
                                    plotOutput("Samples_plot")
                                    #splitLayout(cellWidths = c("50%", "50%"), 
                                    #      plotOutput("Samples_plot"), 
                                    #      plotOutput("Samples_plot2")
                                    #)           
                              )
                          ) #tabsetPanel
                      ) #mainPanel
                   ) #sidebarLayout                
                ), #tabPanel (Samples)
                ##############  Tab2 : Counts  #############
                tabPanel("Counts", fluid = TRUE,
                    sidebarLayout(
                        sidebarPanel(width=3,
                                    fileInput("Counts_file", "Counts_file", accept = ".csv", placeholder = "norm_count_matrix.csv"),
                                    sliderInput(inputId="Counts_slider1", label="Select genes with at least X percentile of variance:", min=0, max=1, value=0.5, step = NULL,round = FALSE,ticks = TRUE),
                                    sliderInput(inputId="Counts_slider2", label="Select genes with at least X samples that are non-zero:", min=0, max=1, value=0.5, round = TRUE,ticks = TRUE),
                                    submitButton("Counts Submit")
                                    #uiOutput("counts_file")
                        ), #sidebarPanel
                        mainPanel(width = 9,
                            tabsetPanel(type = "tabs",
                              tabPanel("Counts_Summary", verbatimTextOutput("Counts_Summary")),
                              tabPanel("Diagnostic Plots",
                                        #plotOutput("Counts_diagnostic_plot"))
                                        splitLayout(cellWidths = c("50%", "50%"), 
                                                  plotOutput("Counts_diagnostic_plot1"), 
                                                  plotOutput("Counts_diagnostic_plot2")
                                        )
                              ),
                              tabPanel("Heatmap", plotOutput("Counts_heatmap")),
                              #tabPanel("PCA Plot", plotOutput("Counts_PCA_plot"))
                              tabPanel("PCA Plot",
                                  sidebarPanel(
                                        titlePanel(h5("PCA Scatter Plot")),
                                        selectInput(inputId = "xcol", label = "X Variable", choices = c("PC1", "PC2"), selected = "PC1"),
                                        selectInput(inputId = "ycol", label = "y Variable", choices = c("PC1", "PC2"), selected = "PC2"),
                                        submitButton("Submit Variables")
                                        
                                  ),
                                  mainPanel(
                                         plotOutput('Counts_PCA_plot')
                                  )
                              )
                            ) #tabsetPanel
                        ) #mainPanel
                     ) #sidebarLayout                
                   ),  #tabPanel("Counts"
                   ##############  Tab3 : DE  #############
                   tabPanel("DE", fluid = TRUE,
                        sidebarLayout(
                            sidebarPanel(width=3,
                                        fileInput("DE_file", "Differential Expression Analysis File", accept = ".csv", placeholder = "norm_count_matrix.csv")
                                        #submitButton("DE_Submit")
                            ), #sidebarPanel
                            mainPanel(width = 9,
                                tabsetPanel(type = "tabs", 
                                  tabPanel("DE Table1", DT::dataTableOutput("DE_Table1")),
                                  #tabPanel("DE_Table2", tableOutput("DE_Table2"))
                                  tabPanel("DE_Plot", 
                                      sidebarPanel(
                                            titlePanel(h3("Volcano Plot")),
                                            h5('X and Y axis variable selection for the Volcano Plot:'),
                                            radioButtons("DE_radio1", "Choose x-axis variable", c("baseMean"="baseMean","log2FoldChange"="log2FoldChange","lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"), selected="log2FoldChange"),
                                            radioButtons("DE_radio2", "Choose y-axis variable", c("baseMean"="baseMean","log2FoldChange"="log2FoldChange","lfcSE"="lfcSE","stat"="stat","pvalue"="pvalue","padj"="padj"), selected="padj"),
                                            colourInput("DE_color1", "Base point color", value = "blue"),
                                            colourInput("DE_color2", "Highlight point color", value = "red"),
                                            sliderInput(inputId="DE_slider", label="Select the magnitude of the p adjusted coloring:", min=-50, max=0, value=-5, step = NULL,round = FALSE,ticks = TRUE),
                                            #actionButton("goPlot", "Start Plotting", width='100%', style = "background-color: lightblue; color: black")
                                            submitButton("Start Plotting")
                                      ), #sidebarPanel
                                      mainPanel(
                                        tabsetPanel(type = "tabs", 
                                            tabPanel("Volcano Plot", plotOutput('DE_Volcano_plot')),
                                            tabPanel("DE Table2", DT::dataTableOutput(outputId = "DE_Table2"))
                                        ) #tabsetPanel
                                      ), #mainPanel
                                      position = c("left", "right"),
                                      fluid = FALSE
                                  )
                                ) #tabsetPanel
                            ) #mainPanel
                          ) #sidebarLayout                
                      ),
                      ##############  Tab4 : Networks  #############
                      #tabPanel("Networks", verbatimTextOutput("summary")),
                      ##############  Tab4 : GSEA  #############
                      tabPanel("GSEA", 
                          sidebarLayout(
                              sidebarPanel(width=3,
                                        fileInput("fgsea_file", "FGSEA File", accept = ".csv", placeholder = "fgsea_file.csv")
                                        #submitButton("Upload FGSEA File")
                              ), #sidebarPanel
                              mainPanel(width = 9,
                                  tabsetPanel(type = "tabs",
                                      tabPanel("Top Results", 
                                            sidebarPanel(width=4,
                                                sliderInput(inputId="GSEA_slider1", label="Top Results by adjusted p-value", min=0, max=100, value=10, step = NULL,round = FALSE,ticks = TRUE),
                                                submitButton("Update Top Pathways")
                                            ), #sidebarPanel
                                            mainPanel(width = 5,
                                                tabsetPanel(#type = "tabs",
                                                    type = "hidden",
                                                    tabPanel("GSEA_Plot1", plotOutput("GSEA_Plot1"))
                                                ) #tabsetPanel
                                            ) #mainPanel
                                       ),
                                       tabPanel("Table", 
                                           sidebarPanel(width=4,
                                                      sliderInput(inputId="GSEA_slider2", label="Filter table by adjusted p-value", min=0, max=1, value=0.1, step = NULL,round = FALSE,ticks = TRUE),
                                                      radioButtons("GSEA_radio2", "Choose all, positive and negative pathways", c("All Pathways"="All","Positive Pathways"="Positive","Negative Pathways"="Negative"), selected="All"),
                                                      submitButton("Update Table View"),
                                                      br(),
                                                      #fileInput("GSEA_file_download", "FGSEA File", accept = ".csv", placeholder = "GSEA_file_download.csv")
                                                      downloadButton("fgsea_NES_Table", "Download Table")
                                           ), #sidebarPanel
                                           mainPanel(width = 5,
                                                tabsetPanel(type = "hidden",
                                                    tabPanel("GSEA_Table2", DT::dataTableOutput("GSEA_Table2"))
                                                  ) #tabsetPanel
                                           ) #mainPanel
                                        ),
                                        tabPanel("Plots",
                                            sidebarPanel(width=4,
                                                  sliderInput(inputId="GSEA_slider3", label="Filter table by adjusted p-value", min=0, max=1, value=0.1, step = NULL,round = FALSE,ticks = TRUE),
                                                  submitButton("Update Padj Threshold")
                                            ), #sidebarPanel
                                            mainPanel(width = 5,
                                                tabsetPanel(type = "hidden",
                                                    tabPanel("GSEA_Plot3", plotOutput("GSEA_Plot3"))
                                                ) #tabsetPanel
                                            ) #mainPanel
                                        ) #tabPanel
                                  ) #tabsetPanel
                              ) #mainPanel
                            ) #sidebarLayout                
                          ) #tabPanel("GSEA")
          ) #tabsetPanel
    
      ) #mainPanel 
  ), #sidebarLayout

  hr(),
  tags$footer(h4("Nested Tab Input Control"))
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  ##############  Tab1 : Samples  #############
  # 12.5.1 Sample Information Exploration
  # The distinct values and distributions of sample information are important to 
  # understand before conducting analysis of corresponding sample data. 
  # This component allows the user to load and examine a sample information matrix.
  # This function takes no arguments. The `reactive({})` bit says "if any of 
  # my inputs (as in, input$...) are changed, run me again". 
  # This is useful when a user clicks a new button or loads a new file. In 
  # our case, look for the uploaded file's datapath argument and load it with 
  # read.csv. Return this data frame in the normal return() style.
  
  ########################## Tab 1: Samples  ############################ 
  load_data <- reactive({
    Samples_file <- input$Samples_file
    if (is.null(Samples_file)) {return(NULL)} 
    Samples_file_df <- as.data.frame(read_tsv(Samples_file$datapath, show_col_types = FALSE))
    return(Samples_file_df)
  })
  
  output$Samples_summary <- renderPrint({
    df <- load_data()
    if (is.null(df)) {return(NULL)} 
    summary(str(df))
    #create a summary table
    describe(df)
  }) 
  
  #tabPanel("Table", tableOutput("Samples_table")),
  output$Samples_table <- DT::renderDataTable({
    dataf <- load_data()
    if(is.null(dataf)){return ()}
    dataf
  }) 

  output$Samples_plot <- renderPlot({
    df <- load_data()
    if(is.null(df)){return ()}
    plot_data <- sample_plot(df, hist=TRUE) 
    plot_data
    
   }, height=300, width=500)
  
  output$Samples_plot2 <- renderPlot({
    return(NULL)
    df <- load_data()
    if(is.null(df)){return ()}
    plot_data <- sample_plot(df, hist=FALSE) 
    plot_data
    
  }, height=300, width=300)
  
  sample_plot <- function(v_counts, hist) {
    #Calculate the number of reads in each sample and plot as a barplot
    #Find the variance of all the genes
    v_count_vars <- apply(v_counts[,-1], 1, var)
    #Now filter out all the genes that have zero variances
    filtered_counts <- v_counts[v_count_vars > 0, ]
    
    sample_reads <- colSums(filtered_counts[-1])
    #Make a barplot showing the number of millions of reads for each sample
    df <- as.data.frame(sample_reads/10^6)
    # make sure the sample names are ordered correctly
    Sample <- factor(rownames(df), levels=rownames(df))
    
    # pivot long so we can plot gene vs the norm_counts
    plot_data <- pivot_longer(filtered_counts, cols=colnames(filtered_counts[-1]), names_to='samplenames', values_to='norm_counts')
   
    ###############################################
    
    if (hist == TRUE) {
      p_plot <- ggplot(df, aes(df[,1]))+
        geom_histogram(color="black", fill="lightblue", bins=80, binwidth=0.5) +
          theme_classic() +
          xlab("Sample Reads") + 
          labs(title = "Histogram of reads (in millions) per sample")
      p_plot  
      
    return(p_plot)
    
    }
    
    bplot <- 
      df %>% ggplot() + geom_bar(aes(x=Sample, y=df[,1], fill=Sample), stat="identity") +
      xlab("Sample") + ylab("# reads") +
      labs(title="Barplot of reads (in millions) for each sample")
    bplot
    return(bplot)
  }
    
  ########################## Tab 2: Counts  ############################
  # Pressing Plot button on web interface will kick off the Plot and Table 
  # rendering.
  observeEvent("Counts Submit", {

    load_Counts_data <- reactive({
      Counts_file <- input$Counts_file
      if (is.null(Counts_file)) {return(NULL)} 
      Counts_file_df <- as.data.frame(read_tsv(Counts_file$datapath, show_col_types = FALSE))
      updateSliderInput(session = session, "Counts_slider2", value = input$Counts_slider2, min = 0, max = ncol(Counts_file_df[-1]), step = 1)
      return(Counts_file_df)
    })
    
  filter_X_var_genes <- function(v_counts, slider1) {
      #Find the variance of all the genes
      v_count_vars <- apply(v_counts[,-1], 1, var)
      v_threshold <- quantile(v_count_vars, slider1)
      filter_X_vars <- v_counts[v_count_vars >= v_threshold, ]
      return(filter_X_vars)
  }
  
  filter_num_of_zeros <- function(v_counts, slider2) {
     thresh <- slider2
     v_counts_zeros <- v_counts
     v_counts_zeros$Zeros<-rowSums(v_counts_zeros==0)
     filter_zeros <- filter(v_counts, (rowSums(v_counts[-1]==0) >= thresh))
    return(filter_zeros)
  } 
  
  filter_counts_data <- function(v_counts, slider1, slider2) {
      filter_var <- filter_X_var_genes(v_counts, slider1)
      thresh <- slider2
      filter_zero_samples <- filter_var[rowSums(filter_var[-1] > 0 ) >= (thresh),]
      filtered <- filter_zero_samples
      return(filtered)
  }
  
  # Shiny Functionalities:
  # 1. Tab with text or a table summarizing the effect of the filtering, including:
  #    1. number of samples
  #    2. total number of genes
  #    3. number and % of genes passing current filter
  #    4. number and % of genes not passing current filter
  output$Counts_Summary <- renderPrint({
      if (is.null(input$Counts_file)) {
        return(NULL)
      }
      Counts <- load_Counts_data()
      #view(Counts)
      if (is.null(Counts)) {return(NULL)} 

      filtered <- filter_counts_data(Counts, input$Counts_slider1, input$Counts_slider2)
      genes_passing <- nrow(filtered)
      percent_genes_passing <- round((100 * genes_passing)/nrow(Counts), 2)
      genes_failing <- nrow(Counts) - nrow(filtered)
      percent_genes_failing <- round((100 * genes_failing)/nrow(Counts), 2)

      print("**************************************************")
      print(paste0("1. number of samples:                          ", ncol(Counts)))      
      print(paste0("2. total number of genes:                      ", nrow(Counts)))
      print(paste0("3. number of genes passing current filter:     ", genes_passing ))
      print(paste0("   % of genes passing current filter:          ", percent_genes_passing ))
      print(paste0("4. number of genes not passing current filter: ", genes_failing))
      print(paste0("   % of genes not passing current filter:      ", percent_genes_failing ))
      print("**************************************************")

      print("summary(str(Counts))")
      summary(str(Counts))
  })
    
  # 2. Tab with diagnostic scatter plots, where genes passing filters are marked in a 
  #    darker color, and genes filtered out are lighter:
  #    1. median count vs variance (consider log scale for plot)
  #    2. median count vs number of zeros
  #tabPanel("Diagnostic Plot", plotOutput("Counts_diagnostic_plot")),
  output$Counts_diagnostic_plot1 <- renderPlot({
    if (is.null(input$Counts_file)) {return(NULL)}
    Counts <- load_Counts_data()
    if (is.null(Counts)) {return(NULL)} 
    
    #rename first col to 'gene'
    names(Counts)[1] <- 'gene'
    filtered <- filter_X_var_genes(Counts, input$Counts_slider1)
    unfiltered <- filter(Counts, !(rownames(Counts) %in% rownames(filtered)))

    # Calculate Median values and Variance values
    Counts_df <- as.matrix(Counts)
    Counts_df <- type.convert(Counts_df, as.is = FALSE)
    sample_median <- rowMedians(Counts_df[,-1])
    sample_median <- rowMedians(as.matrix(Counts[,c(-1)]))
    sample_var <- apply(Counts_df[,-1], 1, var)
    median_var_tibble <- tibble(medians=sample_median, variance=sample_var)
    gene <- Counts$gene
    median_var_tibble <- cbind(median_var_tibble, gene)
    median_var_tibble <- median_var_tibble %>%
      mutate(status = case_when(
        (median_var_tibble$gene %in% filtered$gene) ~ "Passed Filter",
        (median_var_tibble$gene %in% unfiltered$gene) ~ "Failed Filter"
      ))
    
    median_var_tibble$rank <- rank(median_var_tibble$medians)
    view(median_var_tibble)
    
    var_plot <- 
      ggplot(median_var_tibble, aes(x=rank, y=-log10(sample_var), color=status, stat="identity")) + 
      geom_point() + 
      scale_color_manual(values = c('Passed Filter' = "black", 'Failed Filter' = "red")) + 
      labs(title = "Median Vs Variance") +
      xlab("Median") + ylab("-log10(Variance)") 
    var_plot
    
    return(var_plot)
  })
  
  output$Counts_diagnostic_plot2 <- renderPlot({
    if (is.null(input$Counts_file)) {return(NULL)}
    Counts <- load_Counts_data()
    if (is.null(Counts)) {return(NULL)} 
    #rename first col to 'gene'
    names(Counts)[1] <- 'gene'
    filtered <- filter_num_of_zeros(Counts, input$Counts_slider2)
    unfiltered <- filter(Counts, !(rownames(Counts) %in% rownames(filtered)))

    Counts_df <- as.matrix(Counts)
    Counts_df <- type.convert(Counts_df, as.is = FALSE)
    sample_median <- rowMedians(as.matrix(Counts[,c(-1)]))
    sample_var <- apply(Counts_df[,-1], 1, var)
    median_zero_tibble <- tibble(medians=sample_median, variance=sample_var)
    gene <- Counts$gene
    median_zero_tibble <- cbind(median_zero_tibble, gene)
    median_zero_tibble$sample_zero<-rowSums(Counts==0)
    median_zero_tibble <- median_zero_tibble %>%
      mutate(status = case_when(
        (median_zero_tibble$gene %in% filtered$gene) ~ "Passed Filter",
        (median_zero_tibble$gene %in% unfiltered$gene) ~ "Failed Filter"
      ))
    median_zero_tibble$rank <- rank(median_zero_tibble$medians)
    view(median_zero_tibble)
    
    zero_plot <- 
      ggplot(median_zero_tibble, aes(x=rank, y=sample_zero, color=status, stat="identity")) + 
      geom_point() + 
      scale_color_manual(values = c('Passed Filter' = "black", 'Failed Filter' = "red")) + 
      labs(title = "Median Vs Num of Zeros") +
      xlab("Median") + ylab("Num of zeros")
    zero_plot

    return(zero_plot)
  })
  

  # 3. Tab with a clustered heatmap of counts remaining after filtering
  #    1. consider enabling log-transforming counts for visualization
  #    2. be sure to include a color bar in the legend
  output$Counts_heatmap <- renderPlot({
    if (is.null(input$Counts_file)) {return(NULL)}
    Counts <- load_Counts_data()
    if (is.null(Counts)) {return(NULL)} 
    filtered <- filter_counts_data(Counts, input$Counts_slider1, input$Counts_slider2)
    a <- as.matrix(filtered[-1])
    row.names(a) <- a()$Name
    a[is.na(a)] <- 0
    a
    pheatmap(a, kmeans_k = 4)
  })
    
  # 4. Tab with a scatter plot of principal component analysis projections. You may either:
  #    1. allow the user to select which principal components to plot in a scatter plot (e.g. PC1 vs PC2)
  #    2. allow the user to plot the top N principal components as a beeswarm plot
  #    3. be sure to include the % variance explained by each component in the plot labels
  output$Counts_PCA_plot <- renderPlot({
    if (is.null(input$Counts_file)) {return(NULL)}
    Counts <- load_Counts_data()
    if (is.null(Counts)) {return(NULL)} 
    filtered <- filter_counts_data(Counts, input$Counts_slider1, input$Counts_slider2)
    pca_results <- prcomp(scale(t(filtered[-1])), center=FALSE, scale=FALSE)
    pc_variance_explained <- pca_results$sdev^2 / sum((pca_results$sdev)^2)
    pc_variance_explained
    
    #construct a tibble with PCs, variance explained, and cumulative variance explained}
    variance_tibble <- make_variance_tibble(pc_variance_explained, pca_results)
    ################################
    updateSelectInput(session, "xcol", choices = colnames(pca_results$x), selected = "PC1") #input$xcol) #pca_results$x[,1])
    updateSelectInput(session, "ycol", choices = colnames(pca_results$x), selected = "PC2") #input$xcol) #pca_results$x[,2])
    ################################
    pplot <- PCA_plot(pca_results)
    pplot 
    return(pplot)
   }, height=300, width=300)
  
  # Pressing Plot button on web interface will kick off the Plot and Table 
  # rendering.
  observeEvent(
    eventExpr = input[["Submit Variable"]],
    handlerExpr = {
      updateSelectInput(session, "xcol", choices = colnames(pca_results$x), selected = input$xcol) #pca_results$x[,1])
      updateSelectInput(session, "ycol", choices = colnames(pca_results$x), selected = input$ycol) #pca_results$x[,2])
    })  # replace this NULL
      
  
  PCA_plot <- function(pca_results) {
    pca_df <- as.data.frame(pca_results$x)
    #############################
    # Following code converts input$xcol and input$ycol which are in the 'Character' form
    # to actual col names and rearrangs as 1st and 2nd col of a data frame so we can use
    # first 2 cols for plotting.
    ordered_pca <- pca_df %>%  dplyr::arrange(dplyr::desc(!!rlang::sym(input$xcol)),
                                                     dplyr::desc(!!rlang::sym(input$ycol)))
    otherCols <- setdiff(colnames(ordered_pca), unique(c(input$xcol,input$ycol)))
    #Re-arrange columns so that input$xcol and input$ycol are first and second cols
    ordered_pca <- ordered_pca %>%  dplyr::select(!!rlang::sym(input$xcol),
                                                  !!rlang::sym(input$ycol),
                                                  !!otherCols)
    #view(ordered_pca)
    # Calculate % Variance for each selected PCA
    var <- pca_results$sdev^2 / sum(pca_results$sdev^2)
    col1 <- which(colnames(pca_df)==input$xcol)
    col2 <- which(colnames(pca_df)==input$ycol)
    var1_percentage <- round(var[col1] * 100)
    var2_percentage <- round(var[col2] * 100)
    #############################
    pplot <- 
      ggplot(ordered_pca, aes(x=ordered_pca[,1], y=ordered_pca[,2])) + 
      geom_point(size=2.0) + 
      theme_classic() +
      xlab(paste0(input$xcol,": ",var1_percentage,"% Variance")) + 
      ylab(paste0(input$ycol,": ",var2_percentage,"% Variance")) 

    pplot
    return(pplot)
  }
  
  make_variance_tibble <- function(pca_ve, pca_results) {
    pca_cum <- cumsum(pca_results$sdev^2 / sum(pca_results$sdev^2))
    variance_explained = c(pca_ve)
    cumulative = c(pca_cum)
    principal_components = c(colnames(pca_results$rotation))
    var_tibble <- tibble(variance_explained, principal_components, cumulative)
    return(var_tibble)
  }

  })
  ########################## Tab 3: DE  ############################ 
  # 12.5.3 Differential Expression
  # Differential expression identifies which genes, if any, are implicated in a specific 
  # biological comparison. This component allows the user to load and explore a differential 
  # expression dataset.
  # Inputs:
  #   Results of a differential expression analysis in CSV format. If results are already 
  #   made available, you may use those
  #   Otherwise perform a differential expression analysis using DESeq2, limma, or edgeR from 
  #   the provided counts file

  # Shiny Functionalities:
  #   Tab with sortable table displaying differential expression results
  #   Optional: enable gene name search functionality to filter rows of table
  #   Tab with content similar to that described in Assignment 7
 
  load_DE_data <- reactive({
    DE_file <- input$DE_file
    if (is.null(DE_file)) {
      return(NULL)} 
    DE_file_df <- as.data.frame(read_tsv(DE_file$datapath, show_col_types = FALSE))
    return(DE_file_df)
  })
  
  output$DE_Table1 <- DT::renderDataTable({
    de_data <- load_DE_data()
    if(is.null(de_data)){return (NULL)}
    de_data
    return(de_data)
  }) 

  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    if (is.null(x_name) || is.null(y_name)) {return(NULL)}
    if (x_name == y_name) {return(NULL)}
    if (is.null(input$DE_file)) {return(NULL)}
 
    ordered_df <- dataf %>%  dplyr::arrange(dplyr::desc(!!rlang::sym(x_name)),
                                            dplyr::desc(!!rlang::sym(y_name)))
    # Since we do not know which dataf cols are x_name and y_name, we will re-arrange
    # the x_name col as the 1st col and y_name as the 2nd col in the new data frame
    # so that it will be easy to do volcano plot.
    #view(ordered_df)
    otherCols <-setdiff(colnames(ordered_df), unique(c(x_name,y_name)))
    #Re-arrange columns according to the selections
    ordered_df <- ordered_df %>%  dplyr::select(!!rlang::sym(x_name),
                                                !!rlang::sym(y_name),
                                                !!otherCols)
    
    # Add a new col slider_cond whose values will be TRUE, FALSE or NA  and will
    # used to color code the volcano
    # FALSE if y-axis value >= slider_factor, 
    # TRUE if y-axis value < slider_factor, else NA 
    
    slider_factor <- (1 * (10^slider))
    ordered_df <- ordered_df %>%
      mutate(slider_cond = case_when(ordered_df[, 2] < slider_factor ~ "TRUE", 
                                     ordered_df[, 2] >= slider_factor ~ "FALSE", TRUE ~ 'NA'))
    
    x_axis <- colnames(ordered_df[,1])  #x-axis label for the plot
    y_axis <- colnames(ordered_df[,2])  #y-axis label for the plot
    
    # Remove rows with any NA values in plotting cols:
    #ordered_df <- filter(ordered_df, !is.na(ordered_df[,1]) | !is.na(ordered_df[,2]))
    ordered_df <- filter(ordered_df, !is.na(ordered_df[,1]))
    ordered_df <- filter(ordered_df, !is.na(ordered_df[,2]))
    
    # Now to actual plotting
    vol_plot <-  
      ggplot(ordered_df, aes(x=ordered_df[,1], y=-log10(ordered_df[,2]), color=slider_cond)) +
      geom_point() + 
      scale_color_manual(name=paste0(y_name,"< 1 * 10^",slider), values = c('FALSE' = color1, 'TRUE' = color2, 'NA'='black')) + 
      xlab(x_name) + ylab(paste0('-log10(',y_name,')')) +
      labs(title=paste0('Volcano Plot: ', x_name, '  VS  -log10(', y_name, ')' )) +
      theme(legend.position="bottom")
    
    vol_plot
    return(vol_plot)
    ####################################################
  }
  
  #' Draw and filter table
  #'
  #' @param dataf Data frame loaded by load_data()
  #' @param slider Negative number, typically from the slider input.
  #'
  #' @return Data frame filtered to p-adjusted values that are less than 
  #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits 
  #' displayed.
  #' @details Same as above, this function is a standard R function. Tests will 
  #' evaluate it normally. Not only does this function filter the data frame to 
  #' rows that are above the slider magnitude, it should also change the format 
  #' of the p-value columns to display more digits. This is so that it looks 
  #' better when displayed on the web page. I would suggest the function 
  #' `formatC()`
  #'
  #' @examples draw_table(deseq_df, -210)
  #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
  #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
  #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
  draw_table <- function(dataf, slider) {
    data_f <- data.frame(dataf)
    filter_factor <- (1 * (10^slider))
    df <- filter(data_f, (padj < filter_factor))
    df$pvalue <- formatC(df$pvalue, digits=6)
    df$padj <- formatC(df$padj, digits=6)
    
    return(df)
  }
  
  # Pressing Plot button on web interface will kick off the Plot and Table 
  # rendering.
  observeEvent("Start Plotting", {
      output$DE_Volcano_plot <- renderPlot({ 
        if (is.null(input$DE_file)) {return(NULL)}
        dataf <- load_DE_data()
        if (is.null(dataf)) {return(NULL)}
        
        #updateSliderInput(session = session, "DE_slider", value = input$DE_slider, min = 0, max = ncol(Counts_file_df[-1]), step = 1)
        
        #########################################################################
        # Following create_fgsea function is used to create a fgsea file from the
        # de expression. Run it only if we need to create a new fgsea file.
        
        #create_fgsea(dataf)
        #########################################################################
        plot_data <- volcano_plot(dataf, input$DE_radio1, input$DE_radio2, input$DE_slider, input$DE_color1, input$DE_color2)
        plot_data
      },  height = 500, width = 500)  # replace this NULL
      
      output$DE_Table2 <- DT::renderDataTable({
        if (is.null(input$DE_file)) {
          #print("output$DE_Table2 - Null DE File")
          return(NULL)}
        dataf <- load_DE_data()
        if (is.null(dataf)) {
          return(NULL)}
        table_data <- draw_table(dataf, input$DE_slider)
        table_data
      })
  }) #observeEvent("Start Plotting"

  output$DE_Table2 <- DT::renderDataTable({
    if (is.null(input$DE_file)) {
      return(NULL)}
    dataf <- load_DE_data()
    if (is.null(dataf)) {
      return(NULL)}
    table_data <- draw_table(dataf, input$DE_slider)
    table_data
  })
  ########################## Tab4 : GSEA  ############################ 
  # 12.6.1 Gene Set Enrichment Analysis
  # Use your differential gene expression results to compute gene set enrichment analysis 
  # with fgsea. You will need to identify an appropriate gene set database that matches the 
  # organism studied in your dataset.
  
  # Input:
  # a table of fgsea results from the differential expression data.
  # - choose an appropriate ranking metric (log fold change, -log(pvalue), etc) from your 
  #   differential expression results 
  # - run fgsea with appropriate parameters against gene set database of your choice
  # - save the results to a CSV/TSV
  # file upload button
  
  # Shiny Functionalities:
    
  # Tab 1
  # Sidebar:
  #       Slider to adjust number of top pathways to plot by adjusted p-value
  # Main Panel
  #       Barplot of fgsea NES for top pathways selected by slider
  #       Optional: Click on barplot for pathway to display table entry
  # Tab 2
  # Sidebar
  #       Slider to filter table by adjusted p-value (Reactive)
  #       Radio buttons to select all, positive or negative NES pathways
  #       Download button to export current filtered and displayed table results
  # Main panel
  #       Sortable data table displaying the results
  # Tab 3
  # Sidebar
  #       Slider to filter table by adjusted p-value (Reactive)
  # Main panel
  #       Scatter plot of NES on x-axis and -log10 adjusted p-value on y-axis, with gene sets 
  #       below threshold in grey color
  ##############  Tab4 : GSEA  #############

  ensembl_id_to_gene_symbol <- function(values) {
    mart <- useMart('ensembl', dataset='hsapiens_gene_ensembl')
    hgnc_ids <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"), filter="ensembl_gene_id", value=values, mart=mart, uniqueRows=FALSE)
    #hgnc_ids <- getBM(attributes=c("ensembl_gene_id_version","hgnc_symbol"), filter="ensembl_gene_id_version", value=values, mart=mart, uniqueRows=FALSE)
    return(hgnc_ids)
  }
  
  
  #Run fgsea using a ranked list of descending log2FC against the C2 canonical
  #pathways gene set
  #Set minsize to 15 and maxsize to 500, leave the other parameters as defaults
  #fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
  #fgsea_results
  #' Function to run fgsea on DESeq2 results
  #'
  #' @param labeled_results (tibble): the labeled results from DESeq2
  #' @param gmt (str): the path to the GMT file
  #' @param min_size: the threshold for minimum size of the gene set
  #' @param max_size: the threshold for maximum size of the gene set
  #'
  #' @return tibble containing the results from running fgsea using descending
  #' log2foldchange as a ranking metric
  #' @export
  #'
  #' @examples fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
  run_gsea <- function(labeled_results, gmt, min_size, max_size) {
  #run_gsea <- function(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500) {
    # Strip the digit extensions in the gene names in labled_results and pull the 
    # updated/modified gene names to be used in getLDS to retrieve information from
    # human/mouse linked database
    labeled_results_modified <- separate(labeled_results, genes, sep='\\.', into='genes', remove=TRUE)
    # Pull the modified gene names
    #labeled_results_modified <- labeled_results
    genes <- labeled_results_modified$genes
    # Connect to BioMart and retrieve information from human/mouse linked database
    human <- useMart('ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
    mouse <- useMart('ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl', host = "https://dec2021.archive.ensembl.org/")
    
    hgnc_symbols <- getLDS(attributes=c('ensembl_gene_id'), filters='ensembl_gene_id', values=genes, mart=mouse, 
                           attributesL=c('hgnc_symbol'), martL=human, uniqueRows=T)
    
    hgnc_results <- left_join(labeled_results_modified, hgnc_symbols, by=c('genes' = 'Gene.stable.ID'))
    # Create ranked list of log2FoldChange values
    log2FC <- filter(hgnc_results, !is.na(HGNC.symbol) & !is.na(log2FoldChange))
    log2FC <- log2FC %>% arrange(desc(log2FoldChange))
    log2FC_rank <- log2FC[c("HGNC.symbol", "log2FoldChange")]
    log2FoldChange_rank <- deframe(log2FC_rank)
    
    # Get C2 Canonical Pathways gene set collection
    gmt_pathways <- gmtPathways(gmt)
    # Run gsea on ranked log2FoldChange values
    fgsea_results <- fgsea(gmt_pathways, log2FoldChange_rank, minSize=min_size, maxSize=max_size)
    fgsea_results <- as_tibble(fgsea_results)
    #rename first col to 'gene'
    names(fgsea_results)[1] <- 'genes'
    return(fgsea_results)
  }
  
  create_fgsea <- function(DE_file_df) {
    deseq2_res <- DE_file_df
    padj_threshold <- input$DE_slider
    if (is.null(padj_threshold)) {
      padj_threshold <- 0.10
    }
    #rename first col to 'gene'
    names(deseq2_res)[1] <- 'genes'
    
    labeled <- deseq2_res %>%
      #as_tibble(rownames='genes') %>%
      mutate(volc_plot_status = case_when(log2FoldChange > 0 & padj < padj_threshold ~ 'UP', 
                                          log2FoldChange < 0 & padj < padj_threshold ~ 'DOWN', 
                                          TRUE ~ 'NS'))
    #Display the summary of the tibble
    labeled_results <- labeled[order(labeled$padj, decreasing = FALSE), ]
    labeled_results %>% relocate(genes, volc_plot_status, log2FoldChange, padj)
    fgsea_results <- run_gsea(labeled_results, 'c2.cp.v7.5.1.symbols.gmt', 15, 500)
    col_type <- sapply(fgsea_results, class)
    fgsea_res <- fgsea_results[, sapply(fgsea_results, class) != "list"]
    fgsea_res_mat <- as_tibble(fgsea_res)
  }
  
  #' Function to plot top ten positive NES and top ten negative NES pathways
  #' in a barchart
  #'
  #' @param fgsea_results (tibble): the fgsea results in tibble format returned by
  #'   the previous function
  #' @param num_paths (int): the number of pathways for each direction (top or
  #'   down) to include in the plot. Set this at 10.
  #'
  #' @return ggplot with a barchart showing the top twenty pathways ranked by positive
  #' and negative NES
  #' @export
  #'
  #' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
  top_pathways <- function(fgsea_results, num_paths){
    
    #rename first col to 'gene'
    names(fgsea_results)[1] <- 'pathway'
    
    #print("top_pathways")
    #fgsea_results is already sorted in the aescending order of NES
    # make a tibble with only pathways with top 10 positive and top 10 negative NES values
    # gather top 10 negative values
    top_neg <- fgsea_results[1:num_paths,]
    
    # gather top 10 positive values
    start_pos <- nrow(fgsea_results) - (num_paths-1)
    stop_pos  <- nrow(fgsea_results)
    top_pos <- fgsea_results[start_pos:stop_pos,]
    
    # filter out the rows wtih top 10 pos and neg valued pathways
    pos_pathway <- top_pos$pathway
    neg_pathway <- top_neg$pathway
    fgsea_subset <- fgsea_results %>% filter(pathway %in% c(pos_pathway, neg_pathway)) 
    
    # In fgsea_subset create a new col of pathway names (x_axis_names) without '-' 
    # in the name and reorder them for x-axis lables 
    #fgsea_subset <- fgsea_subset %>% mutate(x_axis_names = str_replace_all(pathway, '_', ' ')) %>% 
    fgsea_subset <- fgsea_subset %>% mutate(x_axis_names = str_replace_all(pathway, '_', ' ')) 

    # Create a bar plot of top 10 pos and neg NES pathways
    nes_plot <- 
      ggplot(fgsea_subset) +
      geom_bar(aes(x=x_axis_names, y=NES, fill = NES > 0), stat='identity', show.legend = FALSE) +
      scale_fill_manual(values = c('TRUE' = 'blue', 'FALSE' = 'red')) + 
      #theme_minimal(base_size = 8) +
      theme_bw(base_size = 8) +
      ggtitle('fgsea results for Hallmark MSigDB gene sets') +
      ylab('Normalized Enrichment Score (NES)') +
      xlab('') 
    
    # Since some of the the pathway names are very long let us wrap the names so that we can
    # display the bar plot correctly
    wrap_pathway_name <- function(x) str_wrap(x, width = 50)
    nes_plot <- nes_plot +
      scale_x_discrete(labels = wrap_pathway_name) 
    
    # Filp the plot from vertical to horizontal layout
    nes_plot <- nes_plot + coord_flip()
    nes_plot
    return(nes_plot)
    
  }
  
  load_fgsea_data <- reactive({
    fgsea_file <- input$fgsea_file
    if (is.null(input$fgsea_file)) {
      return(NULL)} 
    fgsea_file_df <- read_tsv(fgsea_file$datapath, show_col_types = FALSE)
    
    if ("log2FoldChange" %in% colnames(fgsea_file_df)) {
      print("This apprears to be a DE file.")
      print("Creating fgsea.csv file now. This will take a while. Please wait.")
      #########################################################################
      # Following create_fgsea function is used to create a fgsea file from the
      # de expression. Run it only if we need to create a new fgsea file.
      fgsea_res_mat <- create_fgsea(fgsea_file_df)
      #########################################################################
      write.csv(fgsea_res_mat, 'fgsea.csv', row.names=FALSE)
      print("Finished creating 'fgsea.csv' file.")
      return(fgsea_res_mat)
    }
    print("This apprears to be FGSEA file. Processing now...")
    fgsea_file_df <- read.csv(fgsea_file$datapath, header = TRUE)
    return(fgsea_file_df)
  })

  output$GSEA_Plot1 <- renderPlot({
      req(input$fgsea_file)
      fgsea_results <- load_fgsea_data()
      view(fgsea_results)
      if(is.null(fgsea_results)){return (NULL)}
      
      #Display the results of the fgsea in a tibble sorted by by NES (ascending)
      filtered_results <- fgsea_results[order(fgsea_results$NES, decreasing = FALSE), ]
      filtered_results
      
      #Plot the top ten pathways with both positive and negative NES (20 total)
      #Color the pathways by the sign of their NES (positive or negative)
      fgsea_plot <- top_pathways(filtered_results, 10)
      fgsea_plot
    })
    
 # }) #observeEvent("Upload FGSEA File"

  observeEvent("Update Top Pathways", {
    val <- input$GSEA_slider1
    updateSliderInput(session, "GSEA_slider1", value = val)
    
    output$GSEA_Plot1 <- renderPlot({
      fgsea_results <- load_fgsea_data()
      if(is.null(fgsea_results)){return (NULL)}
      
      # filter rows with padj < padj_threshold
      pathways_threshold <- input$GSEA_slider1
      #Display the results of the fgsea in a tibble sorted by by NES (ascending)
      #filtered_results <- fgsea_res[order(fgsea_res$NES, decreasing = FALSE), ]
      filtered_results <- fgsea_results[order(fgsea_results$padj, decreasing = FALSE), ]
      filtered_results
      
      #Plot the top ten pathways with both positive and negative NES (20 total)
      #Color the pathways by the sign of their NES (positive or negative)
      fgsea_plot <- top_pathways(filtered_results, pathways_threshold)
      fgsea_plot
    },  height = 400, width = 400)
    
  })#observeEvent("Update Top Pathways")
  

  observeEvent("Update Table", {
    slider2_val <- input$GSEA_slider2
    updateSliderInput(session, "GSEA_slider2", value = slider2_val)
    
    output$GSEA_Table2 <- DT::renderDataTable({
      display_table <- sliderValues()
      return(display_table)
    }) 
  })
  
  output$GSEA_Table2 <- DT::renderDataTable({
    display_table <- sliderValues()
    return(display_table)
  }) 
  # Reactive function
  sliderValues <- reactive({
    fgsea_table <- load_fgsea_data()
    padj_threshold <- input$GSEA_slider2
    fgsea_res <- filter(fgsea_table, padj <= padj_threshold)
    
    pathways_option <- input$GSEA_radio2
    NES_filter <- fgsea_res
    fgsea_res
    if (pathways_option == "Positive") {
      NES_filter <- filter(NES_filter, NES >= 0)
      #print("Positive Pathways Only")
    }
    if (pathways_option == "Negative") {
      NES_filter <- filter(NES_filter, NES < 0)
      #print("Negative Pathways Only")
    }
    NES_filter
    NES_filter_mat <- as_tibble(NES_filter)
    return(NES_filter_mat)
  })
  
  output$fgsea_NES_Table <- downloadHandler(
    filename = function() {
      paste("filtered_fgsea-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(display_table <- sliderValues(), file)
    }
  )
  
  observeEvent("Update Padj Threshold", {
    slider_val <- input$GSEA_slider3
    updateSliderInput(session, "GSEA_slider3", value = slider_val)
   
    output$GSEA_Plot3 <- renderPlot({
      NES_plot <- NES_Padj_Plot()
      return(NES_plot)
    }, height= 400, width=500) 
    
  })#observeEvent("Update Padj Threshold")
  
  NES_Padj_Plot <- reactive({
    fgsea_table <- load_fgsea_data()
    padj_threshold <- input$GSEA_slider3
    fgsea_res_below <- filter(fgsea_table, padj <= padj_threshold)
    fgsea_res_above <- filter(fgsea_table, padj > padj_threshold)
    
    fgsea_res_mod <- fgsea_table %>%
      mutate(slider_status = case_when(padj <= padj_threshold ~ "Below Threshold", 
                                     padj > padj_threshold ~ "Above Threshold"))
    # Now to actual plotting
    nes_padj_plot <-  
      ggplot(fgsea_res_mod, aes(x=NES, y=-log10(padj), color=slider_status)) +
      geom_point() +
      theme_bw() +
      scale_color_manual(name='Padj Threshold', values = c('Below Threshold' = 'gray', 'Above Threshold' = 'red')) + 
      xlab("NES") + ylab(paste0('-log10(padj)')) +
      labs(title= 'Scatter plot of NES vs -log10 adjusted p-value',
           subtitle= paste0('Genes Below Threshold:', nrow(fgsea_res_below),
                            ' Genes Above Threshold:', nrow(fgsea_res_above)))

    nes_padj_plot
    return(nes_padj_plot)
    
  })
  
  NES_Padj_Plot2 <- function(dataf, x_name, y_name, slider, color1, color2) {
    if (is.null(x_name) || is.null(y_name)) {return(NULL)}
    if (x_name == y_name) {return(NULL)}
    if (is.null(input$DE_file)) {return(NULL)}
    
    ordered_df <- dataf %>%  dplyr::arrange(dplyr::desc(!!rlang::sym(x_name)),
                                            dplyr::desc(!!rlang::sym(y_name)))
    # Since we do not know which dataf cols are x_name and y_name, we will re-arrange
    # the x_name col as the 1st col and y_name as the 2nd col in the new data frame
    # so that it will be easy to do volcano plot.
    otherCols <-setdiff(colnames(ordered_df), unique(c(x_name,y_name)))
    #Re-arrange columns according to the selections
    ordered_df <- ordered_df %>%  dplyr::select(!!rlang::sym(x_name),
                                                !!rlang::sym(y_name),
                                                !!otherCols)
    
    # Add a new col slider_cond whose values will be TRUE, FALSE or NA  and will
    # used to color code the volcano
    # FALSE if y-axis value >= slider_factor, 
    # TRUE if y-axis value < slider_factor, else NA 
    
    slider_factor <- (1 * (10^slider))
    ordered_df <- ordered_df %>%
      mutate(slider_cond = case_when(ordered_df[, 2] < slider_factor ~ "TRUE", 
                                     ordered_df[, 2] >= slider_factor ~ "FALSE", TRUE ~ 'NA'))
    
    x_axis <- colnames(ordered_df[,1])  #x-axis label for the plot
    y_axis <- colnames(ordered_df[,2])  #y-axis label for the plot
    
    # Remove rows with any NA values in plotting cols:
    ordered_df <- filter(ordered_df, !is.na(ordered_df[,1]))
    ordered_df <- filter(ordered_df, !is.na(ordered_df[,2]))
    
    # Now to actual plotting
    vol_plot <-  
      ggplot(ordered_df, aes(x=ordered_df[,1], y=-log10(ordered_df[,2]), color=slider_cond)) +
      geom_point() + 
      scale_color_manual(name=paste0(y_name,"< 1 * 10^",slider), values = c('FALSE' = color1, 'TRUE' = color2, 'NA'='black')) + 
      xlab(x_name) + ylab(paste0('-log10(',y_name,')')) +
      labs(title=paste0('Volcano Plot: ', x_name, '  VS  -log10(', y_name, ')' )) +
      theme(legend.position="bottom")
    
    vol_plot
    return(vol_plot)
    ####################################################
  }
  
  ######################################################################

}

# Run the application
shinyApp(ui = ui, server = server)
