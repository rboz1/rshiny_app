# import libraries
library(biomaRt)
library(colourpicker)
library(DESeq2)
library(dplyr)
library(DT)
library(gplots)
library(ggplot2)
library(monocle3)
library(RColorBrewer)
library(Seurat)
library(SeuratObject)
library(shiny)
library(shinycssloaders)
library(tidyverse)


axis_radio <- c("baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")

# ui
ui <- fluidPage(
  tabsetPanel(
    ##### summary tab #####
    tabPanel("Summary",
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId = "sample_csv", label = NULL, placeholder = "sample_info.csv", accept = c(".csv"))
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel(
                     "Summary Table", 
                     tableOutput(outputId = "summary_table") %>% 
                       withSpinner(type = 4, color = "lightgrey")
                   ),
                   tabPanel(
                     "Data Table", 
                     dataTableOutput(outputId = "summary_data_table") %>% 
                       withSpinner(type = 4, color = "lightgrey")
                   ),
                   tabPanel(
                     "Plots", 
                     plotOutput(outputId = "summary_plots", height = "90vh") %>% 
                       withSpinner(type = 4, color = "lightgrey")
                   )
                 )
               )
             )
    ),
    ##### counts matrix tab #####
    tabPanel("Counts",
      sidebarLayout(
        sidebarPanel(
          fileInput("counts_csv", label = NULL, placeholder = "counts_csv.csv", accept = ".csv"),
          sliderInput("var_slider", "Percentile of Variance", min = 0, max = 100, value = 50),
          sliderInput("zero_slider", "Percent of Non-Zero Samples", min = 0, max = 100, value = 50),
          selectInput("pc_x", "Select PC for X-axis", choices = NULL),
          selectInput("pc_y", "Select PC for Y-axis", choices = NULL)
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Summary", tableOutput("counts_summary_table")%>% 
                       withSpinner(type = 4, color = "lightgrey")),
            tabPanel("Diagnostic Plots", plotOutput("counts_scatter_plots", height = "85vh")%>% 
                       withSpinner(type = 4, color = "lightgrey")),
            tabPanel("Heatmap", plotOutput("counts_heatmap")%>% 
                       withSpinner(type = 4, color = "lightgrey"), height = "85vh"),
            tabPanel("PCA", plotOutput("counts_pca_plot", height = "85vh")%>% 
                       withSpinner(type = 4, color = "lightgrey"))
          )
        )
      )
    ),
    ##### differential expression tab #####
    tabPanel("Differential Expression",
             sidebarLayout(
               sidebarPanel(
                 fileInput("upload_file", label = NULL, placeholder = "bulk_counts.csv", accept = c(".csv")),
                 radioButtons("x_axis", "Choose the column for the x-axis", axis_radio, selected = "log2FoldChange"),
                 radioButtons("y_axis", "Choose the column for the y-axis", axis_radio, selected = "padj"),
                 colourInput("base_point_color", "Base point color", "black"),
                 colourInput("highlight_point_color", "Highlight point color", "#A857FF"),
                 sliderInput("slider", "Select the magnitude of the p adjusted coloring:", value = -4, min = -25, max = 0)
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("DE Results", 
                            DT::DTOutput("de_table") %>% withSpinner(type = 4, color = "lightgrey")
                   ),
                   tabPanel("Volcano Plot", 
                            plotOutput("volcano_plot", height = "85vh") %>% withSpinner(type = 4, color = "lightgrey")
                   )
                 )
               )
             )
    ),
    ##### trajectory analysis tab #####
    tabPanel("Trajectory Analysis",
             sidebarLayout(
               sidebarPanel(
                 fileInput("upload_rds", label = NULL, placeholder = "seurat_obj.rds", accept = c(".rds")),
                 selectInput("pathology_input", 
                             "Select Condition(s)", 
                             choices = NULL, 
                             selected = NULL, 
                             multiple = TRUE),
                 selectInput("cluster_input", 
                             "Select ECcluster(s)", 
                             choices = NULL, 
                             selected = NULL, 
                             multiple = TRUE),
                 selectInput("gene_input", 
                             "Select Gene(s)", 
                             choices = NULL, 
                             selected = NULL, 
                             multiple = TRUE),
                 selectInput("color_palette", 
                             "Select Color Palette:", 
                             choices = list("Set1", "Set2", "Dark2", "Spectral"),
                             selected = "Set1"),
                 actionButton("apply_filters", "Apply")
               ),
               
               mainPanel(
                 tabsetPanel(
                   tabPanel("Trajectory Plot", 
                            plotOutput("pseudotime_plot", height = "85vh") %>% 
                              withSpinner(type = 4, color = "lightgrey")),
                   tabPanel("Filtered Data Table", 
                            DTOutput("traj_filtered_table") %>% 
                              withSpinner(type = 4, color = "lightgrey"))
                 )
               )
             )
    )
  )
)


# server
server <- function(input, output, session) {
  # increase max upload size
  options(shiny.maxRequestSize = 3 * 1024^3)
  ##### summary tab ######
  # read in .CSV
  load_sample_data <- reactive({
    req(input$sample_csv)
    
    ext <- tools::file_ext(input$sample_csv$name)
    
    # only accept .CSV files
    data <- switch(ext,
                   csv = read.csv(input$sample_csv$datapath),
                   validate("Invalid file: Please upload a .csv file"))
    return(data)
  })
  
  # make a summary table with general descriptive stats about sample info
  generate_summary_table <- function(sample_info) {
    
    age_stats <- sample_info %>%
      dplyr::select("age") %>%
      pull() %>%
      {sprintf("%.2f ± %.2f", mean(., na.rm = TRUE), sd(., na.rm = TRUE))}
    
    sex_count <- sample_info %>% count(Gender) %>% summarize(output = paste(Gender, "=", n, collapse = ", ")) %>%
      pull(output)
    
    who_grade <- sample_info %>% count(`WHO.grade`) %>% summarize(output = paste(`WHO.grade`, "=", n, collapse = ", ")) %>%
      pull(output)
    
    epilepsy <- sum(grepl("^Y", sample_info$`Lesional.epilepsy..Y.N.`), na.rm = TRUE)
    
    glioblastoma <- sum(grepl("glioblastoma", sample_info$`Integrated.WHO.2016.diagnosis`, ignore.case = TRUE), na.rm = TRUE)
    
    sample_tibble <- tibble(
      age = age_stats,
      sex = sex_count,
      "WHO grade" = who_grade,
      epilepsy = epilepsy,
      glioblastoma = glioblastoma
    )
    
    transposed_tibble <- as_tibble(t(sample_tibble), .name_repair = "unique") %>%
      rename("mean (s.d) or distinct value" = `...1`) %>%
      mutate(column = colnames(sample_tibble),
             type = sapply(`mean (s.d) or distinct value`, typeof)
      ) %>%
      relocate("mean (s.d) or distinct value", .after = "type")
    
    return(transposed_tibble)
  }
  
  # summary table
  output$summary_table <- renderTable({
    generate_summary_table(load_sample_data())
  })
  
  # data table with sortable columns
  output$summary_data_table <- DT::renderDT({
    data <- load_sample_data()
    
    # enable horizontal scrolling
    datatable(data, options = list(
      scrollX = TRUE,  # Enable horizontal scroll
      autoWidth = TRUE  # Automatically adjust column widths
    ))  
  })
  
  # summary plots
  output$summary_plots <- renderPlot({
    data <- load_sample_data()
    
    # violin plot for age
    ggplot(data, aes(x = " ", y = age)) + 
      geom_violin() + 
      geom_jitter(width = 0.1, size = 1.5, alpha = 0.7, color = "black") +
      theme_minimal() +
      labs(title = "Age Distribution", y = "Age (Years)", x = "")
  })
  
  ##### counts matrix tab #####
  
  # reactive value counts matrix
  counts_data <- reactiveVal(NULL)
  
  observeEvent(input$counts_csv, {
    data <- read.csv(input$counts_csv$datapath, row.names = 1)
    counts_data(data)  # Use the full data without redundant column dropping
  })
  
  # filter genes based on variance and nonzero samples
  filtered_data <- reactive({
    req(counts_data())
    
    data <- counts_data()
    
    # cpm normalization
    library_sizes <- colSums(data)
    cpm_data <- t(t(data) / library_sizes * 1e6)
    
    # center data
    cpm_data_centered <- scale(cpm_data, center = TRUE, scale = FALSE)
    
    # calculate variance for each gene
    variances <- apply(cpm_data_centered, 1, var, na.rm = TRUE)
    
    # filter variance percentile
    var_threshold <- quantile(variances, probs = input$var_slider / 100, na.rm = TRUE)
    filtered_by_var <- data[variances >= var_threshold, ]
    
    # calculate non-zero samples for each gene
    non_zero_samples <- rowSums(filtered_by_var != 0)
    
    # filter non-zero samples
    non_zero_threshold <- ncol(filtered_by_var) * (input$zero_slider / 100)
    filtered_by_samples <- filtered_by_var[non_zero_samples >= non_zero_threshold, ]
    
    return(filtered_by_samples)
  })
  
  counts_summary_table <- reactive({
    req(filtered_data())
    
    data <- counts_data()
    filtered <- filtered_data()
    
    total_genes <- nrow(data)
    passing_genes <- nrow(filtered)
    
    summary_df <- data.frame(
      "total_genes" = total_genes,
      "passing_genes" = passing_genes,
      "percent_passing" = round((passing_genes / total_genes) * 100, 2),
      "not_passing" = total_genes - passing_genes,
      "percent_not_passing" = round(((total_genes - passing_genes) / total_genes) * 100, 2)
    )
    
    return(summary_df)
  })
  
  # counts summary table
  output$counts_summary_table <- renderTable({counts_summary_table()})
  
  # scatter plots
  output$counts_scatter_plots <- renderPlot({
    req(filtered_data())
    
    data <- counts_data()
    filtered <- filtered_data()
    
    # calculate variance and number of zeros for each gene
    variances <- apply(data, 1, var, na.rm = TRUE)
    non_zero_samples <- rowSums(data != 0)
    median_counts <- rowMeans(data)
    
    filtered_genes <- rownames(filtered)
    
    plot_data <- data.frame(
      gene = rownames(data),
      median_count = median_counts,
      variance = variances,
      num_zeros = rowSums(data == 0),
      filtered = ifelse(rownames(data) %in% filtered_genes, "Filtered Out", "Remaining")
    )
    
    # median count vs variance
    p1 <- ggplot(plot_data, aes(x = log10(median_count + 1), y = log10(variance + 1), color = filtered)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("Remaining" = "black", "Filtered Out" = "grey")) +
      labs(title = "Median Count vs Variance", x = "log10 Median Count", y = "log10 Variance") +
      theme_minimal()
    
    # median count vs number of zeros
    p2 <- ggplot(plot_data, aes(x = log10(median_count + 1), y = num_zeros, color = filtered)) +
      geom_point(alpha = 0.7) +
      scale_color_manual(values = c("Remaining" = "black", "Filtered Out" = "grey")) +
      labs(title = "Median Count vs Number of Zeros", x = "log10 Median Count", y = "Number of Zeros") +
      theme_minimal()
    
    # plot side by side
    gridExtra::grid.arrange(p1, p2, ncol = 2)
  })
  
  
  # heatmap
  output$counts_heatmap <- renderPlot({
    req(filtered_data())
    
    data <- filtered_data()%>% as.matrix()
    
    # log-transform counts
    log_data <- log10(data + 1)
    
    heatmap.2(log_data,
              col = colorRampPalette(c("blue", "white", "red"))(50),
              trace = "none",
              labRow = FALSE,
              labCol = FALSE,
              key = TRUE)
    
  })
  
  # pca biplot
  output$counts_pca_plot <- renderPlot({
    req(filtered_data())
    
    data <- filtered_data()
    
    # perform pca
    pca_res <- stats::prcomp(data, scale. = TRUE)
    
    pc_x <- as.integer(input$pc_x)
    pc_y <- as.integer(input$pc_y)
    
    # extract user selected pcs
    pca_df <- data.frame(PCx = pca_res$x[, pc_x], PCy = pca_res$x[, pc_y])
    
    # plot
    ggplot(pca_df, aes(x = PCx, y = PCy)) +
      geom_point() +
      labs(
        title = paste("PCA: PC", pc_x, "vs PC", pc_y),
        x = paste("PC", pc_x, "(", round(summary(pca_res)$importance[2, pc_x] * 100, 2), "% variance)", sep=""),
        y = paste("PC", pc_y, "(", round(summary(pca_res)$importance[2, pc_y] * 100, 2), "% variance)", sep="")
      ) +
      theme_minimal()
  })
  
  # update inputs depending on number of pcs available
  observe({
    req(filtered_data())
    
    data <- filtered_data()
    
    # pca
    pca_res <- stats::prcomp(data, scale. = TRUE)
    
    num_pcs <- ncol(pca_res$x)  # Number of PCs from PCA result
    
    updateSelectInput(session, "pc_x", choices = 1:num_pcs, selected = 1) 
    updateSelectInput(session, "pc_y", choices = 1:num_pcs, selected = 2)
  })
  
  ##### differential expression tab ######
  
  # load bulk rna-seq counts
  load_de_data <- reactive({
    req(input$upload_file)
    ext <- tools::file_ext(input$upload_file$name)
    data <- switch(ext,
                   csv = vroom::vroom(input$upload_file$datapath, delim = ","),
                   validate("Invalid file: Please upload a .csv file")
    )
    return(data)
  })
  
  # get hgnc_symbol that corresponds with ensembl gene ids
  get_gene_ids <- reactive({
    
    req(load_de_data())
    data <- load_de_data()
    
    gene_ids_with_names <- tryCatch({
      
      # start ensembl connection
      ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
      
      # get gene information from Ensembl
      gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                         filters = "ensembl_gene_id",
                         values = data$gene,
                         mart = ensembl)
      
      # join with the original data and rearrange columns
      data %>%
        left_join(gene_info, by = c("gene" = "ensembl_gene_id")) %>%
        dplyr::select(-gene) %>%
        dplyr::select(hgnc_symbol, everything()) %>%
        rename(gene = hgnc_symbol)
      
      # if ensembl error then just use ensembl gene ids
    }, error = function(e) {
      return(data)
    })
    return(gene_ids_with_names)
  })
  
  # filter counts, generate coldata, and run deseq2
  run_deseq2 <- reactive({
    
    bulk_counts <- get_gene_ids()
    
    # filter out samples with less than 10 counts
    filtered_counts <- bulk_counts[rowSums(bulk_counts[, -1] >= 10) > 0, ]
    
    gene_names <- filtered_counts[[1]]
    
    # col_data
    col_data <- tibble(sample_name = colnames(filtered_counts[-1])) %>%
      separate(sample_name, c("treatment", "replicate")) %>%
      mutate(
        treatment = as.factor(treatment)
      )
    
    # run deseq2
    deseq_dataset <- DESeqDataSetFromMatrix(countData = as.matrix(filtered_counts[-1]),
                                            colData = col_data,
                                            design = ~ treatment)
    dds <- DESeq(deseq_dataset, test = "Wald")
    res <- results(dds) %>%
      as_tibble() %>%
      mutate(gene = gene_names) %>%
      dplyr::select(gene, everything())
    
    return(res)
  })
  
  # volcano plot generation
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    volc_plt <- dataf %>%
      mutate(diff_expr = padj < 10^slider) %>%
      ggplot(aes(x = !!sym(x_name), y = -log10(!!sym(y_name)), color = diff_expr)) + 
      geom_point() +
      scale_color_manual(values = c("FALSE" = color1, "TRUE" = color2)) +
      labs(color = paste("P-adj <", formatC(10^slider, format = "e", digits = 3))) + 
      theme_bw() +
      theme(legend.position = "bottom")
    
    return(volc_plt)
  }
  
  ## output functions
  
  # differential expression datatable
  output$de_table <- renderDT({
    req(run_deseq2())
    
    df_filtered <- as.data.frame(run_deseq2()) %>%
      filter(padj < 10^input$slider) %>%
      mutate(
        pvalue = as.numeric(pvalue),  
        padj = as.numeric(padj)     
      )
    
    DT::datatable(df_filtered, 
                  options = list(
                    scrollX = TRUE,  
                    autoWidth = TRUE,
                    scrollY = "400px",
                    paging = TRUE))%>%
      formatSignif(columns = c("pvalue", "padj"), digits = 3) 
  })
  
  # volcano plot output
  output$volcano_plot <- renderPlot({
    req(run_deseq2())
    volcano_plot(as.data.frame(run_deseq2()), input$x_axis, input$y_axis, input$slider, input$base_point_color, input$highlight_point_color)})
  
  ##### trajectory analysis tab #####
  # programmatically finds the root node "Large artery" for trajectory analyses
  get_earliest_principal_node <- function(cds, time_bin="Large artery"){
    cell_ids <- which(colData(cds)[, "ECclusters"] == time_bin)
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    root_pr_nodes
  }
  
  # loads seurat .RDS object
  load_traj_data <- reactive({
    
    req(input$upload_rds)
    
    ext <- tools::file_ext(input$upload_rds$name)
    
    # only accept .RDS files
    data <- switch(ext,
                   rds = readRDS(input$upload_rds$datapath),
                   validate("Invalid file: Please upload an .rds file"))
    return(data)
  })
  
  # converts seurat .RDS object into monocle3 cell dataset object
  convert_to_cds <- reactive({
    req(load_traj_data())
    seurat_obj <- load_traj_data()
    expr_matrix <- GetAssayData(seurat_obj, layer = "counts")
    cell_metadata <- seurat_obj@meta.data
    gene_metadata <- data.frame(
      gene_short_name = rownames(expr_matrix),
      row.names = rownames(expr_matrix)
    )
    cds <- new_cell_data_set(
      expression_data = expr_matrix,
      cell_metadata = cell_metadata,
      gene_metadata = gene_metadata
    )
    
    # adds umap coordinates from seurat object to monocle3 object
    cds@int_colData@listData[["reducedDims"]][["UMAP"]] <- seurat_obj@reductions[["umap"]]@cell.embeddings
    
    return(cds)
  })
  
  # cluster and learn trajectory graph
  process_cds <- reactive({
    req(convert_to_cds())
    cds <- convert_to_cds()
    cds <- cluster_cells(cds)
    cds <- learn_graph(cds, verbose = TRUE)
    cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds))
    return(cds)
  })
  
  # find top 5 genes per cluster
  identify_top_markers <- reactive({
    req(process_cds())
    cds <- process_cds()
    marker_test_res <- top_markers(cds, group_cells_by = "partition",
                                   reference_cells = 1000, cores = 2)
    top_specific_markers <- marker_test_res %>%
      filter(fraction_expressing >= 0.10) %>%
      group_by(cell_group) %>%
      top_n(5, pseudo_R2) %>%
      arrange(cell_group, desc(pseudo_R2))
    return(top_specific_markers)
  })
  
  # filter the monocle3 processed object by user input for pathology, cluster, and genes
  filter_user_input <- eventReactive(input$apply_filters, {
    req(identify_top_markers(), process_cds(), input$pathology_input, input$cluster_input, input$gene_input)
    cds <- process_cds()
    pathology_filter <- cds[, colData(cds)$PATH %in% c(input$pathology_input)]
    input_lineage_cds <- pathology_filter[rowData(pathology_filter)$gene_short_name %in% input$gene_input, 
                                          colData(pathology_filter)$ECclusters %in% c(input$cluster_input)]
    input_lineage_cds <- order_cells(input_lineage_cds, root_pr_nodes = get_earliest_principal_node(cds))
    return(input_lineage_cds)
  })
  
  # fill in choices for pathology user input
  observeEvent(convert_to_cds(), {
    cds <- convert_to_cds()
    choices <- unique(colData(cds)$PATH)
    updateSelectInput(inputId = "pathology_input", choices = choices)
  })
  
  # fill in choices for eccluster user input
  observeEvent(convert_to_cds(), {
    cds <- convert_to_cds()
    choices <- unique(colData(cds)$ECclusters)
    updateSelectInput(inputId = "cluster_input", choices = choices)
  })
  
  # fill in choices for gene user input
  observeEvent(identify_top_markers(), {
    markers <- identify_top_markers()
    choices <- unique(markers$gene_short_name)
    updateSelectInput(inputId = "gene_input", choices = choices)
  })
  
  # output trajectory plots
  output$pseudotime_plot <- renderPlot({
    req(filter_user_input(), input$apply_filters)
    input_lineage_cds <- filter_user_input()
    
    # make a new column in the monocle3 object to combine pathology and cluster --> for plot colors
    colData(input_lineage_cds)$combined_path_cluster <- apply(colData(input_lineage_cds), 1, function(x) paste(x['PATH'], x['ECclusters'], sep = " "))
    
    # calculate how many colors needed based on unique pathology + cluster combinations
    n_colors <- length(unique(colData(input_lineage_cds)$combined_path_cluster))
    palette <- brewer.pal(min(n_colors, 9), input$color_palette)
    
    # plot
    plot_genes_in_pseudotime(input_lineage_cds,
                             color_cells_by = "combined_path_cluster",
                             min_expr = 0.5) +
      scale_color_manual(values = palette)
  })
  
  # output filtered data table
  output$traj_filtered_table <- renderDT({
    req(filter_user_input())
    input_lineage_cds <- filter_user_input()
    filtered_data <- as.data.frame(colData(input_lineage_cds))
    DT::datatable(filtered_data, options = list(
      scrollX = TRUE,
      autoWidth = TRUE,
      scrollY = "400px",
      paging = TRUE))
  })
  
}

# run application
shinyApp(ui = ui, server = server)
