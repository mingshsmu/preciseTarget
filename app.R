library(shiny)
library(dplyr)
library(shinythemes)
library(stringr)
library(ggplot2)
library(DT)
library(patchwork)
library(readxl)
library(RColorBrewer)
library(ggsci)
library(reshape2)
library(ggseqlogo)

ui <- navbarPage(title = "PreciseTarget",
                 tabPanel(title = "Introduction", fluidPage(theme = shinytheme("flatly")),
                          
                          icon = icon("dna"),
                          # p("Welcome to ",a("minglab.tech",href="http://minglab.tech/",style="color: #386cb0; font-weight: bold; text-decoration: underline;"),
                          #   "! \nPlease enjoy the ",span("Customized Applications",style="font-weight: bold; font-size:30px; color:#b22222"),"!",style = "font-size:30px; color:#18bc9c"),
                          column(12,
                                 fluidRow(
                                   p(img(src = "preciseTarget.png",width=4308/5,height=1377/5),style="text-align:center;")
                                   ),
                                 fluidRow(p("Copyright © 2024, ",a("CSRCT-SHANGHAI",style="color: #ed0000"),style = "font-size:20px; text-align:center;"))
                          )),
                 
                 tabPanel(title = "Compute",icon = icon("check"),
                          ## key: _crRNA
                          tags$head(
                            tags$style(HTML(".shiny-output-error-validation{color: red;}"))),
                          pageWithSidebar(
                            headerPanel(p("crRNA Design",style = "font-size:25px; color:#18bc9c;")),
                            sidebarPanel(
                              width = 4,
                              fluidRow(
                                column(12,
                                       textAreaInput(inputId = "input_seq_crRNA",label = "Wildtype/reference sequence",
                                                     placeholder = "DNA sequence in plaintext. Direction: 5' -> 3'"),
                                       )
                              ),
                              uiOutput(outputId = "side_crRNA"),
                              fluidRow(column(6,
                                              actionButton(inputId = "calculate_crRNA",label = "Submit & Calculate",icon = icon("calculator"),
                                                           class = "btn-info"))
                              )
                            ),
                            
                            mainPanel(
                              tabsetPanel(
                                type = "tabs",
                                tabPanel(title = "CheckInput",icon = icon("list"),
                                         uiOutput(outputId = "ui_crRNA_checkinput")),
                                tabPanel(title = "Activity", icon = icon("square-poll-vertical"),
                                         uiOutput(outputId = "ui_crRNA_activity")),
                                tabPanel(title = "Discrimination",icon = icon("chart-simple"),
                                         uiOutput(outputId = "ui_crRNA_diff")),
                                tabPanel(title = "Plot (activity)",icon = icon("images"),
                                         uiOutput(outputId = "ui_crRNA_plot_activity_control"),
                                         br(),
                                         uiOutput(outputId = "ui_crRNA_plot_activity_plot")
                                ),
                                tabPanel(title = "Plot (discrimination)",icon = icon("layer-group"),
                                         uiOutput(outputId = "ui_crRNA_plot_diff_control"),
                                         br(),
                                         uiOutput(outputId = "ui_crRNA_plot_diff_plot")
                                ),
                              ),
                            )
                          ),
                 ),  
                 tabPanel(title = "Compute (batch)",icon = icon("list-check"),
                          ## key: _crRNAm
                          tags$head(
                            tags$style(HTML(".shiny-output-error-validation{color: red;}"))),
                          pageWithSidebar(
                            headerPanel(p("crRNA Design (batch)",style = "font-size:25px; color:#18bc9c;")),
                            sidebarPanel(
                              width = 4,
                              fluidRow(
                                column(12,
                                       fileInput(inputId = "file_crRNAm",label = "Upload your file",multiple = F,
                                                 buttonLabel = "Browse...",placeholder = "No file selected")
                                )
                              ),
                              fluidRow(column(6,
                                              actionButton(inputId = "calculate_crRNAm",label = "Submit & Calculate",icon = icon("calculator"),
                                                           class = "btn-info")),
                                       column(6,
                                              downloadButton(outputId = "demo_crRNAm",label = "Demo input"))
                              )
                            ),
                            
                            mainPanel(
                              tabsetPanel(
                                type = "tabs",
                                tabPanel(title = "CheckInput",icon = icon("list"),
                                         uiOutput(outputId = "ui_crRNAm_rawdata")),
                                tabPanel(title = "Activity", icon = icon("square-poll-vertical"),
                                         uiOutput(outputId = "ui_crRNAm_activity")),
                                tabPanel(title = "Discrimination",icon = icon("chart-simple"),
                                         uiOutput(outputId = "ui_crRNAm_diff")),
                                tabPanel(title = "Plot (activity)",icon = icon("images"),
                                         uiOutput(outputId = "ui_crRNAm_plot_activity_control"),
                                         br(),
                                         uiOutput(outputId = "ui_crRNAm_plot_activity_plot")
                                ),
                                tabPanel(title = "Plot (discrimination)",icon = icon("layer-group"),
                                         uiOutput(outputId = "ui_crRNAm_plot_diff_control"),
                                         br(),
                                         uiOutput(outputId = "ui_crRNAm_plot_diff_plot")
                                ),
                              ),
                            )
                          ),
                 ), 
                 tabPanel(title = "PrimerDesign",icon = icon("wand-magic-sparkles"),
                          ## key: _primer
                          tags$head(
                            tags$style(HTML(".shiny-output-error-validation{color: red;}"))),
                          pageWithSidebar(
                            headerPanel(p("Primer Design",style = "font-size:25px; color:#18bc9c;")),
                            sidebarPanel(
                              width = 4,
                              fluidRow(
                                column(12,
                                       fileInput(inputId = "file_primer",label = "Upload your file",
                                                 multiple = FALSE, buttonLabel = "Browse...",
                                                 placeholder = "Upload the result downloaded from the 'Compute' module."))
                              ),
                              fluidRow(
                                column(6,
                                       selectInput(inputId = "cas_protein_primer",label = "CRISPR enzyme",
                                                   choices = "suCas12a2",multiple = F)),
                                column(6,
                                       selectInput(inputId = "trans_method_primer",label = "Transcription",
                                                   choices = c("T7 RNA polymerase")))
                              ),
                              fluidRow(column(6,
                                              actionButton(inputId = "calculate_primer",label = "Design",icon = icon("feather"),
                                                           class = "btn-info"))
                              )
                            ),
                            
                            mainPanel(
                              tabsetPanel(
                                type = "tabs",
                                tabPanel(title = "CheckInput",icon = icon("list"),
                                         uiOutput(outputId = "ui_primer_rawdata")),
                                tabPanel(title = "Primers", icon = icon("square-poll-vertical"),
                                         uiOutput(outputId = "ui_primer_design")),
                              ),
                            )
                          ),
                          
                 ),
                 tabPanel(title = "Palette",icon = icon("palette"),
                          img(src = "00_palette.png",width=6000/4,height=6000/4)
                 ),
                 tabPanel(title = "Developers",icon = icon("user"),
                          
                          
                 ),
                 
                 
)

server <- function(input, output, session) {
  # functions
  source("functions.R")
  # palettes
  palette_in_discrete <- reactive({c("Lancet (9)","NPG (10)", "NEJM (8)","JCO (10)","JAMA (7)",
                                     "AAAS (10)","D3 (10)","Futurama (12)","GSEA (12)","IGV (51)",
                                     "LocusZoom (7)","Rick and Morty (12)","Simpsons (16)",
                                     "Star Trek (7)","Tron Legacy (7)","Chicago (9)","UCSC (26)",
                                     "Set1 (9)","Set2 (8)","Set3 (12)","Pastel1 (9)","Pastel2 (8)",
                                     "Paired (12)","Dark2 (8)","Accent (8)")})
  palette_in_continuous <- reactive({
    div_palette <- paste0("Div_",rownames(RColorBrewer::brewer.pal.info %>% filter(category %in% c("div"))))
    seq_palette <- paste0("Seq_",rownames(RColorBrewer::brewer.pal.info %>% filter(category %in% c("seq"))))
    seq_palette_colorbrewer2 <- paste0("Seq_",
                                       c("BuGn","BuPu","GnBu","OrRd","PuBu","PuBuGn",
                                         "PuRd","RdPu","YlGn","YlGnBu","YlOrBr","YlOrRd",
                                         "Blues","Greens","Greys","Oranges","Purples","Reds"),"3")
    div_palette_colorbrewer2 <- paste0("Div_",
                                       c("BrBG","PiYG","PRGn","PuOr","RdBu",
                                         "RdGy","RdYlBu","RdYlGn","Spectral"),"3")
    viridis_palette <- c("Viridis_magma","Viridis_inferno","Viridis_plasma",
                         "Viridis_viridis","Viridis_cividis","Viridis_rocket",
                         "Viridis_mako","Viridis_turbo")
    palette <- c("navy white firebrick3",div_palette_colorbrewer2,seq_palette_colorbrewer2,
                 viridis_palette,div_palette,seq_palette)
    return(palette)
  })
  palette_word <- reactive({c("black","white","red2","green3","blue","cyan2","magenta3","orange","gray")})
  palette_trans_fun <- reactive({
    function(palette){
      switch (palette,
              "Lancet (9)" = ggsci::pal_lancet()(9),
              "NPG (10)" = ggsci::pal_npg()(10),
              "NEJM (8)" = ggsci::pal_nejm()(8),
              "JCO (10)" = ggsci::pal_jco()(10),
              "JAMA (7)" = ggsci::pal_jama()(7),
              "AAAS (10)" = ggsci::pal_aaas()(10),
              "D3 (10)" = ggsci::pal_d3()(10),
              "Futurama (12)" = ggsci::pal_futurama()(12),
              "GSEA (12)" = ggsci::pal_gsea()(12),
              "IGV (51)" = ggsci::pal_igv()(51),
              "LocusZoom (7)" = ggsci::pal_locuszoom()(7),
              "Rick and Morty (12)" = ggsci::pal_rickandmorty()(12),
              "Simpsons (16)" = ggsci::pal_simpsons()(16),
              "Star Trek (7)" = ggsci::pal_startrek()(7),
              "Tron Legacy (7)" = ggsci::pal_tron()(7),
              "Chicago (9)" = ggsci::pal_uchicago()(9),
              "UCSC (26)" = ggsci::pal_ucscgb()(26),
              "Set1 (9)" = RColorBrewer::brewer.pal(9,"Set1"),
              "Set2 (8)" = RColorBrewer::brewer.pal(8,"Set2"),
              "Set3 (12)" = RColorBrewer::brewer.pal(12,"Set3"),
              "Pastel1 (9)" = RColorBrewer::brewer.pal(9,"Pastel1"),
              "Pastel2 (8)" = RColorBrewer::brewer.pal(8,"Pastel2"),
              "Paired (12)" = RColorBrewer::brewer.pal(12,"Paired"),
              "Dark2 (8)" = RColorBrewer::brewer.pal(8,"Dark2"),
              "Accent (8)" = RColorBrewer::brewer.pal(8,"Accent"),
              "Div_BrBG" = RColorBrewer::brewer.pal(11,"BrBG"),
              "Div_PiYG" = RColorBrewer::brewer.pal(11,"PiYG"),
              "Div_PRGn" = RColorBrewer::brewer.pal(11,"PRGn"),
              "Div_PuOr" = RColorBrewer::brewer.pal(11,"PuOr"),
              "Div_RdBu" = RColorBrewer::brewer.pal(11,"RdBu"),
              "Div_RdGy" = RColorBrewer::brewer.pal(11,"RdGy"),
              "Div_RdYlBu" = RColorBrewer::brewer.pal(11,"RdYlBu"),
              "Div_RdYlGn" = RColorBrewer::brewer.pal(11,"RdYlGn"),
              "Div_Spectral" = RColorBrewer::brewer.pal(11,"Spectral"),
              "Seq_Blues" = RColorBrewer::brewer.pal(9,"Blues"),
              "Seq_BuGn" = RColorBrewer::brewer.pal(9,"BuGn"),
              "Seq_BuPu" = RColorBrewer::brewer.pal(9,"BuPu"),
              "Seq_GnBu" = RColorBrewer::brewer.pal(9,"GnBu"),
              "Seq_Greens" = RColorBrewer::brewer.pal(9,"Greens"),
              "Seq_Greys" = RColorBrewer::brewer.pal(9,"Greys"),
              "Seq_Oranges" = RColorBrewer::brewer.pal(9,"Oranges"),
              "Seq_OrRd" = RColorBrewer::brewer.pal(9,"OrRd"),
              "Seq_PuBu" = RColorBrewer::brewer.pal(9,"PuBu"),
              "Seq_PuBuGn" = RColorBrewer::brewer.pal(9,"PuBuGn"),
              "Seq_PuRd" = RColorBrewer::brewer.pal(9,"PuRd"),
              "Seq_Purples" = RColorBrewer::brewer.pal(9,"Purples"),
              "Seq_RdPu" = RColorBrewer::brewer.pal(9,"RdPu"),
              "Seq_Reds" = RColorBrewer::brewer.pal(9,"Reds"),
              "Seq_YlGn" = RColorBrewer::brewer.pal(9,"YlGn"),
              "Seq_YlGnBu" = RColorBrewer::brewer.pal(9,"YlGnBu"),
              "Seq_YlOrBr" = RColorBrewer::brewer.pal(9,"YlOrBr"),
              "Seq_YlOrRd" = RColorBrewer::brewer.pal(9,"YlOrRd"),
              "Viridis_magma" = viridis::viridis_pal(option = "A")(12),
              "Viridis_inferno" = viridis::viridis_pal(option = "B")(12),
              "Viridis_plasma" = viridis::viridis_pal(option = "C")(12),
              "Viridis_viridis" = viridis::viridis_pal(option = "D")(12),
              "Viridis_cividis" = viridis::viridis_pal(option = "E")(12),
              "Viridis_rocket" = viridis::viridis_pal(option = "F")(12),
              "Viridis_mako" = viridis::viridis_pal(option = "G")(12),
              "Viridis_turbo" = viridis::viridis_pal(option = "H")(12),
              "Div_BrBG3" = c("#d8b365","#f5f5f5","#5ab4ac"),
              "Div_PiYG3" = c("#e9a3c9","#f7f7f7","#a1d76a"),
              "Div_PRGn3" = c("#af8dc3","#f7f7f7","#7fbf7b"),
              "Div_PuOr3" = c("#f1a340","#f7f7f7","#998ec3"),
              "Div_RdBu3" = c("#ef8a62","#f7f7f7","#67a9cf"),
              "Div_RdGy3" = c("#ef8a62","#ffffff","#999999"),
              "Div_RdYlBu3" = c("#fc8d59","#ffffbf","#91bfdb"),
              "Div_RdYlGn3" = c("#fc8d59","#ffffbf","#91cf60"),
              "Div_Spectral3" = c("#fc8d59","#ffffbf","#99d594"),
              "Seq_BuGn3" = c("#e5f5f9","#99d8c9","#2ca25f"),
              "Seq_BuPu3" = c("#e0ecf4","#9ebcda","#8856a7"),
              "Seq_GnBu3" = c("#e0f3db","#a8ddb5","#43a2ca"),
              "Seq_OrRd3" = c("#fee8c8","#fdbb84","#e34a33"),
              "Seq_PuBu3" = c("#ece7f2","#a6bddb","#2b8cbe"),
              "Seq_PuBuGn3" = c("#ece2f0","#a6bddb","#1c9099"),
              "Seq_PuRd3" = c("#e7e1ef","#c994c7","#dd1c77"),
              "Seq_RdPu3" = c("#fde0dd","#fa9fb5","#c51b8a"),
              "Seq_YlGn3" = c("#f7fcb9","#addd8e","#31a354"),
              "Seq_YlGnBu3" = c("#edf8b1","#7fcdbb","#2c7fb8"),
              "Seq_YlOrBr3" = c("#fff7bc","#fec44f","#d95f0e"),
              "Seq_YlOrRd3" = c("#ffeda0","#feb24c","#f03b20"),
              "Seq_Blues3" = c("#deebf7","#9ecae1","#3182bd"),
              "Seq_Greens3" = c("#e5f5e0","#a1d99b","#31a354"),
              "Seq_Greys3" = c("#f0f0f0","#bdbdbd","#636363"),
              "Seq_Oranges3" = c("#fee6ce","#fdae6b","#e6550d"),
              "Seq_Purples3" = c("#efedf5","#bcbddc","#756bb1"),
              "Seq_Reds3" = c("#fee0d2","#fc9272","#de2d26"),
              "navy white firebrick3" = c("navy","white","firebrick3")
      )
    }
  })
  
  # 01. crRNA design +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # 01. key: _crRNA
  
  # ui_side_crRNA
  output$side_crRNA <- renderUI({
    column(12,
           fluidRow(
             column(4,
                    numericInput(inputId = "mutation_site_crRNA",label = "Mutation position",
                                 step = 1,min = 1,value = 1)),
             column(4,
                    selectInput(inputId = "mutate_to_crRNA",label = "Mutate to",
                                choices = c("A","T/U","G","C"),multiple = F)),
             column(4,
                    textInput(inputId = "label_crRNA",label = "Label"))
           ),
           fluidRow(
             column(6,
                    selectInput(inputId = "cas_protein_crRNA",label = "CRISPR enzyme",
                                choices = c("suCas12a2")))
           )
           )
  })
  
  # find_crRNA
  # variables_crRNA
  input_seq_crRNA <- eventReactive(input$calculate_crRNA,{input$input_seq_crRNA})
  mutation_site_crRNA <- eventReactive(input$calculate_crRNA,{input$mutation_site_crRNA})
  mutate_to_crRNA <- eventReactive(input$calculate_crRNA,{input$mutate_to_crRNA})
  label_crRNA <- eventReactive(input$calculate_crRNA,{input$label_crRNA})
  cas_protein_crRNA <- eventReactive(input$calculate_crRNA,{input$cas_protein_crRNA})
  
  crRNAs_list <- eventReactive(input$calculate_crRNA,{
    find_crRNA(input_seq = input_seq_crRNA(),
               mutation_site = mutation_site_crRNA(),
               mutate_to = mutate_to_crRNA(),
               label = label_crRNA(),
               Cas_protein = cas_protein_crRNA())
  })
  
  # plot
  ## plot seq
  plot_seq_crRNA <- reactive({
    seq_wt_mut <- crRNAs_list()$seqs_wt_mut
    return(plot_seq(seq_wt_mut = seq_wt_mut))
  })
  
  output$plot_seq_crRNA <- renderPlot({
    plot_seq_crRNA()
  })
  
  output$ui_crRNA_checkinput <- renderUI({
    eventReactive(input$calculate_crRNA,{
      column(12,
             h3("Input Sequence (WT↑ Mut↓):"),
             downloadButton(outputId = "down_crRNA_plot_seq_png",label = "Download PNG"), 
             downloadButton(outputId = "down_crRNA_plot_seq_pdf",label = "Download PDF"),
             br(),br(),br(),
             plotOutput("plot_seq_crRNA",width = 1000,height = 80)
      )
    })()
  })
  
  # activity
  activity_crRNA <- reactive({
    crRNAs_list()$result
  })
  
  output$activity_crRNA <- renderDataTable({
    dat <- activity_crRNA()
    dat$PFS_activity <- signif(dat$PFS_activity,3)
    dat$guide_activity <- signif(dat$guide_activity,3)
    dat$crRNA_activity <- signif(dat$crRNA_activity,3)
    DT::datatable(data = dat,options = list(pageLength=25))
  })
  
  output$ui_crRNA_activity <- renderUI({
    eventReactive(input$calculate_crRNA,{
      column(12,
             fluidRow(
               downloadButton(outputId = "down_crRNA_activity",
                              label = "Download Table")),
             br(),
             dataTableOutput(outputId = "activity_crRNA",width = "80%"))
    })()
    
  })
  
  # diff
  diff_crRNA <- reactive({
    crRNAs_list()$dat_diff
  })
  
  output$diff_crRNA <- renderDataTable({
    dat <- diff_crRNA()
    dat$activity_input_seq <- signif(dat$activity_input_seq,3)
    dat$activity_mutant_seq <- signif(dat$activity_mutant_seq,3)
    dat$diff_score <- signif(dat$diff_score,3)
    dat <- dat[,-c(8,9,10)]
    DT::datatable(data = dat,options = list(pageLength=25))
  })
  
  output$ui_crRNA_diff <- renderUI({
    eventReactive(input$calculate_crRNA,{
      column(12,
             fluidRow(
               downloadButton(outputId = "down_crRNA_diff",
                              label = "Download Table")),
             br(),
             dataTableOutput(outputId = "diff_crRNA",width = "80%"))
    })()
    
  })
  
  # plot activity
  output$ui_crRNA_plot_activity_control <- renderUI({
    palette_name <- palette_in_discrete()
    eventReactive(input$calculate_crRNA,{
      column(12,
             fluidRow(
               column(3,selectInput(inputId = "palette_crRNA_activity",label = "Palette:",
                                    choices = palette_name, 
                                    selected = "AAAS (10)",multiple = FALSE)),
               column(3,selectInput(inputId = "palette_rev_crRNA_activity",label = "Palette reverse:",
                                    choices = c("No","Yes"), selected = "No",multiple = FALSE)),
               column(3,numericInput(inputId = "alpha_crRNA_activity",label = "Alpha:",
                                     value = 1,min = 0,max = 1,step = 0.1))
             ),
             fluidRow(
               column(4,numericInput(inputId = "width_crRNA_activity",label = "Width (download):",
                                     value = 20,min = 1,step = 1)),
               column(4,numericInput(inputId = "height_crRNA_activity",label = "Height (download):",
                                     value = 5,min = 1,step = 1)),
             ))
    })()
  })
  
  plot_activity_crRNA <- reactive({
    plot_crRNA_activity(res = crRNAs_list()$result,
                        palette = palette_trans_fun()(input$palette_crRNA_activity),
                        palette_rev = input$palette_rev_crRNA_activity,
                        alpha = input$alpha_crRNA_activity)
  })
  output$plot_activity_crRNA <- renderPlot({
    plot_activity_crRNA()
  })
  
  output$ui_crRNA_plot_activity_plot <- renderUI({
    eventReactive(input$calculate_crRNA,{
      column(12,
             downloadButton(outputId = "down_crRNA_plot_activity_png",label = "Download PNG"), 
             downloadButton(outputId = "down_crRNA_plot_activity_pdf",label = "Download PDF"),
             br(),
             plotOutput("plot_activity_crRNA",width = 1200,height = 300)
      )
    })()
  })
  # plot diff
  output$ui_crRNA_plot_diff_control <- renderUI({
    palette_name <- palette_in_discrete()
    eventReactive(input$calculate_crRNA,{
      column(12,
             fluidRow(
               column(3,selectInput(inputId = "palette_crRNA_diff",label = "Palette:",
                                    choices = palette_name, 
                                    selected = "AAAS (10)",multiple = FALSE)),
               column(3,selectInput(inputId = "palette_rev_crRNA_diff",label = "Palette reverse:",
                                    choices = c("No","Yes"), selected = "No",multiple = FALSE)),
               column(3,numericInput(inputId = "alpha_crRNA_diff",label = "Alpha:",
                                     value = 1,min = 0,max = 1,step = 0.1))
             ),
             fluidRow(
               column(4,numericInput(inputId = "width_crRNA_diff",label = "Width (download):",
                                     value = 10,min = 1,step = 1)),
               column(4,numericInput(inputId = "height_crRNA_diff",label = "Height (download):",
                                     value = 5,min = 1,step = 1)),
             ))
    })()
  })
  
  plot_diff_crRNA <- reactive({
    plot_crRNA_diff(crRNAs_list()$dat_diff,
                        palette = palette_trans_fun()(input$palette_crRNA_diff),
                        palette_rev = input$palette_rev_crRNA_diff,
                        alpha = input$alpha_crRNA_diff)
  })
  output$plot_diff_crRNA <- renderPlot({
    plot_diff_crRNA()
  })
  
  output$ui_crRNA_plot_diff_plot <- renderUI({
    eventReactive(input$calculate_crRNA,{
      column(12,
             downloadButton(outputId = "down_crRNA_plot_diff_png",label = "Download PNG"), 
             downloadButton(outputId = "down_crRNA_plot_diff_pdf",label = "Download PDF"),
             br(),
             plotOutput("plot_diff_crRNA",width = 1000,height = 500)
      )
    })()
  })
  
  # Download crRNA ---
  output$down_crRNA_plot_seq_png <- downloadHandler(
    filename = paste0("check_seq_",input$label_crRNA,"_",Sys.time(),".png"),
    content = function(file){
      ggsave(filename = file, device = "png",
             width = 25, height = 2,
             dpi = 300, plot = plot_seq_crRNA())
    }
  )
  output$down_crRNA_plot_seq_pdf <- downloadHandler(
    filename = paste0("check_seq_",input$label_crRNA,"_",Sys.time(),".pdf"),
    content = function(file){
      ggsave(filename = file, device = "pdf",
             width = 25, height = 2,
             plot = plot_seq_crRNA())
    }
  )
  output$down_crRNA_activity <- downloadHandler(
    filename = paste0("crRNA_activity_",input$label_crRNA,"_",Sys.time(),".csv"),
    content = function(file){
      write.csv(activity_crRNA(),file)
    }
  )
  
  output$down_crRNA_diff <- downloadHandler(
    filename = paste0("crRNA_discrimination_",input$label_crRNA,"_",Sys.time(),".csv"),
    content = function(file){
      write.csv(diff_crRNA(),file)
    }
  )
  
  output$down_crRNA_plot_activity_png <- downloadHandler(
    filename = paste0("crRNA_plot_activity_",input$label_crRNA,"_",Sys.time(),".png"),
    content = function(file){
      ggsave(filename = file, device = "png",
             width = input$width_crRNA_activity, height = input$height_crRNA_activity,
             dpi = 300, plot = plot_activity_crRNA())
    }
  )
  output$down_crRNA_plot_activity_pdf <- downloadHandler(
    filename = paste0("crRNA_plot_activity_",input$label_crRNA,"_",Sys.time(),".pdf"),
    content = function(file){
      ggsave(filename = file, device = "pdf",
             width = input$width_crRNA_activity, height = input$height_crRNA_activity,
             plot = plot_activity_crRNA())
    }
  )
  output$down_crRNA_plot_diff_png <- downloadHandler(
    filename = paste0("crRNA_plot_diff_",input$label_crRNA,"_",Sys.time(),".png"),
    content = function(file){
      ggsave(filename = file, device = "png",
             width = input$width_crRNA_diff, height = input$height_crRNA_diff,
             dpi = 300, plot = plot_diff_crRNA())
    }
  )
  output$down_crRNA_plot_diff_pdf <- downloadHandler(
    filename = paste0("crRNA_plot_diff_",input$label_crRNA,"_",Sys.time(),".pdf"),
    content = function(file){
      ggsave(filename = file, device = "pdf",
             width = input$width_crRNA_diff, height = input$height_crRNA_diff,
             plot = plot_diff_crRNA())
    }
  )  
  
  # 02. crRNA multi +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # key: _crRNAm
  infile_crRNAm <- reactive({
    input$file_crRNAm$datapath
  })

  input_paras <- reactive({
    data_read <- function(path1){
      filetype_split <- unlist(strsplit(path1,split = "[.]"))
      filetype <- filetype_split[length(filetype_split)]
      dat <- switch (filetype,
                     "dat" = read.table(path1,header = T),
                     "csv" = read.csv(path1,header = T),
                     "xlsx" = as.data.frame(readxl::read_excel(path1,sheet = 1,col_names = T)),
                     "xls" = as.data.frame(readxl::read_excel(path1,sheet = 1,col_names = T)),
                     "txt" = read.table(path1,header = T)
      )
      dat <- dat[,c("label","mutation_site","mutate_to","input_seq")]
      return(dat)
    }
    data_read(path1 = infile_crRNAm())
  })
  # display parameter df
  output$crRNAm_rawdata <- renderDataTable({
    dat <- data.frame()
    if(!is.null(infile_crRNAm())){
      dat <- input_paras()
    }
    DT::datatable(data = dat,
                  options = list(pageLength = 25))
  })
  
  output$ui_crRNAm_rawdata <- renderUI({
    dataTableOutput(outputId = "crRNAm_rawdata",width = "50%")
  })
  
  # calculate crRNA activity
  res_crRNAm <- eventReactive(input$calculate_crRNAm,{
    find_crRNA_multi(parameters_df = input_paras())
  })
  
  # activity crRNAm
  activity_crRNAm <- reactive({
    res_crRNAm()$result
  })
  
  output$activity_crRNAm <- renderDataTable({
    dat <- activity_crRNAm()
    dat$PFS_activity <- signif(dat$PFS_activity,3)
    dat$guide_activity <- signif(dat$guide_activity,3)
    dat$crRNA_activity <- signif(dat$crRNA_activity,3)
    DT::datatable(data = dat,options = list(pageLength=25))
  })
  
  output$ui_crRNAm_activity <- renderUI({
    eventReactive(input$calculate_crRNAm,{
      column(12,
             fluidRow(
               downloadButton(outputId = "down_crRNAm_activity",
                              label = "Download Table")),
             br(),
             dataTableOutput(outputId = "activity_crRNAm",width = "80%"))
    })()
    
  })
  
  # diff crRNAm
  diff_crRNAm <- reactive({
    res_crRNAm()$dat_diff
  })
  
  output$diff_crRNAm <- renderDataTable({
    dat <- diff_crRNAm()
    dat$activity_input_seq <- signif(dat$activity_input_seq,3)
    dat$activity_mutant_seq <- signif(dat$activity_mutant_seq,3)
    dat$diff_score <- signif(dat$diff_score,3)
    dat <- dat[,-c(8,9,10)]
    DT::datatable(data = dat,options = list(pageLength=25))
  })
  
  output$ui_crRNAm_diff <- renderUI({
    eventReactive(input$calculate_crRNAm,{
      column(12,
             fluidRow(
               downloadButton(outputId = "down_crRNAm_diff",
                              label = "Download Table")),
             br(),
             dataTableOutput(outputId = "diff_crRNAm",width = "80%"))
    })()
  })
  
  # plot activity crRNAm
  output$ui_crRNAm_plot_activity_control <- renderUI({
    palette_name <- palette_in_discrete()
    eventReactive(input$calculate_crRNAm,{
      column(12,
             fluidRow(
               column(3,selectInput(inputId = "palette_crRNAm_activity",label = "Palette:",
                                    choices = palette_name, 
                                    selected = "AAAS (10)",multiple = FALSE)),
               column(3,selectInput(inputId = "palette_rev_crRNAm_activity",label = "Palette reverse:",
                                    choices = c("No","Yes"), selected = "No",multiple = FALSE)),
               column(3,numericInput(inputId = "alpha_crRNAm_activity",label = "Alpha:",
                                     value = 1,min = 0,max = 1,step = 0.1))
             ),
             fluidRow(
               column(4,numericInput(inputId = "width_crRNAm_activity",label = "Width (download):",
                                     value = 20,min = 1,step = 1)),
               column(4,numericInput(inputId = "height_crRNAm_activity",label = "Height (download):",
                                     value = 16,min = 1,step = 1)),
             ))
    })()
  })
  
  plot_activity_crRNAm <- reactive({
    plot_crRNA_activity(res = res_crRNAm()$result,
                        palette = palette_trans_fun()(input$palette_crRNAm_activity),
                        palette_rev = input$palette_rev_crRNAm_activity,
                        alpha = input$alpha_crRNAm_activity)
  })
  output$plot_activity_crRNAm <- renderPlot({
    plot_activity_crRNAm()
  })
  
  output$ui_crRNAm_plot_activity_plot <- renderUI({
    height1 <- length(unique(res_crRNAm()$result$label))*100
    
    eventReactive(input$calculate_crRNAm,{
      column(12,
             downloadButton(outputId = "down_crRNAm_plot_activity_png",label = "Download PNG"), 
             downloadButton(outputId = "down_crRNAm_plot_activity_pdf",label = "Download PDF"),
             br(),
             plotOutput("plot_activity_crRNAm",width = 1200,height = height1)
      )
    })()
  })
  # plot diff crRNAm
  output$ui_crRNAm_plot_diff_control <- renderUI({
    palette_name <- palette_in_discrete()
    eventReactive(input$calculate_crRNAm,{
      column(12,
             fluidRow(
               column(3,selectInput(inputId = "palette_crRNAm_diff",label = "Palette:",
                                    choices = palette_name, 
                                    selected = "AAAS (10)",multiple = FALSE)),
               column(3,selectInput(inputId = "palette_rev_crRNAm_diff",label = "Palette reverse:",
                                    choices = c("No","Yes"), selected = "No",multiple = FALSE)),
               column(3,numericInput(inputId = "alpha_crRNAm_diff",label = "Alpha:",
                                     value = 1,min = 0,max = 1,step = 0.1))
             ),
             fluidRow(
               column(4,numericInput(inputId = "width_crRNAm_diff",label = "Width (download):",
                                     value = 10,min = 1,step = 1)),
               column(4,numericInput(inputId = "height_crRNAm_diff",label = "Height (download):",
                                     value = 10,min = 1,step = 1)),
             ))
    })()
  })
  
  plot_diff_crRNAm <- reactive({
    plot_crRNA_diff_multi(res_crRNAm()$dat_diff,
                    palette = palette_trans_fun()(input$palette_crRNAm_diff),
                    palette_rev = input$palette_rev_crRNAm_diff,
                    alpha = input$alpha_crRNAm_diff)
  })
  output$plot_diff_crRNAm <- renderPlot({
    plot_diff_crRNAm()
  })
  
  output$ui_crRNAm_plot_diff_plot <- renderUI({
    height1 <- ceiling(length(unique(res_crRNAm()$result$label))/3)*250
    
    eventReactive(input$calculate_crRNAm,{
      column(12,
             downloadButton(outputId = "down_crRNAm_plot_diff_png",label = "Download PNG"), 
             downloadButton(outputId = "down_crRNAm_plot_diff_pdf",label = "Download PDF"),
             br(),
             plotOutput("plot_diff_crRNAm",width = 1000,height = height1)
      )
    })()
  })
  
  # Download crRNAm ---
  output$down_crRNAm_activity <- downloadHandler(
    filename = paste0("crRNAm_activity_",Sys.time(),".csv"),
    content = function(file){
      write.csv(activity_crRNAm(),file)
    }
  )
  
  output$demo_crRNAm <- downloadHandler(
    filename = paste0("Demo_compute_batch.csv"),
    content = function(file){
      file.copy(from = "./metadata/05_parameter_df.csv",to=file)
    }
  )
  
  output$down_crRNAm_diff <- downloadHandler(
    filename = paste0("crRNAm_discrimination_",Sys.time(),".csv"),
    content = function(file){
      write.csv(diff_crRNAm(),file)
    }
  )
  
  output$down_crRNAm_plot_activity_png <- downloadHandler(
    filename = paste0("crRNAm_plot_activity_",Sys.time(),".png"),
    content = function(file){
      ggsave(filename = file, device = "png",
             width = input$width_crRNAm_activity, height = input$height_crRNAm_activity,
             dpi = 300, plot = plot_activity_crRNAm())
    }
  )
  output$down_crRNAm_plot_activity_pdf <- downloadHandler(
    filename = paste0("crRNAm_plot_activity_",Sys.time(),".pdf"),
    content = function(file){
      ggsave(filename = file, device = "pdf",
             width = input$width_crRNAm_activity, height = input$height_crRNAm_activity,
             plot = plot_activity_crRNAm())
    }
  )
  output$down_crRNAm_plot_diff_png <- downloadHandler(
    filename = paste0("crRNAm_plot_diff_",Sys.time(),".png"),
    content = function(file){
      ggsave(filename = file, device = "png",
             width = input$width_crRNAm_diff, height = input$height_crRNAm_diff,
             dpi = 300, plot = plot_diff_crRNAm())
    }
  )
  output$down_crRNAm_plot_diff_pdf <- downloadHandler(
    filename = paste0("crRNAm_plot_diff_",Sys.time(),".pdf"),
    content = function(file){
      ggsave(filename = file, device = "pdf",
             width = input$width_crRNAm_diff, height = input$height_crRNAm_diff,
             plot = plot_diff_crRNAm())
    }
  )  
  
  # primer design ==============================================================
  # key: _primer
  infile_primer <- reactive({
    input$file_primer$datapath
  })
  
  input_data_primer <- reactive({
    read.csv(file = infile_primer(),header = T,row.names = 1)
  })
  
  # display parameter df
  output$primer_rawdata <- renderDataTable({
    dat <- data.frame()
    if(!is.null(infile_primer())){
      dat <- input_data_primer()
    }
    DT::datatable(data = dat,
                  options = list(pageLength = 25))
  })
  
  output$ui_primer_rawdata <- renderUI({
    dataTableOutput(outputId = "primer_rawdata",width = "50%")
  })
  
  # primer design
  cas_protein_primer <- reactive({input$cas_protein_primer})
  trans_method_primer <- reactive({input$trans_method_primer})
  designed_primer <- eventReactive(input$calculate_primer,{
    crRNA_primer_design(crRNAs = input_data_primer()$crRNA,trans_method = trans_method_primer(),
                        Cas_protein = cas_protein_primer())
  })
  
  # display designed primers
  output$designed_primer <- renderDataTable({
    DT::datatable(data = designed_primer(),
                  options = list(pageLength = 25))
  })
  
  output$ui_primer_design <- renderUI({
    dataTableOutput(outputId = "designed_primer",width = "50%")
  })
  
}

shinyApp(ui, server)