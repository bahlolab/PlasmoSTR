library(tidyverse)
library(shiny)
library(rsconnect)
library(plotly)
library(shinythemes)
library(ggpubr)
library(DT)
library(dplyr)
library(farver)

#Load data
good_pv <- read.csv("Data/Pv_STR_Top_Table.csv")
good_pv$Heterozygosity <- round(good_pv$Heterozygosity,3)
good_pf <- read.csv("Data/Pf_STR_Top_Table.csv")
good_pf$Heterozygosity <- round(good_pf$Heterozygosity,3)
#Population or Country pairs
compare_Pv_Country <- colnames(good_pv)[13:ncol(good_pv)]
compare_Pf_Country <- colnames(good_pf)[41:ncol(good_pf)]
compare_Pf_Population <- colnames(good_pf)[13:40]
#Plot data
#Pf Population
#Pf population loci plot
geno_pf <- read.csv("Data/Pf_Genotype.csv")
geno_origin_Pf <- geno_pf
Pf_pop_loci <- read.csv("Data/Pf_STR_Top_ShinyPlot.csv")
Pf_pop_loci <- good_pf %>% 
  filter(STRid %in% as.character(Pf_pop_loci$STRid))
Pf_pop_loci <- Pf_pop_loci %>% 
  pivot_longer(cols = CAF.EAF:WSEA.SAM) %>% 
  filter(is.na(value)==FALSE) %>% 
  arrange(value)
vars_pf_population <- setdiff(unique(Pf_pop_loci$value), "STRLoci")
geno_pf_population <- geno_origin_Pf
geno_pf_population <- geno_pf_population %>% 
  filter(STRid %in% as.character(Pf_pop_loci$STRid))
STR_Position <- Pf_pop_loci %>% 
  dplyr::select(STRid,value)
geno_pf_population <- STR_Position %>% 
  left_join(geno_pf_population,by="STRid")
geno_pf_population <- geno_pf_population[,-1]
geno_pf_population$value <- as.character(geno_pf_population$value)
geno_pf_population <- as.data.frame(geno_pf_population)
rownames(geno_pf_population) <- geno_pf_population$value
geno_pf_population <- geno_pf_population[,-1]
geno_pf_population <- t(geno_pf_population)
geno_pf_population <- data.frame(Sample=rownames(geno_pf_population),geno_pf_population)
pop_pf_population <- read.csv("Data/Pf_Sample.csv")
geno_pf_population <- geno_pf_population %>% 
  left_join(pop_pf_population,by="Sample")
vars_pf_population <- gsub(vars_pf_population,pattern="-",replacement=".")
#Pf Country
Pf_country_loci <- read.csv("Data/Pf_STR_Top_ShinyPlot.csv")
Pf_country_loci <- good_pf %>% 
  filter(STRid %in% as.character(Pf_country_loci$STRid))
Pf_country_loci <- Pf_country_loci %>% 
  pivot_longer(cols = Bangladesh.Colombia:VietNam.Thailand) %>% 
  filter(is.na(value)==FALSE) %>% 
  arrange(value)
vars_pf_country <- setdiff(unique(Pf_country_loci$value), "STRLoci")
vars_pf_country[] <- lapply(vars_pf_country, gsub,pattern="Papua New Guinea",replacement="PapuaNewGuinea")
vars_pf_country[] <- lapply(vars_pf_country, gsub,pattern="Congo DR",replacement="CongoDR")
vars_pf_country[] <- lapply(vars_pf_country, gsub,pattern="Ivory Coast",replacement="IvoryCoast")
vars_pf_country[] <- lapply(vars_pf_country, gsub,pattern="Viet Nam",replacement="VietNam")
geno_pf_country <- geno_origin_Pf
geno_pf_country <- geno_pf_country %>% 
  filter(STRid %in% as.character(Pf_country_loci$STRid))
STR_Position_pf_country <- Pf_country_loci %>% 
  dplyr::select(STRid,value)
geno_pf_country <- STR_Position_pf_country %>% 
  left_join(geno_pf_country,by="STRid")
geno_pf_country <- geno_pf_country[,-1]
geno_pf_country$value <- as.character(geno_pf_country$value)
geno_pf_country$value[] <- lapply(geno_pf_country$value, gsub,pattern="Papua New Guinea",replacement="PapuaNewGuinea")
geno_pf_country$value[] <- lapply(geno_pf_country$value[], gsub,pattern="Congo DR",replacement="CongoDR")
geno_pf_country$value[] <- lapply(geno_pf_country$value[], gsub,pattern="Ivory Coast",replacement="IvoryCoast")
geno_pf_country$value[] <- lapply(geno_pf_country$value[], gsub,pattern="Viet Nam",replacement="VietNam")
geno_pf_country <- as.data.frame(geno_pf_country)
rownames(geno_pf_country) <- geno_pf_country$value
geno_pf_country <- geno_pf_country[,-1]
geno_pf_country <- t(geno_pf_country)
geno_pf_country <- data.frame(Sample=rownames(geno_pf_country),geno_pf_country)
pop_pf_country <- read.csv("Data/Pf_Sample.csv")
geno_pf_country <- geno_pf_country %>% 
  left_join(pop_pf_country,by="Sample")
vars_pf_country <- gsub(vars_pf_country,pattern="-",replacement=".")
#Pv country
geno_pv <- read.csv("Data/Pv_Genotype.csv")
geno_origin_Pv <- geno_pv
Pv_country_loci <- read.csv("Data/Pv_STR_Top_ShinyPlot.csv")
Pv_country_loci <- good_pv %>% 
  filter(STRid %in% as.character(Pv_country_loci$STRid))
Pv_country_loci <- Pv_country_loci %>% 
  pivot_longer(cols = Cambodia.Colombia:Thailand.Peru) %>% 
  filter(is.na(value)==FALSE) %>% 
  arrange(value)
vars_pv_country <- setdiff(unique(Pv_country_loci$value), "STRLoci")
geno_pv_country <- geno_origin_Pv
geno_pv_country <- geno_pv_country %>% 
  filter(STRid %in% as.character(Pv_country_loci$STRid))
STR_Position_pv_country <- Pv_country_loci %>% 
  dplyr::select(STRid,value)
geno_pv_country <- STR_Position_pv_country %>% 
  left_join(geno_pv_country,by="STRid")
geno_pv_country <- geno_pv_country[,-1]
geno_pv_country$value <- as.character(geno_pv_country$value)
geno_pv_country <- as.data.frame(geno_pv_country)
rownames(geno_pv_country) <- geno_pv_country$value
geno_pv_country <- geno_pv_country[,-1]
geno_pv_country <- t(geno_pv_country)
geno_pv_country <- data.frame(Sample=rownames(geno_pv_country),geno_pv_country)
pop_pv_country <- read.csv("Data/Pv_Sample.csv")
colnames(pop_pv_country)[1] <- "Sample"
pop_pv_country <- pop_pv_country[,c(1,3)]
pop_pv_country$Country <- as.character(pop_pv_country$Country)
geno_pv_country <- geno_pv_country %>% 
  left_join(pop_pv_country,by="Sample")
vars_pv_country <- gsub(vars_pv_country,pattern="-",replacement=".")
#Pf PCA plot data
Pf_PCA <- read.csv("Data/PCA_Pf.csv")
Pf_PCA$Population <- factor(Pf_PCA$Population,levels = c("CAF","EAF","WAF","SAM","OCE","SAS","ESEA","WSEA"))
#Pv PCA plot data
Pv_PCA <- read.csv("Data/PCA_Pv.csv")

#Define ui
ui <- navbarPage(inverse = TRUE,
                 "PlasmoSTR", 
                 theme = shinytheme("cerulean"),
                 id = "panels_main",
                 tabPanel(icon("home"), 
                          fluidPage(
                            h1("PlasmoSTR: A population-wide database of short tandem repeat variation in ", em("Plasmodium falciparum"), "and ", em("Plasmodium vivax")),
                            br(),
                            h3("About PlasmoSTR"),
                            br(),
                            p(strong("PlasmoSTR contains population and genome-wide information about genetic variation and other characteristics of short tandem repeats (STRs) in ",em("Plasmodium falciparum"), "and ", em("Plasmodium vivax"), " .",style = "font-size:20px;")),
                            br(),
                            p(style="text-align: justify;", "We genotype STRs using bioinformatic tool ", a(strong("HipSTR"),style = "color:#337ab7", href = "https://hipstr-tool.github.io/HipSTR/"), "in more than", 
                              a(strong("3,000", em("P. falciparum ")),"(Data generated by the MalariaGEN Plasmodium falciparum Community Project)",style = "color:#337ab7", href = "https://www.malariagen.net/resource/26"), "and ", 
                              a(strong("174", em("P. vivax")),"(Data generated by the Plasmodium vivax Genome Variation project)", style = "color:#337ab7", href = "https://www.malariagen.net/projects/p-vivax-genome-variation"), "and ", 
                              a("(Data from Hupalo et al 2016).", style = "color:#337ab7", href = "https://www.nature.com/articles/ng.3588"),
                              "We develop a multivariable logistic regression model for the measurement and prediction of the quality of STRs. A set of high-quality STR loci (6,768 from ",
                              em("P. falciparum"), "and 3,496 from ", em("P. vivax"), ") were selected.",style = "font-size: 13pt"), 
                            p(actionLink("link_to_tabpanel_PCA",strong("1. PCA",style = "color:#337ab7")),style = "font-size: 14pt"),
                            p(style="text-align: justify;","Population structure and principal component analysis (PCA) of Plasmodium parasite populations using genome-wide high-quality STRs.",style = "font-size: 13pt"),
                            p(actionLink("link_to_tabpanel_STRdatasets",strong("2. STR datasets",style = "color:#337ab7")),style = "font-size: 14pt"),
                            
                            p(style="text-align: justify;","Genome-wide high-quality STRs.",style = "font-size: 13pt"),
                            p(style="text-align: justify;","Choose the dataset (e.g., Pf or Pv) and the STR genomic location (e.g., All, Coding, Promoter, Intergenic, Intron, or Other) to access the high-quality STR loci information.",style = "font-size: 13pt"),
                            p(style="text-align: justify;","A genomic region of STR locus using the syntax chrom:start-end (e.g., 7:1322349-1322361). PlasmoSTR is currently based on Pf3D7 (v3 PlasmoDB-41) and PvP01 (PlasmoDB release 41) coordinates.",style = "font-size: 13pt"),
                          
                            p(actionLink("link_to_tabpanel_TopSTR",strong("3. TopSTR",style = "color:#337ab7")),style = "font-size: 14pt"),
                            p("The top ten most highly differentiated STRs based on ", em("Jost's D "), "proposed by ", 
                              a(strong("Jost (2008)."),style = "color:#337ab7", href = "https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2008.03887.x"),
                              
                              style= "font-size: 13pt"),
                            p(style="text-align: justify;","Choose the", strong(em("P. falciparum")), "population pair (e.g., CAF.EAF),", strong(em("P. falciparum")), "country pair (e.g., Bangladesh.Colombia), or ", strong(em("P. vivax")), "country pair (e.g., Cambodia.Colombia) to access the top ten most highly differentiated STRs between the population or country pairs.",style = "font-size: 13pt"),
                            p(style="text-align: justify;","Click each row of the datasets can show the samples corresponding motif copies and allele frequency plots. The figures is shown in ", actionLink("link_to_tabpanel_TopPfs",strong("PfPopulation_TopSTRPlots; ",style = "color:#337ab7")), 
                              actionLink("link_to_tabpanel_TopPfsc",strong("PfCountry_TopSTRPlots; ",style = "color:#337ab7")), 
                              actionLink("link_to_tabpanel_TopPvs",strong("PvCountry_TopSTRPlots", ".", style = "color:#337ab7")),
                              style = "font-size: 13pt"),
                            p(strong(em("P. falciparum")), strong("population:"), " Central Africa (CAF), East Africa (EAF), West Africa (WAF), South America (SAM), Oceania (OCE), South Asia (SAS), the eastern part of Southeast Asia (ESEA), the western part of Southeast Asia (WSEA).",style = "font-size: 13pt"),

                            p(strong("4. TopSTRMotifPlots",style = "color:#337ab7"),style = "font-size: 14pt"),
                            p(style="text-align: justify;","Samples corresponding motif copies and allele frequency plots. Plot samples corresponding motif copies and allele frequency distribution of each highly differentiated STR locus (based on ",em("Jost's D"), "value) over different population or country pair.",style = "font-size: 13pt"),
                            p(strong("Access by clicking on each row of each dataset in"),actionLink("link_to_tabpanel_TopSTRSTR",strong("TopSTR",style = "color:#337ab7")), ".",
                              style = "font-size: 13pt"),
                            br(),
                            br(),
                            br(),
                            br(),
                            br(),
                            div(p(strong("PlasmoSTR was developed by Jiru Han with input from Jacob E. Munro, Melanie Bahlo and other "), a("Bahlo Lab members", href = "https://www.wehi.edu.au/people/melanie-bahlo/372/melanie-bahlo-lab-team")), 
                                p(strong("It has been made with Rstudio and Shiny and you may find the code on "), a("github", href = "https://github.com/bahlolab/PlasmoSTR/"), "."),
                                style="text-align: right;")
                          )
                 ),
                 # Sidebar layout with input and output definitions ----
                 #PCA Page
                 tabPanel("PCA",
                          fluidPage(
                            h1("Principal component analysis (PCA)"),
                            br(),
                            br(),
                            mainPanel(
                              p(strong("PCA plot of 6,768 genome-wide high-quality STRs (Pf MalariaGEN 3,047 samples)"),style = "font-size: 17pt"),
                              br(),
                              br(),
                              plotlyOutput('Pf_PCA'),
                              br(),
                              p("Colors represent populations, each sample is represented by a single point."),
                              p("Central Africa (CAF), East Africa (EAF), West Africa (WAF), South America (SAM), Oceania (OCE), South Asia (SAS), the eastern part of Southeast Asia (ESEA), the western part of Southeast Asia (WSEA)"),
                              hr(),
                              br(),
                              br(),
                              p(strong("PCA plot of 3,496 genome-wide high-quality STRs (Pv 175 samples)"),style = "font-size: 17pt"),
                              br(),
                              br(),
                              plotlyOutput('Pv_PCA'),
                              br(),
                              p("Colors represent countries, each sample is represented by a single point.")
                            )
                            
                          )
                 ),
                 
                 #Pf and Pv dataset
                 tabPanel("STR datasets",
                          fluidPage(
                            h1("Genome-wide high-quality STRs"),
                            sidebarLayout(position = "left",
                                          sidebarPanel(
                                            wellPanel(style = "background: white",
                                                      h4(em("Plasmodium falciparum"),": 6,768"),
                                                      h4(em("Plasmodium vivax"),": 3,496"),
                                                      hr(),
                                                      selectInput(inputId = "dataset",
                                                                  label = "Please choose a dataset:",
                                                                  choices = c("","Pf", "Pv")),
                                                      
                                                      selectInput(inputId = "STRcategory",
                                                                  label = "Please choose a genomic location:",
                                                                  choices = c("","All","Coding", "Promoter","Intergenic","Intron", "Other")),
                                                      
                                                      actionButton("View","Generate Table"))
                                          ),
                                          
                                          mainPanel(
                                            br(),
                                            DT::dataTableOutput("view1")
                                          )
                            )
                          )
                 ),
                 #The top ten highly differentiated STRs       
                 tabPanel("TopSTR",
                          h1("Highly differentiated STRs of each pair population or country"),
                          navlistPanel(
                            "Plasmodium falciparum",  
                            tabPanel("1. Pf Population Pair",
                                     fluidRow(column(12,
                                                     wellPanel(#Pf population
                                                       style = "background: white",
                                                       selectInput(inputId = "PfPopulationPair",
                                                                   label = "Please choose a P. falciparum population pair:",
                                                                   choices = c("",compare_Pf_Population)),
                                                       
                                                       actionButton("View2","Generate Table"))),
                                              fluidRow(column(12,
                                                              DT::dataTableOutput("view2"))))), 
                            
                            tabPanel("2. Pf Country Pair",
                                     fluidRow(column(12,
                                                     wellPanel(#Pf population
                                                       style = "background: white",
                                                       #Pf country
                                                       selectInput(inputId = "PfCountrypair",
                                                                   label = "Please choose a P. falciparum country pair:",
                                                                   choices = c("",compare_Pf_Country)),
                                                       actionButton("View3","Generate Table"))),
                                              fluidRow(column(12,
                                                              DT::dataTableOutput("view3"))))), 
                            "Plasmodium vivax",
                            tabPanel("3. Pv Country Pair",
                                     fluidRow(column(12,
                                                     wellPanel(
                                                       style = "background: white",
                                                       #Pv country
                                                       selectInput(inputId = "PvCountrypair",
                                                                   label = "Please choose a P. vivax country pair:",
                                                                   choices = c("",compare_Pv_Country)),
                                                       actionButton("View4","Generate Table"))),
                                              fluidRow(column(12,
                                                              DT::dataTableOutput("view4")))))
                            
                          )
                 ),
                 
                 #Plots Pf population   
                 tabPanel("PfPopulation_TopSTRPlots",
                            tabPanel("Pf Population Pair",
                                     fluidRow(column(6,
                                                     wellPanel(
                                                       style = "background: white",
                                                       selectInput(inputId = 'xcol', 
                                                                   label ='Please choose the plot color for P. falciparum (Population or Country)', 
                                                                   choices = c("Population","Country"))))),
                                     fluidRow(
                                       column(12,
                                              h3(textOutput('textPf')),
                                              splitLayout(style = "border: 1px solid silver:", cellWidths = c("40%","60%"),
                                                          plotlyOutput('Pfx1', height = 500),
                                                          plotlyOutput('Pfx2', height = 500)
                                              ),
                                              h3(textOutput('textPfA')),
                                              splitLayout(style = "border: 1px solid silver:", cellWidths = c("40%","60%"),
                                                          plotlyOutput('Pfx3', height = 500),
                                                          plotlyOutput('Pfx4', height = 500)
                                              )
                                     )))
                          
                 ),
                 #Plots Pf country
                 tabPanel("PfCountry_TopSTRPlots",
                          tabPanel("Pf Country Pair",
                                   fluidRow(column(6,
                                                   wellPanel(
                                                     style = "background: white",
                                                     selectInput(inputId = 'xcol1', 
                                                                 label ='Please choose the plot color for P. falciparum (Population or Country)', 
                                                                 choices = c("Population","Country"))))),
                                   fluidRow(
                                     column(12,
                                            h3(textOutput('textPfc')),
                                            splitLayout(style = "border: 1px solid silver:", cellWidths = c("40%","60%"),
                                                        plotlyOutput('Pfcx1', height = 500),
                                                        plotlyOutput('Pfcx2', height = 500)
                                            ),
                                            h3(textOutput('textPfAc')),
                                            splitLayout(style = "border: 1px solid silver:", cellWidths = c("40%","60%"),
                                                        plotlyOutput('Pfcx3', height = 500),
                                                        plotlyOutput('Pfcx4', height = 500)
                                            )
                                            
                                     )
                                   ))
                          
                 ),
                 #Plots Pv country
                 tabPanel("PvCountry_TopSTRPlots",
                          tabPanel("Pv Country Pair",
                                   fluidRow(
                                     column(12,
                                            h3(textOutput('textPv')),
                                            splitLayout(style = "border: 1px solid silver:", cellWidths = c("40%","60%"),
                                                        plotlyOutput('Pvx1', height = 500),
                                                        plotlyOutput('Pvx2', height = 500)
                                            ),
                                            h3(textOutput('textPvc')),
                                            splitLayout(style = "border: 1px solid silver:", cellWidths = c("40%","60%"),
                                                        plotlyOutput('Pvx3', height = 500),
                                                        plotlyOutput('Pvx4', height = 500)
                                            )
                                            
                                     )
                                   ))
                          
                 )
                
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output,session) {
  
  #Link internal page
  observeEvent(input$link_to_tabpanel_PCA,{
    newvalue <- "PCA"
    updateNavbarPage(session, "panels_main", newvalue)
  }) 
  
  observeEvent(input$link_to_tabpanel_STRdatasets,{
    newvalue <- "STR datasets"
    updateNavbarPage(session, "panels_main", newvalue)
  }) 
  
  
  observeEvent(input$link_to_tabpanel_TopSTR,{
    newvalue <- "TopSTR"
    updateNavbarPage(session, "panels_main", newvalue)
  }) 
  
  observeEvent(input$link_to_tabpanel_TopSTRSTR,{
    newvalue <- "TopSTR"
    updateNavbarPage(session, "panels_main", newvalue)
  }) 
  
  
  
  observeEvent(input$link_to_tabpanel_TopPfs,{
    newvalue <- "PfPopulation_TopSTRPlots"
    updateNavbarPage(session, "panels_main", newvalue)
  }) 
  

  observeEvent(input$link_to_tabpanel_TopPfsc,{
    newvalue <- "PfCountry_TopSTRPlots"
    updateNavbarPage(session, "panels_main", newvalue)
  }) 
  
  
  observeEvent(input$link_to_tabpanel_TopPvs,{
    newvalue <- "PvCountry_TopSTRPlots"
    updateNavbarPage(session, "panels_main", newvalue)
  }) 
 
  
  # Return the requested dataset ----
  Pf = good_pf %>% 
    dplyr::select(STRLocus,Motif,copies,Product.Description,Category,Heterozygosity) %>% 
    dplyr::rename("STR Locus (Pf3D7)" = "STRLocus") %>% 
    dplyr::rename("# copies (Pf3D7)" = "copies")
  Pv = good_pv %>% 
    dplyr::select(STRLocus,Motif,copies,Product.Description,Category,Heterozygosity) %>% 
    dplyr::rename("STR Locus (PvP01)" = "STRLocus") %>% 
    dplyr::rename("# copies (PvP01)" = "copies")

  db <- reactive({
    require(input$View)
    if (input$dataset == "Pf" & input$STRcategory == "All") {
      Pf
    } else if (input$dataset == "Pf" & input$STRcategory == "Coding") {
      Pf %>% 
        filter(Category=="coding")
    } else if (input$dataset == "Pf" & input$STRcategory == "Promoter") {
      Pf %>% 
        filter(Category=="promoter")
    } else if (input$dataset == "Pf" & input$STRcategory == "Intergenic") {
      Pf %>% 
        filter(Category=="intergenic")
    } else if (input$dataset == "Pf" & input$STRcategory == "Intron") {
      Pf %>% 
        filter(Category=="intron")
    } else if (input$dataset == "Pf" & input$STRcategory == "Other") {
      Pf %>% 
        filter(!(Category %in% c("coding", "promoter","intergenic","intron")))
    } else if (input$dataset == "Pv" & input$STRcategory == "All") {
      Pv
    } else if (input$dataset == "Pv" & input$STRcategory == "Coding") {
      Pv %>% 
        filter(Category=="coding")
    } else if (input$dataset == "Pv" & input$STRcategory == "Promoter") {
      Pv %>% 
        filter(Category=="promoter")
    } else if (input$dataset == "Pv" & input$STRcategory == "Intergenic") {
      Pv %>% 
        filter(Category=="intergenic")
    } else if (input$dataset == "Pv" & input$STRcategory == "Intron") {
      Pv %>% 
        filter(Category=="intron")
    } else if (input$dataset == "Pv" & input$STRcategory == "Other") {
      Pv %>% 
        filter(!(Category %in% c("coding", "promoter","intergenic","intron")))
    } 
  })
  
  data <- eventReactive(input$View, {
    db()
  })
  
  output$view1 <- DT::renderDataTable(
    DT::datatable(data(),caption = 'Table: This is genome-wide high-quality STR loci table.',
                  filter = 'top', 
                  extensions = 'Buttons',
                  options = list(dom = 'Blfrtip',
                                 buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All")))
                  )
  )
  
  # Downloadable csv of selected dataset ----
  #Pf population pair
  datasetInput2 <- eventReactive(input$View2,{
    good_pf %>% 
      dplyr::filter(is.na(good_pf[,which(colnames(good_pf)==input$PfPopulationPair)])==FALSE) %>% 
      dplyr::select(STRLocus,Motif,copies,Product.Description,Category,Heterozygosity,input$PfPopulationPair) %>% 
      dplyr::rename("STR Locus (Pf3D7)" = "STRLocus") %>% 
      dplyr::rename("# copies (Pf3D7)" = "copies") %>% 
      arrange(as.numeric(str_extract(.[,7], "\\d+$")))
  })
  
  #Pf country pair
  datasetInput3 <- eventReactive(input$View3,{
    good_pf %>% 
      dplyr::filter(is.na(good_pf[,which(colnames(good_pf)==input$PfCountrypair)])==FALSE) %>% 
      dplyr::select(STRLocus,Motif,copies,Product.Description,Category,Heterozygosity,input$PfCountrypair) %>% 
      dplyr::rename("STR Locus (Pf3D7)" = "STRLocus") %>% 
      dplyr::rename("# copies (Pf3D7)" = "copies") %>% 
      arrange(as.numeric(str_extract(.[,7], "\\d+$")))
  })
  
  #Pv country pair
  datasetInput4 <- eventReactive(input$View4,{
    good_pv %>% 
      dplyr::filter(is.na(good_pv[,which(colnames(good_pv)==input$PvCountrypair)])==FALSE) %>% 
      dplyr::select(STRLocus,Motif,copies,Product.Description,Category,Heterozygosity, input$PvCountrypair) %>% 
      dplyr::rename("STR Locus (PvP01)" = "STRLocus") %>% 
      dplyr::rename("# copies (PvP01)" = "copies") %>% 
      arrange(as.numeric(str_extract(.[,7], "\\d+$")))
  })
  
  #Pf_PCA
  fig1_pf <- Pf_PCA %>% plot_ly(
    type = 'scatter', 
    x = ~PC1, 
    y = ~PC2, 
    color = ~Population, 
    legend = ~Population,
    legendgroup = ~Population,
    text =~Sample
  )
  fig1_pf <- fig1_pf %>% layout(
    xaxis = list(
      showgrid = F
    ),
    yaxis = list(
      showgrid = F
    )
  )
  
  fig2_pf <- Pf_PCA %>% plot_ly(
    type = 'scatter', 
    x = ~PC1, 
    y = ~PC3, 
    color = ~Population,
    legendgroup = ~Population,
    text =~Sample,
    showlegend = F
  )
  
  fig2_pf <- fig2_pf %>% layout(
    xaxis = list(
      showgrid = F
    ),
    yaxis = list(
      showgrid = F
    )
  )
  
  fig_pf <- subplot(fig1_pf, fig2_pf, nrows = 2, shareX = T)
  
  fig_pf
  output$Pf_PCA <- renderPlotly({
    print(
      fig_pf
    )
  })
  
  #Pv_PCA
  fig1_pv <- Pv_PCA %>% plot_ly(
    type = 'scatter', 
    x = ~PC1, 
    y = ~PC2, 
    color = ~Country, 
    legend = ~Country,
    legendgroup = ~Country,
    text =~Sample
  )
  fig1_pv <- fig1_pv %>% layout(
    xaxis = list(
      showgrid = F
    ),
    yaxis = list(
      showgrid = F
    )
  )
  
  fig2_pv <- Pv_PCA %>% plot_ly(
    type = 'scatter', 
    x = ~PC1, 
    y = ~PC3, 
    color = ~Country,
    legendgroup = ~Country,
    text =~Sample,
    showlegend = F
  )
  
  fig2_pv <- fig2_pv %>% layout(
    xaxis = list(
      showgrid = F
    ),
    yaxis = list(
      showgrid = F
    )
  )
  
  fig_pv <- subplot(fig1_pv, fig2_pv, nrows = 2, shareX = T)
  fig_pv
  output$Pv_PCA <- renderPlotly({
    print(
      fig_pv
    )
  })
  
  
  #The top ten highly differentiated STRs P. falciparum population pair
  
  output$view2 <- DT::renderDataTable(
    DT::datatable(datasetInput2(),caption = 'Table: This is top 10 most highly differentitated STR loci table (P. falciparum population pair). Click each row can show the samples corresponding motif copies plots.',
                  filter = 'top', 
                  #options = list(dom = 't'),
                  selection = list(mode="single", target="row"),
                  options = list(paging=FALSE,
                                 dom = 't',
                                 searching=FALSE,
                                 filtering=FALSE,
                                 ordering=FALSE),
                  callback=JS(
                    'table.on("click.dt", "tr", function() {

    tabs = $(".tabbable .nav.nav-tabs li a");
    var data=table.row(this).data();
    $(tabs[1]).click();
    table.row(this).deselect();
    })'                     
                    
                  ))
  )
  
  # A function factory for getting integer y-axis values.
  integer_breaks <- function(n = 5, ...) {
    fxn <- function(x) {
      breaks <- floor(pretty(x, n, ...))
      names(breaks) <- attr(breaks, "labels")
      breaks
    }
    return(fxn)
  }
  
  #Plot pf population
  observeEvent(input$view2_rows_selected,{
    
    output$textPf <- renderText({"STR copies of samples in highly differentiated STR locus"})
    output$textPfA <- renderText({"Allele frequency of samples in highly differentiated STR locus"})
    
    selectedCellsPf <- renderText(as.character(datasetInput2()[input$view2_rows_selected,7]))
    
    
    selectedCellsPfReference <- renderText(datasetInput2()[input$view2_rows_selected,3])
    
    output$Pfx1 = renderPlotly({
      plot1pf <- geno_pf_population %>% 
        filter(is.na(geno_pf_population[,which(colnames(geno_pf_population)==selectedCellsPf())])==FALSE)  %>% 
        filter(Population %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCellsPf()),"\\."))) %>% 
        filter(is.na(selectedCellsPf())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCellsPf()) %>% 
        ggplot(aes_string(x=input$xcol,y="Motif",col=input$xcol,
                          label="Sample",label2=input$xcol,label3="Motif"))+
        geom_jitter(alpha=0.5,height = 0.01)+
        theme_bw()+
        scale_y_continuous(breaks = integer_breaks())+
        labs(
          title = "Pair population (P. falciparum)",
             x="Pairwise population",
             y="STR copies")+
        geom_hline(yintercept = as.numeric(selectedCellsPfReference()),linetype="dashed",size=0.5)+
        annotate("text",0.8,as.numeric(selectedCellsPfReference()),vjust = -1,label = "Pf3D7")
      
      
      ggplotly(plot1pf)
    })

    output$Pfx2 = renderPlotly({
      plot2pf <- geno_pf_population %>% 
        filter(is.na(geno_pf_population[,which(colnames(geno_pf_population)==selectedCellsPf())])==FALSE)  %>% 
        filter(!(Population %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCellsPf()),"\\.")))) %>% 
        filter(is.na(selectedCellsPf())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCellsPf()) %>% 
        ggplot(aes_string(x=input$xcol,y="Motif",col=input$xcol,
                          label="Sample",label2=input$xcol,label3="Motif"))+
        geom_jitter(alpha=0.5,height = 0.01)+
        theme_bw()+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Other populations exclude pair population (P. falciparum)",
             x="Other populations",
             y="STR copies")+
        geom_hline(yintercept = as.numeric(selectedCellsPfReference()),linetype="dashed",size=0.5)+
        annotate("text",0.8,as.numeric(selectedCellsPfReference()),vjust = -1,label = "Pf3D7")
      ggplotly(plot2pf)
    })
    
    output$Pfx3 = renderPlotly({
      plot3pf <- geno_pf_population %>% 
        filter(is.na(geno_pf_population[,which(colnames(geno_pf_population)==selectedCellsPf())])==FALSE) %>% 
        filter(Population %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCellsPf()),"\\."))) %>% 
        filter(is.na(selectedCellsPf())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCellsPf()) %>% 
        group_by(.dots = c("Motif",input$xcol)) %>% 
        summarise(Count=n()) %>% 
        ungroup() %>% 
        ggplot(aes_string(x="Motif",y="Count",fill=input$xcol))+
        geom_col(alpha=0.8) +
        theme_bw()+
        scale_x_continuous(breaks = integer_breaks())+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Pair population (P. falciparum)",
             x="STR copies",
             y="Count")
      ggplotly(plot3pf)
    })
    
    output$Pfx4 = renderPlotly({
      plot4pf <- geno_pf_population %>% 
        filter(is.na(geno_pf_population[,which(colnames(geno_pf_population)==selectedCellsPf())])==FALSE) %>% 
        filter(!(Population %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCellsPf()),"\\.")))) %>% 
        filter(is.na(selectedCellsPf())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCellsPf()) %>% 
        group_by(.dots = c("Motif",input$xcol)) %>% 
        summarise(Count=n()) %>% 
        ungroup() %>% 
        ggplot(aes_string(x="Motif",y="Count",fill=input$xcol))+
        geom_col(alpha=0.8) +
        theme_bw()+
        scale_x_continuous(breaks = integer_breaks())+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Other populations exclude pair population (P. falciparum)",
             x="STR copies",
             y="Count")
      ggplotly(plot4pf)
    })
    
    updateNavbarPage(session, "panels_main", selected = "PfPopulation_TopSTRPlots")
  }) 
  
  #The top ten highly differentiated STRs P. falciparum country pair
  output$view3 <- DT::renderDataTable(
    DT::datatable(datasetInput3(),caption = 'Table: This is top 10 most highly differentitated STR loci table (P. falciparum country pair). Click each row can show the samples corresponding motif copies plots.',
                  filter = 'top', 
                  #options = list(dom = 't'),
                  selection = list(mode="single", target="row"),
                  options = list(paging=FALSE,
                                 dom = 't',
                                 searching=FALSE,
                                 filtering=FALSE,
                                 ordering=FALSE),
                  callback=JS(
                    'table.on("click.dt", "tr", function() {

    tabs = $(".tabbable .nav.nav-tabs li a");
    var data=table.row(this).data();
    $(tabs[1]).click();
    table.row(this).deselect();
    })'                     
                    
                  ))
  )
  
  #plot pf country
  observeEvent(input$view3_rows_selected,{
    
     output$textPfc <- renderText({"STR copies of samples in highly differentiated STR locus"})
     output$textPfAc <- renderText({"Allele frequency of samples in highly differentiated STR locus"})
    
    selectedCellsPfc <- renderText(as.character(datasetInput3()[input$view3_rows_selected,7]))
    selectedCellsPfcReference <- renderText(datasetInput3()[input$view3_rows_selected,3])
    
    output$Pfcx1 = renderPlotly({
      plot1pfc <- geno_pf_country %>% 
        filter(is.na(geno_pf_country[,which(colnames(geno_pf_country)==selectedCellsPfc())])==FALSE) %>% 
        filter(Country %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCellsPfc()),"\\."))) %>% 
        filter(is.na(selectedCellsPfc())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCellsPfc()) %>% 
        ggplot(aes_string(x=input$xcol1,y="Motif",col=input$xcol1,
                          label="Sample",label2=input$xcol1,label3="Motif"))+
        geom_jitter(alpha=0.5,height = 0.01)+
        theme_bw()+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Pair country (P. falciparum)",
             x="Pairwise country",
             y="STR copies")+
        geom_hline(yintercept = as.numeric(selectedCellsPfcReference()),linetype="dashed",size=0.5)+
        annotate("text",0.8,as.numeric(selectedCellsPfcReference()),vjust = -1,label = "Pf3D7")
      
      ggplotly(plot1pfc)
    })
    
    
    output$Pfcx2 = renderPlotly({
      plot2pfc <- geno_pf_country %>% 
        filter(is.na(geno_pf_country[,which(colnames(geno_pf_country)==selectedCellsPfc())])==FALSE) %>% 
        filter(!(Country %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCellsPfc()),"\\.")))) %>% 
        filter(is.na(selectedCellsPfc())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCellsPfc()) %>% 
        ggplot(aes_string(x=input$xcol1,y="Motif",col=input$xcol1,
                          label="Sample",label2=input$xcol1,label3="Motif"))+
        geom_jitter(alpha=0.5,height = 0.01)+
        theme_bw()+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Other countries exclude pair country (P. falciparum)",
             x="Other countries",
             y="STR copies")+
        geom_hline(yintercept = as.numeric(selectedCellsPfcReference()),linetype="dashed",size=0.5)+
        annotate("text",0.8,as.numeric(selectedCellsPfcReference()),vjust = -1,label = "Pf3D7")
      
      ggplotly(plot2pfc)
    })
    
    output$Pfcx3 = renderPlotly({
      plot3pfc <- geno_pf_country %>% 
        filter(is.na(geno_pf_country[,which(colnames(geno_pf_country)==selectedCellsPfc())])==FALSE) %>% 
        filter(Country %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCellsPfc()),"\\."))) %>% 
        filter(is.na(selectedCellsPfc())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCellsPfc()) %>% 
        group_by(.dots = c("Motif",input$xcol1)) %>% 
        summarise(Count=n()) %>% 
        ungroup() %>% 
        ggplot(aes_string(x="Motif",y="Count",fill=input$xcol1))+
        geom_col(alpha=0.8) +
        theme_bw()+
        scale_x_continuous(breaks = integer_breaks())+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Pair country (P. falciparum)",
             x="STR copies",
             y="Count")
      
      ggplotly(plot3pfc)
    })
    
    output$Pfcx4 = renderPlotly({
      plot4pfc <- geno_pf_country %>% 
        filter(is.na(geno_pf_country[,which(colnames(geno_pf_country)==selectedCellsPfc())])==FALSE) %>% 
        filter(!(Country %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCellsPfc()),"\\.")))) %>% 
        filter(is.na(selectedCellsPfc())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCellsPfc()) %>% 
        group_by(.dots = c("Motif",input$xcol1)) %>% 
        summarise(Count=n()) %>% 
        ungroup() %>% 
        ggplot(aes_string(x="Motif",y="Count",fill=input$xcol1))+
        geom_col(alpha=0.8) +
        theme_bw()+
        scale_x_continuous(breaks = integer_breaks())+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Other countries exclude pair country (P. falciparum)",
             x="STR copies",
             y="Count")
      ggplotly(plot4pfc)
    })
    
    updateNavbarPage(session, "panels_main", selected = "PfCountry_TopSTRPlots")
  }) 
  
  
  #The top ten highly differentiated STRs P. vivax country pair
  
  output$view4 <- DT::renderDataTable(
    DT::datatable(datasetInput4(),caption = 'Table: This is top 10 most highly differentitated STR loci table (P. vivax country pair). Click each row can show the samples corresponding motif copies plots.',
                  filter = 'top', 
                  #options = list(dom = 't'),
                  selection = list(mode="single", target="row"),
                  options = list(paging=FALSE,
                                 dom = 't',
                                 searching=FALSE,
                                 filtering=FALSE,
                                 ordering=FALSE),
                  callback=JS(
                    'table.on("click.dt", "tr", function() {

    tabs = $(".tabbable .nav.nav-tabs li a");
    var data=table.row(this).data();
    $(tabs[1]).click();
    table.row(this).deselect();
    })'                     
                    
                  ))
  )
  
  #plot pv country
  observeEvent(input$view4_rows_selected,{
    
      output$textPv <- renderText({"STR copies of samples in highly differentiated STR locus"})
      output$textPvc <- renderText({"Allele frequency of samples in highly differentiated STR locus"})
    
    selectedCells <- renderText(as.character(datasetInput4()[input$view4_rows_selected,7]))
    selectedCellsPvReference <- renderText(datasetInput4()[input$view4_rows_selected,3])
    
    
    output$Pvx1 = renderPlotly({
      plot1 <- geno_pv_country %>% 
        filter(is.na(geno_pv_country[,which(colnames(geno_pv_country)==selectedCells())])==FALSE) %>% 
        filter(Country %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCells()),"\\."))) %>% 
        filter(is.na(selectedCells())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCells()) %>% 
        ggplot(aes_string(x="Country",y="Motif", col="Country",
                          label="Sample",label2="Country",label3="Motif"))+
        geom_jitter(alpha=0.5,height = 0.01)+
        theme_bw()+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Pair country (P. vivax)",
             x="Pairwise country",
             y="STR copies")+
        geom_hline(yintercept = as.numeric(selectedCellsPvReference()),linetype="dashed",size=0.5)+
        annotate("text",0.8,as.numeric(selectedCellsPvReference()),vjust = -1,label = "PvP01")
      
      ggplotly(plot1)
    })
    
    
    output$Pvx2 = renderPlotly({
      plot2 <- geno_pv_country %>% 
        filter(is.na(geno_pv_country[,which(colnames(geno_pv_country)==selectedCells())])==FALSE) %>% 
        filter(!(Country %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCells()),"\\.")))) %>% 
        filter(is.na(selectedCells())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCells()) %>% 
        ggplot(aes_string(x="Country",y="Motif", col="Country",
                          label="Sample",label2="Country",label3="Motif"))+
        geom_jitter(alpha=0.5,height = 0.01)+
        theme_bw()+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Other countries exclude pair country (P. vivax)",
             x="Other countries",
             y="STR copies")+
        geom_hline(yintercept = as.numeric(selectedCellsPvReference()),linetype="dashed",size=0.5)+
        annotate("text",0.8,as.numeric(selectedCellsPvReference()),vjust = -1,label = "PvP01")
      
      ggplotly(plot2)
    })
    
    output$Pvx3 = renderPlotly({
      plot3 <- geno_pv_country %>% 
        filter(is.na(geno_pv_country[,which(colnames(geno_pv_country)==selectedCells())])==FALSE) %>% 
        filter(Country %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCells()),"\\."))) %>% 
        filter(is.na(selectedCells())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCells()) %>% 
        group_by(Motif,Country) %>% 
        summarise(Count=n()) %>% 
        ungroup() %>% 
        ggplot(aes(x=Motif,y=Count,fill=Country))+
        geom_col(alpha=0.8) +
        theme_bw()+
        scale_x_continuous(breaks = integer_breaks())+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Pair country (P. vivax)",
             x="STR copies",
             y="Count")
      ggplotly(plot3)
    })
    
    output$Pvx4 = renderPlotly({
      plot4 <- geno_pv_country %>% 
        filter(is.na(geno_pv_country[,which(colnames(geno_pv_country)==selectedCells())])==FALSE) %>% 
        filter((!Country %in% unlist(strsplit(gsub('[[:digit:]]+', '', selectedCells()),"\\.")))) %>% 
        filter(is.na(selectedCells())==FALSE) %>% 
        dplyr::rename("Motif"=selectedCells()) %>% 
        group_by(Motif,Country) %>% 
        summarise(Count=n()) %>% 
        ungroup() %>% 
        ggplot(aes(x=Motif,y=Count,fill=Country))+
        geom_col(alpha=0.8) +
        theme_bw()+
        scale_x_continuous(breaks = integer_breaks())+
        scale_y_continuous(breaks = integer_breaks())+
        labs(title = "Other countries exclude pair country (P. vivax)",
             x="STR copies",
             y="Count")
      ggplotly(plot4)
    })
    
    updateNavbarPage(session, "panels_main", selected = "PvCountry_TopSTRPlots")
  }) 

}

# Create Shiny app ----
shinyApp(ui, server)



