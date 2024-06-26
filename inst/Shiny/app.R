library(shiny)
library(shinythemes)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggsci)
source("process_data.R")
# Assuming df and other necessary data are loaded here
load("results/COVID_FLU_2016_20231229_res2_0305.RData")
colnames(df) <- c("probabilities", "signals", "logROR","lower","upper","time","AE")
df$time<-format(df$time, "%Y-%m")

ui <- navbarPage(theme = shinytheme("cosmo"),
  title = "SEquentially Adjusted Bayesian Analysis on Safety Surveillance (SEABASS)",
  
  # Page 1: Current COVID vs FLU results
  tabPanel("COVID-19 vs FLU study",
           sidebarLayout(
             sidebarPanel(
               p("Last updated: 2023-12"),
               checkboxInput("includeSignalAEs", "Only show signal AEs", value = FALSE),
               selectInput("showTopAEs", "Select AE number:", choices = c("All","5", "10", "50", "100")),
               uiOutput("dynamicAEs"),  # This will dynamically generate the AE select input
               #uiOutput("monitorTimeSelectors"),
               fluidRow(
                 column(6, selectInput("orderCriterion", "Order by:", choices = c("Risk Probabilities" = "probabilities", "logROR" = "logROR"))),
                 column(6, selectInput("orderTimePoint", "Order at Time:", choices = unique(df$time)))
               ),
               fluidRow(
                 column(6, uiOutput("monitorStartTime")),
                 column(6, uiOutput("monitorEndTime"))
               ),
               numericInput("threshold", "Probability Decision Threshold:", value = 0.975, min = 0, max = 1, step = 0.005)
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Risk Probability", plotOutput("riskPlot"), downloadButton("downloadRiskPlot", "Download Risk Probability Plot")),
                 tabPanel("Reporting Odds Ratio", plotOutput("logRORPlot"), downloadButton("downloadLogRORPlot", "Download logROR Plot"))
               )
             )
           )
  ),
  
  # Page 2: Data Upload and Model Fitting
  tabPanel("Fit New Model",
           tabsetPanel(
             tabPanel("Pre-process data",
                      sidebarLayout(
                        sidebarPanel(
                          textInput("datapath", "Data path (e.g. path/to/VAERS/)"),
                          fluidRow(
                            column(6, selectInput("startYear", "Surveillance Start Year:", choices = 1990:2024, selected = 2016)),
                            column(6, selectInput("endYear", "Surveillance End Year:", choices = 1990:2024, selected = 2023))
                          ),
                          fluidRow(
                            column(6, textInput("targetVaccine", "Target Vaccine Type:", value = "COVID19")),
                            column(6, numericInput("targetVaccineCount", "Filter by count:", value = 100, min = 0))
                          ),
                          fluidRow(
                            column(6, textInput("referenceVaccine", "Reference Vaccine Type:", value = "FLU")),
                            column(6, numericInput("referenceVaccineCount", "Filter by count:", value = 100, min = 0))
                          ),
                          downloadButton("download", "Process & Download")
                        ),
                        mainPanel(
                          # Placeholder for results or outputs related to the analysis
                          h3("Processing..."),
                          #verbatimTextOutput("datahead")
                        )
                      )
             ),
             # You can add more tabs for different sections of Model Fitting here
             tabPanel("Fit model",
             ),
             tabPanel("Update model",
             ),
             tabPanel("Evaluate results",
             ),
             tabPanel("Help",
             ),
           )
  )
  
)


server <- function(input, output, session) {
  
  # Define function to generate select input choices
  monitor_time_choices <- reactive({
    sort(unique(df$time))
  })
  
  # Dynamic UI for AE selection
  output$dynamicAEs <- renderUI({
    req(input$orderCriterion, input$orderTimePoint)
    
    # Filter the data for the selected time point
    timePointData <- df[df$time == input$orderTimePoint, ]
    
    # Filter AEs based on signal selection
    if (input$includeSignalAEs) {
      signal_AEs <- unique(timePointData$AE[timePointData$signals == 1])
    } else {
      signal_AEs <- unique(timePointData$AE)
    }
    
    # Order the AEs based on the selected criterion
    ordered_AEs <- timePointData[order(timePointData[[input$orderCriterion]], decreasing = TRUE), "AE"]
    unique_ordered_AEs <- intersect(ordered_AEs, signal_AEs)  # Ensure AEs are unique and optionally include only signal AEs
    
    if (input$showTopAEs != "All") {
      num_top_AEs <- as.numeric(input$showTopAEs)
      unique_ordered_AEs <- head(unique_ordered_AEs, num_top_AEs)
    }
    
    selectInput("selectedAEs", "Select Adverse Events:",
                choices = unique_ordered_AEs, multiple = TRUE, selectize = TRUE)
  })
  
  observe({
    # Get all unique time points
    unique_time_points <- sort(unique(df$time))
    
    # Update the select input for time points
    updateSelectInput(session, "orderTimePoint", 
                      choices = unique_time_points)
  })
  
  output$monitorStartTime <- renderUI({
      selectInput("monitorStartTime", "Start Date:", choices = monitor_time_choices(), selected = min(monitor_time_choices()))

  })
  
  output$monitorEndTime <- renderUI({
      selectInput("monitorEndTime", "End Date:", choices = monitor_time_choices(), selected = max(monitor_time_choices()))
  })
  
  # Filtered Data based on selected AEs
  filteredData <- reactive({
    req(input$selectedAEs)
    # Adjust the filtering logic according to your actual data structure
    df %>% 
      filter(AE %in% input$selectedAEs, time >= input$monitorStartTime, time <= input$monitorEndTime)
  })
  
  output$riskPlot <- renderPlot({
    req(filteredData())
    data_to_plot <- filteredData()
    threshold <- input$threshold
    
    ggplot(data_to_plot, aes(x = as.Date( paste(time, "01", sep = "-")), y = probabilities, color = AE)) +
      geom_line() +
      geom_point()+
      ylim(0, 1) +
      scale_color_jco() + scale_fill_jco()+
      geom_hline(yintercept = threshold, linetype = "dashed", color = "darkgrey") +
      theme_minimal() +
      theme(legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +  # Adjust text size
      xlab("Date") + ylab("Posterior risk probability") +
      ggtitle("Posterior Risk Probability Over Time")
  })
  
  output$downloadRiskPlot <- downloadHandler(
    filename = function() {
      paste("risk-probability-plot-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(filteredData())
      data_to_plot <- filteredData()
      threshold <- input$threshold
      
      # Generate the plot
      plot <-  ggplot(data_to_plot, aes(x = as.Date( paste(time, "01", sep = "-")), y = probabilities, color = AE)) +
        geom_line() +
        geom_point()+
        ylim(0, 1) +
        geom_hline(yintercept = threshold, linetype = "dashed", color = "darkgrey") +
        theme_minimal() +
        scale_color_jco() + scale_fill_jco()+
        theme(legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +  # Adjust text size
        xlab("Date") + ylab("Posterior risk probability") +
        ggtitle("Posterior Risk Probability Over Time")
      
      # Save the plot to the specified file
      ggsave(file, plot,dpi = 300, width = 10, height = 6)
    }
  )
  
  output$logRORPlot <- renderPlot({
    req(filteredData())  # Make sure filteredData contains the required CI columns
    data_to_plot <- filteredData()
    
    ggplot(filteredData(), aes(x = as.Date( paste(time, "01", sep = "-")), y = logROR, color = AE, fill = AE)) +  # Map fill to AE as well
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, linetype = "dashed", linewidth=0.5) +
      geom_line() +
      geom_point()+
      theme_minimal() +
      scale_color_jco() + scale_fill_jco()+
      theme(legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +  # Adjust text size
      xlab("Date") + ylab("logROR") +
      ggtitle("logROR Over Time")
  })
  
  output$downloadLogRORPlot <- downloadHandler(
    filename = function() {
      paste("logROR-plot-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      # Generate the plot
      plot <-   ggplot(filteredData(), aes(x = as.Date( paste(time, "01", sep = "-")), y = logROR, color = AE, fill = AE)) +  # Map fill to AE as well
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, linetype = "dashed", linewidth=0.5) +
        geom_line() +
        geom_point()+
        theme_minimal() +
        scale_color_jco() + scale_fill_jco()+
        theme(legend.text = element_text(size = 12), axis.text = element_text(size = 12)) +  # Adjust text size
        xlab("Date") + ylab("logROR") +
        ggtitle("logROR Over Time")
      
      # Save the plot to the specified file
      ggsave(file, dpi = 300, plot, width = 10, height = 6)
    }
  )
  
  # Add functionality to update signals based on threshold
  observeEvent(input$threshold, {
    df$signals <- ifelse(df$probabilities > input$threshold, 1, 0)
  })
  
  output$download <- downloadHandler(
    filename = function() {
      paste("processed-data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(input$datapath, input$startYear, input$endYear, input$targetVaccine, input$referenceVaccine, input$targetVaccineCount, input$referenceVaccineCount)
      
      # Call your data processing function with the inputs
      processed_data <- processData(
        datapath = input$datapath, 
        startyear = input$startYear, 
        endyear = input$endYear, 
        case = input$targetVaccine, 
        control = input$referenceVaccine,
        case_cnt = input$targetVaccineCount,
        control_cnt = input$referenceVaccineCount
      )
      
      # Write the processed data to a CSV file
      write.csv(processed_data, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
