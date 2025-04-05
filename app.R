#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# app.R
library(shiny)
library(readxl)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
#install.packages('plotly')
library(plotly)


ui <- fluidPage(
  titlePanel("Tumour Monitoring Dashboard"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload Excel File", accept = ".xlsx"),
      numericInput("threshold", "Humane Endpoint Volume (mm³)", value = 1500),
      actionButton("run", "Run Analysis"),
      br(), br(),
      downloadButton("downloadPlot", "Download Tumour Plot"),
      downloadButton("downloadGroupPlot", "Download Group Plot"),
      downloadButton("downloadAlerts", "Download Alerts Table"),
      downloadButton("downloadWeightPlot", "Download Weight Plot"),
      downloadButton("downloadWeightAlerts", "Download Weight Alerts Table"),
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Tumour Volume",
                 h4("Tumour Progression Per Mouse"),
                 plotlyOutput("mousePlot"),
                 h4("Mice Exceeding Threshold"),
                 tableOutput("alertsTable"),
                 h4("Average Tumour Volume by Genotype"),
                 plotlyOutput("groupPlot")
        ),
        tabPanel("Weight Monitoring",
                 h4("Mouse Weight Over Time"),
                 plotlyOutput("weightPlot"),
                 h4("Mice Below 80% Initial Weight"),
                 tableOutput("weightAlerts")
        )
      )
    )
  )
)

server <- function(input, output) {
  
  output$downloadWeightPlot <- downloadHandler(
    filename = function() { paste0("weight_plot.png") },
    content = function(file) {
      initial_weights <- weight_data() %>%
        arrange(Mouse_ID, Date) %>%
        group_by(Mouse_ID) %>%
        slice(1) %>%
        ungroup() %>%
        mutate(Initial_Weight = Weight)
      
      ggsave(file, plot = ggplot(weight_data(), aes(x = Date, y = Weight, group = Mouse_ID, color = Mouse_ID)) +
               geom_line(linetype = "solid") +
               geom_point() +
               geom_hline(
                 data = initial_weights,
                 aes(yintercept = 0.8 * Initial_Weight, color = Mouse_ID),
                 linetype = "dashed",
                 show.legend = FALSE
               ) +
               labs(x = "Date", y = "Weight (g)", title = "Mouse Weight Over Time") +
               theme_minimal(),
             width = 8, height = 5)
    }
  )
  
  output$downloadWeightAlerts <- downloadHandler(
    filename = function() { paste0("weight_alerts.csv") },
    content = function(file) {
      initial_weights <- weight_data() %>%
        arrange(Mouse_ID, Date) %>%
        group_by(Mouse_ID) %>%
        slice(1) %>%
        ungroup() %>%
        select(Mouse_ID, Initial_Weight = Weight)
      
      weight_w_thresh <- weight_data() %>%
        left_join(initial_weights, by = "Mouse_ID") %>%
        mutate(Below80 = Weight < (0.8 * Initial_Weight)) %>%
        filter(Below80)
      
      write.csv(weight_w_thresh %>%
                  arrange(Mouse_ID, Date) %>%
                  select(Mouse_ID, Date, Weight, Initial_Weight),
                file, row.names = FALSE)
    }
  )
  
  
  
  data_processed <- eventReactive(input$run, {
    req(input$file)
    
    sheets <- excel_sheets(input$file$datapath)
    
    data <- read_excel(input$file$datapath, sheet = sheets[1], col_types = "guess")
    data <- data %>% fill(Cage, Genotype)
    
    d_cols <- grep("^D\\.\\.\\.", colnames(data), value = TRUE)
    d_cols <- d_cols[d_cols != "D"]
    date_cols <- grep("^Date\\.\\.\\.", colnames(data), value = TRUE)
    
    df_list <- list()
    for (i in seq_along(date_cols)) {
      d_col <- d_cols[i]
      d_col_index <- which(colnames(data) == d_col)
      w_col <- colnames(data)[d_col_index + 1]
      v_col <- colnames(data)[d_col_index + 2]
      date_col <- date_cols[i]
      
      df_temp <- data %>%
        select(`Mouse ID`, Cage, Genotype, all_of(d_col), all_of(w_col), all_of(v_col), all_of(date_col)) %>%
        rename(
          Mouse_ID = `Mouse ID`,
          Group = Cage,
          Length = all_of(d_col),
          Width = all_of(w_col),
          Volume = all_of(v_col),
          Date = all_of(date_col)
        ) %>%
        filter(!is.na(Mouse_ID), !is.na(Date))
      df_list[[i]] <- df_temp
    }
    
    bind_rows(df_list) %>%
      mutate(Date = as.Date(Date))  # ensure Date is treated correctly
  })
  
  weight_data <- eventReactive(input$run, {
    req(input$file)
    
    read_excel(input$file$datapath, sheet = "Mice_Weight") %>%
      mutate(Date = as.Date(Date)) %>%       # ✅ convert Date column
      filter(!is.na(Weight), !is.na(Date))
  })
  
  output$mousePlot <- renderPlotly({
    req(data_processed())
    
    p <- ggplot(data_processed(), aes(x = Date, y = Volume, group = Mouse_ID, color = Mouse_ID,
                                      text = paste("Mouse:", Mouse_ID, "<br>Volume:", Volume, "mm³"))) +
      geom_line(size = 0.8) +
      geom_point(size = 2) +
      geom_hline(yintercept = input$threshold, linetype = "dashed", color = "red") +
      labs(x = "Date", y = "Tumour Volume (mm³)") +
      scale_x_date(date_breaks = "2 days", date_labels = "%d %b") +
      theme_minimal()
    
    ggplotly(p, tooltip = "text")
  })
  
  output$groupPlot <- renderPlotly({
    req(data_processed())
    
    summary_data <- data_processed() %>%
      group_by(Date, Genotype) %>%
      summarise(
        Mean_Volume = mean(Volume, na.rm = TRUE),
        SE = sd(Volume, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
    
    p <- ggplot(summary_data, aes(x = Date, y = Mean_Volume, color = Genotype, group = Genotype)) +
      geom_line(aes(text = paste("Genotype:", Genotype,
                                 "<br>Mean Volume:", round(Mean_Volume, 2), "mm³")), size = 1.2) +
      geom_errorbar(aes(ymin = Mean_Volume - SE, ymax = Mean_Volume + SE), width = 0.2) +
      labs(x = "Date", y = "Mean Tumour Volume (mm³)") +
      scale_x_date(date_breaks = "2 days", date_labels = "%d %b") +
      theme_minimal()
    
    ggplotly(p, tooltip = "text")
  })
  
  
  output$alertsTable <- renderTable({
    req(data_processed())
    data_processed() %>%
      filter(Volume > input$threshold) %>%
      distinct(Mouse_ID, .keep_all = TRUE)
  })
  
  output$weightPlot <- renderPlotly({
    req(weight_data())
    
    initial_weights <- weight_data() %>%
      arrange(Mouse_ID, Date) %>%
      group_by(Mouse_ID) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(Initial_Weight = Weight)
    
    p <- ggplot(weight_data(), aes(x = Date, y = Weight, group = Mouse_ID, color = Mouse_ID,
                                   text = paste("Mouse:", Mouse_ID, "<br>Weight:", Weight, "g"))) +
      geom_line() +
      geom_point() +
      geom_hline(
        data = initial_weights,
        aes(yintercept = 0.8 * Initial_Weight, color = Mouse_ID),
        linetype = "dashed",
        show.legend = FALSE
      ) +
      labs(x = "Date", y = "Weight (g)", title = "Mouse Weight Over Time") +
      scale_x_date(date_breaks = "2 days", date_labels = "%d %b") +
      theme_minimal()
    theme_minimal()
    
    ggplotly(p, tooltip = "text")
  })
  
  output$weightAlerts <- renderTable({
    req(weight_data())
    
    initial_weights <- weight_data() %>%
      arrange(Mouse_ID, Date) %>%
      group_by(Mouse_ID) %>%
      slice(1) %>%
      ungroup() %>%
      select(Mouse_ID, Initial_Weight = Weight)
    
    weight_w_thresh <- weight_data() %>%
      left_join(initial_weights, by = "Mouse_ID") %>%
      mutate(Below80 = Weight < (0.8 * Initial_Weight)) %>%
      filter(Below80)
    
    weight_w_thresh %>%
      arrange(Mouse_ID, Date) %>%
      select(Mouse_ID, Date, Weight, Initial_Weight)
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() { paste0("tumour_plot.png") },
    content = function(file) {
      ggsave(file, plot = ggplot(data_processed(), aes(x = Date, y = Volume, group = Mouse_ID, color = Genotype)) +
               geom_line(show.legend = FALSE, size = 0.8) +
               geom_point(size = 2) +
               geom_hline(yintercept = input$threshold, linetype = "dashed", color = "red") +
               labs(x = "Date", y = "Tumour Volume (mm³)") +
               theme_minimal(),
             width = 8, height = 5)
    }
  )
  output$downloadGroupPlot <- downloadHandler(
    filename = function() { paste0("group_plot.png") },
    content = function(file) {
      ggsave(file, plot = ggplot(data_processed(), aes(x = Date, y = Volume, color = Genotype)) +
               stat_summary(fun = mean, geom = "line", size = 1.2) +
               stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
               labs(x = "Date", y = "Mean Tumour Volume (mm³)") +
               theme_minimal(),
             width = 8, height = 5)
    }
  )
  output$downloadAlerts <- downloadHandler(
    filename = function() { paste0("mice_alerts.csv") },
    content = function(file) {
      write.csv(data_processed() %>%
                  filter(Volume > input$threshold) %>%
                  distinct(Mouse_ID, .keep_all = TRUE), file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
