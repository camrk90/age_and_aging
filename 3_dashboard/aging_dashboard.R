library(shinydashboard)
library(shiny)
library(ggplot2)
library(tidyverse)

long_data<- read.table("/scratch/ckelsey4/Cayo_meth/long_data_adjusted.txt")

long_data<- long_data %>%
  group_by(monkey_id) %>%
  mutate(min_age = min(age_at_sampling))

long_data$age_at_sampling<- round(long_data$age_at_sampling, 0)

header<- dashboardHeader(title = "Age vs Aging")

sidebar<- dashboardSidebar(
  sidebarMenu(
    menuItem("Tab 1", tabName = "tab1", icon = icon("chart-line"))
  )
)

body<- dashboardBody(
  tabItem(tabName = "tab1",
          fluidRow(
            box(title = "Plot 1", status = "primary", solidHeader = TRUE,
                plotOutput("plot1", height = 300))
          )
  )
)

ui<- dashboardPage(header, sidebar, body)
                
server <- function(input, output) {
  output$plot1 <- renderPlot({
    
    long_data %>%
      ggplot(aes(x=age_at_sampling, y=reorder(monkey_id, min_age), colour=individual_sex)) +
      geom_path(linewidth = 1.5, alpha = 0.8) +
      geom_point(colour="black") +
      scale_x_continuous(breaks = seq(0, 30, by=2)) +
      scale_colour_manual(values = c("darkolivegreen", "darkmagenta"), name = "Sex") +
      ylab("Individual") +
      xlab("Age") +
      theme_classic(base_size = 24) +
      theme(axis.text.y=element_blank(),
            axis.ticks.y=element_blank())
  })
}
              
shinyApp(ui, server)
