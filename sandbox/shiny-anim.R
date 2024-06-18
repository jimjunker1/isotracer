### * Description

# Script to test the use of shiny for local reactive plots

### * Setup

library(isotracer)
library(tidyverse)
library(shiny)

### * Animation

library(profvis)

p <- profvis({
  
### ** UI

ui <- fluidPage(
  plotOutput("traj"),
  sliderInput("upsilon", label = "upsilon", value = 0.1, min = 0, max = 0.5),
  actionButton("stop_button", label = "Quit")
)

### ** Server

server <- function(input, output, session) {
  # Init
  m <- new_networkModel() %>%
  set_topo("N -> grazer") %>%
  set_steady("N") %>%
  set_init(tibble(comp = c("N", "grazer"),
                  size = c(1, 0.2),
                  prop = c(0.02, 0.00366)),
           comp = "comp", size = "size", prop = "prop") %>%
  set_params(c("eta" = 0.3, "zeta" = 0.3, "lambda_grazer" = 0.05,
               "lambda_N" = 0, "upsilon_N_to_grazer" = 0.02)) %>%
    project(at = seq(0, 10, length.out = 256))
  # Reactive event
  observeEvent(input[["stop_button"]], {
    stopApp()
  })
  color <- reactive({
    invalidateLater(20)
    sample(colors(), 1)
  })
  u <- reactiveVal(0.1)
  # Output
  output[["traj"]] <- renderPlot({
    p <- params(m)
    u(isolate(u()) * 1.05)
    p[["value"]][p[["in_model"]] == "upsilon_N_to_grazer"] <- u()
    m <- set_params(m, deframe(p)) %>%
      project(at = seq(0, 10, length.out = 256))
    z <- m$trajectory[[1]]
    x <- cbind(tibble(time = z[["timepoints"]][[1]]),
               z[["proportions"]][[1]])
    plot(x[["time"]], x[["grazer"]], xlim = c(0, 10), ylim = c(0, 0.02),
         type = "l",  lwd = 3,
         col = color(), pch = 19)
  })
}

### ** Run app

runApp(shinyApp(ui, server))
})

