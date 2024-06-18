### * Description

# Script to test the use of shiny for local reactive plots

### * Setup

library(isotracer)
library(tidyverse)
library(shiny)
library(shinyWidgets)

library(reactlog)
reactlog::reactlog_enable()

### * Prototype

### ** Network

# The user builds and provides a network model
m <- new_networkModel() %>%
  set_topo("NH4 -> algae -> grazer", "algae, grazer -> omni") %>%
  set_steady("NH4") %>%
  set_init(tibble(comp = c("NH4", "algae", "grazer", "omni"),
                  size = c(1, 0.2, 0.4, 0.3),
                  prop = c(0.00366, 0.00366, 0.00366, 0.00366)),
           comp = "comp", size = "size", prop = "prop") %>%
  set_params(tibble::tribble(
                 ~parameter,  ~value,
                      "eta",     0.3,
             "lambda_algae",    0.01,
            "lambda_grazer",    0.01,
               "lambda_NH4",       0,
               "lambda_omni",     0.02,
  "upsilon_algae_to_omni",    0.03,
     "upsilon_grazer_to_omni",    0.02,
  "upsilon_algae_to_grazer",    0.03,
     "upsilon_NH4_to_algae",    0.02,
                     "zeta",     0.4
  ))

pulses <- tibble::tribble(
                    ~ comp, ~ time, ~ unmarked, ~ marked,
                    "NH4",  10, 0, 0.02,
                    "NH4",  20, 0, -0.02
                  )
m <- m %>%
  add_pulse_event(pulses = pulses, comp = "comp", time = "time",
                  unmarked = "unmarked", marked = "marked")

### ** UI

# https://stackoverflow.com/questions/30502870/shiny-slider-on-logarithmic-scale
log10_slider_vals <- 10^seq(-5, 2, length.out = 512)
log10_slider_ticks <- as.character(signif(log10_slider_vals, 3))
log10_slider_ticks_e <- formatC(log10_slider_vals, format = "e", digits = 2)
log10_slider_ticks_e <- gsub("e", "x10^", log10_slider_ticks_e) # If we keep the "e" scientific notation, it is somehow converted back to decimal when displayed by the slider.
log10_slider_ticks_e <- gsub("1.00x10\\^-05", "10^-5", log10_slider_ticks_e)
i <- log10_slider_vals < 0.01
log10_slider_ticks[i] <- log10_slider_ticks_e[i]
names(log10_slider_vals) <- log10_slider_ticks

param_control <- function(p) {
  slider <- sliderTextInput(paste0("param_", p), label = NULL,
                            choices = log10_slider_ticks,
                            selected = "0.1")
  list(fluidRow(column(6, markdown(p)),
                column(6, slider)))
}

network_model_panel <- function(nm) {
  params <- params(nm)
  e <- list(h2("Network model"))
  # TODO Replace this by table with sliders?
  # Add upsilon controls
  e <- c(e, list(h4("Tranfer rates")))
  u_i <- which(grepl("upsilon_", params[["in_model"]]))
  for (i in u_i) {
    p <- params[["in_model"]][i]
    e <- c(e, list(param_control(p)))
  }
  # Add lambda controls
  e <- c(e, list(h4("Loss rates")))
  l_i <- which(grepl("lambda_", params[["in_model"]]))
  for (i in l_i) {
    p <- params[["in_model"]][i]
    e <- c(e, list(param_control(p)))
  }
  # Add eta/zeta controls
  e <- c(e, list(h4("Observations dispersion")))
  d_i <- which(grepl("^[z]*eta", params[["in_model"]]))
  for (i in d_i) {
    p <- params[["in_model"]][i]
    e <- c(e, list(param_control(p)))
  }
  # Build panel
  do.call(wellPanel, e)
}

display_panel <- function() {
  wellPanel(
    h2("Predictions"),
    plotOutput("traj_plot")
  )
}

plot_control_panel <- function() {
  wellPanel(
    h2("Plot controls"),
    "Here go plotting controls"
  )
}

ui <- function(nm) {
  fluidPage(fluidRow(
    column(width = 3, network_model_panel(nm)),
    column(width = 6, display_panel()),
    column(width = 3, plot_control_panel())
  ), theme = bslib::bs_theme(bootswatch = "flatly"))
}

server <- function(nm) {
  function(input, output) {
    params_table <- reactiveVal(params(nm))
    params_names <- params(nm)[["in_model"]]
    lapply(params_names, function(p) {
      slider_id <- paste0("param_", p)
      observeEvent(input[[slider_id]], {
        new_params <- params_table()
        new_params[["value"]][new_params[["in_model"]] == p] <- log10_slider_vals[input[[slider_id]]]
        params_table(new_params)
      })
    })
    output[["traj_plot"]] <- renderPlot({
      m <- m %>% set_params(deframe(params_table())) %>%
        project(at = 0:40)
      plot(m, facet_col = "compartment", facet_row = "type",
           scale = "all")
    })
  }
}

shinyApp(ui(m), server(m))
