### * Description

# A script to try and test code.

### * Setup

library(tidyverse)
library(latex2exp)
library(grid)
devtools::load_all()

### * Trini model

t <- topo(trini_mod)
nodes <- nodes_from_topo(t)
nodes$label <- as.list(nodes$label)
nodes$label[[2]] <- TeX("$\\beta$")
nodes$size <- runif(nrow(nodes), 1, 2)
flows <- flows_from_topo(t)
flows$width <- runif(nrow(flows), 0.2, 2)
save_pdf <- FALSE

if (save_pdf) { pdf("toto.pdf", width = 16, height = 9) }
z <- sankey(t, nodes = nodes, flows = flows, layout = "left2right",
            debug = TRUE, node_f = 1, edge_f = 0.9, edge_n = 32,
            cex_lab = 1.5)
if (save_pdf) { dev.off() }

### * Yuan 2006

y <- new_networkModel() %>%
        set_topo(c("subs -> NH3 -> subs",
                   "NH3 -> Q, E", "E -> Q -> E",
                   "E -> D, M")) %>%
        set_steady("subs") %>%
            set_prop_family("normal_sd")
y <- topo(y)[[1]]
nodes <- nodes_from_topo(y)
nodes$size <- runif(nrow(nodes), 1, 5)
ggtopo(y, edge = "fan")
flows <- flows_from_topo(y)
flows$width <- runif(nrow(flows), 0.2, 5)

grid.newpage()
pushViewport(viewport(layout = grid.layout(ncol = 2)))
pushViewport(viewport(width = 0.5, layout.pos.col = 1))
sankey(y, flows = flows, debug = TRUE, edge_n = 32, edge_f = 0.2, new = FALSE)
grid.rect(gp = gpar(col = "red", lwd = 2))
popViewport()
pushViewport(viewport(width = 0.5, layout.pos.col = 2))
sankey(y, flows = flows, debug = FALSE, edge_f = 0.2, new = FALSE)
grid.rect(gp = gpar(col = "blue", lwd = 2))
popViewport()

z <- sankey(y, nodes = nodes, flows = flows, debug = FALSE, edge_n = 32,
            edge_f = 0.4, node_s = "prop")

### * Radziuk 2001

r <- new_networkModel() %>%
    set_topo("infusion -> plasma -> body -> plasma") %>%
    set_steady(c("infusion", "body"))
r <- topo(r)[[1]]

ggtopo(r, edge = "fan")

sankey(r, debug = TRUE, edge_f = 0.2)
