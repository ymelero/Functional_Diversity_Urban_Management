rao_simplified <- grViz("
digraph SEM_RAO {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50, penwidth = 0.5,
        fontcolor = black, style = rounded]
  edge [color = grey70, penwidth = 0.5, arrowsize = 0.5]

  Tipo [label = 'Space type', shape = box, style = 'rounded,filled',
        colorscheme = gnbu3, fillcolor = 1]

  V1   [label = 'C3-grass', shape = ellipse, style = filled, fillcolor = peachpuff]
  Resp [label = 'RaoQ', shape = box, style = rounded]

  Tipo -> V1  [color = peachpuff, penwidth = 1.2]
  V1   -> Resp [color = peachpuff, penwidth = 1.2]

  graph [nodesep = 0.2, ranksep = 0.3]
}
")

rao_simplified

library(DiagrammeRsvg)
library(rsvg)
svg_code <- export_svg(rao_simplified)
rsvg_pdf(charToRaw(svg_code), file = "Results/Figures/SEM_results_RAO.pdf")

# Bind with boxplot
#install.packages("pdftools")
library(magick)
library(grid)
library(gridExtra)

# Read PDFs
img1 <- image_read_pdf("Results/Figures/SEM_results_RAO.pdf", density = 300)
img2 <- image_read_pdf("Results/Figures/Fig_Predicted_Type_C3grass_RAO.pdf", density = 300)

# Convert to rasterGrob 
rao_g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
rao_g2 <- rasterGrob(as.raster(img2), interpolate = TRUE)

# Combine
grid.arrange(rao_g1, rao_g2, ncol = 2, widths = c(0.45, 3))
grid.text("(a)", x = 0.01, y = 0.95, gp = gpar(fontsize = 10))
grid.text("(b)", x = 0.14, y = 0.95, gp = gpar(fontsize = 10))

