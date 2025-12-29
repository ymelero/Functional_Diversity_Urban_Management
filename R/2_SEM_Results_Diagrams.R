# ---------------------------------------------------------------------------------------------- #
# Resulting SEMs per cluster
# Author: [YM]
# ---------------------------------------------------------------------------------------------- #

library(DiagrammeR)

# 1. SEM C2 
c2 <- grViz("
digraph SEM_B {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50, penwidth = 0.5,
        fontcolor = black, style = rounded]

  # Gris más claro y más fino
  edge [color = '#d9d9d9', penwidth = 0.3, arrowsize = 0.5]

  Area [label = 'Area total', shape = box,style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Conn [label = 'Connectivity (C500)', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Tipo [label = 'Space type', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]

  graph [nodesep = 0.2]

  { rank = same; Area; Conn; Tipo }

  Area -> Conn [style = invis, constraint = false, weight = 50]

  V1 [label = 'C3-grass', shape = ellipse, style = filled, fillcolor = peachpuff]
  V2 [label = 'P-grass',  shape = ellipse, style = filled, colorscheme=gnbu3, fillcolor=1]
  V3 [label = 'Herb-grass', shape = ellipse]
  V4 [label = 'V-grass',  shape = ellipse]

  Resp [label = 'Presence / Abundance', shape = box]

  {rank = same; Area; Conn; Tipo}
  {rank = same; V1; V2; V3; V4}
  {rank = same; Resp}

  # Exógenos -> vegetación
  Area -> V1 [color=none, style=invis, arrowsize=0]
  Area -> V2 [color='#e0f3db', penwidth=0.5]
  Area -> V3 [color=none, style=invis, arrowsize=0]
  Area -> V4

  Tipo -> V1 [color=peachpuff, penwidth=1.7]
  Tipo -> V2 [color='#e0f3db', penwidth=1.1]
  Tipo -> V3
  Tipo -> V4

  # Vegetación -> respuesta
  V1 -> Resp [color=peachpuff, penwidth=1.7]
  V2 -> Resp [color='#e0f3db', penwidth=1.4]
  V3 -> Resp [color=none, style=invis, arrowsize=0]
  V4 -> Resp [color=none, style=invis, arrowsize=0]

  # Directos exógenos -> respuesta
  Tipo -> Resp [color='#e0f3db', penwidth=2]
  Area -> Resp [color=none, style=invis]
  Conn -> Resp [color='#e0f3db', penwidth=0.8]

  Area -> Tipo [dir=both, arrowhead=none, arrowtail=none,
                style=dashed, color='#d9d9d9']
}
")

c2

library(rsvg)
svg_code <- export_svg(c2)
rsvg_pdf(charToRaw(svg_code), file = "Results/Figures/SEM_results_C2.pdf")


# 2. SEM C3
# Resulting SEMs per cluster
c3 <- grViz("
digraph SEM_B {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50, penwidth = 0.5,
        fontcolor = black, style = rounded]

  # Flechas no significativas → gris claro y finas
  edge [color = '#d9d9d9', penwidth = 0.3, arrowsize = 0.5]

  # Nodos
  Area [label = 'Area total', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Conn [label = 'Connectivity (C500)', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Tipo [label = 'Space type', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]

  graph [nodesep = 0.2]
  { rank = same; Area; Conn; Tipo }

  Area -> Conn [style = invis, constraint = false, weight = 50]

  V1 [label = 'C3-grass', shape = ellipse, style = filled, fillcolor = peachpuff]
  V2 [label = 'P-grass',  shape = ellipse]
  V3 [label = 'Herb-grass', shape = ellipse]
  V4 [label = 'V-grass',  shape = ellipse, style = filled, fillcolor = peachpuff]

  Resp [label = 'Presence / Abundance', shape = box]

  {rank = same; Area; Conn; Tipo}
  {rank = same; V1; V2; V3; V4}
  {rank = same; Resp}

  # ---------------------------
  # Exógenos -> Vegetación
  # ---------------------------

  # Area (solo P-grass significativo)
  Area -> V1 [color=none, style=invis]
  Area -> V2
  Area -> V3 [color=none, style=invis]
  Area -> V4 [color=none, style=invis]

  # Space type -> mediadores
  # Nivel 1 (máx C3): Tipo -> V1
  Tipo -> V1 [color=peachpuff, penwidth=1.6]

  # Nivel 3: Tipo -> V4
  Tipo -> V4 [color=peachpuff, penwidth=1.0]

  Tipo -> V2
  Tipo -> V3

  # ---------------------------
  # Vegetación -> Respuesta
  # ---------------------------

  # Nivel 1 (máx C3): C3-grass -> Resp
  V1 -> Resp [color=peachpuff, penwidth=1.6]

  # Nivel 3: V-grass -> Resp
  V4 -> Resp [color=peachpuff, penwidth=1.0]

  # Herb-grass no significativo
  V3 -> Resp [color=none, style=invis]
  V2 -> Resp [color=none, style=invis]

  # ---------------------------
  # Directos -> Respuesta
  # ---------------------------

  # Nivel 4: más fino significativo (Space type → Resp)
  Tipo -> Resp [color='#e0f3db', penwidth=0.7]

  # Nivel 3: Area → Resp
  Area -> Resp [color='#e0f3db', penwidth=1.0]

  # Nivel 2: Connectivity → Resp
  Conn -> Resp [color='#e0f3db', penwidth=1.3]

  # Covarianza
  Area -> Tipo [dir=both, arrowhead=none, arrowtail=none,
                style=dashed, color='#d9d9d9']
}
")

svg_code <- export_svg(c3)
rsvg_pdf(charToRaw(svg_code), file = "Results/Figures/SEM_results_C3.pdf")

# 3. Cluster C4
# Resulting SEMs per cluster
c4 <- grViz("
digraph SEM_B {
  graph [layout = dot, rankdir = TB]

  node [fontname = 'Calibri', fontsize = 10, color = grey50, penwidth = 0.5,
        fontcolor = black, style = rounded]

  # Flechas no significativas → gris claro y finas
  edge [color = '#d9d9d9', penwidth = 0.3, arrowsize = 0.5]

  # Nodos
  Area [label = 'Area total', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Conn [label = 'Connectivity (C500)', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]
  Tipo [label = 'Space type', shape = box, style = 'rounded,filled', colorscheme=gnbu3, fillcolor=1]

  graph [nodesep = 0.2]
  { rank = same; Area; Conn; Tipo }

  Area -> Conn [style = invis, constraint = false, weight = 50]

  V1 [label = 'C3-grass', shape = ellipse, style = filled, fillcolor = peachpuff]
  V2 [label = 'P-grass',  shape = ellipse]
  V3 [label = 'Herb-grass', shape = ellipse, style = filled, fillcolor = peachpuff]
  V4 [label = 'V-grass',  shape = ellipse]

  Resp [label = 'Presence / Abundance', shape = box]

  {rank = same; Area; Conn; Tipo}
  {rank = same; V1; V2; V3; V4}
  {rank = same; Resp}

  # Exógenos -> vegetación

  Area -> V1 [color=none, style=invis]
  Area -> V2
  Area -> V3 [color=none, style=invis]
  Area -> V4 [color=none, style=invis]

  # Space type -> mediadores
  Tipo -> V1 [color=peachpuff, penwidth=1.6]      # Nivel 1
  Tipo -> V3 [color=peachpuff, penwidth=1.0]      # Nivel 2

  Tipo -> V2
  Tipo -> V4

  # Vegetación -> respuesta
  V1 -> Resp [color=peachpuff, penwidth=1.6]            # C3 no sig
  V3 -> Resp [color=peachpuff, penwidth=1.0]       # Nivel 2
  V4 -> Resp [color=none, style=invis]
  V2 -> Resp [color=none, style=invis]

  # Directos -> respuesta
  Area -> Resp [color='#e0f3db', penwidth=0.7]      # Nivel 3
  Conn -> Resp [color=none, style=invis]
  Tipo -> Resp [color=none, style=invis]

  # Covarianza
  Area -> Tipo [dir=both, arrowhead=none, arrowtail=none,
                style=dashed, color='#d9d9d9']
}
")
c4
svg_code <- export_svg(c4)
rsvg_pdf(charToRaw(svg_code), file = "Results/Figures/SEM_results_C4.pdf")


# Bind them

#install.packages("pdftools")
library(magick)
library(grid)
library(gridExtra)

# Read PDFs
img1 <- image_read_pdf("Results/Figures/Working/SEM_results_C2.pdf", density = 300)
img2 <- image_read_pdf("Results/Figures/Working/SEM_results_C3.pdf", density = 300)
img3 <- image_read_pdf("Results/Figures/Working/SEM_results_C4.pdf", density = 300)

# Convert to rasterGrob 
c2_g1 <- rasterGrob(as.raster(img1), interpolate = TRUE)
c3_g3 <- rasterGrob(as.raster(img2), interpolate = TRUE)
c4_g4 <- rasterGrob(as.raster(img3), interpolate = TRUE)

grid.arrange(c2_g1, c3_g3, c4_g4,
             ncol = 1,
             heights = c(1, 1, 1))   # o ajustas proporciones

grid.text("(a)", x = 0.1, y = 0.98, gp = gpar(fontsize = 14))
grid.text("(b)", x = 0.1, y = 0.65, gp = gpar(fontsize = 14))
grid.text("(c)", x = 0.1, y = 0.32, gp = gpar(fontsize = 14))

