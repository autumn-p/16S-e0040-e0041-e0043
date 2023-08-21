library(RColorBrewer)
library(cowplot)

theme_set(theme_cowplot())

DEFAULTS.DARK2 <- brewer.pal(n=8, name="Dark2")
DEFAULTS.SET1 <- brewer.pal(n=8, name="Set1")
DEFAULTS.PAIRED <- brewer.pal(n=8, name="Paired")

DEFAULTS.PALETTE <- DEFAULTS.SET1
DEFAULTS.PALETTE2 <- DEFAULTS.DARK2
DEFAULTS.PALETTE.PAIRED <- c(DEFAULTS.PAIRED[1],DEFAULTS.PAIRED[2],
                             DEFAULTS.PAIRED[5],DEFAULTS.PAIRED[6])
DEFAULTS.PALETTE.IBM.COLORBLINDSAFE <- c("#648FFF", "#DC267F", "#FE6100",
                                         "#FFB000", "#785EF0")
DEFAULTS.PALETTE.TOL.COLORBLINDSAFE <- c("#332288", "#44AA99", "#88CCEE", "#DDCC77",
                                         "#CC6677", "#AA4499", "#882255", "#117733")

DEFAULTS.SHAPES <- c(16,15,17)

# Set common formatting for all plots.
DEFAULTS.THEME_ALL <- theme(
  text=element_text(size=10),
  axis.title=element_text(size=10),
  axis.text=element_text(size=7),
  strip.text.x=element_text(margin = margin(3,0,3,0), size=10),
  strip.text.y=element_text(margin = margin(0,3,0,3), size=10),
  strip.background=element_blank())

DEFAULTS.THEME_PRINT <- theme(
  text=element_text(size=7),
  axis.title=element_text(size=7),
  axis.text=element_text(size=5.5),
  strip.text.x=element_text(margin = margin(3,0,3,0), size=7),
  strip.text.y=element_text(margin = margin(0,3,0,3), size=7),
  strip.background=element_blank(),
  legend.text=element_text(size=6),
  legend.title=element_text(size=7),
  plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))


DEFAULTS.THEME_PRES <- theme(
  text=element_text(size=14),
  axis.title=element_text(size=14, face="bold"),
  axis.text=element_text(size=8),
  strip.text.x=element_text(size=10, margin = margin(3,0,3,0)),
  strip.text.y=element_text(size=10, margin = margin(0,3,0,3)),
  strip.background=element_blank())