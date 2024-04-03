library(hexSticker)
library(ggplot2)
library(icons)

library("ggtree")
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)


p <- ggtree(tree, color="black",layout='circular',size=1.5)+
    theme(
    panel.border = element_rect(colour='transparent',fill=NA),
    panel.background = element_rect(colour='transparent',fill="transparent"),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

plot(p)

ggsave("man/figures/phylo.png")
image <- c("man/figures/CPR.png")

p <- sticker(image, package="CPR",
        p_size=25,
        p_x = 1,
        p_y = 1.5,
        p_color="#00AFAF",
        h_color="#00AFAF",
        h_fill="black",
        s_x=1.0,
        s_y=.8,
        s_width=0.5,
        s_height=0.5,
        h_size=1.2,
        filename="man/figures/logo.png")

plot(p)
ggsave("man/figures/logo.png",width=4.1,height=4.1,dpi=300,units="cm")

img <- c("man/figures/logo.png")
library(usethis)
use_logo(img, geometry = "240x278", retina = TRUE)
