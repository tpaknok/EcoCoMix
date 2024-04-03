library(hexSticker)
library(ggplot2)
library(icons)

library("ggtree")
library(extr)
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)


p <- ggtree(tree, color="#00AFAF",layout='circular',size=1.5)+
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

ggsave("phylo.png")
image <- c("C:/Users/pakno/OneDrive/Desktop/Test/Picture1.png")

sticker(image, package="CPR",
        p_size=30,
        p_x = 1,
        p_y = 1.5,
        p_color="#00AFAF",
        h_color="#00AFAF",
        h_fill="black",
        s_x=1.02,
        s_y=.7,
        s_width=0.8,
        s_height=0.8)

