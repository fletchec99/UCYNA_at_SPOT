#Graphing SPOT on a map 
#Adapted from: https://r-spatial.org/r/2018/10/25/ggplot2-sf-2.html
#09.08.2022

setwd("../SPOT_16S_18S_eLSA/")

install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))

install.packages("devtools")
devtools::install_github("ropensci/rnaturalearthhires")
library("rnaturalearth")

library("ggplot2")
theme_set(theme_bw())

library("ggrepel")

world <- ne_countries(scale = "large", returnclass = "sf")
class(world)

sites <- data.frame(longitude = -118.4, latitude=33.55)
sites

labels <- data.frame(ID=c("Santa Catalina Island", "USC", "SPOT"), 
                     longitude=c(-118.416, -118.285, -118.4), 
                     latitude=c(33.388, 34.022, 33.55))

map=ggplot(data = world) +
  geom_sf() +
  geom_point(data = sites, aes(x = longitude, y = latitude), size = 4, 
             shape = 16, fill = "black") +
  coord_sf(xlim = c(-119, -117.75), ylim = c(33, 34.1), expand = FALSE) +
  geom_text_repel(data = labels, aes(longitude, latitude, label = ID), fontface = "bold", size=4, nudge_x=c(0, 0, 0.5), nudge_y=c(-0.2, 0, 0)) + 
  xlab("Longitude") + ylab("Latitude") + 
  theme(axis.text.x = element_text(face = "bold", color="black", angle=45, size=12, hjust=1), 
        axis.text.y = element_text(colour = "black", face = "bold", size=12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 14, face = "bold")) 
map

ggsave(map, file="SPOT_map_09.08.2022.pdf", height=4, width=4, units="in")


  
