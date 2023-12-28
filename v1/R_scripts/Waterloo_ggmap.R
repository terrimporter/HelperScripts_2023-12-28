setwd("/home/terri/Waterloo_021718/bubbleMaps")

siteTable <- read.csv(file="Sites.csv", header=TRUE, sep=",")
colnames(siteTable) <-c("Site", "Lat", "Lon", "ESVs", "Genera",
                        "Families", "EPTesvs", "EPTgenera", "EPTfamilies")

### Create map plot ###

library(ggmap)

pdf("mapESVs.pdf")
map1 <- get_map(location = c(lon=-80.57074088,lat=43.470726), zoom = 12,  
               maptype="hybrid", scale="auto")
mapPoints1 <- ggmap(map1) +
  geom_point(data=siteTable, aes(x=Lon, y=Lat, size = ESVs), shape=21, 
             colour="yellow", bg="yellow", stroke=1, alpha=0.3) +
  geom_point(data=siteTable, aes(x=Lon, y=Lat, size = EPTesvs), shape=21, 
             colour="yellow", bg="yellow", stroke=1, alpha=1) +
  scale_shape_discrete(solid=F) +
  theme_bw() +
  guides(size = guide_legend(override.aes = list(colour = "black", bg="white", alpha=1)))+
  labs(size="Number of ESVs")+
  scale_size_continuous(range=c(1,10)) +
  scale_x_continuous(limits = c(-80.65, -80.475), expand = c(0, 0)) +
  scale_y_continuous(limits = c(43.45, 43.5), expand = c(0, 0)) 
  #geom_text(data=siteTable, aes(x=Lon, y=Lat, label = Site), color="white", 
   #       size =3, vjust = 0.5, hjust = 0.5)
mapPoints1
dev.off()

pdf("mapGenera.pdf")
map2 <- get_map(location = c(lon=-80.57074088,lat=43.470726), zoom = 12,  
                maptype="hybrid", scale="auto")
mapPoints2 <- ggmap(map2) +
  geom_point(data=siteTable, aes(x=Lon, y=Lat, size = Genera), shape=21, 
             colour="orange", bg="orange", stroke=1, alpha=0.3) +
  geom_point(data=siteTable, aes(x=Lon, y=Lat, size = EPTgenera), shape=21, 
             colour="orange", bg="orange", stroke=1, alpha=1) +
  scale_shape_discrete(solid=F) +
  theme_bw() +
  guides(size = guide_legend(override.aes = list(colour = "black", bg="white", alpha=1)))+
  labs(size="Number of Genera")+
  scale_size_continuous(range=c(1,10)) +
  scale_x_continuous(limits = c(-80.65, -80.475), expand = c(0, 0)) +
  scale_y_continuous(limits = c(43.45, 43.5), expand = c(0, 0)) 
#geom_text(data=siteTable, aes(x=Lon, y=Lat, label = Site), color="white", 
#       size =3, vjust = 0.5, hjust = 0.5)
mapPoints2
dev.off()

pdf("mapFamilies.pdf")
map3 <- get_map(location = c(lon=-80.57074088,lat=43.470726), zoom = 12,  
                maptype="hybrid", scale="auto")
mapPoints3 <- ggmap(map3) +
  geom_point(data=siteTable, aes(x=Lon, y=Lat, size = Families), shape=21, 
             colour="green", bg="green", stroke=1, alpha=0.3) +
  geom_point(data=siteTable, aes(x=Lon, y=Lat, size = EPTfamilies), shape=21, 
             colour="green", bg="green", stroke=1, alpha=1) +
  scale_shape_discrete(solid=F) +
  theme_bw() +
  guides(size = guide_legend(override.aes = list(colour = "black", bg="white", alpha=1)))+
  labs(size="Number of Families")+
  scale_size_continuous(range=c(1,10)) +
  scale_x_continuous(limits = c(-80.65, -80.475), expand = c(0, 0)) +
  scale_y_continuous(limits = c(43.45, 43.5), expand = c(0, 0)) 
#geom_text(data=siteTable, aes(x=Lon, y=Lat, label = Site), color="white", 
#       size =3, vjust = 0.5, hjust = 0.5)
mapPoints3
dev.off()
