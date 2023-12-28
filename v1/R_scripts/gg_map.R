library(tidyverse)

#import records from 2003 to 2016 some of which contain lat/lon and country annotations
#df contains spotty lat/lon and country annot
df<-read.csv('/home/terri/CO1_ncbi/reformatted_taxids/AllEukaryota/gg_latlon.csv', 
             header=T)
names(df)<-c("gb","latitude", "longitude","country_alone")
head(df)

#count number of records per country
#df2 contains one row per country with num/prop records
df2<-as.data.frame(table(df$country_alone))
names(df2)<-c("country_alone","records")
head(df2)

#get proportion of total records
sumrecords<-sum(df2$records)
df2$records_prop<-apply(df2[c("records")], 1, function(x){x/sumrecords*100})

#join original df with records and records_prop
#df3 merges original spotty df with df2 with records, just keeps lines in common
#df3<-merge(df,df2,by="country_alone")
#head(df3)

#get lat and lon to draw world country polygons
wm<-map_data("world")

#mege map data with df2 that contains records and records_prop
wm <- left_join(wm, df2, by=c("region"="country_alone"))

#draw map
mp<-ggplot() +
  theme_void() +
  #color countries by num records with country annot
  geom_polygon(data=wm, aes(x=long, y=lat, group=group, fill=records), color=NA) +
  coord_equal() +
  scale_fill_gradientn(colours = c('#461863','#404E88','#2A8A8C','#7FD157','#F9E53F'),
                       values = scales::rescale(c(25000,50000,100000,200000,400000)),
                       labels = c("< 25,000","25,000 - 49,999","50,000 - 99,999","100,000 - 199,999","200,000 - 399,999"),
                       breaks =c(25000,50000,100000,200000,400000),
                       limits=c(0,400000)) +
  guides(fill = guide_legend(reverse = T))

#add points for records with lat/lon annot from all GenBank CO1 Records
mp <- mp + 
  geom_point(data=df, mapping=aes(x=longitude, y=latitude, shape="."), 
             show.legend=FALSE, color="violetred", size=0.15) 

#add points for BARCODE records

#add legend
mp <- mp +
  labs(fill='Country Annotations',
       x='Longitude',
       y='Latitude',
       shape='Lat/Lon Annotations')

pdf("map.pdf")
mp
dev.off()
