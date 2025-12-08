### Halibut predation story -----

#load libraries ----
library(tidyverse)
library(sf)
library(rnaturalearth)
library(MarConsNetData)
library(ggspatial)
library(grid)
library(lubridate)

s2_as_sf = FALSE

#projections ------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#function for playing with the silouettes

png_to_rgba <- function(img) {
  # img is the array returned by png::readPNG()
  
  # Case 1: grayscale (2D array) → RGB
  if (length(dim(img)) == 2) {
    img <- array(rep(img, 3), dim = c(dim(img), 3))
  }
  
  # Case 2: grayscale + alpha (2 channels) → RGBA
  else if (dim(img)[3] == 2) {
    gray <- img[,,1]
    alpha <- img[,,2]
    rgba <- array(0, dim = c(dim(gray), 4))
    rgba[,,1] <- gray
    rgba[,,2] <- gray
    rgba[,,3] <- gray
    rgba[,,4] <- alpha
    img <- rgba
  }
  
  # Case 3: already RGB or RGBA → unchanged
  # (3 or 4 channels)
  
  return(img)
}

#load the St. Anns Bank MPA shapefiles
sab_zones <- read_sf("data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(CanProj)

sab  <- sab_zones%>%
  st_transform(utm)%>%
  st_buffer(0.0025)%>% #this is small buffer that the vertex coordinate rounding issue
  st_union()%>% #gets rid of the zones
  st_transform(latlong)%>%
  st_as_sf()%>%
  rename(geometry=x)


#load bioregion and get plotting limits----
bioregion <- read_sf("data/Shapefiles/MaritimesPlanningArea.shp")%>%
  st_transform(CanProj)

mar_net <- data_draft_areas()%>%st_transform(latlong)

plot_lim <- bioregion%>% #bounding area for the inset (zoom out)
  st_transform(utm)%>%
  st_buffer(50)%>%
  st_transform(CanProj)%>%
  st_bbox()

#get basemap
global_basemap <- ne_states()%>%st_transform(latlong) #this is a very large scale map
Canada <- ne_states(country = "Canada")%>%st_transform(latlong)

#full tag data
halibut_full <- read.csv("data/Acoustic/Halibut_externaldets_2024_15083_15084.csv")%>%
                filter(!is.na(sensorValue))%>%
                mutate(datetime = mdy_hm(dateCollectedUTC),
                       sensorValue = case_when(sensorType=="pressure" ~ sensorValue*-1,
                                               TRUE~sensorValue),
                       sensorType = ifelse(sensorType == "pressure","Depth (m)","Temperature °C"))


month_grid <- tibble(
  month = seq(
    floor_date(min(halibut_full$datetime), "month"),
    ceiling_date(max(halibut_full$datetime), "month"),
    by = "1 month"
  )
)

quarter_months <- c(3, 6, 9, 12)

season_grid <- month_grid %>%
  mutate(
    season = case_when(
      month(month) %in% c(12, 1, 2) ~ "Winter",
      month(month) %in% c(3, 4, 5) ~ "Spring",
      month(month) %in% c(6, 7, 8) ~ "Summer",
      month(month) %in% c(9,10,11) ~ "Fall"
    ),
    # define rectangles spanning each month
    xmin = month,
    xmax = month + months(1)
  )


season_colors <- c(
  "Winter" = "#cfe0ff",   # cool blue
  "Spring" = "#d9f0c1",   # fresh green
  "Summer" = "#ffe2b3",   # warm orange
  "Fall"   = "#ffd1dc"    # cooling pink
)

season_transitions <- month_grid %>%
  filter(month(month) %in% c(3, 6, 9, 12)) %>%
  pull(month)

# Add the last December if missing
last_dec <- ceiling_date(max(halibut_full$datetime), "year") - months(0)
if (!last_dec %in% season_transitions) {
  season_transitions <- c(season_transitions, last_dec)
}

# Labels like "Mar '23", "Jun '23"
season_labels <- format(season_transitions, "%b '%y")


# Define predation event which occurs between September and December 2024

library(dplyr)
library(lubridate)
library(ggplot2)

# Define Fall 2024 period
fall2024_start <- as.POSIXct("2024-09-01 00:00", tz = "UTC")
fall2024_end   <- as.POSIXct("2024-12-31 23:59", tz = "UTC")

# Filter points in Fall 2024 only
fall2024_points <- halibut_full %>%
  filter(datetime >= fall2024_start & datetime <= fall2024_end)

# Compute xmin/xmax with 15% buffer for x-axis
buffer_fraction <- 0.40

fall2024_rect <- fall2024_points %>%
  group_by(sensorType) %>%
  summarise(
    xmin_box = min(datetime),
    xmax_box = max(datetime),
    ymin_box = min(sensorValue, na.rm = TRUE),
    ymax_box = max(sensorValue, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # buffer x-axis
    x_range = as.numeric(difftime(xmax_box, xmin_box, units = "secs")),
    xmin_box = xmin_box - x_range * buffer_fraction,
    xmax_box = xmax_box + x_range * buffer_fraction,
    # buffer y-axis
    y_range = ymax_box - ymin_box,
    ymin_box = ymin_box - y_range * buffer_fraction,
    ymax_box = ymax_box + y_range * buffer_fraction
  )

# Plot
logger_plot <- ggplot() +
  
  # Seasonal shading
  geom_rect(
    data = season_grid,
    aes(xmin = xmin, xmax = xmax,
        ymin = -Inf, ymax = Inf,
        fill = season),
    alpha = 0.5,
    inherit.aes = FALSE
  ) +
  
  # Vertical dotted lines for seasonal transitions
  geom_vline(
    data = month_grid %>% filter(month(month) %in% quarter_months),
    aes(xintercept = month),
    linetype = "dotted",
    linewidth = 0.4,
    color = "black",
    alpha = 0.7
  ) +
  
  # Predation event
  geom_rect(
    data = fall2024_rect,
    aes(xmin = xmin_box, xmax = xmax_box,
        ymin = ymin_box, ymax = ymax_box),
    inherit.aes = FALSE,
    color = "black",
    fill = season_colors[grepl("Fall",names(season_colors))],
    linetype = 1,
    linewidth = 0.34
  ) +
  
  # Halibut points
  geom_point(
    data = halibut_full,
    aes(x = datetime, y = sensorValue),
    size = 0.8
  ) +
  
  # X-axis with seasonal transitions
  scale_x_datetime(
    breaks = season_transitions[-length(season_transitions)],
    labels = season_labels[-length(season_transitions)],
    expand = expansion(mult = c(0.01, 0.03))
  ) +
  
  scale_fill_manual(values = season_colors) +
  
  facet_wrap(~sensorType, ncol = 2, scales = "free_y") +
  
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  labs(x = "", y = "")


ggsave("output/Acoustic/Predation_event_logger.png",logger_plot,height=4.5,width=7,units="in",dpi=300)

#load the coordinates of post likley being eaten
halibut_tag_data <- read.csv("data/Acoustic/Eaten_Halibut_coords.csv")%>%
                    st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs=latlong)

halibut_release_info <- data.frame(lat=45.98667,long=-59.20944)%>%
                        st_as_sf(coords=c("long","lat"),crs=latlong)

hal_bounds <- rbind(halibut_tag_data%>%dplyr::select(geometry),
                    bioregion%>%dplyr::select(geometry))%>%
              st_buffer(5000)%>%
              st_bbox()

#get the halifax line and st anns bank coordinates

proj_start <- ymd("20190101")
proj_end <- ymd("20251212")
proj_long_upp <- -40.00
proj_long_low <- -70.00
proj_lat_upp <- 60.00
proj_lat_low <- 40.00

geoserver_receivers <- readr::read_csv('https://members.oceantrack.org/geoserver/otn/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=otn:stations_receivers&outputFormat=csv', guess_max = 13579)

otn_stations <- geoserver_receivers %>%
  filter(!is.na(deploy_date)) %>% # Remove rows with missing deploy_date values
  filter((deploy_date > proj_start & deploy_date < proj_end) | # Select stations that fall within the project timeframe defined above
           (recovery_date < proj_end & recovery_date > proj_start) |
           (deploy_date < proj_end & is.na(recovery_date) & deploy_date > proj_start - duration(18, 'months')) |
           (grepl('VR3', model) & deploy_date < proj_end & is.na(recovery_date) & deploy_date > proj_start - duration(4, 'years')) | # Select specific models within certain date ranges
           (grepl('VR4', model) & deploy_date < proj_end & is.na(recovery_date) & deploy_date > proj_start - duration(6, 'years'))) %>%
  filter(stn_lat >= proj_lat_low & stn_lat <= proj_lat_upp &# Filter stations based on latitude and longitude bounds
           stn_long >= proj_long_low & stn_long <= proj_long_upp)%>%
  st_as_sf(coords=c("stn_long","stn_lat"),crs=latlong)%>%
  mutate(year=year(deploy_date),
         array = ifelse(year<2021,"2015-2020","2020-2024"))%>%
  filter(year>2014)#just within the SAB window

hfx_stations <- otn_stations%>%
                filter(collectioncode == "HFX",
                       year%in%c(2020:2025))

sab_array <- otn_stations%>%
             filter(collectioncode=="SABMPA",
                    year%in%c(2022:2024))

#species silouettes
shark_img <- png::readPNG("output/Acoustic/Phylo_White_Shark.png")%>%
             png_to_rgba()%>%
             rasterGrob(.,interpolate = FALSE)

halibut_img <- png::readPNG("output/Acoustic/Phylo_Atlantic_halibut.png")%>%
               png_to_rgba()%>%
               rasterGrob(.,interpolate = FALSE)

#positioning 
library(png)
library(grid)
library(dplyr)

# read + convert + grob
shark_img <- png::readPNG("output/Acoustic/Phylo_White_Shark.png") %>%
  png_to_rgba() %>%
  rasterGrob(interpolate = FALSE)

halibut_img <- png::readPNG("output/Acoustic/Phylo_Atlantic_halibut.png") %>%
  png_to_rgba() %>%
  rasterGrob(interpolate = FALSE)

#positioning for the silouettes (from Chat)
# Extract x/y bounds
xmin <- hal_bounds[1]
xmax <- hal_bounds[3]
ymin <- hal_bounds[2]
ymax <- hal_bounds[4]

# placement positions
y_mid <- ymin + 0.25 * (ymax - ymin)

# vertical size of images (15% of map height)
img_height <- 0.18 * (ymax - ymin)

# recompute widths
hal_width   <- img_height * hal_ratio
shark_width <- img_height * shark_ratio

# --- horizontal placement ---

hal_xmin <- xmin + 0.25 * (xmax - xmin) + 0.15 * (xmax - xmin)
hal_xmax <- hal_xmin + hal_width

shark_xmax <- xmax - 0.05 * (xmax - xmin) + 0.10 * (xmax - xmin)
shark_xmin <- shark_xmax - shark_width

# vertical placement unchanged
ymin_img <- y_mid - img_height / 2
ymax_img <- y_mid + img_height / 2


#make map

pred_map <- ggplot()+
  geom_sf(data=global_basemap)+
  geom_sf(data=Canada,fill="grey50")+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=mar_net)+
  geom_sf(data=sab,fill="cornflowerblue")+
  geom_sf(data=hfx_stations)+
  geom_sf(data=sab_array)+
  geom_sf(data=halibut_release_info,pch=21,fill="white",size=2.5)+
  geom_sf(data=halibut_tag_data,pch=21,fill="white",size=2.5)+
  coord_sf(xlim=hal_bounds[c(1,3)],ylim=hal_bounds[c(2,4)])+
  theme_bw()+
  annotation_scale(location="br")+

annotation_custom(   #Atlantic halibut silouette
  halibut_img,
  xmin = hal_xmin, xmax = hal_xmax,
  ymin = ymin_img, ymax = ymax_img
) +
  
annotation_custom( #White shark silouette
  shark_img,
  xmin = shark_xmin, xmax = shark_xmax,
  ymin = ymin_img, ymax = ymax_img
)

ggsave("output/Acoustic/Acoustic_Predation_event.png",pred_map,height=5,width=6,units="in",dpi=300)

