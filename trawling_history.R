# used packages
library(sf)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(marmap)
library(tidyterra)
library(ggspatial)
library(terra)

# load data
todas <- readRDS("todas_final.RDS") |> 
  mutate(yr = substr(fecha, 1, 4))

prod <- read.csv("prod_final.csv")

# trawling per year
spatraster <- rast(extent = ext(st_bbox(todas)), resolution = 1/3600)
rasterize(vect(todas), spatraster, fun = "count", by = todas$yr) -> rasterizado_fr

# PRODUCTION ###########
fechas <- unique(prod$fecha)
# dates with production info
todas |> 
  filter(fecha %in% fechas) |> 
  mutate(fecha = substr(fecha, 1, 10)) |> 
  group_by(fecha, fk_id_marea) |> 
  mutate(geometry = st_centroid(geometry)) |> 
  summarise(geometry = st_combine(geometry), .groups = "drop") |> 
  mutate(geometry = st_centroid(geometry)) -> fechas_con_prod

prod |> 
  left_join(fechas_con_prod, by=c("fk_id_marea", "fecha")) -> prod_x_marea

# daily production
prod_x_marea |> 
  filter(!st_is_empty(geometry)) |> st_as_sf() |> 
  mutate(yr = as.numeric(substr(fecha, 1, 4)))-> d

# raster with trawled areas
rasterize(vect(todas), spatraster, fun = "min", by = todas$yr) -> rasterizado_fr_pa
# calculate the intensity of previous trawling on areas where production is available
out<-NULL
for(i in unique(d$yr)){
  chosen_year <- i
  end_index <- chosen_year - 1995 + 1
  r_sum <- sum(rasterizado_fr_pa[[1:end_index]], na.rm = T)
  bind_cols(d |> data.frame() |> filter(yr == i) |> select(-geometry), intensidad = extract(r_sum, d |>  filter(yr == i))$sum) -> tmp
  out<-bind_rows(out, tmp)
}

# figure
ggplot(out |> 
         mutate(g = ifelse(intensidad > 1, "Yes", "No")) |> 
         filter(!is.na(g)), aes(g, tons))+
  geom_hline(yintercept = 5, color = "lightgray", alpha = 0.5) +
  geom_boxplot(coef = 9999, aes(color = as.factor(intensidad)))+
  # geom_point()+
  stat_summary(shape = 21, alpha = 0.75)+
  facet_wrap(.~yr)+
  theme_cowplot()+
  labs(x = "Previously trawled?", y = "Daily production (t)")+
  scale_color_discrete(name = "Number of years trawled")+
  theme(legend.position = "top")+ 
  guides(colour = guide_legend(nrow = 1))

# t-tests
out |> 
  mutate(g = ifelse(intensidad > 1, "Yes", "No")) |> 
  filter(!is.na(g)) -> tt

testst<-NULL
for(i in unique(tt$yr)){
  prueba <- t.test(tons~g, data = tt |> filter(yr == i))
  data.frame(yr = i, t = round(prueba$statistic, digits = 3), df = round(prueba$parameter, digits = 2), p = round(prueba$p.value, digits = 3)) -> tmp
  bind_rows(testst, tmp) -> testst
}

# table
library(gt)
library(gtExtras)

testst %>%
  mutate(year = yr) |> 
  select(year, t, df, p) |> 
  gt() %>%
  gt_theme_dot_matrix() |> 
  gtsave("testst.png")

# supplementary figure
out |> 
  group_by(fk_id_marea, yr) |> 
  mutate(g = ifelse(is.na(intensidad), 0, 1)) |> 
  group_by(fk_id_marea, g, yr) |> 
  summarise(p = mean(tons)) |> 
  filter(fk_id_marea %in% cods) |> 
  ggplot(aes(as.factor(g), p))+
  geom_boxplot(coef = 9)+
  geom_point(position = position_jitter(width = 0.1), aes(fill = as.factor(yr)), shape = 21 ,alpha = 0.6, size = 2)+
  # stat_summary()+
  theme_cowplot()+
  labs(x = "", y = "Daily production (t)")+
  scale_x_discrete(labels = c("Unreliable position", "Precise position"))+
  # guides(colour = guide_legend(nrow = 1))+
  scale_fill_discrete(name = "")


# OVERALL TRAWLING EFFORT ##########
# total trawling
rasterize(vect(todas), spatraster, fun = "count") -> rasterizado

# area trawled by number of visits
expanse(rasterizado, byValue = T) -> expandido

# sum of trawled area
sum(expandido$area*expandido$value) -> total_raster

# percent seafloor trawled only once
expandido$area[1]/sum(expandido$area)*100

# sefloor area trawled
sum(expandido$area)

# HOW MUCH EFFORT ##########
# trawled surface by year
spatraster_año <- rast(extent = ext(st_bbox(todas)), resolution = 1/3600)
sup_x_año <- NULL
for(i in sort(unique(todas$yr))){
  print(i)
  todas |> filter(yr == i) -> año
  rasterize(vect(año), spatraster_año, fun = "count") -> rasterizado_año
  expanse(rasterizado_año, byValue = T) -> expandido_año
  data.frame(value = expandido_año$value, area = expandido_año$area, yr = i) -> tmp
  sup_x_año <- bind_rows(sup_x_año, tmp)
}

# figure
sup_x_año |> 
  mutate(tot = area*value) |> 
  group_by(yr) |> 
  summarise(total = sum(tot)) |> 
  mutate(acumulado = cumsum(total)) |> 
  ggplot(aes(as.numeric(yr), total / 1000000*20))+
  geom_bar(fill = "gray70", aes(x = as.numeric(yr), y = acumulado / 1000000), stat = "identity", inherit.aes = F)+
  geom_point(size = 3)+
  geom_line() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), sec.axis = sec_axis(transform=~./20, name=expression(paste("Annual surface (", km^2, ")")))) +
  scale_x_continuous(breaks = seq(1995, 2023, by = 2))+
  labs(x = "", y = expression(paste("Accumulated surface (", km^2, ")"))) +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# WHICH SPATIAL DISTRIBUTION ##########
# estimation of overlapping and non-overlapping areas per year
sup_x_año %>%
  pivot_wider(names_from = value, values_from = area) %>% 
  replace(is.na(.), 0) %>%
  mutate(no_solap = rowSums(across(`1`:ncol(.))), 
         solap =rowSums(across(`2`:ncol(.))),
         yr_i = as.numeric(yr),
         yr_j = yr) |> 
  select(yr_i, yr_j, no_solap, solap) -> intra

# cumulative seafloor area trawled across years
spatraster_año <- rast(extent = ext(st_bbox(todas)), resolution = 1/3600)
salida_sol <-NULL
for(i in sort(unique(todas$yr))){
  print(i)
  todas |> filter(yr %in% 1995:as.numeric(i)) -> año
  rasterize(vect(año), spatraster_año, fun = "count") -> rasterizado_acum
  data.frame(yr = i, ar_sol = expanse(rasterizado_acum)$area) -> dt
  salida_sol <- bind_rows(salida_sol, dt)
}

# figure
left_join(intra |> rename(yr = yr_i) |> select(-yr_j), salida_sol |> mutate(yr = as.numeric(yr))) |> arrange(yr) |> 
  mutate(ar = no_solap/1000000, # area not overlapped per year
         ar_sol = ar_sol/1000000, # cumulative seafloor surface
         sup_nueva=ar_sol-lag(ar_sol, default = 0), # new seafloor surface per year
         prop_nueva = sup_nueva/ar*100, # same as %
         yr = as.numeric(yr)) |> 
  left_join(sup_x_año |> 
              mutate(yr = as.numeric(yr),
                     tot = area*value) |> 
              group_by(yr) |> 
              summarise(total = sum(tot)) |> 
              mutate(acumulado = cumsum(total))) |>
  mutate(prop_solap = ((acumulado/1000000)- ar_sol)/(acumulado/1000000),
         prop_nuevo = ar_sol *1000000 / acumulado) |>
  ggplot()+
  geom_line(aes(yr, prop_nueva*500), size = 1, color = "#DC322F")+
  geom_line(aes(yr, acumulado/1000000), size = 1, color = "#073642")+ # barrido totales con sf
  geom_line(aes(yr, ar_sol), size = 1, color = "#268BD2")+
  geom_point(aes(yr, acumulado/1000000), size = 3, fill = "#073642", shape = 21)+
  geom_point(aes(yr, ar_sol), size = 3, fill = "#268BD2", shape = 21)+
  theme_cowplot()+
  geom_ribbon(aes(x = yr, ymin=ar_sol,ymax=acumulado/1000000), fill="blue", alpha=0.15) +
  geom_label(aes(2018, max(acumulado/1000000), label = "All trawls"), size = 6, color = "#073642")+
  geom_label(aes(2014.25, 25000, label = "Overlap"), size = 6, color = "mediumpurple3")+
  geom_label(aes(2019, max(ar_sol)-12000, label = "Accounting\ntrawl overlap"), size = 6, color = "#268BD2")+
  geom_label(aes(2002, max(ar_sol)+32000, label = "New seafloor"), size = 6, color = "#DC322F")+
  scale_y_continuous(sec.axis = sec_axis(~ . / 500, name = "New seafloor (%)"))+
  labs(x = "", y = bquote("Surface"~(km^2)))+
  theme(
    axis.title.y.right = element_text(colour = "#DC322F"),
    axis.line.y.right = element_line(color = "#DC322F"),
    axis.ticks.y.right = element_line(color = "#DC322F"),
    axis.text.y.right = element_text(color = "#DC322F"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)) -> f5

# 614 x 452

ggplot(out |> 
         mutate(g = ifelse(intensidad > 1, "Yes", "No")) |> 
         filter(!is.na(g)), aes(g, tons, color = g))+
  geom_boxplot(coef = 9999)+
  geom_point(alpha = 0.25, shape = 21, position = position_jitter(0.25))+
  theme_cowplot()+
  scale_color_manual(values = c("red", "gray30"))+
  labs(x = "Previously trawled?", y = "Daily production (t)")+
  theme(legend.position = "none", axis.title = element_text(size = 18),
        axis.text = element_text(size = 14))+ 
  guides(colour = guide_legend(nrow = 1)) ->f5b

plot_grid(f5, f5b, ncol = 1, labels = "AUTO", align = "vh")

# test 
out |> 
  mutate(g = ifelse(intensidad > 1, "Yes", "No"),
         fk_id_marea= as.factor(fk_id_marea)) |> 
  filter(!is.na(g)) -> prods

library(DHARMa)
library(glmmTMB)
mod <- glmmTMB(tons ~g, prods, ziformula = ~.)
mod2 <- glmmTMB(tons ~1, prods, ziformula = ~.)
anova(mod, mod2)
plot(simulateResiduals(mod2))


# FREQUENCY OF TRAWLING ###########
# main fishing ground area
bbox_coords <- st_bbox(c(
  xmin = -61,  # Minimum longitude
  ymin = -45.5, # Minimum latitude
  xmax = -54.3, # Maximum longitude
  ymax = -36.7  # Maximum latitude
), crs = st_crs(4326)) # CRS 4326 (WGS 84)

spatraster_año <- rast(extent = ext(bbox_coords), resolution = 1/3600)
raster_extent <- as.polygons(ext(spatraster_año), crs = crs(spatraster_año))

# keep only trawls within the main fishing ground
dentro <- todas[st_within(todas, st_as_sf(raster_extent), sparse = F), ]

100-nrow(dentro)/nrow(todas)*100 # % trawls lost due to this subset

# new raster
rasterize(vect(dentro), spatraster_año, fun = "count") -> rasterizado_punto3
gc() # free memory to prevent crash

# data frame with cell id (to detect new visits)
total_punto3 <- as.data.frame(rasterizado_punto3, cells = TRUE, na.rm = TRUE)
gc() # free memory to prevent crash

# number of cells visited
total_celdas_punto3 <- nrow(total_punto3)
length(unique(total_punto3$cell)) -> con_visita
gc()

# cell visits per year
revisitas_x_año<- list()
for(i in sort(unique(dentro$yr))){
  print(i)
  dentro |> filter(yr == i) -> año
  rasterize(vect(año), spatraster_año, fun = "min") -> rasterizado_año
  gc()
  # Combine all tile results into a single data frame
  final_results_tmp <- as.data.frame(rasterizado_año, cells = TRUE, na.rm = TRUE)
  rm(rasterizado_año)
  gc()
  data.frame(cell = final_results_tmp$cell, yr = i) -> revisitas_x_año[[as.numeric(i) - 1994]]
  rm(final_results_tmp)
  gc()
  # saveRDS(revisitas_x_año, "revisitas_x_año.rds")
  # gc()
}

# join loop results
combined_data <- do.call(rbind, revisitas_x_año)

# Step 2: Sort by id and year
combined_data <- combined_data[order(combined_data$cell, combined_data$yr), ]
gc()

# add time lag to those visited more than one year
combined_data |> 
  group_by(cell) |> 
  mutate(yr = as.numeric(yr),
         lag = lag(yr),
         lag_yr = yr-lag) |> 
  filter(!is.na(lag_yr)) -> punto3


# cells visited more than once
length(unique(punto3$cell)) -> con_visitas

# number of years visited per cell
combined_data |> 
  group_by(cell) |> 
  summarise(n = n()) -> pepe

# panel A
pepe |> 
  ggplot()+
  geom_bar(aes(x = as.character(n), y = after_stat(count)/con_visita*100, fill = ifelse(n == "1", "highlight", "normal")))+
  scale_x_discrete(limits = factor(1:16))+
  scale_fill_manual(values = c("highlight" = "red", "normal" = "gray30")) +
  labs(x = "Number of years trawled", y = "Percent of total area trawled")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_cowplot() +
  theme(legend.position = "none")-> externa

# panel B
pepe |> 
  filter(n == 1) |> 
  left_join(combined_data) |> 
  group_by(yr) |> 
  summarise(n = n()) |> 
  ggplot(aes(x = as.numeric(yr), y = n/(con_visita-con_visitas)*100))+
  geom_line(color = "red")+
  geom_hline(linetype = 2, aes(yintercept = mean(n/(con_visita-con_visitas)*100))) + 
  geom_point(size = 3, color = "red")+
  theme_cowplot() + 
  labs(x = "", y = "Percent of single-year visits") -> inset

# upper panel
ggdraw(externa) +
  draw_plot(inset, .35, .35, .6, .6) +
  draw_plot_label(
    c("A", "B"),
    c(0, 0.45),
    c(1, 0.95),
    size = 12
  ) -> paneles_ab

# mean number of visits per cell
nrow(punto3)/length(unique(punto3$cell))
100-con_visitas/con_visita*100 # cells visited only once

# lower panel
punto3 |> 
  ggplot()+
  geom_bar(aes(x = as.character(lag_yr), y = after_stat(count)/total_celdas_punto3*100), fill = "gray30") +
  scale_x_discrete(limits = factor(1:30))+
  labs(x = "Time lag between trawlings (yrs)", y = "Percent of trawling on repeated locations")+
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))+
  theme_cowplot() -> panel_c

# figure
plot_grid(paneles_ab, panel_c, ncol = 1, labels = c("", "C"))

# 736 x 869

# proportion of visits with time lag > 4
punto3 |> 
  group_by(lag_yr) |> 
  summarise(n = n()) |> 
  filter(lag_yr > 4) |> 
  pull(n) |> sum() / nrow(punto3)

# MAP #############

bat <- getNOAA.bathy(-70, -54, -57, -32, res = 1, keep = TRUE)
bat_xyz <- as.xyz(bat)

rasterizado <- rast("rasterizado.tif")

ggplot()+
  stat_contour(data = bat_xyz, aes(x = V1, y = V2, z = V3, color= after_stat(level)), breaks = c(75, 150, 300, 600, 1200, 2400, 4800)*-1) +
  scale_color_continuous(name = "Depth (m)")+
  annotation_map(map_data("world")) +
  geom_sf(data = udms, alpha =0.5) + 
  geom_spatraster(data = crop(rasterizado, raster::extent(-61, -54.5, -45.25, -37)), maxcell = 50000000, aes(fill = after_stat(value)))+  
  scale_fill_viridis_c(option = "inferno", direction = -1, na.value = NA, name = "Fishing tows")+
  coord_sf(xlim=c(-61,-54.5),ylim=c(-45.25,-37))+
  ggOceanMaps::theme_map(grid.size = 0.1, grid.col = "gray")+
  theme(legend.position = "right")+
  geom_rect(aes(xmin=-59,xmax = -58,ymin=-42.5,ymax = -41.5), fill = NA, linewidth = 1, color = scales::alpha("red", 0.5))+
  scale_x_continuous(breaks = c(-61,-59,-57, -55))+
  geom_rect(aes(xmin=-61,xmax = -54.5,ymin=-45.25,ymax = -37), fill = NA, linewidth = 1, color = scales::alpha("black", 0.5))+
  labs(x = "Longitude", y = "Latitude") -> resultado


ggplot()+
  annotation_map(map_data("world")) +
  stat_contour(data = bat_xyz, aes(x = V1, y = V2, z = V3, color= after_stat(level)), breaks = c(75, 150, 300, 600, 1200, 2400, 4800)*-1) +
  geom_sf(data = udms, alpha =0.5) + 
  scale_color_continuous(name = "Depth (m)")+
  geom_spatraster(data = rasterizado, maxcell = 50000000, aes(fill = after_stat(value)))+  coord_sf(xlim=c(-70,-54.5),ylim=c(-57,-32))+
  scale_fill_viridis_c(option = "inferno", direction = -1, na.value = NA, name = "Fishing tows")+
  ggOceanMaps::theme_map(grid.size = 0.1, grid.col = "gray")+
  labs(x = "Longitude", y = "Latitude")+
  geom_rect(aes(xmin=-61,xmax = -54.5,ymin=-45.25,ymax = -37), fill = NA, linewidth = 1, color = scales::alpha("black", 0.5))+
  scale_x_continuous(breaks = c(-69,-62,-55))+
  theme(legend.position = "none")-> referencia


ggplot()+
  annotation_map(map_data("world")) +
  stat_contour(data = bat_xyz, aes(x = V1, y = V2, z = V3, color= after_stat(level)), breaks = c(75, 150, 300, 600, 1200, 2400, 4800)*-1) +
  geom_sf(data = udms, alpha =0.5) + 
  scale_color_continuous(name = "Depth (m)")+
  geom_spatraster(data = crop(rasterizado, raster::extent(-59, -58, -42.5, -41.5)), maxcell = 50000000, aes(fill = after_stat(value)))+  
  coord_sf(xlim=c(-59,-58),ylim=c(-42.5,-41.5))+
  scale_fill_viridis_c(option = "inferno", direction = -1, na.value = NA, name = "Fishing tows")+
  ggOceanMaps::theme_map(grid.size = 0.1, grid.col = "gray")+
  labs(x = "Longitude", y = "Latitude")+
  theme(legend.position = "none")+
  geom_rect(aes(xmin=-59,xmax = -58,ymin=-42.5,ymax = -41.5), fill = NA, linewidth = 1, color = scales::alpha("red", 0.5))+
  scale_x_continuous(breaks = c(-59,-58.5,-58))+
  scale_y_continuous(breaks = c(-42.4,-42,-41.6))-> zoom


plot_grid(plot_grid(referencia, zoom, ncol = 1, labels = c("A", "C")), resultado, nrow = 1, rel_widths = c(0.4,0.6), labels = c("", "B"))
# 1051*818
