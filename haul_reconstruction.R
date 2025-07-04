# companion functions
# function to draw parallel lines
segment.shift <- function(x, y, d){
  # calculate vector
  v <- c(x[2] - x[1],y[2] - y[1])
  # normalize vector
  v <- v/sqrt((v[1]**2 + v[2]**2))
  # perpendicular unit vector
  vnp <- c( -v[2], v[1] )
  return(list(x = c( x[1] + d*vnp[1], x[2] + d*vnp[1]),
              y = c( y[1] + d*vnp[2], y[2] + d*vnp[2])))
}

# function needed to extend line
st_ends_heading <- function(line)
{
  M <- line
  i1 <- c(2, nrow(M) - 1)
  j1 <- c(1, -1)
  headings <- mapply(i1, j1, FUN = function(i1, j1) {
    Ax <- M[i1-j1,1]
    Ay <- M[i1-j1,2]
    Bx <- M[i1,1]
    By <- M[i1,2]
    unname(atan2(Ay-By, Ax-Bx))
  })
  return(headings)
}

# function to extend line
st_extend_line <- function(line, distance, end = "BOTH")
{
  if (!(end %in % c("BOTH", "HEAD", "TAIL")) | length(end) != 1) stop("'end' must be 'BOTH', 'HEAD' or 'TAIL'")
  M <- line
  keep <- !(end == c("TAIL", "HEAD"))
  ends <- c(1, nrow(M))[keep]
  headings <- st_ends_heading(line)[keep]
  distances <- if (length(distance) == 1) rep(distance, 2) else rev(distance[1:2])
  M[ends,] <- M[ends,] + distances[keep] * c(cos(headings), sin(headings))
  newline <- M
  return(newline)
}

# load packages
library(terra)
library(tidyterra)
library(sf)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(data.table)
library(lwgeom)

# projection for manipulations in meters
utm = "+proj=utm +zone=21 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# projecting linestrings in a sf object to EPSG:32721
linestrings |> 
  mutate(geometry = st_transform(geometry, crs = utm)) -> linesUTM

# extracts initial and final position of each linestring
coords <- (linesUTM %>% st_coordinates())[,1:2]

# First, move path backwards ##########
# adds separation between center of the vessel and center of the net (dis),
# and also moves the linestring 280 m back (distance between vessel and nets)
linesUTM2 <- linesUTM %>% 
  filter(!st_is_empty(geometry)) |> # filters out empty geometries
  data.frame() |> 
  dplyr::select(-geometry) |> 
  mutate(dis = case_when(barco=="Capesante" ~ 10.5,
                         barco=="ASI" ~ 10.47,
                         barco =="ASIII" ~ 10.72,
                         barco =="MT" ~ 12.25,
                         barco =="EB" ~ 11.65,
                         barco =="EBII" ~ 11.75)) %>% 
  bind_cols(data.frame(matrix(coords[,1], ncol = 2, byrow = T),
                       matrix(coords[,2], ncol = 2, byrow = T)) %>% 
              rename_with(~c("loni", "lonf", "lati", "latf"))) |> 
  mutate(Dy = latf-lati,
         Dx = lonf-loni,
         alfa = abs(atan(Dy/Dx)),
         ynewi = lati + sign(Dy)*280*sin(alfa),
         ynewf = latf + sign(Dy)*280*sin(alfa),
         xnewi = loni + sign(Dx)*280*cos(alfa),
         xnewf = lonf + sign(Dx)*280*cos(alfa))



# Second, separate nets from vessel ##########

# empty table to fill with info of each haul
offsets <- matrix(NA, nrow = nrow(linesUTM2), ncol = 8) %>%
  data.frame() %>%
  rename_with(~c("xii", "xif", "yii", "yif", "xdi", "xdf", "ydi", "ydf"))

# fill in the new coordinates of each net path
for ( i in 1:nrow(linesUTM2) ) {
  print(i)
  xs <- c(linesUTM2$ynewi[i], linesUTM2$ynewf[i])
  ys <- c(linesUTM2$xnewi[i], linesUTM2$xnewf[i])
  d <- linesUTM2$dis[i]
  new.si <- segment.shift( xs, ys, d )
  new.sd <- segment.shift( xs, ys, -d )
  offsets$yii[i] <- new.si$x[1] ; offsets$yif[i] <- new.si$x[2]
  offsets$xii[i] <- new.si$y[1] ; offsets$xif[i] <- new.si$y[2]
  offsets$ydi[i] <- new.sd$x[1] ; offsets$ydf[i] <- new.sd$x[2]
  offsets$xdi[i] <- new.sd$y[1] ; offsets$xdf[i] <- new.sd$y[2]
}
offsets2 <- offsets

# Third, add extra length to one net ##########

# define extra length for one net when both used
extra <- 30

# extends central line of fishing (with some conditionals due to differences depending on vessel)
for ( j in 1:nrow(linesUTM2) ) {
  print(j)
  if(linesUTM2$bolsa_b[j] <= 1 & linesUTM2$bolsa_e[j] <= 1){ #if both nets used
    line <- if(linesUTM2$barco[j] == "ASIII"){ # define which side to extend
      matrix(c(offsets[j,1],offsets[j,3], offsets[j,2],offsets[j,4]), ncol=2, byrow = T)} else{
      matrix(c(offsets[j,5],offsets[j,7], offsets[j,6],offsets[j,8]), ncol=2, byrow = T)
    }
    a <- st_extend_line(line, extra) # extend line
    if(linesUTM2$barco[j] == "ASIII"){offsets2[j, 1:4] <- c(a[,1], a[,2])} else {
      offsets2[j, 5:8] <- c(a[,1], a[,2])
    }
  }
}

# creates a linestring between initial and final positions of the new lines
dti <- as.data.table(data.frame(offsets2[,1:4], id=linesUTM2$id))
dti1 <- dti[, .(id, lon = xii, lat = yii)]
dti2 <- dti[, .(id, lon = xif, lat = yif)]
dti1[, seq := 1L ]
dti2[, seq := 2L ]
dti <- rbindlist(list(dti1, dti2), use.names = TRUE)
setorder(dti, id, seq)
sfi <- sfheaders::sf_linestring(
  obj = dti
  , x = "lon"
  , y = "lat"
  , linestring_id = "id"
)
st_crs(sfi)<-utm

dtd <- as.data.table(data.frame(offsets2[,5:8], id=linesUTM2$id))
dtd1 <- dtd[, .(id, lon = xdi, lat = ydi)]
dtd2 <- dtd[, .(id, lon = xdf, lat = ydf)]
dtd1[, seq := 1L ]
dtd2[, seq := 2L ]
dtd <- rbindlist(list(dtd1, dtd2), use.names = TRUE)
setorder(dtd, id, seq)
sfd <- sfheaders::sf_linestring(
  obj = dtd
  , x = "lon"
  , y = "lat"
  , linestring_id = "id"
)
st_crs(sfd)<-utm

# filter out lines where there was no fishing
sfd <- sfd %>% filter(!id %in%linesUTM2$id[which(linesUTM2$bolsa_e == 9 | (linesUTM2$bolsa_e == 999 & linesUTM2$bolsa_b != 999 & linesUTM2$Simple..Doble == "Simple"))])
sfi <- sfi %>% filter(!id %in%linesUTM2$id[which(linesUTM2$bolsa_b == 9 | (linesUTM2$bolsa_b == 999 & linesUTM2$bolsa_e != 999 & linesUTM2$Simple..Doble == "Simple"))])

# combine lines of the same fishing event
veremos<-bind_rows(sfi,sfd)
v2 <- veremos %>% group_by(id) %>% summarise(geometry= st_combine(geometry))

# Fourth, add buffer to turn the linestring into a polygon ##########
buff <- st_buffer(v2, ifelse(linesUTM2$barco=="EBII", 8.8, 7), endCapStyle = "FLAT") 

# get back to geographic coordinates
st_transform(buff, 4326)->b2
