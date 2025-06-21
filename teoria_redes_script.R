setwd("tu_directorio_de_trabajo")



###### librerias ########
library(dplyr) ##para manejar datos
library(ggplot2) ##para graficar
library(sf) ## para manejar datos espaciales
library(ggspatial) ##para graficar
library(hrbrthemes) ##para graficar
library(ggthemes) ##para graficar
library(gtools)
library(terra) ##manejar rasters
library(velociraptr) ##descargar datos fosiles desde paleobiology database
library(rgplates) ### reconstruir coordenadas paleogeograficas
library(purrr) ##### calcular riqueza de especies
library(smoothr) #### calcular riqueza de especies
library(vegan) ##manejar datos ecologicos
library(glue)
library(viridis) ##colores
library(igraph) ##teoria de redes


## vamos a descargar datos de caballos fosiles (familia Equidae) y sus relativos del Mioceno
## a nivel de Genero (23.03 a 5.33 Ma atras)
caballo_miocene=downloadPBDB(Taxa="Equidae",StartInterval="Miocene",StopInterval="Miocene")

## quitamos ocurrencias que no llegan a genero
caballo_miocene=subset(caballo_miocene,!caballo_miocene$genus=="")
caballo_miocene<-caballo_miocene[,-c(16:18)]

### veamos cuantos generos tenemos
collections_miocene=tapply(caballo_miocene$collection_no,list(caballo_miocene$collection_no,caballo_miocene$genus),length)

collections_miocene[is.na(collections_miocene)]=0
collections_miocene[collections_miocene>1]=1
ngenus_miocene<-ncol(collections_miocene)

##tenemos 32 generos de caballos para el Mioceno
print(ngenus_miocene)

## veamos la completitud de nuestro ensamble (aproximacion tipo presencia-ausencia por coleccion)


specpool(collections_miocene) ##basado en CHAO 2 abrian 47 generos totales para el Mioceno

print(glue("en total tenemos {specpool(collections_miocene)[1]} generos para el Mioceno,
           y una completitud de ensambles del {round((specpool(collections_miocene)[1]/specpool(collections_miocene)[2])*100,2)}%"))

### vamos acalcular las paleolatitudes de nuestras pcurrencias ####
web<-getgws() ##llink a gplates

setgws(web) ##conectar a gplates

miocene_coordinates<-data.frame(lon=caballo_miocene$lng,lat=caballo_miocene$lat)

### ponemos la reconstruccion en 20 millones de años
### si decides cambiar el codigo a distintas edades del Mioceno
### puede que te de resultados ligeramente diferentes a los mios
### no es lo mas adecuado poner esta edad ya que el abanico de edades
### de las ocurrencias es muy grande, pero por fines de experimentacion
### y demostrativos plotearemos todo a los 20 millones de anos

miocene_coordinates<-as.data.frame(reconstruct(miocene_coordinates,age=20))

## juntamos datos
caballo_miocene=cbind(caballo_miocene[,c(8)],miocene_coordinates)

caballo_miocene=na.omit(caballo_miocene) ##quitamos ocurrencias que no se pudieron reconstruir

colnames(caballo_miocene)[1]<-c("Genus") ##arreglar nombre de columnas

####cargamos paleomapas
miocene_poly = st_read("paleo_map/20Ma_CM_v7.shp")

####mapa de ocurrencias
miocene_points <- ggplot() + 
  geom_sf(data = miocene_poly, fill = "grey70", col = 'black', lwd = 0.2)+ 
  geom_point(data = caballo_miocene, aes(x = paleolong, y = paleolat),col="#8B2323")+
  labs(x = '', y = '') + ##Nombres de los ejes x e y
  ggtitle(label = 'Sampling Points') + ##titulo del mapa
  theme_ipsum_es() + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5), 
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom', 
        legend.key.width = unit(3.5, 'line'),
        panel.background = element_rect(fill = "white"),
        plot.subtitle = element_text(face = 'italic', hjust = 0.5)) +
  annotation_scale(location = "bl")

miocene_points
####pasamos nuestros puntos a geometrias
mio.caballo_sf <- st_as_sf(caballo_miocene, coords = c("paleolong", "paleolat"), crs = 4326)

### creamos un buffer
mio.buff <- sf::st_buffer(mio.caballo_sf, dist = 500000) |> 
  sf::st_union() |> 
  sf::st_sf() |>
  sf::st_transform(crs = terra::crs(mio.caballo_sf))

### creamos poligonos para cada genero
mio.caballo_poly <- mio.caballo_sf %>%
  group_by(Genus) %>%
  summarize(geometry = st_union(geometry)) %>%
  st_convex_hull()

### creamos una grilla espacial donde vamos a calcular la riqueza de generos por celda
### de acuerdo a la spuerposicion espacial de los generos de caballos
### en este caso es una grilla de 2 grado x 2 grado de latitud
### tambien puedes modificar esto e ir experimentando

mio.grid <- st_make_grid(mio.caballo_sf, cellsize = c(2, 2))

# Convertir la grilla en un objeto sf con IDs únicos
mio.grid <- st_sf(geometry = mio.grid)
mio.grid$cell_id <- paste0("Cell_", 1:nrow(mio.grid))

# Filtrar celdas que intersectan con el buffer
mio.buff <- st_make_valid(mio.buff)
buffer_cells <- st_intersects(mio.grid, mio.buff, sparse = FALSE)
mio.grid <- mio.grid[rowSums(buffer_cells) > 0, ]

##filtrar por ambiente terrestre
terrestrial_cells <- st_intersects(mio.grid, miocene_poly, sparse = FALSE)
mio.grid <- mio.grid[rowSums(terrestrial_cells) > 0, ]


# Verificar el número de celdas
cat("Número de celdas:", nrow(mio.grid), "\n")

##### ver nuestra grilla #######

###ver nuestra grilla espacial, buffer y ocurrencias
ggplot() +
  geom_sf(data = miocene_poly, fill = "grey70", col = "black", lwd = 0.2) +
  geom_sf(data = mio.grid, fill = NA, col = "blue", lwd = 0.5) +
  geom_sf(data = mio.buff, fill = NA, col = "forestgreen", lwd = 0.5, alpha = 0.3) +
  geom_point(data = caballo_miocene, aes(x = paleolong, y = paleolat), col = "brown") +
  labs(x = "", y = "") +
  ggtitle("Grilla restringida al buffer") +
  theme_ipsum_es() +
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white")
  ) +
  annotation_scale(location = "bl")

# Crear matriz presencia-ausencia
intersections <- st_intersects(mio.grid, mio.caballo_poly, sparse = FALSE)
rownames(intersections) <- mio.grid$cell_id
colnames(intersections) <- mio.caballo_poly$Genus
incidence_matrix <- t(intersections) * 1  # Géneros como filas, celdas como columnas

###red biogeografica
g <- graph_from_biadjacency_matrix(incidence_matrix, directed = FALSE)
V(g)$type <- c(rep(FALSE, nrow(incidence_matrix)), rep(TRUE, ncol(incidence_matrix)))
V(g)$name <- c(rownames(incidence_matrix), colnames(incidence_matrix))

# Calcular modularidad
communities <- cluster_infomap(g)
mod <- modularity(g, membership(communities))
cat("Modularidad de la red biogeográfica:", mod, "\n")
V(g)$community <- membership(communities)

# Asignar comunidades a la grilla
cell_communities <- data.frame(
  cell_id = V(g)$name[V(g)$type == TRUE],
  community = V(g)$community[V(g)$type == TRUE]
)
mio.grid <- left_join(mio.grid, cell_communities, by = "cell_id")

# Visualizar unidades biogeográficas
mapa_units<-ggplot() +
  geom_sf(data = miocene_poly, fill = "grey70", col = "black", lwd = 0.2) +
  #geom_sf(data = mio.buff, fill = NA, col = "green", lwd = 0.5, alpha = 0.3) +
  geom_sf(data = mio.grid, aes(fill = as.factor(community)), alpha = 0.5) +
  #geom_point(data = caballo_miocene, aes(x = paleolong, y = paleolat), col = "#8B2323") +
  scale_fill_viridis_d(name = "Unidad Biogeográfica") +
  labs(x = "", y = "") +
  ggtitle("Unidades biogeográficas de los caballos del Mioceno en 2D (infomap)") +
  theme_ipsum_es() +
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white")
  ) +
  annotation_scale(location = "bl")
ggsave(plot = mapa_units, filename = glue('miocene_units.jpg'), units = 'in', width = 10, height = 7, dpi = 300)


############# tambien podemos plotear la red biogeografica ##########
plot(g,
     vertex.color=V(g)$community,
     vertex.label=NA,
     layout=layout_with_fr,
     vertex.size=7)

plot(g,
     vertex.color=V(g)$community,
     vertex.label=NA,
     layout=layout_nicely,
     vertex.size=7)

#####podemos calcular otras metricas#########
# Calcular betweenness#########
V(g)$betweenness <- betweenness(g)
V(g)$betweenness
cell_nodes <- which(V(g)$type == TRUE)
centrality_data <- data.frame(
  cell_id = V(g)$name[cell_nodes],
  betweenness = V(g)$betweenness[cell_nodes]
)

# Unir centralidad a la grilla
mio.grid <- left_join(mio.grid, centrality_data, by = "cell_id")

ggplot() +
  geom_sf(data = miocene_poly, fill = "grey70", col = "black", lwd = 0.2) +
  geom_sf(data = mio.grid, aes(fill = betweenness), alpha = 0.7) +
  scale_fill_viridis(option="viridis") +
  labs(x = "", y = "") +
  ggtitle("Betweeness") +
  theme_ipsum_es() +
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white")
  ) +
  annotation_scale(location = "bl") 

##closeness######
V(g)$closeness <- closeness(g)

closeness_data <- data.frame(
  cell_id = V(g)$name[cell_nodes],
  closeness = V(g)$closeness[cell_nodes]
)

# Unir a la grilla
mio.grid$closeness<-closeness_data$closeness/max(closeness_data$closeness)

max(mio.grid$closeness)
ggplot() +
  geom_sf(data = miocene_poly, fill = "grey70", col = "black", lwd = 0.2) +
  geom_sf(data = mio.grid, aes(fill = closeness), alpha = 0.7) +
  scale_fill_viridis(option="viridis") +
  labs(x = "", y = "") +
  ggtitle("Closeness normalizado") +
  theme_ipsum_es() +
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white")
  ) +
  annotation_scale(location = "bl") 

###eigen centrality#####
V(g)$eigen_centrality <- eigen_centrality(g)$vector

eigen_centrality_data <- data.frame(
  cell_id = V(g)$name[cell_nodes],
  eigen_centrality = V(g)$eigen_centrality[cell_nodes]
)

# Unir a la grilla

mio.grid$eigen_centrality <- eigen_centrality_data$eigen_centrality/max(eigen_centrality_data$eigen_centrality)
ggplot() +
  geom_sf(data = miocene_poly, fill = "grey70", col = "black", lwd = 0.2) +
  geom_sf(data = mio.grid, aes(fill = eigen_centrality), alpha = 0.7) +
  scale_fill_viridis(option="viridis") +
  labs(x = "", y = "") +
  ggtitle("eigen centrality normalizado") +
  theme_ipsum_es() +
  theme(
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.background = element_rect(fill = "white")
  ) +
  annotation_scale(location = "bl") 
