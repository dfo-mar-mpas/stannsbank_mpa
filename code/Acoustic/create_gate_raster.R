create_gate_raster <- function(df_gate, spacing_m = 1000, cell_width = 1000, crs_proj = 32620) { 
  
  require(sf)
  require(dplyr)
  require(terra)
  require(purrr)
  
  # 1️⃣ Project points to meters
  df_gate <- st_transform(df_gate, crs_proj)
  
  # 2️⃣ Order points along station_id and create line
  line <- df_gate %>%
    arrange(station_id) %>%
    st_union() %>%
    st_cast("LINESTRING")
  
  # 3️⃣ Sample points along the line side-by-side (start at 0)
  line_length <- as.numeric(st_length(line))
  n_points <- ceiling(line_length / spacing_m)
  
  # fractions along the line (start at 0)
  sample_fracs <- (0:(n_points-1)) * spacing_m / line_length
  sample_fracs <- pmin(sample_fracs, 1)
  sampled_pts <- st_line_sample(line, sample = sample_fracs)
  sampled_pts <- st_cast(sampled_pts, "POINT")
  
  # 4️⃣ Create perpendicular rectangles for each sampled point
  rectangles <- map(sampled_pts, function(pt) {
    pt_coords <- st_coordinates(pt)
    line_coords <- st_coordinates(line)
    
    # Find nearest segment for tangent
    dists <- sqrt((line_coords[,1]-pt_coords[1])^2 + (line_coords[,2]-pt_coords[2])^2)
    idx <- which.min(dists)
    if(idx == nrow(line_coords)) idx <- idx - 1
    
    dx <- line_coords[idx+1,1] - line_coords[idx,1]
    dy <- line_coords[idx+1,2] - line_coords[idx,2]
    theta <- atan2(dy, dx)
    perp <- theta + pi/2
    
    half_len <- spacing_m / 2
    half_width <- cell_width / 2
    
    # Corners of polygon, closed by repeating first point
    corners <- matrix(
      c(
        pt_coords[1] - half_len*cos(theta) - half_width*cos(perp),
        pt_coords[2] - half_len*sin(theta) - half_width*sin(perp),
        
        pt_coords[1] + half_len*cos(theta) - half_width*cos(perp),
        pt_coords[2] + half_len*sin(theta) - half_width*sin(perp),
        
        pt_coords[1] + half_len*cos(theta) + half_width*cos(perp),
        pt_coords[2] + half_len*sin(theta) + half_width*sin(perp),
        
        pt_coords[1] - half_len*cos(theta) + half_width*cos(perp),
        pt_coords[2] - half_len*sin(theta) + half_width*sin(perp),
        
        # Close polygon
        pt_coords[1] - half_len*cos(theta) - half_width*cos(perp),
        pt_coords[2] - half_len*sin(theta) - half_width*sin(perp)
      ),
      ncol = 2,
      byrow = TRUE
    )
    
    st_polygon(list(corners))
  })
  
  # 5️⃣ Convert to sf and assign CRS
  rectangles_sf <- st_sf(geometry = st_sfc(rectangles, crs = crs_proj))
  
  # 6️⃣ Create raster over polygon extent
  r <- rast(ext(rectangles_sf),
            resolution = c(spacing_m, spacing_m),
            crs = st_crs(rectangles_sf)$proj4string)
  
  # 7️⃣ Rasterize polygons
  r <- rasterize(vect(rectangles_sf), r, field = 1)
  
  # 8️⃣ Return both raster and polygons
  return(list(raster = r, cells = rectangles_sf))
}