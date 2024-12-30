library(tidyverse)

sample_within_cell <- function(centroid_coords, dx, dy, num_of_samples){
  
  # Function to draw samples uniformly within a cell
  # - centroid_coords: 2d vector with x and y coordinates of the cell centroid
  # - dx: cell width in the x direction
  # - dy: cell width in the y direction
  # - num_of_samples: Number of samples to be drawn
  # returns: a dataframe with columns x and y, representing the x and y coordinates
  # of the samples
 
  x <- centroid_coords[1]$x
  y <- centroid_coords[2]$y
  
  random_x <- runif(num_of_samples, -dx, dx)
  random_y <- runif(num_of_samples, -dy, dy)
  
  x_coords <- x + 1 / 2 * random_x #Factor 1/2 because things are relative to the
  #centroid
  y_coords <- y + 1 / 2 *random_y
  
  coords <- data.frame(x = x_coords, y = y_coords)
  
  return(coords)
}

sample_locations <- function(sampled_cells, exp_grid, num_of_grid_points=200){
  # Function that finds the coordinates of the sampled data locations  
  #   -sampled_cells: Vector with indices of sampled cells
  #   -exp_grid: Grid that realization is defined over
  # 
  # returns: Coordinates of data points, and the amount of data points
  # that are sampled within each cell (including zero if there are no
  #                                    data points within a cell)
  
  cell_counts<-data.frame(Cell=sampled_cells) %>%
    count(Cell) %>%
    complete(Cell = 1:(num_of_grid_points^2), fill = list(n = 0)) %>%
    rename(Count = n)
  
  sampled_locations <- data.frame(x = numeric(), y = numeric())
  dx<-diff(exp_grid[, 1])[1] #Assuming quadratic grid cells
  
  #Note that the following for loop follows the ordering of the cells
  for(i in 1:(num_of_grid_points^2)){
    times_cell_sampled <- cell_counts[i, 2]$Count
    if(times_cell_sampled>0){
      sampled_location<- sample_within_cell(exp_grid[i, ],dx, dx, 
                          num_of_samples = times_cell_sampled)
      
      sampled_locations <-rbind(sampled_locations, sampled_location)
    }
  }
  return(list(sampled_locations=sampled_locations,
              cell_counts = cell_counts))
}

generate_data <- function(realization, sampling_density, exp_grid,
                          num_of_data_points = 2000,
                          num_of_grid_points = 200){
  # Function that generates a discretize dataset, sampled according
  # to a sampling density, given a realization over some grid
  #   -realization: The realization to sample from
  #   -sampling_density: Vector of sampling density values evaluated at the
  # grid points 
  #   -exp_grid: The grid over which the realization is generated
  #   -num_of_data_points: The amount of datapoints in the discretized dataset
  # 
  # returns: The disretized dataset (a num_of_data_points by 3) dataframe with
  # columns for the x, y and z values
  
  #This function has been updated to take in a sampling density instead
  #of the sampled_cells, because we will have to sample the cells anyhow.
  
  sampled_cells <- sample(1:num_of_grid_points^2, size = num_of_data_points, 
                          replace = TRUE, prob = sampling_density)
  
  sample_results <- sample_locations(sampled_cells, exp_grid)
  
  sampled_locations <- sample_results[[1]]
  
  cell_counts <- sample_results[[2]]
  
  df <- exp_grid
  df$z <- as.vector(realization)
  
  sampled_locations$z <- rep(df$z, cell_counts$Count) #Populates the
  #z values of the data points
  
  return(sampled_locations)
}
