


clean.data <- function(matrix, total_counts = NULL, sim_doublet_ratio = 2.0, n_neighbors = NULL, expected_doublet_rate = 0.075, random_state = 0L) {
  tmp <- scrublet(matrix = matrix, 
                  total_counts = total_counts, 
                  sim_doublet_ratio = sim_doublet_ratio, 
                  n_neighbors = n_neighbours, 
                  expected_doublet_rate = expected_doublet_rate, 
                  stdev_doublet_rate = stdev_doublet_rate, 
                  random_state = random_state)
  tmp <- decontaminate(matrix = tmp)
  return(tmp)
}
