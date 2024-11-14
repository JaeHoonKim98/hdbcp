# Function to post hoc set the minimum segment length in Inspect
inspect_prune <- function(changepoints, minseglen) {
  if (is.null(changepoints)) {
    return(integer(0))
  }
  sorted_changepoints <- changepoints[order(-changepoints[, "max.proj.cusum"]), ]

  keep <- rep(TRUE, nrow(sorted_changepoints))

  for (i in 1:nrow(sorted_changepoints)) {
    if (keep[i]) {
      current_loc <- sorted_changepoints[i, "location"]
      for (j in i:nrow(sorted_changepoints)) {
        if (current_loc == sorted_changepoints[j, "location"]) {
          next
        }
        if (abs(current_loc - sorted_changepoints[j, "location"]) <= minseglen) {
          keep[j] <- FALSE
        }
      }
    }
  }
  return(sort(sorted_changepoints[keep, "location"]))
}

# Function to apply the reconcile method within the GeomCP algorithm
geomcp_reconcile <- function(distance_points, angle_points, xi=nws[1]) {
  final_points <- sort(unique(angle_points))

  for (point in distance_points) {
    interval <- (point - xi + 1):(point + xi - 1)

    if (!any(angle_points %in% interval)) {
      final_points <- c(final_points, point)
    }
  }

  final_points <- sort(unique(final_points))

  return(final_points)
}
