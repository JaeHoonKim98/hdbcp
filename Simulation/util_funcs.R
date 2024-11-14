`%notin%` <- Negate(`%in%`)

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

# Function to apply the CROPS algorithm to GeomCP results
geom_CROPS <- function(data, num_distance, num_angle, xi=nws[1]) {
  cps_distance <- cpt_distance@cpts.full[which.min(abs(apply(cpt_distance@cpts.full, 1, function(x) sum(!is.na(x))) - num_distance)),][!is.na(cpt_distance@cpts.full[which.min(abs(apply(cpt_distance@cpts.full, 1, function(x) sum(!is.na(x))) - num_distance)),])]
  cps_angle <- cpt_angle@cpts.full[which.min(abs(apply(cpt_angle@cpts.full, 1, function(x) sum(!is.na(x))) - num_angle)),][!is.na(cpt_angle@cpts.full[which.min(abs(apply(cpt_angle@cpts.full, 1, function(x) sum(!is.na(x))) - num_angle)),])]
  return(geomcp_reconcile(cps_distance, cps_angle, xi))
}

# Function to convert detected change points into segmented data
cps_to_partition <- function(values, n) {
  if (length(values) == 0) {
    return(list(c(1, n)))
  }
  sorted_values <- sort(values)
  starts <- c(1, sorted_values + 1)
  ends <- c(sorted_values, n)
  partitions <- lapply(seq_along(starts), function(i) {
    c(starts[i], ends[i])
  })
  return(partitions)
}

# Function to calculate the Jaccard index for two intervals (using only start and end values of intervals)
jaccard_index_intervals <- function(interval, interval_prime) {
  start1 <- interval[1]
  end1 <- interval[2]
  start2 <- interval_prime[1]
  end2 <- interval_prime[2]

  intersection_start <- max(start1, start2)
  intersection_end <- min(end1, end2)
  intersection_length <- max(0, intersection_end - intersection_start + 1)

  length1 <- end1 - start1 + 1
  length2 <- end2 - start2 + 1
  union_length <- length1 + length2 - intersection_length

  return(intersection_length / union_length)
}

# Function to calculate the covering metric (using only start and end values of intervals)
calculate_covering_metric <- function(partition_G, partition_G_prime, n) {
  covering_metric <- sum(sapply(partition_G, function(A) {
    max_jaccard <- max(sapply(partition_G_prime, function(A_prime) jaccard_index_intervals(A, A_prime)))
    (diff(A) + 1) / n * max_jaccard
  }))
  return(covering_metric)
}

# Function to calculate the proportion of detected points that are close to points from another method
calc_prop_within_range <- function(cp_to_calc, cp_ref, margin) {
  within_range <- sapply(cp_to_calc, function(x) {
    any(abs(x - cp_ref) <= margin)
  })
  return(mean(within_range))
}
