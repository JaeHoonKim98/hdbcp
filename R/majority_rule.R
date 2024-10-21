majority_rule_mxPBF <- function(result_mxPBFs, nws) {
  post_process <- function(detected_points_list, nw) {
    W <- length(detected_points_list)
    groups <- list()
    for (w in 1:W) {
      points <- detected_points_list[[w]]
      points <- points[points %notin% unlist(groups)]
      while (length(points)>0) {
        grouping <- sapply(seq_along(points), function(i){
          point <- points[i]
          interval <- (point - nw + 1):(point + nw - 1)
          interval <- interval[interval %notin% unlist(groups)]
          points_in_group <- unlist(sapply(detected_points_list[w:W], function(x) x[x %in% interval]))
          return(list(len = length(points_in_group), points = points_in_group))
        })
        if (max(unlist(grouping[1,])) <= 1) {
          for (i in 1:ncol(grouping)) {
            groups <- append(groups, grouping[2,i])
          }
          break
        }
        most_interval <- which(unlist(grouping[1,]) == max(unlist(grouping[1,])))
        if (length(most_interval)>1) {
          variances <- sapply(most_interval, function(i) var(unlist(grouping[2,i])))
          most_interval <- most_interval[which.min(variances)]
        }
        groups <- append(groups,grouping[2,most_interval])
        points <- points[points %notin% unlist(groups)]
      }
    }
    return(groups)
  }
  pre_groups <- lapply(result_mxPBFs, function(res) {
    res$Change_points
  })
  post_groups <- post_process(pre_groups, nws[1])
  major_groups <- post_groups[unlist(lapply(post_groups, function(x) length(x) >= length(pre_groups)/2))]
  minor_groups <- post_groups[unlist(lapply(post_groups, function(x) length(x) < length(pre_groups)/2))]
  for (i in seq_along(minor_groups)) {
    loc <- mean(unlist(minor_groups[[i]]))
    criteria <- length(pre_groups)
    while(loc < nws[criteria] || loc >= n - nws[criteria]) {
      if (criteria == 1) {
        break
      }
      criteria = criteria - 1
    }
    if (length(unlist(minor_groups[[i]])) >= criteria / 2) {
      major_groups <- append(major_groups, minor_groups[[i]])
    }
  }
  if (length(major_groups) == 0) {
    return(integer(0))
  }
  return(unlist(lapply(major_groups, function(x) round(mean(x)))))
}
