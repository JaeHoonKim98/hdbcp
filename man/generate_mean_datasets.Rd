% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_mean_datasets.R
\name{generate_mean_datasets}
\alias{generate_mean_datasets}
\title{Generate Simulated Datasets with Change Points in Mean Vector}
\usage{
generate_mean_datasets(
  n = 500,
  p = 200,
  signal_size = 1,
  pre_proportion = 0.4,
  pre_value = 0.3,
  single_point = round(n/2),
  multiple_points = c(round(n/4), round(2 * n/4), round(3 * n/4)),
  type = c(1, 2, 3, 4, 5)
)
}
\arguments{
\item{n}{Number of observations to generate.}

\item{p}{Number of features or dimensions for each observation.}

\item{signal_size}{Magnitude of the signal to apply at change points.}

\item{pre_proportion}{Proportion of the covariance matrix's off-diagonal elements to be set to a pre-defined value (default is 0.4).}

\item{pre_value}{Value assigned to selected off-diagonal elements of the covariance matrix (default is 0.3).}

\item{single_point}{Location of a single change point in the dataset (default is n/2).}

\item{multiple_points}{Locations of multiple change points within the dataset (default is quartiles of n).}

\item{type}{Integer specifying the type of dataset to return. Options are as follows:
- 1: No change points (H0 data)
- 2: Single change point with rare signals
- 3: Single change point with many signals
- 4: Multiple change points with rare signals
- 5: Multiple change points with many signals
The default options are 1, 2, 3, 4, and 5.}
}
\value{
A 3D array containing the generated datasets. Each slice represents a different dataset type.
}
\description{
This function generates simulated datasets that include change points in the mean vector for change point detection.
Users can specify various parameters to control the dataset size, dimension, size of signal, and change point locations.
The generated datasets include datasets with and without change points, allowing for comparisons in simulation studies.
}
\examples{
# Generate a default dataset
datasets <- generate_mean_datasets(100, 50, 1)

null_data <- datasets[,,1]
single_many_data <- datasets[,,3]

}
