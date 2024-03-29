% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{read_pat_cross_df_to_list}
\alias{read_pat_cross_df_to_list}
\title{Read a df of pat cross functions and convert to a list of spotwise matrices}
\usage{
read_pat_cross_df_to_list(
  file,
  NETWORK_PATH,
  distance = 20,
  value = "gx",
  split_on = "spot"
)
}
\arguments{
\item{file}{filename to be read from NETWORK_PATH subfolder}

\item{NETWORK_PATH}{path where networks are stored}

\item{distance}{distance at which to evaluate pat cross function}

\item{value}{name of method used}

\item{split_on}{variable to split on}
}
\value{
list of "split_on"-wise matrices
}
\description{
Read a df of pat cross functions and convert to a list of spotwise matrices
}
