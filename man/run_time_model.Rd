% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{run_time_model}
\alias{run_time_model}
\title{Fit or load model along with fitting time}
\usage{
run_time_model(fit.model, file, ...)
}
\arguments{
\item{fit.model}{model fitting function}

\item{file}{filename}

\item{...}{arguments to pass to fit.model}
}
\value{
fitted model, from disk if already exists, otherwise will fit
}
\description{
Fit or load model along with fitting time
}
