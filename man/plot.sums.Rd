\name{plot.sums}
\alias{plot.sums}
\title{Summarize SOD data by plot within a subsite}
\usage{
  plot.sums(data, subsite, year = "2002",
    sporulators = c("UMCA", "LIDE"), size.sp = "LIDE",
    sizeclass.breaks = c(0, 2, 10, 30, Inf))
}
\arguments{
  \item{data}{Tree census data such as
  \code{data(SOD.data)}}

  \item{subsite}{The subsite for which the data should be
  summarized}

  \item{year}{the year for which the summary is wanted}

  \item{sporulators}{a character vector of species which
  are considered sporulators all other species are
  classified as "OTH" with no disease}

  \item{size.sp}{Species which should be broken into size
  classes.  All other species will only have a single
  sizeclass (1)}

  \item{sizeclass.breaks}{The break points, in cm DBH,
  between size classes. Should include lower (0) and upper
  (Inf) extrema.}
}
\description{
  Summarize SOD data by plot within a subsite
}

