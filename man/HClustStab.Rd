% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bubble.R
\name{HClustStab}
\alias{HClustStab}
\title{Calculate the average snip-rebuild distance using hclust()}
\usage{
HClustStab(data, method = "complete")
}
\arguments{
\item{data}{A numerical matrix}

\item{method}{Desired methodology for joining clusters in hclust() and as such should be  one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).}
}
\value{
Returns the average snip-rebuild distance for each sample in \code{data} using hclust()
}
\description{
Calculate the average snip-rebuild distance using hclust()
}
\examples{
BubStab(myData, "complete")
BubStab(myData, method= "ward.D")

}
