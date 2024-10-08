% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.cosimmrSTAN_pred_out.R
\name{summary.cosimmrSTAN_pred_out}
\alias{summary.cosimmrSTAN_pred_out}
\title{Summarises the output created with \code{\link{predict.cosimmrSTAN_output}}}
\usage{
\method{summary}{cosimmrSTAN_pred_out}(
  object,
  type = c("quantiles", "statistics", "correlations"),
  obs = 1,
  ...
)
}
\arguments{
\item{object}{An object of class \code{cosimmr_pred_output} produced by the
function  \code{\link{predict.cosimmr_output}}}

\item{type}{The type of output required. At least none of quantiles',
'statistics', or 'correlations'.}

\item{obs}{The observation to generate a summary for. Defaults to 1.}

\item{...}{Not used}
}
\value{
A list containing the following components:
\item{quantiles }{The quantiles of each parameter from the posterior
distribution} \item{statistics }{The means and standard
deviations of each parameter} \item{correlations }{The posterior
correlations between the parameters} Note that this object is reported
silently so will be discarded unless the function is called with an object
as in the example below.
}
\description{
Produces textual summaries and convergence diagnostics for an object created
with  \code{\link{predict.cosimmrSTAN_output}}. The different
options are: 'quantiles' which produces credible intervals
for the parameters, 'statistics' which produces means and standard
deviations, and 'correlations' which produces correlations between the
parameters.
}
\details{
The quantile output allows easy calculation of 95 per cent credible
intervals of the posterior dietary proportions. The correlations allow the
user to judge which sources are non-identifiable.
}
\examples{
\donttest{
# A simple example with 10 observations, 2 tracers and 4 sources

# The data
data(geese_data_day1)
cosimmr_1 <- with(
  geese_data_day1,
  cosimmrSTAN_load(
    formula = mixtures ~ c(1,2,3,3,2,3,1,2,1),
    source_names = source_names,
    source_means = source_means,
    source_sds = source_sds,
    correction_means = correction_means,
    correction_sds = correction_sds,
    concentration_means = concentration_means
  )
)

# Plot
plot(cosimmr_1)


# FFVB run
cosimmr_1_out <- cosimmr_ffvb(cosimmr_1)

# Summarise
summary(cosimmr_1_out) # This outputs all the summaries
summary(cosimmr_1_out, type = "quantiles") # Just the diagnostics
# Store the output in ans
ans <- summary(cosimmr_1_out,
  type = c("quantiles", "statistics")
)
}
}
\seealso{
See \code{\link{cosimmr_stan}}for creating objects suitable for
this function, and many more examples.
See also \code{\link{cosimmr_stan}} for creating cosimmr objects,
\code{\link{plot.cosimmr_input}} for creating isospace plots,
\code{\link{plot.cosimmr_output}} for plotting output.
}
\author{
Emma Govan <emmagovan@gmail.com> Andrew Parnell
}
