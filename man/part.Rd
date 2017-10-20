\name{part}
\alias{part}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Fits a linear model based on spline terms with additional support for other independent variables.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The user will input a data frame, then designate the variable that is the outcome. Then the spline term is selected along with any other independent variables. Finally, a number K partitions is chosen for the algorithm to search for potential cubic spline knots based on the spline term and partition.
}
\usage{
part(d, outcomeVariable, splineTerm, additionalVars = NULL, K)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{
%%     ~~Describe \code{x} here~~
A data frame data set with column names.
}
\item{outcomeVariable}{
The variable from 'd' that is the outcome.
}
\item{splineTerm}{
The spline term, inherited from 'd'.
}
\item{additionalVars}{
A vector of additional variables to be included in the model.
}
\item{K}{
The number of evenly spaced partitions to be searched.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{fits}{
The fitted values of the linear model.
}
\item{xhat}{
The entire feature matrix.
}
\item{coefs}{
The significant coefficients of the model.
}
\item{adjr2}{
The adjusted R^2 value.
}
}
\author{
%%  ~~who you are~~
Eric Golinko
}

\examples{

## for simple spline model.
data(LakeHuron)
d <- data.frame(seq(1875, 1972, 1), LakeHuron)
names(d) <- c('date', 'lh')
fit <- part(d = d, outcomeVariable = 'lh', splineTerm = 'date', K = 20)
fit
plot(d$date, d$lh)
lines(d$date, fit$fits, col = 'red')


## multivariate
data(freeny)
freeny$time <- as.numeric(rownames(freeny))
fit <- part(d = freeny, outcomeVariable = 'y',
    splineTerm = 'time', additionalVars = c('market.potential', 'income.level'), K =2)
fit$coefs
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ regression }
% \keyword{ nonlinear }% __ONLY ONE__ keyword per line
