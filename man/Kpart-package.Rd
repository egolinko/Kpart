\name{Kpart-package}
\alias{Kpart-package}
\alias{Kpart}
\docType{package}
\title{
Kpart
}
\description{
Cubic spline regression using the absolute maximum deviate to determine potential knots. This version also includes support for addidtional independednt variables to be included in the model.
}
\details{
\tabular{ll}{
Package: \tab Kpart\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2012-08-02\cr
License: \tab Open Source\cr
}
~~ This package is intended for use with non-linearly associated data. The function part firsts selects points for cubic spline knots using an algorithm to find the absolute maximum deviate from the partition mean, then fits a best fitting model by using the best subset method and maximum adjR2. The function returns the values selected as knots in the model. The function part(d, outcomeVariable, splineTerm, additionalVars = NULL, K) takes five arguments. K is a positive integer that indicates how many equally spaced partitions the user would like to search.~~


-- Recent update includes support for additional variables, 2016-07-23.  --
}
\author{
Eric Golinko

Maintainer: egolinko@gmail.com
}
\references{
Based upon the thesis 'A min/max algorithm for cubic splines over k-partitions'
http://fau.catalog.fcla.edu/fa.jsp?ix=pm&V=D&I=0&st=FA029952462
}


\keyword{ package }
