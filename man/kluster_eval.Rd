\name{kluster_eval}
\alias{kluster_eval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{function to perform evaluation analysis on kluster implementations.
%%  ~~function to perform evaluation analysis on kluster implementations.  ~~
}
\description{If a specific algorithm is not specified by the user, it will perform the kluster implementations of all cluster number approximation algorithms will provide data for evaluation of the best algorithms as well as the processing time. The actual number of clusters needs to be provided for the function to calculate approximation error.
%%  ~~ If a specific algorithm is not specified by the user, it will perform the kluster implementations of all cluster number approximation algorithms will provide data for evaluation of the best algorithms as well as the processing time. The actual number of clusters needs to be provided for the function to calculate approximation error. ~~
}
\usage{
kluster_eval(data, clusters, iter_sim = 1, iter_klust, smpl, algorithm = "Default", cluster = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{clusters}{number of clusters, as we know, for calculating error. This is a requirement for this function. If you don't know the number of clusters, user `kluster` function instead.
%%     ~~Describe \code{clusters} here~~
}
  \item{iter_sim}{number of simulation iterations, default at 1
%%     ~~Describe \code{iter_sim} here~~
}
  \item{iter_klust}{number of iterations for clustering with sample_n size x
%%     ~~Describe \code{iter_klust} here~~
}
  \item{smpl}{size of the sample_n to be taken with replacement out of data
%%     ~~Describe \code{smpl} here~~
}
  \item{algorithm}{select analysis algorithm from BIC, PAMK, CAL, and AP. "Default" returns results from all available algorithms.
%%     ~~Describe \code{algorithm} here~~
}
  \item{cluster}{if TURE it'll do clustering which will take a lot longer! Not available for now...
%%     ~~Describe \code{cluster} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{returns the following values:
\item{sim }{For the selected algorithm, returns both the most frequent and the average approximated number of clusters produced by kluster procedure, processing time, and error}
\item{m_bic_k,m_cal_k,m_ap_k,m_pam_k }{the average approximated number of cluster for each selected algorithm}
\item{f_bic_k,f_cal_k,f_ap_k,f_pam_k }{the most frequent approximated number of cluster for each selected algorithm}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
dat = read.csv("data/Breast_Cancer_Wisconsin.csv")
##returning kluster's most frequent product using the BIC algorithm:
k = kluster_eval(data = dat[,c("area_mean","texture_mean")], clusters = 2, iter_sim = 1, iter_klust = 100, smpl = 100)$sim
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
