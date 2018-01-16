\name{kluster}
\alias{kluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{a function to scalably approximate the number of clusters in a dataset.
%%  ~~function to do scalably approximate the number of clusters in a dataset ~~
}
\description{kluster is the main kluster function. If an algorithm is not pre-defined, it will use the best implementation of kluster (most frequent product on BIC) for the production purpose. If a sample size is not pre-defined, it will use the recommended sample size (if n> 3000, sample size = 500, otherwise, sample size = 100) as default. If an iteration is not pre-set, it will iterate 100 times, as recommended through our simulation analyses.
%%  ~~ kluster is the main kluster function. If an algorithm is not pre-defined, it will use the best implementation of kluster (most frequent product on BIC) for the production purpose. If a sample size is not pre-defined, it will use the recommended sample size (if n> 3000, sample size = 500, otherwise, sample size = 100) as default. If an iteration is not pre-set, it will iterate 100 times, as recommended through our simulation analyses. ~~
}
\usage{
kluster(data, iter_klust = 100, smpl = 100, algorithm = "BIC")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{iter_klust}{number of iterations for clustering with sample_n size x -- present to 100
%%     ~~~~
}
  \item{smpl}{size of the sample_n to be taken with replacement out of data -- preset to 100
%%     ~~~~
}
  \item{algorithm}{select analysis algorithm from BIC, PAMK, CAL, and AP -- preset to BIC
%%     ~~~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{returns the following values:
\item{sim }{returns both the most frequent and the average approximated number of clusters for the selected algorithm}
\item{m_bic_k,m_cal_k,m_ap_k,m_pam_k }{the average approximated number of cluster for each selected algorithm}
\item{f_bic_k,f_cal_k,f_ap_k,f_pam_k }{the most frequent approximated number of cluster for each selected algorithm}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~Hossein Estiri~~
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
k = kluster(data = dat[,c("area_mean","texture_mean")],iter_klust = 100, smpl=100, algorithm = "BIC")$f_BIC_k

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
