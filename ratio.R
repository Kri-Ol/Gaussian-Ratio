library(ggplot2)

m.x <- 5; s.x <- 4
m.y <- 4; s.y <- 1

# Implementation of ratio of two gaussian PDFs, as computed by D.Hinkley,
# see https://en.wikipedia.org/wiki/Ratio_distribution#Gaussian_ratio_distribution
gauss_ratio_PDF <- function(x, m.x = 0.0, s.x = 1.0, m.y = 0.0, s.y = 1.0) {

    a <- function(x) {
        sqrt( (x/s.x)^2 + (1.0/s.y)^2 )
    }

    b <- function(x) {
        (m.x*x)/s.x^2 + m.y/s.y^2
    }

    c <- (m.x/s.x)^2 + (m.y/s.y)^2

    d <- function(x) {
        u <- b(x)^2 - c*a(x)^2
        l <- 2.0*a(x)^2
        exp( u / l )
    }

    # PDF for the ratio of the two different gaussians, X/Y
    r <- b(x)/a(x)
    q <- pnorm(r) - pnorm(-r)

    (r*d(x)/a(x)^2) * (1.0/(sqrt(2.0*pi)*s.x*s.y)) * q + exp(-0.5*c)/(pi*s.x*s.y*a(x)^2)
}

# normalization
nn <- integrate(function(x) gauss_ratio_PDF(x, m.x=m.x, s.x=s.x, m.y=m.y, s.y=s.y), -Inf, Inf)
nn <- nn[["value"]]
print(nn) # should be 1

# plot PDF
p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
p <- p + stat_function(fun = function(x) gauss_ratio_PDF(x, m.x=m.x, s.x=s.x, m.y=m.y, s.y=s.y)/nn) + xlim(-2.0, 6.0)
print(p)


# first momentum
m1 <- integrate(function(x) x*gauss_ratio_PDF(x, m.x=m.x, s.x=s.x, m.y=m.y, s.y=s.y), -Inf, Inf)
m1 <- m1[["value"]]
print(m1)

# mean
m.xy <- m1/nn
print(m.xy)

# sampling case
set.seed(32345)
n <- 10^6L

x <- rnorm(n, mean = m.x, sd = s.x); y <- rnorm(n, mean = m.y, sd = s.y)
print(mean(x/y))
print(sd(x/y))

# second momentum
m2 <- integrate(function(x) x*x*gauss_ratio_PDF(x, m.x=m.x, s.x=s.x, m.y=m.y, s.y=s.y), -Inf, Inf)
