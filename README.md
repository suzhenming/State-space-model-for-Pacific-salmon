# State-space-model-for-Pacific-salmon

The relationship between spawning stock abundance (S) and
subsequent recruitment (R) is fundamental for fish population
dynamics (Quinn and Deriso, 1999). Parameter estimates for
such relationships provide the basis for setting key values of
variables used in fisheries management. However, traditional
stock-recruitment analyses that use standard regression and treat
S as an independent variable and R as a dependent variable are
fraught with estimation problems (Hilborn and Walters, 1992). One
major issue with these analyses is that they rarely consider error
caused by mismeasuring S (measurement or observation error).
However, mismeasurement of S is quite common in fisheries science
(Hilborn and Walters, 1992; Needle, 2001). Analyses that
ignore this error may be subject to severe bias in parameter estimates,
which may result in serious management problems (Walters
and Ludwig, 1981). Furthermore, large natural variation in survival
rate from eggs to recruitment (McGurk, 1986) and dynamic feedbacks
in stock-recruitment systems also exacerbate the difficulty of 
estimating stock-recruitment parameters (Hilborn and Walters,
1992).

To deal with the drawbacks of traditional stock-recruitment
analysis, several researchers have developed other methods. One
of them is an “errors-in-variables” (EV) approach (Ludwig and
Walters, 1981; Schnute, 1994; Schnute and Kronlund, 2002). The
EV method explicitly represents both natural variability and measurement
errors inherent in the data. Unfortunately, the EV method
requires that an assumption be made about what proportion of the
total variation is due to measurement error; this step is necessary to
avoid a model identification problem, i.e., a situation where model
parameters cannot be uniquely estimated (Schnute, 1994). However,
we rarely have data on that proportion, thus, many authors
assume it is 0.5 (e.g., Ludwig and Walters, 1981). Nevertheless,
there is little basis for this assumption, and this proportion’s value
can potentially greatly influence the resulting parameter estimates
of the stock-recruitment function, as we show later.
Another more advanced method, state-space modeling (Harvey,
1989), has been applied to stock-recruitment analysis by several
researchers, including Meyer and Millar (2001), Rivot et al.
(2001), Schnute and Kronlund (2002), and Walters and Martell
(2004). Here, we also used a state-space model applied to stockrecruitment
data on Pacific salmon (Oncorhynchus spp.), but
we went beyond what other authors have done by developing a version of it that contains three sources of uncertainty, i.e.,
measurement error in spawning stock size, process stochasticity
in recruitment, and the effect of time series bias. The latter
bias arises due to non-random scatter of data points caused by
linkage between a given year’s recruitment anomaly and the subsequent
spawner abundance (Walters, 1985). For the fitting process,
we develop a Bayesian approach via Markov chain Monte Carlo
(MCMC) sampling to make inferences about the state-space model
(Clark, 2007). The Bayesian approach is the most general method
for fitting non-linear state-space models (Millar and Meyer, 2000),
whereas other methods, such as the EV and extended Kalman filter
methods (Schnute and Kronlund, 2002), have more limitations
(Rivot et al., 2004).
