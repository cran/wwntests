## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
library(wwntests)
set.seed(1234)
b <- brown_motion(N = 200, J = 100)
f <- far_1_S(N = 200, J = 100, S = 0.75)

## ------------------------------------------------------------------------
fport_test(f_data = b, test = 'single-lag', lag = 1, suppress_raw_output = TRUE)

## ------------------------------------------------------------------------
fport_test(f_data = f, test = 'single-lag', lag = 1, suppress_raw_output = TRUE)

## ------------------------------------------------------------------------
fport_test(f_data = b, test = 'multi-lag', lag = 10, suppress_raw_output = TRUE)

## ------------------------------------------------------------------------
fport_test(f_data = f, test = 'multi-lag', lag = 10, suppress_raw_output = TRUE)

## ---- fig.width=6, fig.height=4.25---------------------------------------
autocorrelation_coeff_plot(f_data = b, K = 20)

## ---- fig.width=6, fig.height=4.25---------------------------------------
autocorrelation_coeff_plot(f_data = f, K = 20)

## ------------------------------------------------------------------------
fport_test(b, test='spectral', bandwidth = 'static', suppress_raw_output = TRUE)

## ------------------------------------------------------------------------
fport_test(b, test='spectral', kernel = 'Bartlett', bandwidth = 3, suppress_raw_output = TRUE)

## ------------------------------------------------------------------------
fport_test(f, test='spectral', kernel = 'Parzen', bandwidth = 10, suppress_raw_output = TRUE)

## ------------------------------------------------------------------------
fport_test(f, test='spectral', bandwidth = 'adaptive', alpha = 0.01, suppress_raw_output = TRUE)

## ------------------------------------------------------------------------
fport_test(b, test = 'independence', components = 3, lag = 3, suppress_raw_output = TRUE)

## ------------------------------------------------------------------------
fport_test(f, test = 'independence', components = 16, lag = 10, suppress_raw_output = TRUE)

