# Simulated data with a known boundary at x = 0.45 to illustrate boundary detection

study.area <- data.frame(
  x = c(0, 0.8, 1, 0.6, 0, 0),
  y = c(0, 0, 0.4, 1.1, 1.1, 0)
)

set.seed(852)
# Simulate 100 occurrence coordinates
N <- 100
xy.dat_arg <- data.frame(
  x = rpois(N, 4) / 10,
  y = rpois(N, 4) / 10
)

xy.dat_arg <- xy.dat_arg[order(xy.dat_arg$x), ]

# Simulate a boundary at x = 0.45
n_left <- sum(xy.dat_arg$x < 0.45)

set.seed(1)
common_data <- letters[1:5]

left_data <- sample(
  c(common_data, LETTERS[1:10]),
  size = n_left,
  replace = TRUE,
  # common-endemic probability ratio 3:4
  prob = c(rep(3, 5), rep(4, 10))
)

right_data <- sample(
  c(common_data, LETTERS[11:20]),
  size = N - n_left,
  replace = TRUE,
  # common-endemic probability ratio 3:4
  prob = c(rep(3, 5), rep(4, 10))
)

data_arg <- c(left_data, right_data)

# Apply hespdiv
example_hespdiv <- hespdiv(
  data = data_arg,
  xy.dat = xy.dat_arg,
  n.split.pts = 6,  # small value used here for illustration
  method = "sor",   # subdivision minimizing Sorensen-Dice similarity
  S.crit = 0.3,     # minimum area size is 30% of study area
  study.pol = study.area,
  use.chull = FALSE
)
usethis::use_data(example_hespdiv, overwrite = TRUE, compress = "xz")
