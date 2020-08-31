# Load and describe dataset ----
dat <- read.csv(paste0(data.location, "mpi.csv"))
colnames(dat)

# Columns 7-12 are living standards data
describe(dat)

# Make data long format ----
X <- reshape(
  dat, direction = "long", varying = 7:12, v.names = "response", timevar = "item")
# Add item names
X$item_n <- factor(X$item, levels = 1:6, labels = gsub(".", " ", colnames(dat)[7:12], fixed = TRUE))
X$item_t <- paste0("x", X$item)  # Add item id: x1:x6
X$id <- rep(1:nrow(dat), 6) # Country id
X <- X[order(X$id), ]  # Sort by country
X$response <- X$response / 100  # Rescale percentages to 0-1 interval
X <- X[!is.na(X$response), ]  # Drop missing data, 3 cases
X$response.01 <- X$response  # Duplicate responses
# Check minimum and maximum values
aggregate(response.01 ~ item_n, X, function (x) c(min(x), max(x)))
# Rescale cooking fuel and electricity using transformation from
# Smithson & Verkuilen, 2006
X$response.01[X$item_n == "Cooking fuel"] <-
  (X$response.01[X$item_n == "Cooking fuel"] * (sum(X$item_n == "Cooking fuel") - 1) + .5) /
  sum(X$item_n == "Cooking fuel")
X$response.01[X$item_n == "Electricity"] <-
  (X$response.01[X$item_n == "Electricity"] * (sum(X$item_n == "Electricity") - 1) + .5) /
  sum(X$item_n == "Electricity")
# Check minimum and maximum values again
aggregate(response.01 ~ item_n, X, function (x) c(min(x), max(x)))

# Extend original dataset with 0-1 interval data ----
dat <- cbind(dat, dat[, 7:12] / 100)
colnames(dat)[15:20] <- paste0("x", 1:6)
