library(testthat)

test_that("dataset with fewer than 10 observations throws error", {
    expect_error(ks_test(list(Small = rnorm(5), Normal = rnorm(100))), "Each dataset must have at least 10 observations.")
})

test_that("dataset containing NAs throws error", {
    bad_data <- rnorm(50)
    bad_data[10] <- NA
    good_data <- rnorm(100)

    expect_error(ks_test(list(Good = good_data, Bad = bad_data)), "Dataset contains NA values.")
})

test_that("less than 2 datasets throws error", {
    expect_error(ks_test(list(A = rnorm(50))), "At least 2 datasets must be provided.")
})


test_that("more than 10 datasets throws error", {
    too_many <- lapply(1:11, function(i) rnorm(50))
    names(too_many) <- paste0("D", 1:11)

    expect_error(ks_test(too_many), "A maximum of 10 datasets can be compared at once.")
})

test_that("non-numeric data throws error", {
    x <- rnorm(100)
    x2 <- c(x, "test")
    y <- rnorm(100)

    expect_error(ks_test(list(X = x2, Y = y)), "All datasets must be numeric")
})

test_that("ks_test summary method works as expected", {
    set.seed(123)   # ensures reproducibility

    x <- rnorm(100, mean = 0)
    y <- rnorm(100, mean = 1)

    res <- ks_test(list(X = x, Y = y))
    res.summ <- summary(res)

    expect_output(print(res.summ),"Reject H0")
})

test_that("D-statistic is symmetric across input order", {
    x <- rnorm(100)
    y <- a <- rexp(200)
    res1 <- ks_test(list(A = x, B = y))
    res2 <- ks_test(list(B = y, A = x))
    d1 <- res1$results$D.Statistic[1]
    d2 <- res2$results$D.Statistic[1]
    expect_equal(d1, d2)
})


test_that("identical datasets return p-value = 1", {
    x <- rnorm(50)
    res <- ks_test(list(A = x, B = x))
    expect_equal(res$results$p.value[1], 1)
})
