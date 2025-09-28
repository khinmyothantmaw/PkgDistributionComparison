#' Pairwise Two-sample Kolmogorov–Smirnov Test
#'
#' Perform pairwise two-sample Kolmogorov–Smirnov (KS) tests to compare whether
#' multiple numeric datasets come from the same distribution.
#'
#' @param data_list A named list of numeric vectors. Minimum 2 datasets, maximum 10.
#' Each vector must have at least 10 observations and contain no NAs.
#'
#' @details
#' The Kolmogorov–Smirnov test statistic \eqn{D} is defined as the maximum absolute
#' difference between the empirical cumulative distribution functions (ECDFs) of two samples.
#' The null hypothesis is that the two samples are drawn from the same distribution.
#'
#' This function performs **pairwise comparisons** for all datasets in `data_list`.
#'
#' @return An object of class \code{"ks_test"} containing:
#' \item{results}{A data.frame with columns \code{Sample1}, \code{Sample2}, \code{D.Statistic}, \code{p.value},
#' and \code{conclusion}}
#' \item{data_list}{The original list of datasets.}
#' \item{method}{Description of the test.}
#'
#' @examples
#' set.seed(1)
#' a <- rnorm(100)
#' b <- rnorm(100, 0.5)
#' c <- rnorm(100, 1)
#'
#' res <- ks_test(list(A = a, B = b, C = c))
#'
#' # Print results
#' print(res)
#'
#' # Show summary of Kolmogorov–Smirnov test results
#' summary(res)
#'
#' # Show summary including base R ks.test results
#' summary(res, base_ks = TRUE)
#'
#' # Plot ECDFs
#' ecdfplot(res)
#'
#' # Plot ECDFs and annotate maximum D-statistic for each pair
#' ecdfplot(res, show_pairwise_D = TRUE)
#'
#' # Plot Heatmap of D-statistic and p-values
#' heatmap(res) # Shows D-statistic
#' heatmap(res, type = "p") # Show p-values
#'
#'
#' @export
ks_test <- function(data_list) {
    res <- pairwise_ks_test(data_list, n_terms = 100)
    class(res) <- "ks_test"
    return(res)
}

#' @export
print.ks_test <- function(x, ...) {
    cat(x$method, "\n\n")
    print(x$results)
    invisible(x)
}

#' @export
summary.ks_test <- function(object, base_ks = FALSE, ...) {

    res <- object$results
    res$Conclusion <- ifelse(res$p.value < 0.05, "Reject H0 (The two samples come from different distributions)", "Fail to reject H0 (The two samples come from the same distribution)")

    if (base_ks) {
        # Compute base R ks.test for each pair
        base_results <- lapply(seq_len(nrow(res)), function(i) {
            xdat <- object$dataList[[res$Sample1[i]]]
            ydat <- object$dataList[[res$Sample2[i]]]
            ks.test(xdat, ydat, alternative = "two.sided")
        })
        res$base_ks_D.Statistic <- sapply(base_results, function(z) z$statistic)
        res$base_ks_p.value <- sapply(base_results, function(z) z$p.value)
    } else {
        res$base_ks_D.Statistic <- NA
        res$base_ks_p.value <- NA
    }

    out <- list(
        results = res,
        method = object$method,
        n_datasets = length(object$data_list),
        dataset_names = names(object$data_list),
        base_ks = base_ks
    )
    class(out) <- "summary.ks_test"
    return(out)
}

#' @export
print.summary.ks_test <- function(x, ...) {
    cat("Pairwise Two-sample Kolmogorov–Smirnov Test Summary\n")
    cat("Number of datasets:", x$n_datasets, "\n")
    cat("Datasets:", paste(x$dataset_names, collapse = ", "), "\n\n")

    df <- x$results
    if (x$base_ks) {
        print(df[, c("Sample1", "Sample2", "D.Statistic", "p.value", "base_ks_D.Statistic", "base_ks_p.value", "Conclusion")],
              row.names = FALSE)
    } else {
        print(df[, c("Sample1", "Sample2", "D.Statistic", "p.value", "Conclusion")], row.names = FALSE)
    }
    invisible(x)
}

#' @export
ecdfplot <- function(x, ...) {
    UseMethod("ecdfplot")
}

#' Generic heatmap
#' @export
heatmap <- function(x, ...) {
    UseMethod("heatmap")
}

#' ECDF Plot for KS Test Objects
#'
#' Plot empirical cumulative distribution functions (ECDFs) for datasets
#' included in a \code{ks_test} object.
#'
#' @param x An object of class \code{ks_test}.
#' @param show_pairwise_D Logical. If \code{TRUE}, annotate the plot with
#'   the pairwise maximum D-statistics between datasets.
#' @param ... Additional graphical parameters passed to \code{plot}.
#' @return A base R plot is produced. No return value.
#' @export
ecdfplot.ks_test <- function(x, show_pairwise_D = FALSE, ...) {

    data_list <- x$dataList
    n <- length(data_list)

    if (n == 0) stop("No datasets to plot")

    # Assign names if missing
    if (is.null(names(data_list))) {
        names(data_list) <- paste0("Sample", seq_along(data_list))
    }

    if (n > 10) warning("More than 10 datasets may make the plot hard to read")

    # Assign colors automatically
    colors <- grDevices::rainbow(n)

    # Determine overall range for x-axis
    all_values <- unlist(data_list)
    rng <- range(all_values, finite = TRUE)

    # Plot first ECDF
    plot(ecdf(data_list[[1]]), col = colors[1], lwd = 2,
         xlim = rng, ylim = c(0, 1), xlab = "Value", ylab = "ECDF",
         main = "Empirical CDFs of Datasets", ...)

    # Add remaining ECDFs
    if (n > 1) {
        for (i in 2:n) {
            lines(ecdf(data_list[[i]]), col = colors[i], lwd = 2)
        }
    }

    # Add legend
    legend("bottomright", legend = names(data_list), col = colors, lwd = 2, bty = "n")

    # Optional: show pairwise D statistics
    if (show_pairwise_D && n > 1) {
        pairs <- combn(n, 2)
        for (k in 1:ncol(pairs)) {
            i <- pairs[1, k]
            j <- pairs[2, k]
            Fx <- ecdf(data_list[[i]])
            Fy <- ecdf(data_list[[j]])
            pooled <- sort(unique(c(data_list[[i]], data_list[[j]])))
            diffs <- abs(Fx(pooled) - Fy(pooled))
            D <- max(diffs)
            t_star <- pooled[which.max(diffs)]

            # Draw vertical segment at D
            segments(t_star, Fx(t_star), t_star, Fy(t_star),
                     col = "darkgreen", lwd = 2, lty = 2)
            # Annotate D value
            text(t_star, mean(c(Fx(t_star), Fy(t_star))),
                 labels = paste0("D=", round(D, 3)),
                 pos = 3, col = "darkgreen", cex = 0.8)
        }
    }
}

#' Heatmap for KS Test Objects
#'
#' Plot a heatmap of pairwise KS statistics or p-values for a \code{ks_test} object.
#'
#' @param x An object of class \code{ks_test}.
#' @param type Character, either \code{"D"} to show D-statistics or \code{"p"} to show p-values.
#' @param ... Additional graphical parameters passed to \code{plot}.
#' @return A base R heatmap is produced. No return value.
#' @export
heatmap.ks_test <- function(x, type = c("D", "p"), ...) {
    # type = "D" for D-statistic, "p" for p-value
    type <- match.arg(type)
    heatmap_type <- ifelse(type == "D", "D-statistic", "p-value")


    df <- x$results
    datasets <- unique(c(df$Sample1, df$Sample2))
    n <- length(datasets)

    if(n < 3) stop("Heatmap requires at least 3 datasets for meaningful visualization.")

    # Build matrix of values
    mat <- matrix(NA, nrow = n, ncol = n)
    rownames(mat) <- datasets
    colnames(mat) <- datasets

    for(i in 1:nrow(df)) {
        s1 <- df$Sample1[i]
        s2 <- df$Sample2[i]
        val <- if(type == "D") df$D.Statistic[i] else df$p.value[i]
        mat[s1, s2] <- val
        mat[s2, s1] <- val
    }
    diag(mat) <- 0  # optional self-comparison

    # Prepare color palette
    vals <- mat[lower.tri(mat) | upper.tri(mat)]
    vals <- vals[!is.na(vals)]
    col_fun <- colorRampPalette(c("white", "red"))
    ncol <- 100
    cols <- col_fun(ncol)
    breaks <- seq(min(vals), max(vals), length.out = ncol+1)

    # Plot empty
    plot(1:n, 1:n, type = "n", xlim = c(0.5, n+0.5), ylim = c(0.5, n+0.5),
         xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = paste("Pairwise KS Test Heatmap:", heatmap_type), ...)
    axis(1, 1:n, labels = datasets)
    axis(2, 1:n, labels = datasets, las = 2)

    # Draw rectangles and values
    for(i in 1:n) {
        for(j in 1:n) {
            if(i == j) {
                # Diagonal: grey
                rect(j-0.5, n-i+0.5, j+0.5, n-i+1.5, col = "grey", border = "grey")
                text(j, n-i+1, labels = "", col = "black")  # optional: leave blank
            } else if(!is.na(mat[i,j])) {
                val <- mat[i,j]
                # Determine color
                col_idx <- findInterval(val, breaks, all.inside = TRUE)
                rect(j-0.5, n-i+0.5, j+0.5, n-i+1.5, col = cols[col_idx], border = "grey")
                # Choose text color based on background
                txt_col <- if(col_idx > ncol/2) "white" else "black"
                text(j, n-i+1, labels = round(val, 3), col = txt_col)
            }
        }
    }
}
