#' Pairwise Hellinger Distance via Kernel Density Estimation
#'
#' Computes pairwise Hellinger distances between numeric datasets using
#' Gaussian kernel density estimation (KDE).
#'
#' @param data_list A named list of numeric vectors. Each element represents a dataset.
#'   Must contain at least two datasets, each with at least 10 observations and no \code{NA}s.
#' @param grid_points Number of grid points used to approximate the densities (default = 200).
#' @param bandwidth Bandwidth for the Gaussian kernel. If not provided, it is estimated
#'   using Silverman's rule of thumb from the pooled data.
#'
#' @details
#' The Hellinger distance between two probability density functions \(p(x)\) and \(q(x)\) is
#'
#' \deqn{ H(p, q) = \sqrt{\frac{1}{2} \int \left(\sqrt{p(x)} - \sqrt{q(x)}\right)^2 dx}. }
#'
#' Kernel density estimates are computed for each dataset on a common grid,
#' normalized, and pairwise Hellinger distances are calculated.
#'
#' @return A list of class \code{"hellinger_result"} with the following elements:
#' \itemize{
#'   \item \code{distance_matrix} - Symmetric matrix of pairwise Hellinger distances.
#'   \item \code{pairs} - Data frame of pairwise distances with columns
#'         \code{dataset1}, \code{dataset2}, and \code{hellinger}.
#'   \item \code{data_names} - Names of the datasets.
#'   \item \code{sizes} - Sample sizes of each dataset.
#'   \item \code{bandwidth} - Bandwidth used for KDE.
#'   \item \code{grid_points} - Number of grid points used.
#'   \item \code{n_datasets} - Number of datasets compared.
#' }
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100)
#' y <- rnorm(100, mean = 1)
#' z <- rnorm(100, mean = 2)
#'
#' res <- hellinger_test(list(X = x, Y = y, Z = z))
#' print(res)
#' summary(res)
#' heatmap(res)
#' @export
hellinger_test <- function(data_list) {
    res <- pairwise_hellinger_test(data_list)
    class(res) <- "hellinger_test"
    attr(res, "original_data") <- data_list
    return(res)
}

#' @export
print.hellinger_test <- function(x, ...) {
    cat("Hellinger Distance Result\n")
    cat("-------------------------\n")
    cat("Number of datasets:", x$n_datasets, "\n")

    cat("Dataset Sizes:\n")
    sizes <- setNames(x$sizes, x$data_names)
    print(sizes)
    cat("\n")

    cat("Pairwise Hellinger Distances:\n")
    print(x$pairs)

    invisible(x)
}

#' @export
summary.hellinger_test<- function(object, ...) {
    pairs <- object$pairs
    distances <- pairs$hellinger

    summary_list <- list(
        n_datasets = object$n_datasets,
        bins = object$bins,
        sizes = setNames(object$sizes, object$data_names),
        min_distance = min(distances),
        max_distance = max(distances),
        mean_distance = mean(distances),
        median_distance = median(distances),
        pairwise = pairs,
        distance_matrix = object$distance_matrix
    )

    class(summary_list) <- "summary.hellinger_test"
    summary_list
}

#' @export
print.summary.hellinger_test <- function(x, ...) {
    cat("Hellinger Distance Summary\n")
    cat("--------------------------\n")
    cat("Number of datasets:", x$n_datasets, "\n")
    cat("Histogram bins:", x$bins, "\n\n")

    cat("Dataset Sizes:\n")
    print(x$sizes)
    cat("\n")

    cat("Minimum distance:", round(x$min_distance, 4), "\n")
    cat("Maximum distance:", round(x$max_distance, 4), "\n")
    cat("Mean distance   :", round(x$mean_distance, 4), "\n")
    cat("Median distance :", round(x$median_distance, 4), "\n")
    cat("\nPairwise Distances:\n")
    print(x$pairwise)
}

#' @export
heatmap <- function(x, ...) {
    UseMethod("heatmap")
}

#' Heatmap for Hellinger Test Objects
#'
#' Plot a heatmap of pairwise Hellinger distance for a \code{hellinger_test} object.
#'
#' @param x An object of class \code{hellinger_test}.
#' @param ... Additional graphical parameters passed to \code{plot}.
#' @return A base R heatmap is produced. No return value.
#' @export
heatmap.hellinger_test <- function(x, ...) {
    mat <- x$distance_matrix
    n <- nrow(mat)
    labels <- x$data_names

    # Color scale: light pink <0.1, pink->red above
    max_val <- max(mat, na.rm = TRUE)
    breaks <- c(seq(0, 0.1, length.out = 10), seq(0.1001, max_val, length.out = 41))
    cols <- c(rep("linen", 9),
              colorRampPalette(c("wheat", "tan", "peru"))(41))

    par(mar = c(6, 6, 3, 3))
    image(1:n, 1:n, t(mat[n:1, ]), col = cols, breaks = breaks,
          xaxt = "n", yaxt = "n", xlab = "", ylab = "",
          main = "Hellinger Distance Heatmap", ...)

    # Overlay numbers and grey main diagonal
    for (i in 1:n) {
        for (j in 1:n) {
            val <- mat[i, j]
            y <- n - j + 1  # flip y-axis

            # Grey main diagonal
            if (i == j) {
                rect(i-0.5, y-0.5, i+0.5, y+0.5, col = "grey80", border = NA)
            } else {
                # Text color: dark for <0.1, white for higher

                text(i, y, sprintf("%.3f", val), col = "black", cex = 0.8)
            }
        }
    }

    # Rotated x-axis labels
    x_pos <- 1:n
    y_pos <- par("usr")[3] - 0.5
    text(x = x_pos, y = y_pos, labels = labels, srt = 45, adj = 1, xpd = TRUE)

    # Y-axis labels
    axis(2, at = 1:n, labels = rev(labels), las = 1)

    box()
}

#' @export
plot_kde <- function(x, ...) {
    UseMethod("plot_kde")
}

# S3 method
#' Plot KDE for Hellinger test objects
#'
#' @param x A hellinger_test object
#' @param n_points Number of grid points for KDE
#' @param colors Vector of colors
#' @param ... Additional graphical parameters
#' @export
plot_kde.hellinger_test <- function(x, n_points = 200, colors = NULL, ...) {
    data_list <- attr(x, "original_data")
    if (is.null(data_list)) stop("Original data not stored in object.")

    n <- length(data_list)
    if (is.null(colors)) colors <- rainbow(n)

    # Determine global range
    all_values <- unlist(data_list)
    global_min <- min(all_values)
    global_max <- max(all_values)
    grid <- seq(global_min, global_max, length.out = n_points)

    # Plot setup
    plot(NULL, xlim = c(global_min, global_max), ylim = c(0, max(sapply(data_list, function(d) {
        d <- density(d, from = global_min, to = global_max, n = n_points)
        max(d$y)
    }))), xlab = "Value", ylab = "Density", main = "KDE of datasets", ...)

    # Overlay densities
    for (i in seq_along(data_list)) {
        d <- density(data_list[[i]], from = global_min, to = global_max, n = n_points)
        lines(d$x, d$y, col = colors[i], lwd = 2)
    }

    # Add legend
    legend("topright", legend = names(data_list), col = colors, lwd = 2)
}
