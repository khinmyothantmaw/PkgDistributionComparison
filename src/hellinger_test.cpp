#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Gaussian kernel density estimatio
std::vector<double> kde_estimate(const std::vector<double>& data,
                                  const std::vector<double>& grid,
                                  double bandwidth) {
    int n = data.size();
    int m = grid.size();
    std::vector<double> dens(m, 0.0);
    double norm_const = 1.0 / (n * bandwidth * std::sqrt(2.0 * M_PI));

    for (int j = 0; j < m; j++) {
        double sum_val = 0.0;
        for (int i = 0; i < n; i++) {
            double u = (grid[j] - data[i]) / bandwidth;
            sum_val += std::exp(-0.5 * u * u);
        }
        dens[j] = norm_const * sum_val;
    }
    return dens;
}

// Compute mean and standard deviation
double mean_cpp(const std::vector<double>& data) {
    double s = 0.0;
    for (double v : data) s += v;
    return s / data.size();
}


double sd_cpp(const std::vector<double>& data, double data_mean) {
    double s = 0.0;
    for (double v : data) s += (v - data_mean) * (v - data_mean);
    return std::sqrt(s / (data.size() - 1));
}

// [[Rcpp::export]]
List pairwise_hellinger_test(List data_list,
                             int grid_points = 200,
                             double bandwidth = NA_REAL) {

    int n = data_list.size();
    if (n < 2) stop("At least 2 datasets must be provided.");
    if (n > 10) stop("A maximum of 10 datasets can be compared at once.");

    std::vector<int> sizes(n);
    double global_min = R_PosInf;
    double global_max = R_NegInf;

    // Convert R list to vector of vectors and validate
    std::vector< std::vector<double> > datasets(n);
    for (int i = 0; i < n; i++) {
        SEXP elem = data_list[i];
        if (!(Rf_isInteger(elem) || Rf_isReal(elem))) {
            stop("All datasets must be numeric (integer or double).");
        }
        NumericVector x = as<NumericVector>(elem);
        if (x.size() < 10) stop("Each dataset must have at least 10 observations.");
        if (is_true(any(is_na(x)))) stop("Dataset contains NA values.");

        sizes[i] = x.size();
        datasets[i] = std::vector<double>(x.begin(), x.end());

        double min_x = *std::min_element(datasets[i].begin(), datasets[i].end());
        double max_x = *std::max_element(datasets[i].begin(), datasets[i].end());
        if (min_x < global_min) global_min = min_x;
        if (max_x > global_max) global_max = max_x;
    }

    // Bandwidth: Silverman's rule if not provided
    if (NumericVector::is_na(bandwidth)) {
        std::vector<double> all_data;
        for (auto& d : datasets) all_data.insert(all_data.end(), d.begin(), d.end());

        double m = mean_cpp(all_data);
        double s = sd_cpp(all_data, m);
        int N = all_data.size();
        bandwidth = 1.06 * s * std::pow(N, -1.0 / 5.0);
    }

    // Create grid
    std::vector<double> grid(grid_points);
    double step = (global_max - global_min) / (grid_points - 1);
    for (int i = 0; i < grid_points; i++) grid[i] = global_min + i * step;
    double dx = step;  // grid spacing for integration

    // Distance matrix
    NumericMatrix out(n, n);
    CharacterVector data_names = data_list.names();
    if (data_names.size() == 0) {
        data_names = CharacterVector(n);
        for (int i = 0; i < n; i++) data_names[i] = "D" + std::to_string(i + 1);
    }

    std::vector<std::string> pair1, pair2;
    std::vector<double> hellinger_vals;

    // Compute pairwise Hellinger distances
    for (int i = 0; i < n; i++) {
        std::vector<double> dens_i = kde_estimate(datasets[i], grid, bandwidth);

        // Normalize densities (integrate to 1)
        double sum_i = 0.0;
        for (double v : dens_i) sum_i += v;
        sum_i *= dx;
        for (auto &v : dens_i) v /= sum_i;

        for (int j = i; j < n; j++) {
            std::vector<double> dens_j = kde_estimate(datasets[j], grid, bandwidth);

            double sum_j = 0.0;
            for (double v : dens_j) sum_j += v;
            sum_j *= dx;
            for (auto &v : dens_j) v /= sum_j;

            // Compute ∫ (sqrt(p) - sqrt(q))^2 dx
            double integral = 0.0;
            for (int k = 1; k < grid_points; k++) {
                double f1 = std::sqrt(dens_i[k-1]) - std::sqrt(dens_j[k-1]);
                double f2 = std::sqrt(dens_i[k]) - std::sqrt(dens_j[k]);
                integral += 0.5 * (f1*f1 + f2*f2) * dx;
            }

            double h = std::sqrt(integral / 2.0);

            if (h > 1.0) h = 1.0; // numerical stability
            if (h < 0.0) h = 0.0;

            out(i, j) = h;
            out(j, i) = h;

            if (i < j) {
                pair1.push_back(std::string(data_names[i]));
                pair2.push_back(std::string(data_names[j]));
                hellinger_vals.push_back(h);
            }
        }
    }

    DataFrame pairwise_table = DataFrame::create(
        Named("dataset1") = pair1,
        Named("dataset2") = pair2,
        Named("hellinger") = hellinger_vals
    );

    return List::create(
        Named("distance_matrix") = out,
        Named("pairs") = pairwise_table,
        Named("data_names") = data_names,
        Named("original_data") = data_list,
        Named("sizes") = sizes,
        Named("bandwidth") = bandwidth,
        Named("grid_points") = grid_points,
        Named("n_datasets") = n,
        Named("class") = "hellinger_result"
    );
}
