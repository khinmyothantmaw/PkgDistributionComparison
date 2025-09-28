#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List pairwise_ks_test(List data_list, int n_terms = 100) {
    int n = data_list.size();
    if (n < 2) stop("At least 2 datasets must be provided.");
    if (n > 10) stop("A maximum of 10 datasets can be compared at once.");

    CharacterVector names = data_list.names();
    std::vector<std::string> sample_names(n);
    for (int i = 0; i < n; i++) {
        sample_names[i] = as<std::string>(names[i]);
    }

    // Prepare result storage
    std::vector<std::string> Sample1, Sample2;
    std::vector<double> Dstats, Pvals;
    std::vector<int> n1s, n2s;

    // Loop over all pairwise combinations
    for (int i = 0; i < n - 1; i++) {
        SEXP elem = data_list[i];
        if (!(Rf_isInteger(elem) || Rf_isReal(elem))) {
            stop("All datasets must be numeric (integer or double).");
        }
        NumericVector x = as<NumericVector>(elem);
        std::sort(x.begin(), x.end());
        int nx = x.size();
        if (nx < 10) stop("Each dataset must have at least 10 observations.");
        if (is_true(any(is_na(x)))) stop("Dataset contains NA values.");

        for (int j = i + 1; j < n; j++) {
            SEXP elem2 = data_list[j];
            if (!(Rf_isInteger(elem2) || Rf_isReal(elem2))) {
                stop("All datasets must be numeric (integer or double).");
            }
            NumericVector y = as<NumericVector>(elem2);
            std::sort(y.begin(), y.end());
            int ny = y.size();

            if (ny < 10) stop("Each dataset must have at least 10 observations.");
            if (is_true(any(is_na(y)))) stop("Dataset contains NA values.");

            // Two-sample KS statistic
            int i1 = 0, j1 = 0;
            double cdf_x = 0.0, cdf_y = 0.0, D = 0.0;

            while (i1 < nx && j1 < ny) {
                if (x[i1] <= y[j1]) {
                    cdf_x = (i1 + 1.0) / nx;
                    i1++;
                } else {
                    cdf_y = (j1 + 1.0) / ny;
                    j1++;
                }
                D = std::max(D, std::fabs(cdf_x - cdf_y));
            }
            while (i1 < nx) {
                cdf_x = (i1 + 1.0) / nx;
                i1++;
                D = std::max(D, std::fabs(cdf_x - cdf_y));
            }
            while (j1 < ny) {
                cdf_y = (j1 + 1.0) / ny;
                j1++;
                D = std::max(D, std::fabs(cdf_x - cdf_y));
            }

            // Effective sample size
            double n_eff = (double)nx * ny / (nx + ny);
            long double lambda = (sqrt(n_eff) + 0.12 + 0.11 / sqrt(n_eff)) * D;

            // Asymptotic p-value
            long double pval = 0.0L;
            for (int t = 1; t <= n_terms; t++) {
                long double sign = (t % 2 == 1) ? 1.0L : -1.0L;
                pval += sign * expl(-2.0L * t * t * lambda * lambda);

            }
            double final_pval = std::min(std::max((double)(2 * pval), 0.0), 1.0);

            // Store results
            Sample1.push_back(sample_names[i]);
            Sample2.push_back(sample_names[j]);
            Dstats.push_back(D);
            Pvals.push_back(final_pval);
            n1s.push_back(nx);
            n2s.push_back(ny);
        }
    }

    DataFrame results = DataFrame::create(
        _["Sample1"] = Sample1,
        _["Sample2"] = Sample2,
        _["D.Statistic"] = Dstats,
        _["p.value"] = Pvals,
        _["n1"] = n1s,
        _["n2"] = n2s,
        _["stringsAsFactors"] = false
    );

    return List::create(
        _["results"] = results,
        _["dataList"] = data_list,
        _["method"] = "Pairwise Two-sample Kolmogorov–Smirnov test (Rcpp)"
    );
}
