# PkgDistributionCompare 

R package for comparing numeric datasets using **Kolmogorov–Smirnov test** and **Hellinger distance**.

## Features
- Pairwise Kolmogorov–Smirnov (KS) test across multiple datasets  
- Pairwise Hellinger distance between distributions  
- Summary statistics for better understanding of the results
- Visualization tools:
  - `ecdfplot()` – empirical CDF plots with KS statistics  
  - `heatmap()` – matrix visualization of KS test results or Hellinger distances  
- Designed for numeric continuous datasets 

## Installation

```{r installation, eval=FALSE}
# Install devtools if needed
install.packages("devtools")

# Load devtools
library(devtools)

# Install the package from the local directory
install("path/to/your/package")

# Or, to load all functions without installing
load_all("path/to/your/package")

# Generate documentation
document()
```

## Usage Examples

### 1. Create Example Data

```{r}
set.seed(1)
a <- rnorm(100)
b <- rnorm(100, 0.5)
c <- rnorm(80, 1)
```

### 2. KS Test

The KS test compares the empirical cumulative distribution functions (ECDFs) of two datasets.  

In this package, KS tests are applied pairwise to all datasets in your list.  
- **Null hypothesis (H0):** each pair of samples comes from the same distribution.  
- The D-statistic is the maximum vertical difference between the two ECDFs.  
- A small D (with a large p-value) suggests the distributions are similar.  
- A large D (with a small p-value) indicates significant differences.   


```{r}
res <- ks_test(list(A = a, B = b, C = c))
res
```

**Print Results**

```{r}
print(res)
```

**Summary**

```{r}
# Summary of KS test result
summary(res)

# Summary including base R ks.test results
summary(res, base_ks = TRUE)
```

### Interpretation

- **Reject null hypothesis** if p-value < 0.05  
  > The two samples come from different distributions.

- **Fail to reject null hypothesis** if p-value ≥ 0.05  
  > The two samples come from the same distribution.


**Plot ECDFs**

```{r}
# Plot ECDFs
ecdfplot(res)

# Plot ECDFs and annotate maximum D-statistic for each pair
ecdfplot(res, show_pairwise_D = TRUE)
```

**Heatmap**

```{r}
# Heatmap for D-statistic
heatmap(res)

# Heatmap for p-values
heatmap(res, type = 'p')
```

### 3. Hellinger Distance

The Hellinger distance measures similarity between two probability distributions.

Pairwise Hellinger distances are computed across all datasets.  
- Values range from **0 to 1**.  
- Lower values → distributions are closer or more similar.  
- Higher values → greater dissimilarity.  

This makes it a natural complement to the KS test:  
- If datasets come from the same distribution (KS test does not reject H₀), they usually show low Hellinger distances.  
- If KS test rejects H₀, Hellinger distances tend to be higher, reflecting separation.  


```{r}
res_h <- hellinger_test(list(A = a, B = b, C = c))
res_h
```

**Summary**

```{r}
# Summary of Hellinger distance
summary(res_h)
```

**Heatmap**

```{r}
heatmap(res_h)
```

