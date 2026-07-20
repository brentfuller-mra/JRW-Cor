# Johnson Relative Weights - Correlation matrix as input

An R function to calculate **Johnson's Relative Weights** (Relative Importance Analysis) directly from a correlation matrix rather than raw respondent-level data. 

Used for identifying the relative impact/influence of different drivers of an outcome (e.g., brand health, overall satisfaction) when independent variables are highly collinear.

Standard multiple regression breaks down when drivers are highly correlated, leading to unstable or flipped beta coefficients. Relative weights decompose the $R^2$ to show the contribution of each driver.

For general RW and related analyses, just use the relaimpo, domir, leaps, glmnet, pls, gbm, randomForest, or other packages that address multicolinearity in different ways.

## Why use the correlation matrix as an input?

1. You can account for missing values via pairwise deletion in a correlation matrix.
2. You can feed the function polychoric (or tetrachoric) correlation matrices, which in some circle as seen as more mathematically rigorous for the ordinal Likert scales typically found in consumer surveys.

## Usage and output

- **Matrix Input:** Accepts a square, symmetric correlation matrix or data frame.
- **Basic Output:** Returns a `tibble` with the raw weights and relative importance for each predictor variable.

## Installation

Source this function directly:

```R
source("https://raw.githubusercontent.com/brentfuller-mra/R-dp-helpers/main/jrelwgts_cor.R")
