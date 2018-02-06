# spatial_regression
Stata code to do Conley standard errors in OLS and 2SLS. I am indebted to Solomon Hsiang whose Stata code is the basis for my code. http://www.fight-entropy.com/2010/06/standard-error-adjustment-ols-for.html. Unfortunately Hsiang's code appears to no longer be available on his website.

One innovation of this code byond Hsiang is that this code allows for an optional prestep where distances between observations are calculated prior to the regression being run. These distances are saved as a distance matrix, where each element of the matrix gives the distance between the row observation and the column observation. Creating this distance matrix ex-ante reduces computational time when running multiple regressions, as the matrix only needs to be read in once rather than computed each time the regression is run. An additional benefit of this method is that because I use Mata to construct the distance matrix, the number of observations is not constrained by Stata's matrix size limitations (11,000 observations). These regression run fairly quickly. For example, in Lewis (2018) I have about 12,500 observations and find that regressions typically take about 5 seconds to run.

Another innovation of this code is that it allows for 2SLS. Both OLS and 2SLS are run using the same command. The 2SLS uses the parantheses syntax that is also used in Stata's ivregress.

Files include:
sample_code.do --> shows how each of these programs is used. Demonstrates different options for spatial_reg.
create_sample_data.do --> creates a sample (fake) dataset.
create_distance_matrix.ado --> program that creates a distance matrix using Mata; saves as a .txt file.
read_distance_matrix.ado --> program that reads in the .txt file and saves the distance matrix in Mata.
spatial_reg.ado --> program that does the regression itself.

This method is based on Conley (1999) "GMM Estimation with Cross Sectional Dependence." Journal of Econoometrics, Vol. 92 Issue 1 (September 1999) pp. 1-45.

If you use this code please cite Lewis (2018) "Patchwork Policies, Spillovers, and the Search for Oil and Gas". I also use this methodology in a paper with Paul Brehm entitled "To Trade or Not to Trade: Oil Leases, Information Asymmetry, and Coase"

Check out my research at https://sites.google.com/site/erickylelewis
