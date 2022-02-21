#include <Rcpp.h>
using namespace Rcpp;

//' Checks if there is a $2K_2$ between two vertices in graph `g`
//' by calculating $|N(v1) cap N(v2)| / min(|N(v1)|, |N(v2)|)$.
//' This version uses only the neighborhoods of said vertices, optionally
//' excluding `v1` and `v2` from each other's neighborhoods if `v1` or
//' `v2` are non-negative.
//' @param n1 The neighborhood of `v1`
//' @param n2 The neighborhood of `v2`
//' @param v1 (Optional) The first vertex of the pair
//' @param v2 (Optional) The second vertex of the pair
//' @return A boolean indicating whether there is a 2K_2 between `v1` and `v2`
// [[Rcpp::export]]
bool is_2k2_neigh(NumericVector n1, NumericVector n2, int v1 = -1, int v2 = -1)
{
    if (v1 >= 0)
    {
        n1 = setdiff(n1, NumericVector::create(v2));
    }
    if (v2 >= 0)
    {
        n2 = setdiff(n2, NumericVector::create(v1));
    }

    int intersect_length = intersect(n1, n2).length();

    return !(intersect_length > 0 && (intersect_length == n1.length() || intersect_length == n2.length()));
}

//' Checks if there is a $2K_2$ between two vertices in graph `g`
//' by calculating $|N(v1) cap N(v2)| / min(|N(v1)|, |N(v2)|)$.
//' This version uses only the neighborhoods of said vertices, optionally
//' excluding `v1` and `v2` from each other's neighborhoods if `v1` or
//' `v2` are non-negative.
//' @param n1 The neighborhood of `v1`
//' @param n2 The neighborhood of `v2`
//' @param v1 (Optional) The first vertex of the pair
//' @param v2 (Optional) The second vertex of the pair
//' @param treshold A value (between 0 and 1) that is still accepted as nested.
//' @return A boolean indicating whether there is a 2K_2 between `v1` and `v2`
// [[Rcpp::export]]
bool is_2k2_neigh_treshold(NumericVector n1, NumericVector n2, int v1 = -1, int v2 = -1, float treshold = 1.0)
{
    if (v1 >= 0)
    {
        n1 = setdiff(n1, NumericVector::create(v2));
    }
    if (v2 >= 0)
    {
        n2 = setdiff(n2, NumericVector::create(v1));
    }

    int minlen = n1.length() < n2.length() ? n1.length() : n2.length();

    if (minlen < 1)
    {
        return false;
    }

    return ((float)intersect(n1, n2).length() / minlen) < treshold;
}

//' Checks if there is a $2K_2$ between two vertices in graph `g`
//' by calculating $|N(v1) cap N(v2)| / min(|N(v1)|, |N(v2)|)$.
//' This version uses only the neighborhoods of said vertices, optionally
//' excluding `v1` and `v2` from each other's neighborhoods if `v1` or
//' `v2` are non-negative.
//' @param n1 The neighborhood of `v1`
//' @param n2 The neighborhood of `v2`
//' @param v1 (Optional) The first vertex of the pair
//' @param v2 (Optional) The second vertex of the pair
//' @param treshold A value (between 0 and 1) that is still accepted as nested.
//' @return A boolean indicating whether there is a 2K_2 between `v1` and `v2`
//' @export
// [[Rcpp::export]]
float is_2k2_neigh_value(NumericVector n1, NumericVector n2, int v1 = -1, int v2 = -1)
{
    if (v1 >= 0)
    {
        n1 = setdiff(n1, NumericVector::create(v2));
    }
    if (v2 >= 0)
    {
        n2 = setdiff(n2, NumericVector::create(v1));
    }

    int minlen = n1.length() < n2.length() ? n1.length() : n2.length();

    if (minlen < 1)
    {
        return 0;
    }

    return (float)(intersect(n1, n2).length()) / minlen;
}