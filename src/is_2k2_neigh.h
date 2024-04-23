#include <Rcpp.h>

bool is_2k2_neigh_cpp(std::vector<int>& n1, std::vector<int>& n2, int v1, int v2);

/// @brief Calculates the direction of the nested relationship between two nodes
/// @param n1 The neighborhood of vertex #1
/// @param n2 The neighborhood of vertex #2
/// @param v1 Vertex #1, if >= 0, will be removed from `n2`
/// @param v2 Vertex #2, if >= 0, will be removed from `n1`
/// @return 1 if `n1` contains `n2`, -1 if `n2` contains `n1` and 0 if the two vertices are not nested.
int nested_direction_cpp(const std::vector<int> &n1, const std::vector<int> &n2, int v1, int v2);

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
//' @export
// [[Rcpp::export]]
bool is_2k2_neigh(Rcpp::NumericVector n1, Rcpp::NumericVector n2, int v1, int v2);

//' Checks if there is a $2K_2$ between two vertices in graph `g`
//' by calculating $|N(v1) cap N(v2)| / min(|N(v1)|, |N(v2)|)$.
//' This version uses only the neighborhoods of said vertices, optionally
//' excluding `v1` and `v2` from each other's neighborhoods if `v1` or
//' `v2` are non-negative.
//' @param n1 The neighborhood of `v1`
//' @param n2 The neighborhood of `v2`
//' @param v1 (Optional) The first vertex of the pair
//' @param v2 (Optional) The second vertex of the pair
//' @return The direction of nestedness between `v1` and `v2`. Returns 1 if `n1` contains `n2`, -1 if `n2` contains `n1` and 0 if the two nodes are not nested.
//' @export
// [[Rcpp::export]]
int nested_direction(Rcpp::NumericVector n1, Rcpp::NumericVector n2, int v1, int v2);

//' @export
// [[Rcpp::export]]
double nestedness_value(Rcpp::NumericVector n1, Rcpp::NumericVector n2, int v1, int v2);