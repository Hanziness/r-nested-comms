#include <Rcpp.h>
#include "is_2k2_neigh.h"
using namespace Rcpp;

bool is_2k2_neigh_cpp(std::vector<int>& n1, std::vector<int>& n2, int v1 = -1, int v2 = -1) {
    std::vector<int> n1_vec, n2_vec;

    n1_vec.reserve(n1.size());
    n2_vec.reserve(n2.size());

    std::copy_if(n1.begin(), n1.end(), std::back_inserter(n1_vec), [v2](int x)
                 { return x != v2; });
    std::copy_if(n2.begin(), n2.end(), std::back_inserter(n2_vec), [v1](int x)
                 { return x != v1; });

    std::vector<int> intersection;

    std::set_intersection(
        n1_vec.begin(),
        n1_vec.end(),
        n2_vec.begin(),
        n2_vec.end(),
        std::back_inserter(intersection));

    auto intersect_length = intersection.size();

    if (intersect_length > 0) {
        return !(intersect_length == n1_vec.size() || intersect_length == n2_vec.size());
    } else if (intersect_length == 0 && v1 >= 0 && v2 >= 0) {
        bool n1_match = std::count(n1_vec.begin(), n1_vec.end(), v2);
        bool n2_match = std::count(n2_vec.begin(), n2_vec.end(), v1);

        return !(n1_match || n2_match);
    }

    return true;
}

int nested_direction_cpp(const std::vector<int>& n1, const std::vector<int>& n2, int v1 = -1, int v2 = -1) {
    std::vector<int> n1_vec, n2_vec;

    n1_vec.reserve(n1.size());
    n2_vec.reserve(n2.size());

    std::copy_if(n1.begin(), n1.end(), std::back_inserter(n1_vec), [v2](int x)
                 { return x != v2; });
    std::copy_if(n2.begin(), n2.end(), std::back_inserter(n2_vec), [v1](int x)
                 { return x != v1; });

    std::vector<int> intersection;

    std::set_intersection(
        n1_vec.begin(),
        n1_vec.end(),
        n2_vec.begin(),
        n2_vec.end(),
        std::back_inserter(intersection));

    auto intersect_length = intersection.size();

    if (intersect_length > 0) {
        if (intersect_length == n1_vec.size()) {
            return 1;
        }
        else if (intersect_length == n2_vec.size())
        {
            return -1;
        }
        else
        {
            return 0;
        }
    } else if (intersect_length == 0 && v1 >= 0 && v2 >= 0) {
        // the intersection is empty but n1 contains v2 or n2 contains v1
        bool n1_match = std::count(n1_vec.begin(), n1_vec.end(), v2);
        bool n2_match = std::count(n2_vec.begin(), n2_vec.end(), v1);

        if (n1_match) {
            return 1;
        } else if (n2_match) {
            return -1;
        } else {
            return 0;
        }
    }

    return 0;
}

bool is_2k2_neigh(NumericVector n1, NumericVector n2, int v1 = -1, int v2 = -1) {
    std::vector<int> n1_vec, n2_vec;
    n1_vec.reserve(n1.size());
    n2_vec.reserve(n2.size());

    std::copy(n1.begin(), n1.end(), std::back_inserter(n1_vec));
    std::copy(n2.begin(), n2.end(), std::back_inserter(n2_vec));

    return nested_direction_cpp(n1_vec, n2_vec, v1, v2) == 0;
}

int nested_direction(NumericVector n1, NumericVector n2, int v1 = -1, int v2 = -1) {
    std::vector<int> n1_vec, n2_vec;
    n1_vec.reserve(n1.size());
    n2_vec.reserve(n2.size());

    std::copy(n1.begin(), n1.end(), std::back_inserter(n1_vec));
    std::copy(n2.begin(), n2.end(), std::back_inserter(n2_vec));

    return nested_direction_cpp(n1_vec, n2_vec, v1, v2);
}

double nestedness_value(NumericVector n1, NumericVector n2, int v1 = -1, int v2 = -1) {
    std::vector<int> n1_vec, n2_vec;
    n1_vec.reserve(n1.size());
    n2_vec.reserve(n2.size());

    std::copy_if(n1.begin(), n1.end(), std::back_inserter(n1_vec), [v2](int x)
                 { return x != v2; });
    std::copy_if(n2.begin(), n2.end(), std::back_inserter(n2_vec), [v1](int x)
                 { return x != v1; });

    if (n1_vec.size() == 0 || n2_vec.size() == 0) {
        return 0;
    }

    std::vector<int> intersection;

    std::set_intersection(
        n1_vec.begin(),
        n1_vec.end(),
        n2_vec.begin(),
        n2_vec.end(),
        std::back_inserter(intersection));

    auto intersect_length = intersection.size();

    return (double)intersect_length / (double)std::min(n1_vec.size(), n2_vec.size());
}