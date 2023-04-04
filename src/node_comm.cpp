#include <Rcpp.h>
#include <stack>
#include "is_2k2_neigh.h"
using namespace Rcpp;

using LKSet_inner = std::vector<int>;
using LKSet = std::vector<LKSet_inner>;

LKSet list_to_lkset(List &lk)
{
    LKSet result;
    result.reserve(lk.size());

    for (int i = 0; i < lk.size(); i++)
    {
        NumericVector current_item = lk[i];
        std::vector<int> current_row;
        current_row.reserve(current_item.size());
        current_item = current_item - 1; // subtract 1 to convert them to indices
        std::copy(current_item.begin(), current_item.end(), std::back_inserter(current_row));
        result.push_back(current_row);
    }

    return result;
}

List lkset_to_list(const LKSet &lkset)
{
    List result;

    for (std::size_t i = 0; i < lkset.size(); i++)
    {
        IntegerVector current_item;
        for (const auto &el : lkset[i]) {
            current_item.push_back(el + 1);
        }

        result.push_back(current_item);
    }

    return result;
}

LKSet get_lk2_from_neighs(const LKSet &neighs)
{
    LKSet lk2;

    for (std::size_t i = 0; i < neighs.size(); i++)
    {
        // NOTE we subtract 1 at indexing due to C++-R indexing differences
        for (int neigh : neighs[i])
        {
            for (int neigh2 : neighs[neigh])
            {
                if (i >= neigh2)
                {
                    continue;
                }

                int direction = nested_direction_cpp(neighs[i], neighs[neigh2], i, neigh2);
                if (direction == 1)
                {
                    lk2.push_back(std::vector<int>({(int)i, (int)neigh2}));
                }
                else if (direction == -1)
                {
                    lk2.push_back(std::vector<int>({(int)neigh2, (int)i}));
                }
            }
        }
    }

    // deduplicate list
    std::set<std::vector<int>> deduplicated_list(lk2.begin(), lk2.end());
    lk2 = std::vector<LKSet_inner>(deduplicated_list.begin(), deduplicated_list.end());

    return lk2;
}

/// @brief Container structure for a topsort level: it contains the current level's iterator and whether there was a further push from this level.
struct topsort_level
{
    unsigned long vertex;
    unsigned long neigh_index;
    bool pushed;

    topsort_level(unsigned long vertex) : vertex(vertex), pushed(false), neigh_index(0) {}
};

LKSet longest_path_topsort(const LKSet &l2)
{
    // the paths that will be returned
    LKSet paths;

    // find all roots
    // 1. find root nodes
    std::set<int> all_nodes;
    std::set<int> non_roots;
    std::deque<int> roots;

    LKSet neighs;
    int max_size = 0;

    for (const auto &row : l2)
    {
        all_nodes.insert(row[0]);
        all_nodes.insert(row[1]);
        non_roots.insert(row[1]);

        if (max_size < row[0]) {
            max_size = row[0];
        }

        if (max_size < row[1]) {
            max_size = row[1];
        }
    }
    std::set_difference(all_nodes.begin(), all_nodes.end(), non_roots.begin(), non_roots.end(), std::back_inserter(roots));

    neighs.resize(max_size + 1);
    for (auto &row : l2)
    {
        neighs[row[0]].push_back(row[1]);
    }

    // 2. create level stack

    std::stack<topsort_level> levels;
    LKSet_inner current_path;

    while (!roots.empty())
    {
        Rcpp::checkUserInterrupt();

        // pick next root
        int root = roots.back();
        roots.pop_back();

        levels.push(topsort_level(root));

        current_path.push_back(root);
        while (!levels.empty())
        {
            Rcpp::checkUserInterrupt();
            auto &top = levels.top();
            const auto path_last_item = current_path.back();
            auto &current_neighs = neighs[top.vertex];

            // if all neighbors have been visited
            if (top.neigh_index >= current_neighs.size())
            {
                // if we haven't branched off of this level, we can add the path to the list
                if (!top.pushed)
                {
                    paths.push_back(current_path);
                }

                current_path.pop_back();
                levels.pop();
                continue;
            }

            // get the next neighbor
            const auto next_neigh = current_neighs[top.neigh_index];
            top.neigh_index++;

           
            // push the next neighbor into the path
            current_path.push_back(next_neigh);

            top.pushed = true;
            levels.push(topsort_level(next_neigh));
        }
    }

    return paths;
}

/// @brief Returns the transitive reduction of a DAG
/// @param l2 The edgelist of the graph
/// @return The reduced edgelist
LKSet transitive_reduction(const LKSet &l2) {
    LKSet l2_reduced = l2;
    int current_node;
    int second_node;

    auto first_it = l2_reduced.begin();
    while (first_it != l2_reduced.end())
    {
        bool was_deletion = false;
        current_node = (*first_it)[0];
        second_node = (*first_it)[1];

        auto second_it = l2.begin();

        while (second_it != l2.end())
        {
            // we can break here because the vector is ordered (ascending) by the first item
            if ((*second_it)[0] > second_node) {
                break;
            } else if ((*second_it)[0] == second_node) {
                auto found_position = std::find_if(l2_reduced.begin(), l2_reduced.end(), [current_node, second_it](std::vector<int> pair)
                                                { return pair[0] == current_node && pair[1] == (*second_it)[1]; });
                if (found_position != l2_reduced.end()) {
                    auto after_delete_it = l2_reduced.erase(found_position);
                    if (after_delete_it <= first_it) {
                        first_it = after_delete_it;
                    }

                    was_deletion = true;
                    continue;
                }
            }

            second_it++;
        }

        if (!was_deletion) {
            first_it++;
        }
    }

    return l2_reduced;
}

/// Return a map where vertices with the same neighborhood are mapped to the same index
std::map<int, int> get_similarity_map(const LKSet &neighs) {
    const std::size_t n = neighs.size();
    std::map<int, int> vertex_map;
    std::vector<std::set<int>> similar_vertices;
    for (std::size_t i = 0; i < n; i++)
    {
        for (std::size_t j = i + 1; j < n; j++)
        {
            auto neighs_i = neighs[i];
            auto neighs_j = neighs[j];

            auto pos_i_in_j = std::find(neighs_j.begin(), neighs_j.end(), i);

            // remove i and j from each other's neighborhoods
            if (pos_i_in_j != neighs_j.end()) {
                neighs_j.erase(pos_i_in_j);
            }

            auto pos_j_in_i = std::find(neighs_i.begin(), neighs_i.end(), j);
            if (pos_j_in_i != neighs_i.end()) {
                neighs_i.erase(pos_j_in_i);
            }

            if (neighs_i.size() > 0 && neighs_i == neighs_j) {
                bool found = false;

                // check if there's already a group for vertex i
                for (auto &group : similar_vertices)
                {
                    if (std::find(group.begin(), group.end(), i) != group.end()) {
                        found = true;
                        group.insert(i);
                        group.insert(j);
                        break;
                    }
                }

                // create new group if there's no group for them
                if (!found) {
                    similar_vertices.push_back(std::set<int>({(int)i, (int)j}));
                }
            }
        }
    }

    // convert the list of sets to a map
    int min_item;
    for (const auto &group : similar_vertices) {
        min_item = *std::min_element(group.begin(), group.end());
        for (const auto &element : group) {
            if (element == min_item) {
                continue;
            }
            vertex_map.insert(std::make_pair(element, min_item));
        }
    }

    return vertex_map;
}

//' @export
// [[Rcpp::export]]
List get_lk_all_topsort(List neighs)
{
    LKSet neighs_lkset = list_to_lkset(neighs);

    std::map<int, int> similarity_map = get_similarity_map(neighs_lkset);
    // replace vertices in neighs_lkset
    for (const auto &pair : similarity_map) {
        // remove the neighborhood of the mapped vertex
        neighs_lkset[pair.first].clear();

        // remove edges to vertices with empty neighborhood
        for (auto &neigh_row : neighs_lkset) {
            const auto find_index = std::find(neigh_row.begin(), neigh_row.end(), pair.first);
            if (find_index != neigh_row.end()) {
                neigh_row.erase(find_index);
            }
        }
    }
    LKSet l2 = get_lk2_from_neighs(neighs_lkset);
    l2 = transitive_reduction(l2);

    LKSet paths = longest_path_topsort(l2);

    // Add nodes with no entries in L2
    // List of all nodes
    std::vector<int> full_nodes(neighs.size());
    std::iota(full_nodes.begin(), full_nodes.end(), 0); // 0..(n-1)
    // Set of found nodes in L2
    std::set<int> found_nodes;
    // The resulting list of nodes not in L2
    std::vector<int> detached_nodes;
    for (auto &row : l2)
    {
        found_nodes.insert(row[0]);
        found_nodes.insert(row[1]);
    }
    std::set_difference(full_nodes.begin(), full_nodes.end(), found_nodes.begin(), found_nodes.end(), std::back_inserter(detached_nodes));

    // Create a separate community for each detached node
    for (const auto &el : detached_nodes) {
        // do not add these for collapsed nodes!
        if (similarity_map.find(el) == similarity_map.end()) {
            paths.push_back(std::vector<int>({el}));
        }
    }

    // restore compressed vertices
    std::map<int, std::vector<int>> restore_map;
    for (const auto &pair : similarity_map) {
        if (restore_map.count(pair.second) == 0) {
            restore_map[pair.second] = std::vector<int>({pair.first});
        } else {
            restore_map[pair.second].push_back(pair.first);
        }
    }

    for (auto &restore_data : restore_map)
    {
        std::sort(restore_data.second.begin(), restore_data.second.end());

        // find compressed vertices
        for (auto &path : paths) {
            const auto find_index = std::find(path.begin(), path.end(), restore_data.first);

            // insert mapped vertex into path
            if (find_index != path.end()) {
                path.insert(find_index + 1, restore_data.second.begin(), restore_data.second.end());
            }
        }
    }

    std::sort(paths.begin(), paths.end());

    return lkset_to_list(paths);
}