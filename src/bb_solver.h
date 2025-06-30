#pragma once
#include <queue>
#include <random>
#include <chrono>
#include "bb_instance.h"

class BabNode {
    using SolutionSet = std::vector<int>;
public:
    BabNode(SolutionSet lb, SolutionSet ub, size_t id) : lb(std::move(lb)), ub(std::move(ub)), id(id) {}
    bool operator<(const BabNode &a) const { return id < a.id; }

    SolutionSet lb;
    SolutionSet ub;
    size_t id;
};

class BabSolver {
public:
    using SolutionSet = std::vector<int>;

    /* constructor */
    explicit BabSolver(ProbInstance& problem_data_);

    /* methods */
    void branch_and_bound(std::chrono::time_point<std::chrono::high_resolution_clock> start, double time_limit);

    [[nodiscard]] size_t get_node_count() const { return num_nodes;}
    [[nodiscard]] size_t get_initial_node_count() const { return num_initial_nodes; }
    [[nodiscard]] size_t get_remained_node_count() const { return nodes.size(); }
    [[nodiscard]] bool is_time_limit() const { return time_limit_reached; }

    [[nodiscard]] double get_best_solution_value() const { return best_solution_value; }
    [[nodiscard]] SolutionSet get_best_solution() const { return best_solution; }

    [[nodiscard]] bool check_resource_constraints(const SolutionSet &x) const;

private:
    ProbInstance& problem;
    size_t nbSubSystems;

    static bool has_time_limit(const double value) {
        if (static_cast<int>(value) > 0) return true;
        return false;
    }

    SolutionSet lb_initial;
    SolutionSet ub_initial;

    SolutionSet best_solution;
    double best_solution_value;

    size_t num_nodes{0u};
    size_t num_initial_nodes{0u};

    bool time_limit_reached{false};

    std::priority_queue<BabNode> nodes;

    std::vector<std::vector<size_t>> homo_subsystems_in_subsystem;

    /* methods */
    [[nodiscard]] SolutionSet find_local_optimum(const SolutionSet& x, const SolutionSet& ub) const;
    [[nodiscard]] SolutionSet compute_upper_limit(const BabSolver::SolutionSet &x, const BabSolver::SolutionSet &ub) const;
    void branch(const BabNode& current_node);
    [[nodiscard]] static SolutionSet compute_lb_x(const SolutionSet& x, const SolutionSet& lb, int j);
    [[nodiscard]] static SolutionSet compute_ub_minus(const SolutionSet& x, const SolutionSet& ub, int j);

    [[nodiscard]] double compute_solution_value(const SolutionSet& x) const;
    [[nodiscard]] double get_system_reliability(const std::vector<int> &x) const;

    [[nodiscard]] bool is_feasible(const SolutionSet& x) const;

    static bool is_proper_dominated(const SolutionSet& x1, const SolutionSet& x2);

    static bool is_strict_dominated(const SolutionSet& x1, const SolutionSet& x2);

    static bool is_strong_dominated(const SolutionSet& x1, const SolutionSet& x2);
};