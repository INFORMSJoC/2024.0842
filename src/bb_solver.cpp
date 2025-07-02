#include <iostream>
#include "bb_utility.h"
#include "bb_solver.h"

BabSolver::BabSolver(ProbInstance &problem_data_) : problem(problem_data_) {
    nbSubSystems = problem.num_subsystems * problem.num_component_types;

    homo_subsystems_in_subsystem.resize(problem.num_subsystems);
    for (int j = 0; j < problem.num_subsystems; ++j) {
        homo_subsystems_in_subsystem.reserve(problem.num_component_types);
        for (int h = 0; h < problem.num_component_types; h++) {
            homo_subsystems_in_subsystem[j].push_back(j * problem.num_component_types + h);
        }
    }

    /* initialize the initial bounds */
    for (const auto &row: problem.component_lb_at_subsystem) {
        lb_initial.insert(lb_initial.end(), row.begin(), row.end());
    }

    ub_initial.reserve(nbSubSystems);
    for (const auto &row: problem.component_ub_at_subsystem) {
        ub_initial.insert(ub_initial.end(), row.begin(), row.end());
    }

    best_solution_value = 0.0;
}

bool BabSolver::check_resource_constraints(const SolutionSet &x) const {
    // Example 1
    if (problem.test_type == TEST_TYPE::EXAMPLE_ONE) {
        {
            auto e1_c1_lhs = 10.0 * std::exp(x[0]/2.0) * x[1] + 20 * x[2] + 3.0 * std::pow(x[3], 2) + 8.0 * x[4];
            auto e1_c1_rhs = 200.0;
            if (e1_c1_lhs > e1_c1_rhs) return false;
        }

        {
            auto e1_c2_lhs = 10.0 * std::exp(x[0]/2.0) + 4.0 * std::exp(x[1]) + 2.0 * std::pow(x[2], 3)
                  + 6.0 * (std::pow(x[3], 2.0) + std::exp(x[3] / 4.0)) + 7.0 * std::exp(x[4] / 4.0);
            auto e1_c2_rhs = 310.0;
            if (e1_c2_lhs > e1_c2_rhs) return false;
        }

        {
            auto e1_c3_lhs = 12.0 * (std::pow(x[1], 2.0) + std::exp(x[1])) + 5.0 * x[2] * std::exp(x[2] / 4.0)
                  + 3.0 * x[0] * std::pow(x[3], 2.0) + 2.0 * std::pow(x[4], 3.0);
            auto e1_c3_rhs = 520.0;
            if (e1_c3_lhs > e1_c3_rhs) return false;
        }
    }
    // Example 2
    else if (problem.test_type == TEST_TYPE::EXAMPLE_TWO) {
        {
            auto e2_c1_lhs = 8.0 * std::exp(x[0] / 2.0) * x[1]
                     + 4.0 * std::exp(x[2] / 2.0)
                     + 2.0 * x[3]
                     + 2.0 * (x[4] + std::exp(x[4] / 4.0))
                     + 6.0 * std::pow(x[5], 2.0) * x[6]
                     + 2.0 * x[7]
                     + 8.0 * std::pow(x[8], 3.0) * std::exp(x[9] / 2.0);
            auto e2_c1_rhs = 120.0;

            if (e2_c1_lhs > e2_c1_rhs) return false;
        }

        {
            auto e2_c2_lhs = 16.0 *  std::pow(x[0], 2.0) * x[1]
                   + 6.0 * std::exp(x[2] * x[3] / 6.0)
                   + 7.0 * x[4] * std::exp(x[5] / 4.0)
                   + 12.0 * x[6] * std::pow(x[7], 3.0)
                   + 1.0 * (x[8] + std::exp(x[8] / 2.0))
                   + 9.0 * x[1] * std::exp(x[9] / 4.0);
            auto e2_c2_rhs = 300.0;

            if (e2_c2_lhs > e2_c2_rhs) return false;
        }
    }
    else if (problem.test_type == TEST_TYPE::TEST_ONE) {
        {// the 1st constraint
            auto t1_c1_lhs = 0.0;
            for (int j = 0; j < problem.num_subsystems; ++j) {
                for (auto h = 0u; h < problem.num_component_types; ++h) {
                    t1_c1_lhs += problem.resource_for_component_at_subsystem[0u][j][h] * x[j * problem.num_component_types + h];
                }
            }
            auto t1_c1_rhs = problem.amount_of_each_resources[0u];
            if (t1_c1_lhs > t1_c1_rhs) return false;
        }
        {// the 2nd constraint
            auto t1_c2_lhs = 0.0;
            for (int j = 0; j < problem.num_subsystems; ++j) {
                for (auto h = 0u; h < problem.num_component_types; ++h) {
                    t1_c2_lhs += problem.resource_for_component_at_subsystem[1u][j][h] *
                        (x[j * problem.num_component_types + h] + std::exp(x[j * problem.num_component_types + h] / 4.0) );
                }
            }
            auto t1_c2_rhs = problem.amount_of_each_resources[1u];
            if (t1_c2_lhs > t1_c2_rhs) return false;
        }
        {// the 3rd constraint
            auto t1_c3_lhs = 0.0;
            for (int j = 0; j < problem.num_subsystems; ++j) {
                for (auto h = 0u; h < problem.num_component_types; ++h) {
                    t1_c3_lhs += problem.resource_for_component_at_subsystem[2u][j][h] *
                        std::pow(x[j * problem.num_component_types + h], 3.0);
                }
            }
            auto t1_c3_rhs = problem.amount_of_each_resources[2u];
            if (t1_c3_lhs > t1_c3_rhs) return false;
        }
    }
    else if (problem.test_type == TEST_TYPE::TEST_TWO) {
        {// the 1st constraint
            auto t2_c1_lhs = 0.0;
            for (int j = 0; j < problem.num_subsystems; ++j) {
                for (auto h = 0u; h < problem.num_component_types; ++h) {
                    t2_c1_lhs += problem.resource_for_component_at_subsystem[0u][j][h] * x[j * problem.num_component_types + h];
                }
            }
            auto t2_c1_rhs = problem.amount_of_each_resources[0u];

            if (t2_c1_lhs > t2_c1_rhs) return false;
        }
        {// the 2nd constraint
            auto t2_c2_lhs = 0.0;
            for (int j = 0; j < problem.num_subsystems; ++j) {
                for (auto h = 0u; h < problem.num_component_types; ++h) {
                    t2_c2_lhs += problem.resource_for_component_at_subsystem[1u][j][h] *
                        std::pow(x[j * problem.num_component_types + h], 2.0);
                }
            }
            auto t2_c2_rhs = problem.amount_of_each_resources[1u];

            if (t2_c2_lhs > t2_c2_rhs) return false;
        }

        {// the 3rd constraint
            auto t2_c3_lhs = 0.0;
            for (int j = 0; j < problem.num_subsystems; ++j) {
                for (auto h = 0u; h < problem.num_component_types; ++h) {
                    t2_c3_lhs += problem.resource_for_component_at_subsystem[2u][j][h] * x[j * problem.num_component_types + h]
                       * std::exp(x[j * problem.num_component_types + h] / 4.0);
                }
            }
            auto t2_c3_rhs = problem.amount_of_each_resources[2u];

            if (t2_c3_lhs > t2_c3_rhs) return false;
        }
    }
    else if (problem.test_type == TEST_TYPE::MIXED) {
        for (auto i = 0u; i < problem.num_resource_types; ++i)
        {
            auto mixed_c_lhs = 0.0;

            for (auto j = 0u; j < problem.num_subsystems; ++j) {
                for (auto h = 0u; h < problem.num_component_types; ++h) {
                    mixed_c_lhs += problem.resource_for_component_at_subsystem[i][j][h] * x[j * problem.num_component_types + h];
                }
            }

            auto mixed_rhs = problem.amount_of_each_resources[i];

            if (mixed_c_lhs > mixed_rhs) return false;
        }
    }

    return true;
}

void BabSolver::branch_and_bound(std::chrono::time_point<std::chrono::high_resolution_clock> start, const double time_limit) {
    auto last_output_time = 0.0;

    // initialize
    best_solution = lb_initial;
    best_solution_value = compute_solution_value(best_solution);

    // std::cout << "best_solution: " << get_best_solution_value() << std::endl;

    const auto lb = lb_initial;
    const auto ub = ub_initial;

    // generate the root nodes
    const auto num_subsystems = static_cast<int>(homo_subsystems_in_subsystem.size());
    std::vector<int> indices(num_subsystems, 0);
    while (true) {
        // Print the current combination
        auto ll = lb;
        auto uu = ub;
        for (int j = 0; j < num_subsystems; ++j) {
            // std::cout << homo_subsystems_in_subsystem[j][indices[j]] << (j == num_subsystems - 1 ? "\n" : " ");
            ll[homo_subsystems_in_subsystem[j][indices[j]]] = 1;
            for (auto k = homo_subsystems_in_subsystem[j][0]; k < homo_subsystems_in_subsystem[j][indices[j]]; ++k) {
                uu[k] = 0;
            }
        }
        uu = compute_upper_limit(ll, uu);
        auto node = BabNode(ll, uu, ++num_nodes);
        nodes.push(node);

        // find the rightmost vector that can be incremented
        int i = num_subsystems - 1;
        while (i >= 0) {
            indices[i]++;
            if (indices[i] < homo_subsystems_in_subsystem[i].size()) {
                break; // successfully incremented, move to the next combination
            }
            indices[i] = 0; // reset to the first element of this vector
            i--; // move to the next left vector
        }

        // if all indices are reset to zero, we've printed all combinations
        if (i < 0) break;

        // check time limit
        const auto end = std::chrono::high_resolution_clock::now();
        const auto used_time = std::chrono::duration<double>(end - start).count();

        if (used_time - last_output_time >= 60) {
            std::cout << "used time is " << used_time << " seconds\n";
            last_output_time = used_time;
            std::cout << "remained nodes: " << get_remained_node_count() << "\n";
        }

        if (has_time_limit(time_limit) && used_time >= time_limit) {
            std::cout << "time limit (" << time_limit << " seconds) is reached." << std::endl;
            time_limit_reached = true;
            return;
        }
    }

    // check time limit
    {
        const auto end = std::chrono::high_resolution_clock::now();
        const auto used_time = std::chrono::duration<double>(end - start).count();
        if (used_time - last_output_time >= 60) {
            std::cout << "used time is " << used_time << " seconds\n";
            last_output_time = used_time;
            std::cout << "remained nodes: " << get_remained_node_count() << "\n";
        }
        if (has_time_limit(time_limit) && used_time >= time_limit) {
            std::cout << "time limit (" << time_limit << " seconds) is reached." << std::endl;
            time_limit_reached = true;
            return;
        }
    }

    // the number of initial nodes
    num_initial_nodes = nodes.size();

    while(!nodes.empty()) {
        auto current_node = nodes.top();
        nodes.pop();

        // std::cout << "current_node: " << compute_solution_value(current_node.ub) << std::endl;
        if (!is_feasible(current_node.lb) || compute_solution_value(current_node.ub) < best_solution_value) {
            continue;
        }
        const auto end = std::chrono::high_resolution_clock::now();
        const auto used_time = std::chrono::duration<double>(end - start).count();
        if (used_time - last_output_time >= 60) {
            std::cout << "used time is " << used_time << " seconds\n";
            last_output_time = used_time;
            std::cout << "remained nodes: " << get_remained_node_count() << "\n";
        }
        if (has_time_limit(time_limit) && used_time >= time_limit) {
            std::cout << "time limit (" << time_limit << " seconds) is reached." << std::endl;
            time_limit_reached = true;
            return;
        }

        branch(current_node);
    }
}

BabSolver::SolutionSet BabSolver::find_local_optimum(const BabSolver::SolutionSet& x, const BabSolver::SolutionSet &ub) const {
    auto local_optimum = x;

    while (is_feasible(local_optimum)) {
        auto max_obj = 0.0;
        int k = 0;
        bool found = false;
        for (int j = 0; j < x.size(); j++) {
            if (local_optimum[j] >= ub[j]) continue;

            auto x_up = local_optimum;
            ++x_up[j];

            if (is_feasible(x_up))
            {
                if (const auto obj = compute_solution_value(x_up); obj > max_obj) {
                    max_obj = obj;
                    k = j;
                    found = true;
                }
            }
        }

        if (found) { ++local_optimum[k]; }
        else { break; }
    }

    return local_optimum;
}


BabSolver::SolutionSet BabSolver::compute_upper_limit(const BabSolver::SolutionSet &x, const BabSolver::SolutionSet &ub) const {
    auto x_x = x;

    for (int j = 0; j < x.size(); ++j)
    {
        if (x[j] >= ub[j]) { continue; }

        auto x_up = x;
        ++x_up[j];
        while (is_feasible(x_up) && x_up[j] <= ub[j]) {
            x_x[j] = x_up[j];
            ++x_up[j];
        }
    }

    return x_x;
}

void BabSolver::branch(const BabNode& current_node) {
#ifndef NDEBUG
    std::cout << "Branching..." << std::endl;
#endif

    const auto lb = current_node.lb;
    const auto ub = current_node.ub;

#ifndef NDEBUG
    std::cout << "LB:\t";
    for (const auto m : lb) {
        std::cout << m << "\t";
    }
    std::cout << "\n" << "UB:\t";
    for (const auto m : ub) {
        std::cout << m << "\t";
    }
    std::cout << std::endl;
#endif

    const auto local_optimum = find_local_optimum(lb, ub);
#ifndef NDEBUG
    std::cout << "local_optimum:\t";
    for (const auto m : local_optimum) {
        std::cout << m << "\t";
    }
    std::cout << std::endl;
#endif

    if (const auto local_optimum_value = compute_solution_value(local_optimum);
        local_optimum_value > best_solution_value)
    {
        best_solution = local_optimum;
        best_solution_value = local_optimum_value;
    }

    for (int j = 0; j < nbSubSystems; ++j)
    {
        auto lb_new = compute_lb_x(local_optimum, lb, j);
        auto ub_new = compute_ub_minus(local_optimum, ub, j);

        if (!is_strict_dominated(lb_new, ub_new) || !is_feasible(lb_new) ||
            compute_solution_value(ub_new) < best_solution_value)
        {
            // prune this sub-problem
            continue;
        }

        if (auto ub_computed = compute_upper_limit(lb_new, ub_new);
            !is_strict_dominated(lb_new, ub_computed) || compute_solution_value(ub_computed) < best_solution_value)
        {
            // prune this sub-problem
            continue;
        }

        // add this new sub-problem to the queue
        auto node = BabNode(lb_new, ub_new, ++num_nodes);
        nodes.emplace(node);
    }
}

BabSolver::SolutionSet BabSolver::compute_lb_x(const BabSolver::SolutionSet &x, const BabSolver::SolutionSet &lb, int j) {
    auto x_new = SolutionSet(x.size());
    for(int i = 0; i < j; ++i)
    {
        x_new[i] = x[i];
    }
    for(int i = j; i < x_new.size(); ++i)
    {
        x_new[i] = lb[i];
    }

    return x_new;
}

BabSolver::SolutionSet BabSolver::compute_ub_minus(const BabSolver::SolutionSet &x, const BabSolver::SolutionSet &ub, int j) {
    auto x_new = ub;

    x_new[j] = x[j] - 1;

    return x_new;
}

double BabSolver::compute_solution_value(const BabSolver::SolutionSet &x) const {
    return get_system_reliability(x);
}

bool BabSolver::is_feasible(const BabSolver::SolutionSet &x) const {
    /* check variable (upper) bounds */
    for (int j = 0; j < x.size(); ++j) {
        if (x[j] > ub_initial[j]) return false;
    }

    /* check resource constraint if any */
    const auto result = check_resource_constraints(x);

    return result;
}

double BabSolver::get_system_reliability(const std::vector<int> &x) const {
    if (problem.test_type == TEST_TYPE::EXAMPLE_ONE) {
        std::vector<double> r1{0.000, 0.800, 0.850, 0.900, 0.925, 0.950, 0.975};
        auto R1 = r1[x[0]];
        auto R2 = parallel_reliability(x[1], 0.75);
        auto R3 = k_out_of_n_reliability(x[2] + 1, 2, 0.88);
        auto R4 = parallel_reliability(x[3], 0.70);
        auto R5 = parallel_reliability(x[4], 0.85);

        auto Q1 = 1.0 - R1;
        auto Q2 = 1.0 - R2;
        auto Q3 = 1.0 - R3;
        auto Q4 = 1.0 - R4;
        auto Q5 = 1.0 - R5;

        return R5 * (1 - Q1 * Q3) * (1 - Q2 * Q4) + Q5 * (1.0 - (1.0 - R1 * R2) * (1.0 - R3 * R4));
    }
    if (problem.test_type == TEST_TYPE::EXAMPLE_TWO) {
        auto R1 = parallel_reliability(x[0], 0.83);
        auto R2 = parallel_reliability(x[1], 0.89);
        auto R3 = parallel_reliability(x[2], 0.92);
        auto R4 = parallel_reliability(x[3], 0.85);
        auto R5 = parallel_reliability(x[4], 0.89);
        auto R6 = parallel_reliability(x[5], 0.93);
        auto R7 = parallel_reliability(x[6], 0.83);
        auto R8 = parallel_reliability(x[7], 0.94);
        auto R9 = parallel_reliability(x[8], 0.82);
        auto R10 = parallel_reliability(x[9], 0.91);

        auto Q3 = 1 - R3;
        auto Q7 = 1 - R7;
        auto Q8 = 1 - R8;
        auto Q9 = 1 - R9;

        return (1.0 - (1.0 - (1.0 - Q3 * (1.0 - R1 * R2)) * R4) * (1.0 - R5 * R6)) * (1.0 - Q7 * Q8 * Q9) * R10;
    }
    if (problem.test_type == TEST_TYPE::TEST_ONE) {
        std::vector<double> R(problem.num_subsystems);
        for (int j = 0; j < problem.num_subsystems; ++j) {
            std::vector<double> r(problem.num_component_types);
            for (int h = 0; h < problem.num_component_types; ++h) {
                r[h] = parallel_reliability(x[j * problem.num_component_types + h],
                                            problem.component_reliability_at_subsystem[j][h]);
            }
            R[j] = parallel_reliability(r);
        }

        std::vector<double> Q(problem.num_subsystems);
        for (int j = 0; j < problem.num_subsystems; ++j) {
            Q[j] = 1.0 - R[j];
        }

        constexpr auto num_blocks = 5;
        auto num_sub_per_block = problem.num_subsystems / num_blocks;
        std::vector<double> RR(num_blocks);

        for (auto j = 0u; j < num_blocks; ++j) {
            RR[j] = 1.0;
            for (auto i = 0u; i < num_sub_per_block; ++i) {
                RR[j] *= R[num_blocks * i + j];
            }
        }
        std::vector<double> QQ(num_blocks);
        for (auto j = 0u; j < num_blocks; ++j) {
            QQ[j] = 1.0 - RR[j];
        }

        const auto S = RR[4] * (1 - QQ[0] * QQ[2]) * (1 - QQ[1] * QQ[3]) +
                       QQ[4] * (1 - (1 - RR[0] * RR[1]) * (1 - RR[2] * RR[3]));
        // std::cout << "S: " << S << std::endl;
        return S;
    }
    if (problem.test_type == TEST_TYPE::TEST_TWO) {
        std::vector<double> R(problem.num_subsystems);
        for (int j = 0; j < problem.num_subsystems; ++j) {
            R[j] = parallel_reliability(x[j], problem.component_reliability_at_subsystem[j][0]);
        }

        std::vector<double> Q(problem.num_subsystems);
        for (int j = 0; j < problem.num_subsystems; ++j) {
            Q[j] = 1.0 - R[j];
        }

        constexpr auto num_sub_per_hsp = 5;
        const auto num_hsp = problem.num_subsystems / num_sub_per_hsp;
        double S = 1.0;
        for (auto i = 0u; i < num_hsp; ++i) {
            S *= 1 - (1 - Q[i * num_sub_per_hsp + 0u] * Q[i * num_sub_per_hsp + 1u] * R[i * num_sub_per_hsp + 2u]) *
                         (1 - R[i * num_sub_per_hsp + 3u] * R[i * num_sub_per_hsp + 4u]);
        }
        return S;
    }
    if (problem.test_type == TEST_TYPE::MIXED) {
        std::vector<double> R(problem.num_subsystems);
        for (int j = 0; j < problem.num_subsystems; ++j) {
            std::vector<double> r(problem.num_component_types);
            for (int h = 0; h < problem.num_component_types; ++h) {
                r[h] = parallel_reliability(x[j * problem.num_component_types + h], problem.component_reliability_at_subsystem[j][h]);
            }
            R[j] = parallel_reliability(r);
        }

        std::vector<double> Q(problem.num_subsystems);
        for (int j = 0; j < problem.num_subsystems; ++j) {
            Q[j] = 1.0 - R[j];
        }

        switch (problem.system_index) {
            case 1: // 5 subsystems
            {
                const auto S1 = R[4] * (1 - Q[0] * Q[2]) * (1 - Q[1] * Q[3]) +
                          Q[4] * (1 - (1 - R[0] * R[1]) * (1 - R[2] * R[3]));
                return S1;
            }
            case 2: // 5 subsystems
            {
                const auto S2 = R[4] * (1 - Q[1] * Q[3]) + Q[4] * (1 - (1 - R[0] * R[1]) * (1 - R[2] * R[3]));
                return S2;
            }
            case 3: // 6 subsystems
            {
                const auto S3 = R[5] * (1 - Q[0] * Q[2] * Q[3]) * (1 - Q[1] * Q[4]) +
                          Q[5] * (1 - (1 - R[0] * R[1]) * (1 - (1 - Q[2] * Q[3]) * R[4]));
                return S3;
            }
            case 4: // 7 subsystems
            {
                const auto Q234 = Q[2] + R[2] * Q[3] * Q[4];
                const auto R234 = 1 - Q234;
                const auto S4 = Q[6] * (1 - (1 - R[0] * R[1]) * (1 - R234 * R[5])) +
                          R[6] * (1 - Q[0] * Q234) * (1 - Q[1] * Q[5]);
                return S4;
            }
            case 5: { // 7 subsystems
                const auto S5 = Q[2] * Q[4] * (1 - (1 - R[0] * R[3] * R[5]) * (1 - R[1] * R[6]))
                          + R[2] * R[4] * (1 - Q[0] * Q[1]) * (1 - Q[5] * Q[6])
                          + R[2] * Q[4] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * R[5]) * Q[6])
                          + Q[2] * R[4] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[5] * Q[6]);
                return S5;
            }
            case 6: { // 8 subsystems
                const auto S6 = Q[3] * Q[5] * (1 - (1 - R[0] * R[2] * R[6]) * (1 - R[1] * R[4] * R[7]))
                          + R[3] * R[5] * (1 - Q[0] * (1 - R[1] * R[4])) * (1 - Q[6] * Q[7])
                          + R[3] * Q[5] * (1 - Q[0] * (1 - R[1] * R[4])) * (1 - (1 - R[2] * R[6]) * Q[7])
                          + Q[3] * R[5] * (1 - (1 - R[0] * R[2]) * (1 - R[1] * R[4])) * (1 - Q[6] * Q[7]);
                return S6;
            }
            case 7: { // 8 subsystems
                const auto S7 = Q[3] * Q[4] * (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[6]) * (1 - R[2] * R[7]))
                          + R[3] * R[4] * (1 - Q[0] * Q[1] * Q[2]) * (1 - Q[5] * Q[6] * Q[7])
                          + R[3] * Q[4] * (1 - (1 - R[2] * R[7]) * (1 - (1 - Q[0] * Q[1]) * (1 - Q[5] * Q[6])))
                          + Q[3] * R[4] * (1 - (1 - R[0] * R[5]) * (1 - (1 - Q[1] * Q[2]) * (1 - Q[6] * Q[7])));
                return S7;
            }
            case 8: { // 9 subsystems
                const auto S8 = Q[4] * Q[2] * Q[6] * (1 - (1 - R[0] * R[3] * R[7]) * (1 - R[1] * R[5] * R[8]))
                          + Q[4] * R[2] * Q[6] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * R[7]) * (1 - R[5] * R[8]))
                          + Q[4] * Q[2] * R[6] * (1 - (1 - R[0] * R[3]) * (1 - R[1] * R[5])) * (1 - Q[7] * Q[8])
                          + Q[4] * R[2] * R[6] * (1 - Q[0] * Q[1]) * (1 - Q[3] * Q[5]) * (1 - Q[7] * Q[8])
                          + R[4] * Q[2] * Q[6] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - (1 - R[5] * R[8]) * Q[7])
                          + R[4] * R[2] * Q[6] * (1 - Q[0] * Q[1]) * (1 - (1 - R[5] * R[8]) * Q[7])
                          + R[4] * Q[2] * R[6] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[7] * Q[8])
                          + R[4] * R[2] * R[6] * (1 - Q[0] * Q[1]) * (1 - Q[7] * Q[8]);
                return S8;
            }
            case 9: // 10 subsystems
            {
                const auto S9 = (1 - (1 - (1 - Q[2] * (1 - R[0] * R[1])) * R[3]) * (1 - R[4] * R[5])) * (1 - Q[6] * Q[7] * Q[8]) *
                          R[9];
                return S9;
            }
            case 10: { // 10 subsystems
                constexpr auto num_sub_per_hsp = 5;
                const auto num_hsp = problem.num_subsystems / num_sub_per_hsp;
                double S10 = 1.0;
                for (auto i = 0u; i < num_hsp; ++i) {
                    S10 = S10 * (1 -
                            (1 - Q[i * num_sub_per_hsp + 0u] * Q[i * num_sub_per_hsp + 1u] * R[i * num_sub_per_hsp + 2u]) *
                                 (1 - R[i * num_sub_per_hsp + 3u] * R[i * num_sub_per_hsp + 4u]));
                }
                return S10;
            }
            case 11: { // 12 subsystems
                const auto S11 = Q[2] * Q[4] * Q[9] *
                           (1 - (1 - R[0] * R[3] * (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[1] * R[5] * R[11]))
                           + Q[2] * Q[4] * R[9] * Q[7] * Q[8] *
                             (1 - (1 - R[0] * R[3] * R[6] * R[10]) * (1 - R[1] * R[5] * R[11]))
                           + Q[2] * Q[4] * R[9] * Q[7] * R[8] * (1 - (1 - R[0] * R[3] * R[6]) * (1 - R[1] * R[5])) *
                             (1 - Q[10] * Q[11])
                           + Q[2] * Q[4] * R[9] * R[7] * Q[8] * (1 - (1 - R[0] * R[3]) * (1 - R[1] * R[5])) *
                             (1 - (1 - R[6] * R[10]) * Q[11])
                           + Q[2] * Q[4] * R[9] * R[7] * R[8] * (1 - (1 - R[0] * R[3]) * (1 - R[1] * R[5])) *
                             (1 - Q[10] * Q[11])
                           + Q[2] * R[4] * Q[9] * (1 - (1 - R[0] * R[3]) * Q[1]) *
                             (1 - (1 - (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[5] * R[11]))
                           + Q[2] * R[4] * R[9] * Q[7] * Q[8] * (1 - (1 - R[0] * R[3]) * Q[1]) *
                             (1 - (1 - R[6] * R[10]) * (1 - R[5] * R[11]))
                           + Q[2] * R[4] * R[9] * Q[7] * R[8] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[5] * Q[6]) *
                             (1 - Q[10] * Q[11])
                           + Q[2] * R[4] * R[9] * R[7] * Q[8] * (1 - (1 - R[0] * R[3]) * Q[1]) *
                             (1 - (1 - R[6] * R[10]) * Q[11])
                           + Q[2] * R[4] * R[9] * R[7] * R[8] * (1 - (1 - R[0] * R[3]) * Q[1]) * (1 - Q[10] * Q[11])
                           + R[2] * Q[4] * Q[9] * (1 - Q[0] * Q[1]) *
                             (1 - (1 - R[3] * (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * Q[4] * R[9] * Q[7] * Q[8] * (1 - Q[0] * Q[1]) *
                             (1 - (1 - R[3] * R[6] * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * Q[4] * R[9] * Q[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - (1 - R[3] * R[6]) * Q[5]) *
                             (1 - Q[10] * Q[11])
                           + R[2] * Q[4] * R[9] * R[7] * Q[8] * (1 - Q[0] * Q[1]) * (1 - Q[3] * Q[5]) *
                             (1 - (1 - R[6] * R[10]) * Q[11])
                           + R[2] * Q[4] * R[9] * R[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - Q[3] * Q[5]) * (1 - Q[10] * Q[11])
                           + R[2] * R[4] * Q[9] * (1 - Q[0] * Q[1]) *
                             (1 - (1 - (1 - Q[6] * (1 - R[7] * R[8])) * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * R[4] * R[9] * Q[7] * Q[8] * (1 - Q[0] * Q[1]) *
                             (1 - (1 - R[6] * R[10]) * (1 - R[5] * R[11]))
                           + R[2] * R[4] * R[9] * Q[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - Q[5] * Q[6]) * (1 - Q[10] * Q[11])
                           + R[2] * R[4] * R[9] * R[7] * Q[8] * (1 - Q[0] * Q[1]) * (1 - (1 - R[6] * R[10]) * Q[11])
                           + R[2] * R[4] * R[9] * R[7] * R[8] * (1 - Q[0] * Q[1]) * (1 - Q[10] * Q[11]);
                return S11;
            }
            case 12: { // 12 subsystems
                const auto S12 = (1 - R[6] * R[7]) * (1 - R[3] * R[4]) * Q[9] *
                           (1 - (1 - R[0] * R[5] * R[10]) * (1 - R[1] * R[2] * R[8] * R[11]))
                           + (1 - R[6] * R[7]) * R[3] * R[4] * Q[9] * (1 - Q[0] * (1 - R[1] * R[2])) *
                             (1 - (1 - R[5] * R[10]) * (1 - R[8] * R[11]))
                           + (1 - R[6] * R[7]) * (1 - R[3] * R[4]) * R[9] *
                             (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[2] * R[8])) * (1 - Q[10] * Q[11])
                           + (1 - R[6] * R[7]) * R[3] * R[4] * R[9] * (1 - Q[0] * (1 - R[1] * R[2])) * (1 - Q[5] * Q[8]) *
                             (1 - Q[10] * Q[11])
                           + R[6] * R[7] * (1 - R[3] * R[4]) * Q[9] * (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[2])) *
                             (1 - (1 - R[8] * R[11]) * Q[10])
                           + R[6] * R[7] * R[3] * R[4] * Q[9] * (1 - Q[0] * (1 - R[1] * R[2])) *
                             (1 - (1 - R[8] * R[11]) * Q[10])
                           + R[6] * R[7] * (1 - R[3] * R[4]) * R[9] * (1 - (1 - R[0] * R[5]) * (1 - R[1] * R[2])) *
                             (1 - Q[10] * Q[11])
                           + R[6] * R[7] * R[3] * R[4] * R[9] * (1 - Q[0] * (1 - R[1] * R[2])) * (1 - Q[10] * Q[11]);
                return S12;
            }
            default: {
                std::cout << "ERROR: invalid number of subsystems: " << problem.system_index << std::endl;
                return 0.0;
            }
        }
    }

    std::cout << "ERROR: invalid number of subsystems: " << problem.system_index << std::endl;
    return 0.0;
}


bool BabSolver::is_proper_dominated(const BabSolver::SolutionSet &x1, const BabSolver::SolutionSet &x2) {
    if(x1.size() != x2.size())
    {
        return false;
    }

    for (size_t i = 0; i < x1.size(); ++i){
        if (x1[i] > x2[i]) return false;
    }

    return true;
}

bool BabSolver::is_strict_dominated(const BabSolver::SolutionSet &x1, const BabSolver::SolutionSet &x2) {
    if(!is_proper_dominated(x1, x2))
    {
        return false;
    }

    for (size_t i = 0; i < x1.size(); ++i){
        if (x1[i] < x2[i]) return true;
    }

    return true;
}

bool BabSolver::is_strong_dominated(const BabSolver::SolutionSet &x1, const BabSolver::SolutionSet &x2) {
    if(x1.size() != x2.size())
    {
        return false;
    }

    for (size_t i = 0; i < x1.size(); ++i){
        if (x1[i] >= x2[i]) return false;
    }

    return true;
}
