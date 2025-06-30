#include "mm_reliability.h"
#include "mm_solution.h"
#include "mm_solver.h"
#include <algorithm>
#include <cmath>
#include <iterator>
#include <optional>
#include <chrono>
#include "gurobi_c++.h"

namespace rrap
{
    SOLUTION_METHOD get_solution_method_code(const int method_identifier) {
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::Q_BASIC)) {
            return SOLUTION_METHOD::Q_BASIC;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::Q_BASIC_LOG_X)) {
            return SOLUTION_METHOD::Q_BASIC_LOG_X;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::Q_TREE)) {
            return SOLUTION_METHOD::Q_TREE;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::Q_TREE_LOG_X)) {
            return SOLUTION_METHOD::Q_TREE_LOG_X;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::Q_OB_TREE)) {
            return SOLUTION_METHOD::Q_OB_TREE;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::Q_OB_TREE_LOG_X)) {
            return SOLUTION_METHOD::Q_OB_TREE_LOG_X;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::GENERAL_BASIC)) {
            return SOLUTION_METHOD::GENERAL_BASIC;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::GENERAL_LOG_X)) {
            return SOLUTION_METHOD::GENERAL_LOG_X;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::GENERAL_TREE)) {
            return SOLUTION_METHOD::GENERAL_TREE;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::GENERAL_TREE_LOG_X)) {
            return SOLUTION_METHOD::GENERAL_TREE_LOG_X;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::GENERAL_OB_TREE)) {
            return SOLUTION_METHOD::GENERAL_OB_TREE;
        }
        if (method_identifier == static_cast<int>(SOLUTION_METHOD::GENERAL_OB_TREE_LOG_X)) {
            return SOLUTION_METHOD::GENERAL_OB_TREE_LOG_X;
        }
        return SOLUTION_METHOD::ERROR_METHOD;
    }

    std::string get_solution_method_name(const SOLUTION_METHOD method_code) {
        switch (method_code) {
            case SOLUTION_METHOD::Q_BASIC: {
                return "Q_BASIC";
            }
            case SOLUTION_METHOD::Q_BASIC_LOG_X: {
                return "Q_BASIC_LOG_X";
            }
            case SOLUTION_METHOD::Q_TREE: {
                return "Q_TREE";
            }
            case SOLUTION_METHOD::Q_TREE_LOG_X: {
                return "Q_TREE_LOG_X";
            }
            case SOLUTION_METHOD::Q_OB_TREE: {
                return "Q_OB_TREE";
            }
            case SOLUTION_METHOD::Q_OB_TREE_LOG_X: {
                return "Q_OB_TREE_LOG_X";
            }
            case SOLUTION_METHOD::GENERAL_BASIC: {
                return "GENERAL_BASIC";
            }
            case SOLUTION_METHOD::GENERAL_LOG_X: {
                return "GENERAL_BASIC_LOG_X";
            }
            case SOLUTION_METHOD::GENERAL_TREE: {
                return "GENERAL_TREE";
            }
            case SOLUTION_METHOD::GENERAL_TREE_LOG_X: {
                return "GENERAL_TREE_LOG_X";
            }
            case SOLUTION_METHOD::GENERAL_OB_TREE: {
                return "GENERAL_OB_TREE";
            }
            case SOLUTION_METHOD::GENERAL_OB_TREE_LOG_X: {
                return "GENERAL_OB_TREE_LOG_X";
            }
            default: {
                return "error_method";
            }
        }
    }

    std::string get_solution_method_name(const int method_identifier) {
        const auto method_code = get_solution_method_code(method_identifier);
        const auto method_name = get_solution_method_name(method_code);
        return method_name;
    }

    std::vector<std::vector<int>> compute_num_x(const ProbInstance& inst) {
        auto num_x = std::vector<std::vector<int>>(inst.num_subsystems);
        for (auto j = 0u; j < inst.num_subsystems; ++j) {
            num_x[j] = std::vector<int>(inst.num_component_types);
            for (auto h = 0u; h < inst.num_component_types; ++h) {
                num_x[j][h] = static_cast<int>(std::ceil(std::log2(inst.component_ub_at_subsystem[j][h] + 1)));
            }
        }

        return num_x;
    }

    std::vector<std::vector<std::vector<double>>> compute_num_k(const ProbInstance& inst,
        const std::vector<std::vector<int>>& num_x)
    {
        auto num_k = std::vector<std::vector<std::vector<double>>>(inst.num_subsystems);
        for (auto j = 0u; j < inst.num_subsystems; ++j) {
            num_k[j] = std::vector<std::vector<double>>(inst.num_subsystems);
            for (auto h = 0u; h < inst.num_component_types; ++h) {
                num_k[j][h] = std::vector<double>(num_x[j][h]);
                auto sum = 0.0;

                for (auto k = 0u; k < num_x[j][h] - 1u; ++k) {
                    num_k[j][h][k] = std::pow(2, k);
                    sum += num_k[j][h][k];
                }

                num_k[j][h][num_x[j][h] - 1u] = static_cast<double>(inst.component_ub_at_subsystem[j][h]) - sum;
            }
        }

        return num_k;
    }

    void add_mixed_log_x_constraints(
        const ProbInstance& inst,
        GRBModel& model,
        const std::vector<std::vector<std::vector<GRBVar>>>& x,
        const std::vector<std::vector<int>>& num_x,
        const std::vector<std::vector<std::vector<double>>>& num_k)
    {
        for (auto i = 0u; i < inst.num_resource_types; ++i)
        {
            GRBLinExpr expr = 0;

            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    for (auto k = 0u; k < num_x[j][h]; ++k)
                    {
                        expr += inst.resource_for_component_at_subsystem[i][j][h] * num_k[j][h][k] * x[j][h][k];
                    }
                }
            }

            model.addConstr(expr <= inst.amount_of_each_resources[i]);
        }
    }

    void add_mixed_constraints(const ProbInstance& inst, GRBModel& model,
        const std::vector<std::vector<std::vector<GRBVar>>>& x)
    {
        for (auto i = 0u; i < inst.num_resource_types; ++i) {
            GRBLinExpr expr = 0;

            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += inst.resource_for_component_at_subsystem[i][j][h] * static_cast<double>(k) *
                                x[j][h][k];
                    }
                }
            }

            model.addConstr(expr <= inst.amount_of_each_resources[i]);
        }
    }

    void add_test_one_constraints(const ProbInstance& inst, GRBModel& model,
        const std::vector<std::vector<std::vector<GRBVar>>>& x)
    {
        {// 1st constraint
            GRBLinExpr expr = 0;

            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += inst.resource_for_component_at_subsystem[0u][j][h] * static_cast<double>(k) *
                                x[j][h][k];
                    }
                }
            }

            model.addConstr(expr <= inst.amount_of_each_resources[0u]);
        }
        {// 2nd constraint
            GRBLinExpr expr = 0;

            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += inst.resource_for_component_at_subsystem[1u][j][h] *
                                (static_cast<double>(k) + std::exp(static_cast<double>(k) / 4.0)) *
                                x[j][h][k];
                    }
                }
            }

            model.addConstr(expr <= inst.amount_of_each_resources[1u]);
        }
        {// 3rd constraint
            GRBLinExpr expr = 0;

            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += inst.resource_for_component_at_subsystem[2u][j][h] *
                                std::pow(static_cast<double>(k), 3.0) *
                                x[j][h][k];
                                std::pow(static_cast<double>(k), 3.0);
                        std::cout << inst.resource_for_component_at_subsystem[2u][j][h] *
                                std::pow(static_cast<double>(k), 3.0) << std::endl;
                    }
                }
            }

            model.addConstr(expr <= inst.amount_of_each_resources[2u]);
        }
    }

    void add_test_two_constraints(const ProbInstance& inst, GRBModel& model,
        const std::vector<std::vector<std::vector<GRBVar>>>& x) {
        {// 1st constraint
            GRBLinExpr expr = 0;

            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += inst.resource_for_component_at_subsystem[0u][j][h] * static_cast<double>(k) *
                                x[j][h][k];
                    }
                }
            }

            model.addConstr(expr <= inst.amount_of_each_resources[0u]);
        }
        {// 2nd constraint
            GRBLinExpr expr = 0;

            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += inst.resource_for_component_at_subsystem[1u][j][h] *
                                std::pow(static_cast<double>(k), 2.0) *
                                x[j][h][k];
                    }
                }
            }

            model.addConstr(expr <= inst.amount_of_each_resources[1u]);
        }
        {// 3rd constraint
            GRBLinExpr expr = 0;

            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += inst.resource_for_component_at_subsystem[2u][j][h] *
                                static_cast<double>(k) * std::exp(static_cast<double>(k) / 4.0) *
                                x[j][h][k];
                    }
                }
            }

            model.addConstr(expr <= inst.amount_of_each_resources[2u]);
        }
    }

    auto ProbSolver::solve(const int method) const -> std::optional<ProbSolution> {
        if (method == 110) { return solve_Q_basic_model(); }
        if (method == 111) { return solve_Q_basic_logx_model(); }
        if (method == 120) { return solve_Q_tree_model(inst.logic_system_tree.systems); }
        if (method == 121) { return solve_Q_tree_logx_model(inst.logic_system_tree.systems); }
        if (method == 130) { return solve_Q_tree_model(inst.sorted_logic_system_tree.systems); }
        if (method == 131) { return solve_Q_tree_logx_model(inst.sorted_logic_system_tree.systems); }
        if (method == 210) { return solve_general_basic_model(); }
        if (method == 211) { return solve_general_basic_logx_model(); }
        if (method == 220) { return solve_general_tree_model(inst.qr_system_tree.systems); }
        if (method == 221) { return solve_general_tree_logx_model(inst.qr_system_tree.systems); }
        if (method == 230) { return solve_general_tree_model(inst.qr_sorted_system_tree.systems); }
        if (method == 231) { return solve_general_tree_logx_model(inst.qr_sorted_system_tree.systems); }
        return std::nullopt;
    }

    /* Q-form */
    std::optional<ProbSolution> ProbSolver::solve_Q_basic_model() const
    {
        /* the number of logical subsystem */
        const auto num_logic_subs { inst.system_reliability.m_other_terms.size() };
#ifndef  NDEBUG
        std::cout << "NUM_LOGIC_SUBSYSTEMS: " << num_logic_subs << std::endl;
#endif

        /* save logical-subsystems in a vector for following usage */
        std::vector<std::pair<std::set<int>, int>> logic_subsystems;  /* all logical subsystems */
        std::vector<size_t> num_subs_in_logic;     /* number of subsystems in each logic subsystem */
        logic_subsystems.reserve(num_logic_subs);
        num_subs_in_logic.reserve(num_logic_subs);
        for (auto& [term, coeff] : inst.system_reliability.m_other_terms){
            logic_subsystems.emplace_back(term, coeff);
            num_subs_in_logic.push_back(term.size());
        }
#ifndef NDEBUG
        std::cout << ">>>>> SUM = " << num_logic_subs << std::endl;
        for (auto l = 0; l < num_logic_subs; ++l)
        {
            std::cout << l << " =====> " << num_subs_in_logic[l] << "\t" << logic_subsystems[l].second << std::endl;
        }
#endif

        /* the subsystems in each logic subsystem */
        std::vector<std::vector<size_t>> subsystems(num_logic_subs);
        for (auto l = 0; l < num_logic_subs; ++l)
        {
            std::copy(logic_subsystems[l].first.begin(), logic_subsystems[l].first.end(),
                      std::back_inserter(subsystems[l]));
        }
#ifndef NDEBUG
        for (auto l = 0; l < num_logic_subs; ++l)
        {
            std::cout << "logic subsystem [" << l << "]: " << std::endl;
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s)
            {
                std::cout << "\tsubsystem [" << s << "]: " << subsystems[l][s] << "\n";
            }
            std::cout << std::endl;
        }
#endif

        /* find the component type that has the minimum reliability at each subsystem */
        std::vector<std::vector<double>> min_r(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            min_r[l] = std::vector<double>(num_subs_in_logic[l]);
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                auto j = subsystems[l][s];
                min_r[l][s] = *std::min_element(inst.component_reliability_at_subsystem[j].begin(),
                                                inst.component_reliability_at_subsystem[j].end());
                // std::cout << "=> min_r[" << l << "][" << s << "] = " <<  min_r[l][s] << std::endl;
            }
        }

        bool combination_model_use_rho_ub = true;
        bool combination_model_use_valid_inequalities = true;
        bool combination_model_use_aggregated_expression = false;
        bool combination_model_use_dominated_components = true;
        bool combination_model_use_ub_one = false;

        /* calculate the upper bounds of decision variables \varphi_{ls} */
        std::vector<std::vector<double>> varphi_ub_(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            varphi_ub_[l] = std::vector<double>(num_subs_in_logic[l], 1.0);
        }
        for (auto l = 0u; l < num_logic_subs; ++l) {
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                varphi_ub_[l][s] = (s == 0 ? 1.0 : varphi_ub_[l][s - 1u]) * (1.0 - min_r[l][s]);
            }
        }

        /* calculate the upper bounds of decision variables \rho_{ls*} */
        std::vector<std::vector<std::vector<double>>> rho_ub(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            rho_ub[l] = std::vector<std::vector<double>>(num_subs_in_logic[l]);
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                rho_ub[l][s] = std::vector<double>(inst.num_component_types, 1.0);
            }
        }

        if (combination_model_use_rho_ub) {
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    for (int s2 = 0; s2 < s; ++s2) {
                        if (s != 0u) {
                            rho_ub[l][s] = std::vector<double>(inst.num_component_types, varphi_ub_[l][s - 1u]);
                        }

                        const auto h_last = inst.num_component_types - 1u;
                        rho_ub[l][s][h_last] *= 1 - min_r[l][s];
                    }
                }
            }
#ifndef  NDEBUG
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    for (auto h = 0u; h < inst.num_component_types; ++h) {
                        std::cout <<"rho_ub[" << l << "][" << s << "]["<< h << "] = " << rho_ub[l][s][h] << "\n";
                    }
                }
            }
#endif
        }

        /* sort the subsystems according to the number of their occurrence in logic subsystems */
#if 0
        std::map<int, int> num_occurrence;
        for (auto &subsystem: logic_subsystems) {
            for (auto &term: subsystem.first) {
                num_occurrence[term] += 1;
            }
        }

        std::vector<int> sorted_subsystems(inst.num_subsystems);
        std::iota(sorted_subsystems.begin(), sorted_subsystems.end(), 0);
        std::sort(sorted_subsystems.begin(), sorted_subsystems.end(),
                  [&num_occurrence](int a, int b) { return num_occurrence[a] < num_occurrence[b]; });

        std::vector<int> branch_priorities(inst.num_subsystems);
        for (int s = 0; s < inst.num_subsystems; ++s){
            branch_priorities[sorted_subsystems[s]] = s * 10;
        }
#ifndef  NDEBUG
        for (int x = 0; x < inst.num_subsystems; ++x) {
            std::cout << "sorted_subsystems[" << x << "]=" << sorted_subsystems[x] << std::endl;
        }
#endif
#ifndef NDEBUG
        for (int x = 0; x < inst.num_subsystems; ++x) {
            std::cout << "num_occurrence[" << x << "]=" << num_occurrence[x] << std::endl;
        }
#endif
#endif

        /* model and solve */

        try
        {
            GRBEnv env{};
            GRBModel model{env};

            // variables

            auto rho = std::vector<std::vector<std::vector<GRBVar>>>(num_logic_subs);
            for (auto l = 0u; l < num_logic_subs; ++l) {
                rho[l] = std::vector<std::vector<GRBVar>>(num_subs_in_logic[l]);
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    rho[l][s] = std::vector<GRBVar>(inst.num_component_types);
                    for (auto h = 0u; h < inst.num_component_types; ++h) {
                        rho[l][s][h] = model.addVar(0.0, rho_ub[l][s][h], 0.0, GRB_CONTINUOUS);
                    }
                }
            }

            auto x = std::vector<std::vector<std::vector<GRBVar>>>(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                x[j] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    x[j][h] = std::vector<GRBVar>(inst.component_ub_at_subsystem[j][h] + 1u);
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k)
                    {
                        x[j][h][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                        // x[j][h][k].set(GRB_IntAttr_BranchPriority, branch_priorities[j]);
                    }
                }
            }

            /* fixing variable according to domination between different component types */
            if (combination_model_use_dominated_components) {
                for(auto [j, h] : inst.dominated_components) {
                    x[j][h][0].set(GRB_DoubleAttr_LB, 1.0);
                    for (auto k = 1u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        x[j][h][k].set(GRB_DoubleAttr_UB, 0.0);
                    }
                }
            }


            // Constraints:

            /*
             * boundary condition:
             * for each logic subsystem l, calculate the reliability of the first subsystem (s = 0) after the first type
             * of components (h = 0) is configured.
             * */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                const auto j = subsystems[l][0u];
                GRBLinExpr expr = 0;

                for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                    expr += compute_jhk(j, 0u, k) * x[j][0u][k];
                }

                model.addConstr(rho[l][0u][0u] == expr);
            }

            /* recursive relation between different subsystems of each logic subsystem */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 1u; s < num_subs_in_logic[l]; ++s) {
                    const auto j = subsystems[l][s]; // j = \delta(s)
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                        const auto h_last = inst.num_component_types - 1u;
                        // model.addGenConstrIndicator(x[j][0u][k], 1, rho[l][s][0u] == rho[l][s - 1u][h_last] * std::pow(1.0 - inst.component_reliability_at_subsystem[j][0u], static_cast<double>(k)));
                        if (logic_subsystems[l].second > 0) {
                            const auto big_m = rho_ub[l][s - 1u][h_last] * (1 - compute_jhk(j,0u,k));//compute_comb_big_m(subsystems, l, s, 0u, k);
                            model.addConstr(rho[l][s][0u] <= rho[l][s - 1u][h_last] * compute_jhk(j,0u,k) + big_m * (1 - x[j][0u][k]));
                        } else {
                            const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][0u]);
                            const auto big_m = rho_ub[l][s - 1u][h_last] * (compute_jhk(j,0u,k) - compute_jhk(j,0u,K)); //compute_comb_big_m(subsystems, l, s, 0u, k);
                            model.addConstr(rho[l][s][0u] >= rho[l][s - 1u][h_last] * compute_jhk(j,0u,k) - big_m * (1 - x[j][0u][k]));
                        }
                    }
                }
            }

            // recursive relation within the same subsystem
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    const auto j = subsystems[l][s]; // j = \delta(s)
                    for (auto h = 1u; h < inst.num_component_types; ++h) {
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                            // model.addGenConstrIndicator(x[j][h][k], 1, rho[l][s][h] == rho[l][s][h - 1u] * std::pow(1.0 - inst.component_reliability_at_subsystem[j][h], static_cast<double>(k)));

                            if (logic_subsystems[l].second > 0) {
                                auto big_m = rho_ub[l][s][h - 1u] * (1 - compute_jhk(j,h,k)); // compute_comb_big_m_positive(subsystems, l, s, h, k);
                                model.addConstr(rho[l][s][h] <= rho[l][s][h - 1u] * compute_jhk(j,h,k) + big_m * (1 - x[j][h][k]));
                            } else {
                                const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][h]);
                                const auto big_m = rho_ub[l][s][h - 1u] * (compute_jhk(j,h,k) - compute_jhk(j,h,K)); // compute_comb_big_m_negative(subsystems, l, s, h, k);
                                model.addConstr(rho[l][s][h] >= rho[l][s][h - 1u] * compute_jhk(j,h,k) - big_m * (1 - x[j][h][k]));
                            }
                        }
                    }
                }
            }

            // sum_{k=0}^{u_jh} = 1, \forall j, h: exactly a particular number of components of type h is used
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    GRBLinExpr expr = 0;

                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += x[j][h][k];
                    }

                    model.addConstr(expr == 1);
                }
            }

            // at least one component is used
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                GRBLinExpr expr = 0;

                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 1u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += x[j][h][k];
                    }
                }

                model.addConstr(expr >= 1);
            }

            // resource constraints
            if (inst.test_type == TEST_TYPE::MIXED) {
                add_mixed_constraints(inst, model, x);
            }
            else if (inst.test_type == TEST_TYPE::TEST_ONE) {
                add_test_one_constraints(inst, model, x);
            }
            else if (inst.test_type == TEST_TYPE::TEST_TWO) {
                add_test_two_constraints(inst, model, x);
            }

            /* valid inequalities */
            if (combination_model_use_valid_inequalities) {
                // the value of \rho is non-increasing from component to component within each subsystem
                // for (auto l = 0u; l < num_logic_subs; ++l) {
                //     for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                //         for (auto h = 1u; h < inst.num_component_types; ++h) {
                //             model.addConstr(rho[l][s][h] <= rho[l][s][h - 1u]);
                //         }
                //     }
                // }

                // and a similar relationship exists among the transition between different subsystems
                // for (auto l = 0u; l < num_logic_subs; ++l) {
                //     for (auto s = 1u; s < num_subs_in_logic[l]; ++s) {
                //         const auto last_com = inst.num_component_types - 1u;
                //         model.addConstr(rho[l][s][0u] <= rho[l][s - 1u][last_com]);
                //     }
                // }

                // upper bound on the reliability of each subsystem after the components are configured
                for (auto l = 0u; l < num_logic_subs; ++l) {
                    for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                        const auto j = subsystems[l][s];
                        for (auto h = 1u; h < inst.num_component_types; ++h) {
                            GRBLinExpr expr = 0;
                            for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                                expr += compute_jhk(j, h, k) * x[j][h][k];
                            }

                            model.addConstr(rho[l][s][h] <= rho_ub[l][s][0u] * expr);
                        }
                    }
                }

                for (auto l1 = 0u; l1 < num_logic_subs; ++l1) {
                    std::set<size_t> a;
                    for (auto s1 = 0u; s1 < num_subs_in_logic[l1]; ++s1) {
                        a.insert(subsystems[l1][s1]);
                    }
                    for (auto l2 = 0u; l2 < num_logic_subs; ++l2) {
                        if (l1 == l2) continue;

                        std::set<size_t> b;
                        for (auto s2 = 0u; s2 < num_subs_in_logic[l2]; ++s2) {
                            b.insert(subsystems[l2][s2]);

                            // check if 'a' is a subset of 'b'
                            if (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                                const auto h_last = inst.num_component_types - 1u;
                                model.addConstr(rho[l2][s2][h_last] <= rho[l1][num_subs_in_logic[l1] - 1u][h_last]);
                                break;
                            }
                        }
                    }
                }
                //                for (auto l1 = 0u; l1 < num_logic_subs; ++l1) {
                //                    std::set<size_t> a;
                //                    for (auto s1 = 0u; s1 < num_subs_in_logic[l1]; ++s1) {
                //                        a.insert(subsystems[l1][s1]);
                //                        for (auto l2 = 0u; l2 < num_logic_subs; ++l2) {
                //                            if (l1 == l2) continue;
                //
                //                            std::set<size_t> b;
                //                            for (auto s2 = 0u; s2 < num_subs_in_logic[l2]; ++s2) {
                //                                b.insert(subsystems[l2][s2]);
                //
                //                                // check if 'a' is a subset of 'b'
                //                                if (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                //                                    const auto h_last = inst.num_component_types - 1u;
                //                                    model.addConstr(rho[l2][s2][h_last] <= rho[l1][s1][h_last]);
                //                                    break;
                //                                }
                //                            }
                //                        }
                //                    }
                //                }
            }


            /* Objective functions */

            {
                GRBLinExpr expr = inst.system_reliability.const_term;
                for (auto l = 0u; l < num_logic_subs; ++l)
                {
                    const auto last_sub = num_subs_in_logic[l] - 1u;
                    const auto last_com = inst.num_component_types - 1u;
                    expr += logic_subsystems[l].second * rho[l][last_sub][last_com];
                }

                if (combination_model_use_ub_one) {
                    model.addConstr(expr <= 1.0);
                }

                model.setObjective(expr, GRB_MAXIMIZE);
            }

            // Solve:

            // Time limit: 5 mins
            if (has_time_limit) {
                model.set(GRB_DoubleParam_TimeLimit, time_limit);
            }

            // Gap tolerance: 0.1%
            // model.set(GRB_DoubleParam_MIPGap, 0.001);

            // model.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
            // model.set(GRB_IntParam_NodeMethod, GRB_METHOD_BARRIER);

            model.write("model.lp");
            // model.write("model.mps");

            auto start_time = std::chrono::high_resolution_clock::now();

            model.optimize();

            auto end_time = std::chrono::high_resolution_clock::now();
            auto cpu_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            const auto status = model.get(GRB_IntAttr_Status);

            if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
            {
                std::cerr << "\t-- Model infeasible or unbounded!\n";
                return std::nullopt;
            }

            if (status != GRB_OPTIMAL)
            {
                std::cerr << "\t-- Warning: Solution not optimal!\n";
            }

            std::cout << "\t-- OBJ: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;

            // Extract solution:

            auto s = ProbSolution{inst};
            s.obj = model.get(GRB_DoubleAttr_ObjVal);
            s.cpu_time = cpu_time;

            if (status == GRB_TIME_LIMIT) {
                s.time_limit = 1;
                s.mip_gap = model.get(GRB_DoubleAttr_MIPGap);
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }
            else {
                s.time_limit = 0;
                s.mip_gap = 0.0;
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }

            s.num_components_at_subsystem.resize(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                s.num_components_at_subsystem[j].resize(inst.num_component_types);

                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    for (auto k = 0u; k < inst.component_ub_at_subsystem[j][h]; ++k)
                    {
                        if (x[j][h][k].get(GRB_DoubleAttr_X) > 0.5)
                        {
                            s.num_components_at_subsystem[j][h] = k;
                        }
                    }
                }
            }

            s.system_reliability = compute_system_reliability(s.num_components_at_subsystem);
            std::cout << "\t-- OBJ: " << s.system_reliability << std::endl;

            s.num_bin_vars = model.get(GRB_IntAttr_NumBinVars);
            s.num_vars = model.get(GRB_IntAttr_NumConstrs);
            s.num_cons = model.get(GRB_IntAttr_NumVars);
            s.num_nzs = model.get(GRB_IntAttr_NumNZs);
            s.num_bb_nodes = model.get(GRB_DoubleAttr_NodeCount);

            return s;
        }
        catch (GRBException &e)
        {
            std::cerr << "GUROBI ERROR:" << "\n\t" << e.getErrorCode() << "\n\t" << e.getMessage() << "\n";
        }
        catch (...)
        {
            std::cerr << "OTHER_ERRORS" << "\n";
        }

        return std::nullopt;
    }

    std::optional<ProbSolution> ProbSolver::solve_Q_basic_logx_model() const
    {
        /* the number of product terms */
        const auto num_logic_subs { inst.system_reliability.m_other_terms.size() };
        // std::cout << "NUM_LOGIC_SUBSYSTEMS: " << num_logic_subs << std::endl;

        /* save product terms in a vector for the following usage */
        std::vector<std::pair<std::set<int>, int>> logic_subsystems;  /* all product terms */
        std::vector<size_t> num_subs_in_logic;  /* number of subsystems in each product term */
        logic_subsystems.reserve(num_logic_subs);
        num_subs_in_logic.reserve(num_logic_subs);
        for (auto& [term, coeff] : inst.system_reliability.m_other_terms){
            logic_subsystems.emplace_back(term, coeff);
            num_subs_in_logic.push_back(term.size());
        }

        /* the subsystems in each product term */
        std::vector<std::vector<size_t>> subsystems(num_logic_subs);
        for (auto l = 0; l < num_logic_subs; ++l) {
            std::copy(logic_subsystems[l].first.begin(), logic_subsystems[l].first.end(),
                      std::back_inserter(subsystems[l]));
        }


        /* calculate the upper bounds of decision variables \rho_{ls*} */
        std::vector<std::vector<double>> upper_rho(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            upper_rho[l] = std::vector<double>(num_subs_in_logic[l], 1.0);
        }

        for (auto l = 0u; l < num_logic_subs; ++l) {
            for (auto s = 1u; s < num_subs_in_logic[l]; ++s) {
                for (int s2 = 0; s2 < s; ++s2) {
                    auto j = subsystems[l][s2];
                    auto pos = std::min_element(inst.component_reliability_at_subsystem[j].begin(),
                                                inst.component_reliability_at_subsystem[j].end());
                    auto h_min = static_cast<size_t>(pos - inst.component_reliability_at_subsystem[j].begin());
                    upper_rho[l][s] *= 1 - inst.component_reliability_at_subsystem[j][h_min];
                }
            }
        }

        auto num_x = compute_num_x(inst);
        auto num_k = compute_num_k(inst, num_x);


        /* model and solve */
        try
        {
            GRBEnv env{};
            GRBModel model{env};

            // variables
            auto rho = std::vector<std::vector<std::vector<std::vector<GRBVar>>>>(num_logic_subs);
            for (auto l = 0u; l < num_logic_subs; ++l) {
                rho[l] = std::vector<std::vector<std::vector<GRBVar>>>(num_subs_in_logic[l]);

                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    const auto j = subsystems[l][s]; // j = \delta(s)
                    rho[l][s] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                    for (auto h = 0u; h < inst.num_component_types; ++h) {
                        rho[l][s][h] = std::vector<GRBVar>(num_x[j][h]);
                        for (auto k = 0u; k < num_x[j][h]; ++k) {
                            rho[l][s][h][k] = model.addVar(0.0, upper_rho[l][s], 0.0, GRB_CONTINUOUS);
                        }
                    }
                }
            }

            auto x = std::vector<std::vector<std::vector<GRBVar>>>(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                x[j] = std::vector<std::vector<GRBVar>>(inst.num_component_types);

                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    x[j][h] = std::vector<GRBVar>(num_x[j][h]);

                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        x[j][h][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    }
                }
            }


            // Constraints

            /*
             * boundary condition
             * */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                const auto j = subsystems[l][0u];

                model.addConstr(rho[l][0u][0u][0u] == 1 -
                                                      (1 -
                                                       std::pow(1.0 - inst.component_reliability_at_subsystem[j][0u],
                                                                num_k[j][0u][0u])) * x[j][0u][0u]);
            }

            /* recursive relation between different numbers of each component-type of each subsystem of each term */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    auto j = subsystems[l][s];

                    for (auto h = 0u; h < inst.num_component_types; ++h) {
                        for (auto k = 1u; k < num_x[j][h]; ++k) {
                            auto big_m = 1.0;

                            if (logic_subsystems[l].second > 0) {
                                model.addConstr(rho[l][s][h][k] <= rho[l][s][h][k - 1u]);
                                model.addConstr(rho[l][s][h][k] <= rho[l][s][h][k - 1u] * std::pow(
                                        1.0 - inst.component_reliability_at_subsystem[j][h], num_k[j][h][k])
                                                                   + big_m * (1 - x[j][h][k]));
                            }
                            else {
                                model.addConstr(rho[l][s][h][k] >= rho[l][s][h][k - 1u] * std::pow(
                                        1.0 - inst.component_reliability_at_subsystem[j][h], num_k[j][h][k]));
                                model.addConstr(rho[l][s][h][k] >= rho[l][s][h][k - 1u] - big_m * x[j][h][k]);
                            }
                        }
                    }
                }
            }

            /* recursive relation between different components of each subsystem of each logical subsystem */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    auto j = subsystems[l][s];

                    for (auto h = 1u; h < inst.num_component_types; ++h) {
                        auto big_m = 1.0;

                        if (logic_subsystems[l].second > 0) {
                            model.addConstr(rho[l][s][h][0u] <= rho[l][s][h - 1u][num_x[j][h - 1u] - 1u]);
                            model.addConstr(rho[l][s][h][0u] <= rho[l][s][h - 1u][num_x[j][h - 1u] - 1u] * std::pow(
                                    1.0 - inst.component_reliability_at_subsystem[j][h], num_k[j][h][0u])
                                                                + big_m * (1 - x[j][h][0u]));
                        }
                        else {
                            model.addConstr(rho[l][s][h][0u] >= rho[l][s][h - 1u][num_x[j][h - 1u] - 1u] * std::pow(
                                    1.0 - inst.component_reliability_at_subsystem[j][h], num_k[j][h][0u]));
                            model.addConstr(rho[l][s][h][0u] >= rho[l][s][h - 1u][num_x[j][h - 1u] - 1u]
                                                                - big_m * x[j][h][0u]);
                        }
                    }
                }
            }

            /* recursive relation between different logic subsystems */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 1u; s < num_subs_in_logic[l]; ++s) {
                    auto j = subsystems[l][s];
                    auto jm = subsystems[l][s - 1u];
                    const auto h_last = inst.num_component_types - 1u;
                    auto big_m = 1.0;

                    if (logic_subsystems[l].second > 0) {
                        model.addConstr(rho[l][s][0u][0u] <= rho[l][s - 1u][h_last][num_x[jm][h_last] - 1u]);
                        model.addConstr(rho[l][s][0u][0u] <= rho[l][s - 1u][h_last][num_x[jm][h_last] - 1u] *
                                                             std::pow(1.0 - inst.component_reliability_at_subsystem[j][0u], num_k[j][0u][0u])
                                                             + big_m * (1 - x[j][0u][0u]));
                    }
                    else {
                        model.addConstr(rho[l][s][0u][0u] >= rho[l][s - 1u][h_last][num_x[jm][h_last] - 1u] * std::pow(
                                1.0 - inst.component_reliability_at_subsystem[j][0u], num_k[j][0u][0u]));
                        model.addConstr(rho[l][s][0u][0u] >= rho[l][s - 1u][h_last][num_x[jm][h_last] - 1u]
                                                             - big_m * x[j][0u][0u]);
                    }
                }
            }


            // at least one component is used
            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                GRBLinExpr expr = 0;

                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    for (auto k = 0u; k < num_x[j][h]; ++k)
                    {
                        expr += x[j][h][k];
                    }
                }

                model.addConstr(expr >= 1);
            }

            // resource constraints
            add_mixed_log_x_constraints(inst, model, x, num_x, num_k);


            /* Objective functions */
            {
                GRBLinExpr expr = inst.system_reliability.const_term;
                for (auto l = 0u; l < num_logic_subs; ++l)
                {
                    const auto last_sub = num_subs_in_logic[l] - 1u;
                    const auto last_com = inst.num_component_types - 1u;
                    const auto j = subsystems[l][last_sub];

                    expr += logic_subsystems[l].second * rho[l][last_sub][last_com][num_x[j][last_com] - 1u];
                }

                model.setObjective(expr, GRB_MAXIMIZE);
            }

            // solve

            if (has_time_limit) {
                model.set(GRB_DoubleParam_TimeLimit, time_limit);
            }

            auto start_time = std::chrono::high_resolution_clock::now();

            model.optimize();

            auto end_time = std::chrono::high_resolution_clock::now();
            auto cpu_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            const auto status = model.get(GRB_IntAttr_Status);

            if (status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
            {
                std::cerr << "\t-- Model infeasible or unbounded!\n";
                return std::nullopt;
            }

            if (status != GRB_OPTIMAL)
            {
                std::cerr << "\t-- Warning: Solution not optimal!\n";
            }

            std::cout << "\t-- OBJ: " << model.get(GRB_DoubleAttr_ObjVal) << std::endl;

            // extract solution

            auto s = ProbSolution{inst};
            s.obj = model.get(GRB_DoubleAttr_ObjVal);
            s.cpu_time = cpu_time;

            if (status == GRB_TIME_LIMIT) {
                s.time_limit = 1;
                s.mip_gap = model.get(GRB_DoubleAttr_MIPGap);
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }
            else {
                s.time_limit = 0;
                s.mip_gap = 0.0;
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }

            s.num_components_at_subsystem.resize(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                s.num_components_at_subsystem[j].resize(inst.num_component_types);

                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    s.num_components_at_subsystem[j][h] = 0;

                    for (auto k = 0u; k < num_x[j][h]; ++k)
                    {
                        if (x[j][h][k].get(GRB_DoubleAttr_X) > 0.5)
                        {
                            s.num_components_at_subsystem[j][h] += num_k[j][h][k];
                        }
                    }
                }
            }

            s.system_reliability = compute_system_reliability(s.num_components_at_subsystem);
            std::cout << "\t-- OBJ: " << s.system_reliability << std::endl;

            s.num_bin_vars = model.get(GRB_IntAttr_NumBinVars);
            s.num_vars = model.get(GRB_IntAttr_NumConstrs);
            s.num_cons = model.get(GRB_IntAttr_NumVars);
            s.num_nzs = model.get(GRB_IntAttr_NumNZs);
            s.num_bb_nodes = model.get(GRB_DoubleAttr_NodeCount);

            return s;
        }
        catch (GRBException &e)
        {
            std::cerr << "GUROBI ERROR:" << "\n";
            std::cerr << "\t" << e.getErrorCode() << "\n";
            std::cerr << "\t" << e.getMessage() << "\n";
        }
        catch (...)
        {
            std::cerr << "OTHER_ERRORS" << "\n";
        }

        return std::nullopt;
    }

    template<typename T>
    std::optional<ProbSolution> ProbSolver::solve_Q_tree_model(const std::vector<T>&  systems) const
    {
        auto subsystems = std::vector<size_t>(systems.size());
        for (auto l = 1u; l < systems.size(); ++l) {
            subsystems[l] = systems[l]->m_subsystem_id.value();
        }

        auto fathers = std::vector<size_t>(systems.size());
        for (auto l = 1u; l < systems.size(); ++l) {
            fathers[l] = systems[l]->m_father->m_id;
        }

        const bool tree_model_use_rho_ub = true;
        const bool tree_model_use_valid_inequalities = true;
        const bool tree_model_use_aggregated_expression = false;
        const bool tree_model_use_dominated_components = true;
        const bool tree_model_use_ub_one = false;


        /* calculate the upper bounds */

        auto rho_ub = std::vector<std::vector<double>>(systems.size());
        for (auto l = 0u; l < systems.size(); ++l) {
            rho_ub[l] = std::vector<double>(inst.num_component_types, 1.0);
        }
        if (tree_model_use_rho_ub) {
            for (auto l = 0u; l < systems.size(); ++l) {
                if (systems[l]->m_father) {
                    const auto f = fathers[l];
                    const auto j = subsystems[l];
                    const auto h_last = inst.num_component_types - 1u;

                    const auto min = *std::min_element(inst.component_reliability_at_subsystem[j].begin(),
                                                inst.component_reliability_at_subsystem[j].end());

                    rho_ub[l] = std::vector<double>(inst.num_component_types, rho_ub[f][h_last]);
                    rho_ub[l][h_last] *= 1 - min;
                }
            }
        }


        // model and solve
        try {
            GRBEnv env{};
            GRBModel model{env};

            /* variables */

            auto rho = std::vector<std::vector<GRBVar>>(systems.size());
            for (auto l = 0u; l < systems.size(); ++l) {
                rho[l] = std::vector<GRBVar>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    rho[l][h] = model.addVar(0.0, rho_ub[l][h], 0.0, GRB_CONTINUOUS);
                }
            }

            // x_{jh}^{k} = 1 iif k type-h components are configured at the j-th subsystem
            auto x = std::vector<std::vector<std::vector<GRBVar>>>(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                x[j] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    x[j][h] = std::vector<GRBVar>(inst.component_ub_at_subsystem[j][h] + 1u);
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        x[j][h][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    }
                }
            }

            // fixing variable according to domination between different component types
            if (tree_model_use_dominated_components) {
                for(auto [j, h] : inst.dominated_components) {
                    x[j][h][0].set(GRB_DoubleAttr_LB, 1.0);
                    for (auto k = 1u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        x[j][h][k].set(GRB_DoubleAttr_UB, 0.0);
                    }
                }
            }

            /* constraints */

            // root node [l = 0u]
            for (auto h = 0u; h < inst.num_component_types; ++h) {
                rho[0u][h].set(GRB_DoubleAttr_LB, 1.0);
            }

            // other nodes, first component type [l\ne 0u, h = 0u]
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto f = fathers[l];
                const auto j = subsystems[l];
                const auto h_last = inst.num_component_types - 1u;

                if (tree_model_use_aggregated_expression) {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                        GRBLinExpr expr = 0;
                        for (auto k2 = k; k2 <= inst.component_ub_at_subsystem[j][0u]; ++k2) {
                            expr += x[j][0u][k2];
                        }

                        const auto big_m = rho_ub[f][h_last] * (1.0 - compute_jhk(j, 0u, k));
                        model.addConstr(rho[l][0u] <= rho[f][h_last] * compute_jhk(j, 0u, k) + big_m * (1 - expr));
                    }
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                        GRBLinExpr expr = 0;
                        for (auto k2 = 0; k2 <= k; ++k2) {
                            expr += x[j][0u][k2];
                        }

                        const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][0u]);
                        const auto big_m = rho_ub[f][h_last] * (compute_jhk(j, 0u, k) - compute_jhk(j, 0u,K));
                        model.addConstr(rho[l][0u] >= rho[f][h_last] * compute_jhk(j, 0u, k) - big_m * (1 - expr));
                    }
                } else {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                        {
                            const auto big_m = rho_ub[f][h_last] * (1.0 - compute_jhk(j, 0u, k));
                            model.addConstr(rho[l][0u] <= rho[f][h_last] * compute_jhk(j, 0u, k) + big_m * (1 - x[j][0u][k]));
                        }
                        {
                            const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][0u]);
                            const auto big_m = rho_ub[f][h_last] * (compute_jhk(j, 0u, k) - compute_jhk(j, 0u, K));
                            model.addConstr(rho[l][0u] >= rho[f][h_last] * compute_jhk(j, 0u, k) - big_m * (1 - x[j][0u][k]));
                        }
                    }
                }
            }

            // the propagation within each node
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = subsystems[l];
                for (auto h = 1u; h < inst.num_component_types; ++h) {
                    const auto f = fathers[l];
                    const auto h_last = inst.num_component_types - 1u;

                    if (tree_model_use_aggregated_expression) {
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                            GRBLinExpr expr = 0;
                            for (auto k2 = k; k2 <= inst.component_ub_at_subsystem[j][h]; ++k2) {
                                expr += x[j][h][k2];
                            }

                            const auto big_m = rho_ub[f][h_last] * (1.0 - compute_jhk(j, h, k));
                            model.addConstr(rho[l][h] <= rho[l][h - 1u] * compute_jhk(j, h, k) + big_m * (1 - expr));
                        }
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                            GRBLinExpr expr = 0;
                            for (auto k2 = 0; k2 <= k; ++k2) {
                                expr += x[j][h][k2];
                            }

                            const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][h]);
                            const auto big_m = rho_ub[f][h_last] * (compute_jhk(j, h, k) - compute_jhk(j, h, K));
                            model.addConstr(rho[l][h] >= rho[l][h - 1u] * compute_jhk(j, h, k) - big_m * (1 - expr));
                        }
                    } else {
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                            // model.addGenConstrIndicator(x[j][h][k], 1, rho[l][h] == rho[l][h - 1u] * std::pow(
                            //         1.0 - inst.component_reliability_at_subsystem[j][h], static_cast<double>(k)));
                            {
                                const auto big_m = rho_ub[f][h_last] * (1.0 - compute_jhk(j, h, k));
                                model.addConstr(rho[l][h] <= rho[l][h - 1u] * compute_jhk(j, h, k) + big_m * (1 - x[j][h][k]));
                            }
                            {
                                const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][h]);
                                const auto big_m = rho_ub[f][h_last] * (compute_jhk(j, h, k) - compute_jhk(j, h, K));
                                model.addConstr(rho[l][h] >= rho[l][h - 1u] * compute_jhk(j, h, k) - big_m * (1 - x[j][h][k]));
                            }
                        }
                    }
                }
            }


            // valid inequalities
            if (tree_model_use_valid_inequalities) {
                for (auto l = 1u; l < systems.size(); ++l) {
                    const auto j = subsystems[l];
                    for (auto h = 0u; h < inst.num_component_types; ++h) {
                        GRBLinExpr expr = 0;
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                            expr += compute_jhk(j, h, k) * x[j][h][k];
                        }
                        const auto f = systems[l]->m_father->m_id;
                        const auto h_last = inst.num_component_types - 1u;
                        model.addConstr(rho[l][h] <= rho_ub[f][h_last] * expr);
                    }
                }

                /* relationships between different logic subsystems */
                for (auto l1 = 1u; l1 < systems.size(); ++l1) {
                    for (auto l2 = 1u; l2 < systems.size(); ++l2) {
                        if (l1 == l2) continue;

                        // we check starting from leaf nodes
                        if (!systems[l2]->m_children.empty()) continue;

                        // note: the nodes are stored in such a vector that the index of a father node
                        // is always less that of a child node
                        bool on_same_route = false;
                        {
                            size_t lmin = std::min(l1, l2);
                            size_t lmax = std::max(l1, l2);
                            auto f = systems[lmax]->m_father;
                            while (f) {
                                if (f->m_id == lmin) {
                                    on_same_route = true;
                                    break;
                                }
                                f = f->m_father;
                            }
                        }
                        if (on_same_route) continue;

                        std::vector<size_t> aa{l1};
                        std::vector<size_t> bb{l2};

                        auto f1 = systems[l1]->m_father;
                        while (f1 && f1->m_father) {
                            aa.push_back(f1->m_id);
                            f1 = f1->m_father;
                        }

                        auto f2 = systems[l2]->m_father;
                        while (f2 && f2->m_father) {
                            bb.push_back(f2->m_id);
                            f2 = f2->m_father;
                        }

                        std::reverse(bb.begin(), bb.end());

                        std::vector<size_t> a;
                        std::vector<size_t> b;
                        for (auto l: aa) {
                            a.push_back(systems[l]->m_subsystem_id.value());
                        }
                        for (auto l: bb) {
                            b.push_back(systems[l]->m_subsystem_id.value());
                        }

                        if (a.empty() || b.empty()) continue;

                        if (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                            auto last_id = bb.back();
                            const auto h_last = inst.num_component_types - 1u;
                            model.addConstr(rho[last_id][h_last] <= rho[l1][h_last]);
                        }
                    }
                }
            }


            // exactly a particular number of components of type h is used
            // - \sum_{k=0}^{u_jh} = 1, \forall j\in J, h\in H_j.
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    GRBLinExpr expr = 0;

                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += x[j][h][k];
                    }

                    model.addConstr(expr == 1);
                }
            }

            // at least one component is used
            // \sum_{h\in H_j} \sum_{k\in K_jh} x_{jh}^{k} >= 1, \forall j\in J.
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                GRBLinExpr expr = 0;

                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 1u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += x[j][h][k];
                    }
                }

                model.addConstr(expr >= 1);
            }

            // resource constraints
            if (inst.test_type == TEST_TYPE::MIXED) {
                add_mixed_constraints(inst, model, x);
            }
            else if (inst.test_type == TEST_TYPE::TEST_ONE) {
                add_test_one_constraints(inst, model, x);
            }
            else if (inst.test_type == TEST_TYPE::TEST_TWO) {
                add_test_two_constraints(inst, model, x);
            }


            /* Objective functions */

            {
                GRBLinExpr expr = 0;
                for (auto l = 0u; l < systems.size(); ++l) {
                    if (systems[l]->m_logic) {
                        const auto h_last = inst.num_component_types - 1u;
                        expr += systems[l]->m_coeff * rho[l][h_last];
                    }
                }

                if (tree_model_use_ub_one) {
                    model.addConstr(expr <= 1.0);
                }

                model.setObjective(expr, GRB_MAXIMIZE);
            }

            // Solve:

            if (has_time_limit) {
                model.set(GRB_DoubleParam_TimeLimit, time_limit);
            }


            auto start_time = std::chrono::high_resolution_clock::now();

            model.optimize();

            auto end_time = std::chrono::high_resolution_clock::now();
            auto cpu_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();


            const auto status = model.get(GRB_IntAttr_Status);

            if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
            {
                std::cerr << "\t-- Model infeasible or unbounded!\n";
                return std::nullopt;
            }

            if (status != GRB_OPTIMAL)
            {
                std::cerr << "\t-- Warning: Solution not optimal!\n";
            }

            std::cout << "\t-- OBJ: " << model.get(GRB_DoubleAttr_ObjVal) << "\n";

            // extract solution

            auto s = ProbSolution{inst};
            s.obj = model.get(GRB_DoubleAttr_ObjVal);
            s.cpu_time = cpu_time;

            if (status == GRB_TIME_LIMIT) {
                s.time_limit = 1;
                s.mip_gap = model.get(GRB_DoubleAttr_MIPGap);
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }
            else {
                s.time_limit = 0;
                s.mip_gap = 0.0;
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }

            s.num_components_at_subsystem.resize(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                s.num_components_at_subsystem[j].resize(inst.num_component_types);

                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    for (auto k = 0u; k < inst.component_ub_at_subsystem[j][h]; ++k)
                    {
                        if (x[j][h][k].get(GRB_DoubleAttr_X) > 0.5)
                        {
                            s.num_components_at_subsystem[j][h] = k;
                        }
                    }
                }
            }

            s.system_reliability = compute_system_reliability(s.num_components_at_subsystem);
            std::cout << "\t-- OBJ: " << s.system_reliability << "\n";

            s.num_bin_vars = model.get(GRB_IntAttr_NumBinVars);
            s.num_vars = model.get(GRB_IntAttr_NumConstrs);
            s.num_cons = model.get(GRB_IntAttr_NumVars);
            s.num_nzs = model.get(GRB_IntAttr_NumNZs);
            s.num_bb_nodes = model.get(GRB_DoubleAttr_NodeCount);

            return s;
        }
        catch (GRBException &e)
        {
            std::cerr << "GUROBI ERROR:\n";
            std::cerr << "\t" << e.getErrorCode() << "\n";
            std::cerr << "\t" << e.getMessage() << "\n";
        }
        catch (...)
        {
            std::cerr << "OTHER_ERRORS\n";
        }

        return std::nullopt;
    }

    template<typename T>
    std::optional<ProbSolution> ProbSolver::solve_Q_tree_logx_model(const std::vector<T>&  systems) const
    {
        const bool sorted_tree_logx_model_use_rho_ub = true;
        const bool sorted_tree_logx_model_use_valid_inequalities = true;
        const bool sorted_tree_logx_model_use_dominated_components = true;
        const bool sorted_tree_logx_model_use_ub_one = false;

        auto num_x = compute_num_x(inst);
        auto num_k = compute_num_k(inst, num_x);

        /* calculate the upper bounds */
        auto rho_ub = std::vector<std::vector<double>>(systems.size());
        for (auto l = 0u; l < systems.size(); ++l) {
            rho_ub[l] = std::vector<double>(inst.num_component_types, 1.0);
        }
        if (sorted_tree_logx_model_use_rho_ub) {
            for (auto l = 0u; l < systems.size(); ++l) {
                if (systems[l]->m_father) {
                    const auto f = systems[l]->m_father->m_id;
                    const auto j = systems[l]->m_subsystem_id.value();
                    const auto h_last = inst.num_component_types - 1u;

                    auto min = *std::min_element(inst.component_reliability_at_subsystem[j].begin(),
                                                 inst.component_reliability_at_subsystem[j].end());

                    rho_ub[l] = std::vector<double>(inst.num_component_types, rho_ub[f][h_last]);
                    rho_ub[l][h_last] *= 1 - min;
                }
            }
        }


        // model and solve
        try
        {
            GRBEnv env{};
            GRBModel model{env};

            /* variables */

            auto rho = std::vector<std::vector<std::vector<GRBVar>>>(systems.size());
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = systems[l]->m_subsystem_id.value();
                rho[l] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    rho[l][h] = std::vector<GRBVar>(num_x[j][h]);
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        rho[l][h][k] = model.addVar(0.0, rho_ub[l][h], 0.0, GRB_CONTINUOUS);
                    }
                }
            }

            auto x = std::vector<std::vector<std::vector<GRBVar>>>(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                x[j] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    x[j][h] = std::vector<GRBVar>(num_x[j][h]);
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        x[j][h][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    }
                }
            }

            /* fixing variable according to domination between different component types */
            if (sorted_tree_logx_model_use_dominated_components) {
                for (auto [j, h]: inst.dominated_components) {
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        x[j][h][k].set(GRB_DoubleAttr_UB, 0.0);
                    }
                }
            }


            // Constraints:

            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = systems[l]->m_subsystem_id.value();
                const auto f = systems[l]->m_father->m_id;

                if (f == 0u) {//root node
                    model.addConstr(rho[l][0u][0u] == 1 - (1 - compute_jhk(j, 0u, num_k[j][0u][0u])) * x[j][0u][0u]);
                } else {
                    const auto f_j = systems[l]->m_father->m_subsystem_id.value();
                    const auto f_h_last = inst.num_component_types - 1u;
                    const auto f_k_last = num_x[f_j][f_h_last] - 1u;

                    const auto big_m = rho_ub[f][f_h_last];
                    model.addConstr(rho[l][0u][0u] <= rho[f][f_h_last][f_k_last]);
                    model.addConstr(rho[l][0u][0u] <= rho[f][f_h_last][f_k_last] * compute_jhk(j, 0u, num_k[j][0u][0u]) + big_m * (1 - x[j][0u][0u]));
                    model.addConstr(rho[l][0u][0u] >= rho[f][f_h_last][f_k_last] * compute_jhk(j, 0u, num_k[j][0u][0u]));
                    model.addConstr(rho[l][0u][0u] >= rho[f][f_h_last][f_k_last] - big_m * x[j][0u][0u]);
                }
            }


            // other nodes [l\ne 0u, h = 0u]
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = systems[l]->m_subsystem_id.value();
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 1u; k < num_x[j][h]; ++k) {
                        const auto big_m = rho_ub[l][h];
                        model.addConstr(rho[l][h][k] <= rho[l][h][k - 1u]);
                        model.addConstr(rho[l][h][k] <= rho[l][h][k - 1u] * compute_jhk(j, h, num_k[j][h][k]) + big_m * (1 - x[j][h][k]));
                        model.addConstr(rho[l][h][k] >= rho[l][h][k - 1u] * compute_jhk(j, h, num_k[j][h][k]));
                        model.addConstr(rho[l][h][k] >= rho[l][h][k - 1u] - big_m * x[j][h][k]);
                    }
                }
            }

            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = systems[l]->m_subsystem_id.value();
                for (auto h = 1u; h < inst.num_component_types; ++h) {
                    const auto last_k = num_x[j][h - 1u] - 1u;
                    const auto big_m = rho_ub[l][h - 1u];
                    model.addConstr(rho[l][h][0u] <= rho[l][h - 1u][last_k]);
                    model.addConstr(rho[l][h][0u] <= rho[l][h - 1u][last_k] * compute_jhk(j, h, num_k[j][h][0u]) + big_m * (1 - x[j][h][0u]));
                    model.addConstr(rho[l][h][0u] >= rho[l][h - 1u][last_k] * compute_jhk(j, h, num_k[j][h][0u]));
                    model.addConstr(rho[l][h][0u] >= rho[l][h - 1u][last_k] - big_m * x[j][h][0u]);
                }
            }


            // exactly a particular number of components of type h is used
            // - \sum_{k=0}^{u_jh} x_{jh}^{k} >= 1, \forall j\in J, h\in H_j.
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                GRBLinExpr expr = 0;

                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        expr += x[j][h][k];
                    }
                }

                model.addConstr(expr >= 1);
            }

            // resource constraints
            add_mixed_log_x_constraints(inst, model, x, num_x, num_k);


            // valid inequalities
            if (sorted_tree_logx_model_use_valid_inequalities) {
                /* relationships between different logic subsystems */
                for (auto l1 = 1u; l1 < systems.size(); ++l1) {
                    for (auto l2 = 1u; l2 < systems.size(); ++l2) {
                        if (l1 == l2) continue;

                        // we check starting from leaf nodes
                        if (!systems[l2]->m_children.empty()) continue;

                        // note: the nodes are stored in such a vector that the index of a father node
                        // is always less that of a child node
                        bool on_same_route = false;
                        {
                            size_t lmin = std::min(l1, l2);
                            size_t lmax = std::max(l1, l2);
                            auto f = systems[lmax]->m_father;
                            while (f) {
                                if (f->m_id == lmin) {
                                    on_same_route = true;
                                    break;
                                }
                                f = f->m_father;
                            }
                        }
                        if (on_same_route) continue;

                        std::vector<size_t> aa{l1};
                        std::vector<size_t> bb{l2};

                        auto f1 = systems[l1]->m_father;
                        while (f1 && f1->m_father) {
                            aa.push_back(f1->m_id);
                            f1 = f1->m_father;
                        }

                        auto f2 = systems[l2]->m_father;
                        while (f2 && f2->m_father) {
                            bb.push_back(f2->m_id);
                            f2 = f2->m_father;
                        }

                        std::reverse(bb.begin(), bb.end());

                        std::vector<size_t> a;
                        std::vector<size_t> b;
                        for (auto l: aa) {
                            a.push_back(systems[l]->m_subsystem_id.value());
                        }
                        for (auto l: bb) {
                            b.push_back(systems[l]->m_subsystem_id.value());
                        }

                        if (a.empty() || b.empty()) continue;

                        if (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                            auto last_id = bb.back();
                            b.pop_back();
                            bb.pop_back();
                            while (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                                last_id = bb.back();
                                b.pop_back();
                                bb.pop_back();
                            }

                            const auto h_last = inst.num_component_types - 1u;
                            const auto last_j = systems[last_id]->m_subsystem_id.value();
                            const auto last_k_last = num_x[last_j][h_last] - 1u;
                            const auto l1_j = systems[l1]->m_subsystem_id.value();
                            const auto l1_k_last = num_x[l1_j][h_last] - 1u;
                            model.addConstr(rho[last_id][h_last][last_k_last] <= rho[l1][h_last][l1_k_last]);
                        }
                    }
                }
            }


            /* Objective functions */

            {
                GRBLinExpr expr = systems[0u]->m_coeff;

                for (auto l = 1u; l < systems.size(); ++l) {
                    if (systems[l]->m_logic) {
                        const auto j = systems[l]->m_subsystem_id.value();
                        const auto h_last = inst.num_component_types - 1u;
                        const auto k_last = num_x[j][h_last] - 1u;

                        expr += systems[l]->m_coeff * rho[l][h_last][k_last];
                    }
                }

                if (sorted_tree_logx_model_use_ub_one) {
                    model.addConstr(expr <= 1.0);
                }

                model.setObjective(expr, GRB_MAXIMIZE);
            }

            // solve

            if (has_time_limit) {
                model.set(GRB_DoubleParam_TimeLimit, time_limit);
            }

            auto start_time = std::chrono::high_resolution_clock::now();

            model.optimize();

            auto end_time = std::chrono::high_resolution_clock::now();
            auto cpu_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            const auto status = model.get(GRB_IntAttr_Status);

            if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
            {
                std::cerr << "\t-- Model infeasible or unbounded!\n";
                return std::nullopt;
            }

            if (status != GRB_OPTIMAL)
            {
                std::cerr << "\t-- Warning: Solution not optimal!\n";
            }

            std::cout << "\t-- OBJ: " << model.get(GRB_DoubleAttr_ObjVal) << "\n";

            // Extract solution:

            auto s = ProbSolution{inst};
            s.obj = model.get(GRB_DoubleAttr_ObjVal);
            s.cpu_time = cpu_time;

            if (status == GRB_TIME_LIMIT) {
                s.time_limit = 1;
                s.mip_gap = model.get(GRB_DoubleAttr_MIPGap);
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }
            else {
                s.time_limit = 0;
                s.mip_gap = 0.0;
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }

            s.num_components_at_subsystem.resize(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                s.num_components_at_subsystem[j].resize(inst.num_component_types, 0.0);

                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    for (auto k = 0u; k < num_x[j][h]; ++k)
                    {
                        if (x[j][h][k].get(GRB_DoubleAttr_X) > 0.5)
                        {
                            s.num_components_at_subsystem[j][h] += num_k[j][h][k];
                        }
                    }
                }
            }

            s.system_reliability = compute_system_reliability(s.num_components_at_subsystem);
            std::cout << "\t-- OBJ: " << s.system_reliability << "\n";

            s.num_bin_vars = model.get(GRB_IntAttr_NumBinVars);
            s.num_vars = model.get(GRB_IntAttr_NumConstrs);
            s.num_cons = model.get(GRB_IntAttr_NumVars);
            s.num_nzs = model.get(GRB_IntAttr_NumNZs);
            s.num_bb_nodes = model.get(GRB_DoubleAttr_NodeCount);

            return s;
        }
        catch (GRBException &e)
        {
            std::cerr << "GUROBI ERROR:\n";
            std::cerr << "\t" << e.getErrorCode() << "\n";
            std::cerr << "\t" << e.getMessage() << "\n";
        }
        catch (...)
        {
            std::cerr << "OTHER_ERRORS\n";
        }

        return std::nullopt;
    }

    /*  general-form */
    std::optional<ProbSolution> ProbSolver::solve_general_basic_model() const
    {
        /* the number of product terms */
        const auto num_logic_subs{inst.qr_system_reliability.m_other_terms.size()};

        /* save product terms in a vector for the following usage */
        std::vector<std::pair<std::set<std::pair<int, bool>>, int>> logic_subsystems;
        logic_subsystems.reserve(num_logic_subs);
        for (auto &[term, coeff]: inst.qr_system_reliability.m_other_terms) {
            logic_subsystems.emplace_back(term, coeff);
        }

        /* the subsystems in each logic subsystem */
        std::vector<std::vector<std::pair<int, bool>>> subsystems(num_logic_subs);
        for (auto l = 0; l < num_logic_subs; ++l) {
            std::copy(logic_subsystems[l].first.begin(), logic_subsystems[l].first.end(),
                      std::back_inserter(subsystems[l]));
        }

        // pre-process
        for (auto l = 0u; l < num_logic_subs; ++l) {
            auto x = inst.component_reliability_at_subsystem;
            std::sort(subsystems[l].begin(), subsystems[l].end(),
                      [&x](auto &a, auto &b) {
                          auto j1 = *std::min_element(x[a.first].begin(), x[a.first].end());
                          auto j2 = *std::min_element(x[b.first].begin(), x[b.first].end());
                          if (a.second && b.second)
                              return j1 > j2;
                          if (!a.second && !b.second)
                              return j1 < j2;
                          if (a.second)
                              return false;
                          return true;
                      });
        }

        /* number of subsystems in each logic subsystem */
        std::vector<size_t> num_subs_in_logic;
        num_subs_in_logic.reserve(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            num_subs_in_logic.push_back(subsystems[l].size());
        }

        /* find the component type that has the minimum reliability at each subsystem */
        std::vector<std::vector<double>> min_r(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            min_r[l] = std::vector<double>(num_subs_in_logic[l]);
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                auto j = subsystems[l][s].first;
                min_r[l][s] = *std::min_element(inst.component_reliability_at_subsystem[j].begin(),
                                                inst.component_reliability_at_subsystem[j].end());
            }
        }


        const bool general_model_use_varphi_ub = true;
        const bool general_model_use_rho_ub = true;
        const bool general_model_use_valid_inequalities = true;
        const bool general_model_use_aggregated_expression = false;
        const bool general_model_use_dominated_components = true;
        const bool general_model_use_ub_one = false;

        /* calculate the upper bounds of decision variables \varphi_{ls} */
        std::vector<std::vector<double>> varphi_ub_(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            varphi_ub_[l] = std::vector<double>(num_subs_in_logic[l], 1.0);
        }
        for (auto l = 0u; l < num_logic_subs; ++l) {
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                if (subsystems[l][s].second) { // R
                    varphi_ub_[l][s] = (s == 0 ? 1.0 : varphi_ub_[l][s - 1u]);
                } else { // Q
                    varphi_ub_[l][s] = (s == 0 ? 1.0 : varphi_ub_[l][s - 1u]) * (1.0 - min_r[l][s]);
                }
            }
        }

        /* calculate the upper bounds of decision variables \rho_{lsh} */
        std::vector<std::vector<std::vector<double>>> rho_ub(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            rho_ub[l] = std::vector<std::vector<double>>(num_subs_in_logic[l]);
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                rho_ub[l][s] = std::vector<double>(inst.num_component_types, 1.0);
            }
        }

        if (general_model_use_rho_ub) {
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    if (s != 0u) {
                        rho_ub[l][s] = std::vector<double>(inst.num_component_types, varphi_ub_[l][s - 1u]);
                    }

                    const auto h_last = inst.num_component_types - 1u;
                    rho_ub[l][s][h_last] *= 1 - min_r[l][s];
                }
            }
        }

        /* should we use the upper bounds of varphi */
        std::vector<std::vector<double>> varphi_ub(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            varphi_ub[l] = std::vector<double>(num_subs_in_logic[l], 1.0);
        }
        if (general_model_use_varphi_ub) {
            varphi_ub = std::move(varphi_ub_);
        }


        /* model and solve */

        try
        {
            GRBEnv env;
            GRBModel model{env};

            // variables rho, varphi, x

            auto rho = std::vector<std::vector<std::vector<GRBVar>>>(num_logic_subs);
            for (auto l = 0u; l < num_logic_subs; ++l) {
                rho[l] = std::vector<std::vector<GRBVar>>(num_subs_in_logic[l]);
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    rho[l][s] = std::vector<GRBVar>(inst.num_component_types);
                    for (auto h = 0u; h < inst.num_component_types; ++h) {
                        rho[l][s][h] = model.addVar(0.0, rho_ub[l][s][h], 0.0, GRB_CONTINUOUS);
                    }
                }
            }

            auto varphi = std::vector<std::vector<GRBVar>>(num_logic_subs);
            for (auto l = 0u; l < num_logic_subs; ++l) {
                varphi[l] = std::vector<GRBVar>(num_subs_in_logic[l]);
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    varphi[l][s] = model.addVar(0.0, varphi_ub[l][s], 0.0, GRB_CONTINUOUS);
                }
            }

            auto x = std::vector<std::vector<std::vector<GRBVar>>>(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                x[j] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    x[j][h] = std::vector<GRBVar>(inst.component_ub_at_subsystem[j][h] + 1u);
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        x[j][h][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    }
                }
            }

            if (general_model_use_dominated_components) {
                for(auto [j, h] : inst.dominated_components) {
                    x[j][h][0].set(GRB_DoubleAttr_LB, 1.0);
                    for (auto k = 1u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        x[j][h][k].set(GRB_DoubleAttr_UB, 0.0);
                    }
                }
            }

            // constraints

            /*
             * boundary condition
             * */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                const auto j = subsystems[l][0u].first;
                GRBLinExpr expr = 0;

                for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                    expr += compute_jhk(j,0u,k) * x[j][0u][k];
                }

                model.addConstr(rho[l][0u][0u] == expr);
            }

            /* recursive relation between different subsystems of each logic subsystem */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    const auto h_last = inst.num_component_types - 1u;

                    if (subsystems[l][s].second) { // R
                        if (s == 0u) {
                            model.addConstr(varphi[l][s] == 1.0 - rho[l][s][h_last]);
                        } else {
                            model.addConstr(varphi[l][s] == varphi[l][s - 1u] - rho[l][s][h_last]);
                        }
                    } else { // Q
                        model.addConstr(varphi[l][s] == rho[l][s][h_last]);
                    }
                }
            }

            /* recursive relation between different subsystems of each logic subsystem */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 1u; s < num_subs_in_logic[l]; ++s) {
                    auto j = subsystems[l][s].first; // j = \delta(s)
                    if (general_model_use_aggregated_expression) {
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                            {
                                GRBLinExpr expr = 0;
                                for (auto k2 = k; k2 <= inst.component_ub_at_subsystem[j][0u]; ++k2) {
                                    expr += x[j][0u][k2];
                                }

                                const auto big_m = varphi_ub[l][s - 1u] * (1 - compute_jhk(j, 0u, k));
                                model.addConstr(rho[l][s][0u] <= varphi[l][s - 1u] * compute_jhk(j, 0u, k) + big_m * (1 - expr));
                            }
                            {
                                GRBLinExpr expr = 0;
                                for (auto k2 = 0u; k2 <= k; ++k2) {
                                    expr += x[j][0u][k2];
                                }

                                const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][0u]);
                                const auto big_m =  varphi_ub[l][s - 1u] * (compute_jhk(j,0u,k) - compute_jhk(j,0u,K));

                                model.addConstr(rho[l][s][0u] >= varphi[l][s - 1u] * compute_jhk(j,0u,k) - big_m * (1 - expr));
                            }
                        }
                    }
                    else {
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                            // model.addGenConstrIndicator(x[j][0u][k], 1, rho[l][s][0u] == varphi[l][s - 1u] * std::pow(1.0 - inst.component_reliability_at_subsystem[j][0u], static_cast<double>(k)));
                            {
                                const auto big_m = varphi_ub[l][s - 1u] * (1 - compute_jhk(j,0u,k));
                                model.addConstr(rho[l][s][0u] <= varphi[l][s - 1u] * compute_jhk(j,0u,k) + big_m * (1 - x[j][0u][k]));
                            }
                            {
                                const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][0u]);
                                const auto big_m = varphi_ub[l][s - 1u] * (compute_jhk(j,0u,k) - compute_jhk(j,0u,K));
                                model.addConstr(rho[l][s][0u] >= varphi[l][s - 1u] * compute_jhk(j,0u,k) - big_m * (1 - x[j][0u][k]));
                            }
                        }
                    }
                }
            }

            // recursive relation within the same subsystem
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    auto j = subsystems[l][s].first; // j = \delta(s)
                    for (auto h = 1u; h < inst.num_component_types; ++h) {
                        if (general_model_use_aggregated_expression) {
                            for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                                GRBLinExpr expr = 0;
                                for (auto k2 = k; k2 <= inst.component_ub_at_subsystem[j][h]; ++k2) {
                                    expr += x[j][h][k2];
                                }

                                const auto big_m =  rho_ub[l][s][h - 1u] * (1 - compute_jhk(j,h,k));
                                model.addConstr(rho[l][s][h] <= rho[l][s][h - 1u] * compute_jhk(j,h,k) + big_m * (1 - expr));
                            }

                            for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                                GRBLinExpr expr = 0;
                                for (auto k2 = 0u; k2 <= k; ++k2) {
                                    expr += x[j][h][k2];
                                }

                                const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][h]);
                                const auto big_m =  rho_ub[l][s][h - 1u] * (compute_jhk(j,h,k) - compute_jhk(j,h,K));

                                model.addConstr(rho[l][s][h] >= rho[l][s][h - 1u] * compute_jhk(j,h,k) - big_m * (1 - expr));
                            }
                        } else {
                            for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                                // indicator version
                                // model.addGenConstrIndicator(x[j][h][k], 1, rho[l][s][h] == rho[l][s][h - 1u] * std::pow(1.0 - inst.component_reliability_at_subsystem[j][h], static_cast<double>(k)));
                                {
                                    const auto big_m = rho_ub[l][s][h - 1u] * (1 - compute_jhk(j,h,k));
                                    model.addConstr(rho[l][s][h] <= rho[l][s][h - 1u] * compute_jhk(j,h,k) + big_m * (1 - x[j][h][k]));
                                }
                                {
                                    const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][h]);
                                    const auto big_m = rho_ub[l][s][h - 1u] * (compute_jhk(j,h,k) - compute_jhk(j,h,K));
                                    model.addConstr(rho[l][s][h] >= rho[l][s][h - 1u] * compute_jhk(j,h,k) - big_m * (1 - x[j][h][k]));
                                }
                            }
                        }
                    }
                }
            }

            // sum_{k=0}^{u_jh} = 1, \forall j, h: exactly a particular number of components of type h is used
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    GRBLinExpr expr = 0;

                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += x[j][h][k];
                    }

                    model.addConstr(expr == 1);
                }
            }

            // at least one component is used
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                GRBLinExpr expr = 0;

                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 1u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += x[j][h][k];
                    }
                }

                model.addConstr(expr >= 1);
            }

            // resource constraints
            if (inst.test_type == TEST_TYPE::MIXED) {
                add_mixed_constraints(inst, model, x);
            }
            else if (inst.test_type == TEST_TYPE::TEST_ONE) {
                add_test_one_constraints(inst, model, x);
            }
            else if (inst.test_type == TEST_TYPE::TEST_TWO) {
                add_test_two_constraints(inst, model, x);
            }


            /* valid inequalities */
            if (general_model_use_valid_inequalities) {
                // upper bound on the reliability of each subsystem after the components are configured
                for (auto l = 0u; l < num_logic_subs; ++l) {
                    for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                        const auto j = subsystems[l][s].first;
                        for (auto h = 1u; h < inst.num_component_types; ++h) {
                            GRBLinExpr expr = 0;
                            for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                                expr += compute_jhk(j, h, k) * x[j][h][k];
                            }

                            model.addConstr(rho[l][s][h] <= rho_ub[l][s][0u] * expr);
                        }
                    }
                }

                for (auto l1 = 0u; l1 < num_logic_subs; ++l1) {
                    std::set<std::pair<int, bool>> a;
                    for (auto s1 = 0u; s1 < num_subs_in_logic[l1]; ++s1) {
                        a.insert(subsystems[l1][s1]);
                        for (auto l2 = 0u; l2 < num_logic_subs; ++l2) {
                            if (l1 == l2) continue;

                            std::set<std::pair<int, bool>> b;
                            for (auto s2 = 0u; s2 < num_subs_in_logic[l2]; ++s2) {
                                b.insert(subsystems[l2][s2]);

                                // check if 'a' is a subset of 'b'
                                if (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                                    model.addConstr(varphi[l2][s2] <= varphi[l1][s1]);
                                    break;
                                }
                            }
                        }
                    }
                }
            }


            /* Objective functions */

            {
                GRBLinExpr expr = inst.qr_system_reliability.const_term;
                // std::cout << "system_reliability.const_term: " << inst.system_reliability.const_term << "\n";

                for (auto l = 0u; l < num_logic_subs; ++l) {
                    const auto last_sub = num_subs_in_logic[l] - 1u;

                    expr += logic_subsystems[l].second * varphi[l][last_sub];
                }

                if (general_model_use_ub_one) {
                    model.addConstr(expr <= 1.0);
                }

                model.setObjective(expr, GRB_MAXIMIZE);
            }

            // Solve:

            if (has_time_limit) {
                model.set(GRB_DoubleParam_TimeLimit, time_limit);
            }

            auto start_time = std::chrono::high_resolution_clock::now();

            model.optimize();

            auto end_time = std::chrono::high_resolution_clock::now();
            auto cpu_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            const auto status = model.get(GRB_IntAttr_Status);

            if (status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
            {
                std::cout << "\t-- Model infeasible or unbounded!\n";
                return std::nullopt;
            }

            if (status != GRB_OPTIMAL)
            {
                std::cout << "\t-- Warning: Solution not optimal!\n";
            }

            std::cout << "\t-- OBJ: " << model.get(GRB_DoubleAttr_ObjVal) << "\n";

            // extract solution

            auto s = ProbSolution{inst};
            s.obj = model.get(GRB_DoubleAttr_ObjVal);
            s.cpu_time = cpu_time;

            if (status == GRB_TIME_LIMIT) {
                s.time_limit = 1;
                s.mip_gap = model.get(GRB_DoubleAttr_MIPGap);
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }
            else {
                s.time_limit = 0;
                s.mip_gap = 0.0;
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }

            s.num_components_at_subsystem.resize(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                s.num_components_at_subsystem[j].resize(inst.num_component_types);

                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    for (auto k = 0u; k < inst.component_ub_at_subsystem[j][h]; ++k)
                    {
                        if (x[j][h][k].get(GRB_DoubleAttr_X) > 0.5)
                        {
                            s.num_components_at_subsystem[j][h] = k;
                        }
                    }
                }
            }

            s.system_reliability = compute_system_reliability(s.num_components_at_subsystem);
            std::cout << "\t-- OBJ: " << s.system_reliability << "\n";

            s.num_bin_vars = model.get(GRB_IntAttr_NumBinVars);
            s.num_vars = model.get(GRB_IntAttr_NumConstrs);
            s.num_cons = model.get(GRB_IntAttr_NumVars);
            s.num_nzs = model.get(GRB_IntAttr_NumNZs);
            s.num_bb_nodes = model.get(GRB_DoubleAttr_NodeCount);

            return s;
        }
        catch (GRBException &e)
        {
            std::cerr << "GUROBI ERROR:\n" << "\t" << e.getErrorCode() << "\n\t" << e.getMessage() << "\n";
        }
        catch (...)
        {
            std::cerr << "OTHER_ERRORS" << "\n";
        }

        return std::nullopt;
    }

    std::optional<ProbSolution> ProbSolver::solve_general_basic_logx_model() const
    {
        /* the number of product terms */
        const auto num_logic_subs { inst.qr_system_reliability.m_other_terms.size() };

        /* save product terms in a vector for the following usage */
        std::vector<std::pair<std::set<std::pair<int, bool>>, int>> logic_subsystems;
        logic_subsystems.reserve(num_logic_subs);
        for (auto& [term, coeff] : inst.qr_system_reliability.m_other_terms){
            logic_subsystems.emplace_back(term, coeff);
        }

        /* the subsystems in each logic subsystem */
        std::vector<std::vector<std::pair<size_t, bool>>> subsystems(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            std::copy(logic_subsystems[l].first.begin(), logic_subsystems[l].first.end(),
                      std::back_inserter(subsystems[l]));
        }

        auto num_x = compute_num_x(inst);
        auto num_k = compute_num_k(inst, num_x);

        /* number of subsystems in each logic subsystem */
        std::vector<size_t> num_subs_in_logic;
        num_subs_in_logic.reserve(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            num_subs_in_logic.push_back(subsystems[l].size());
        }

        /* find the component type that has the minimum reliability at each subsystem */
        std::vector<std::vector<double>> min_r(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            min_r[l] = std::vector<double>(num_subs_in_logic[l]);
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                auto j = subsystems[l][s].first;
                min_r[l][s] = *std::min_element(inst.component_reliability_at_subsystem[j].begin(),
                                                inst.component_reliability_at_subsystem[j].end());
            }
        }

        const bool general_logx_model_use_varphi_ub = true;
        const bool general_logx_model_use_rho_ub = true;
        const bool general_logx_model_use_valid_inequalities = false;
        const bool general_logx_model_use_dominated_components = true;
        const bool general_logx_model_use_ub_one = false;

        /* calculate the upper bounds of decision variables \varphi_{ls} */
        std::vector<std::vector<double>> varphi_ub_(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            varphi_ub_[l] = std::vector<double>(num_subs_in_logic[l], 1.0);
        }
        for (auto l = 0u; l < num_logic_subs; ++l) {
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                if (subsystems[l][s].second) { // R
                    varphi_ub_[l][s] = (s == 0 ? 1.0 : varphi_ub_[l][s - 1u]);
                } else { // Q
                    varphi_ub_[l][s] = (s == 0 ? 1.0 : varphi_ub_[l][s - 1u]) * (1.0 - min_r[l][s]);
                }
            }
        }

        /* calculate the upper bounds of decision variables \rho_{ls*} */
        std::vector<std::vector<std::vector<double>>> rho_ub(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            rho_ub[l] = std::vector<std::vector<double>>(num_subs_in_logic[l]);
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                rho_ub[l][s] = std::vector<double>(inst.num_component_types, 1.0);
            }
        }
        if (general_logx_model_use_rho_ub) {
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    if (s != 0u) {
                        rho_ub[l][s] = std::vector<double>(inst.num_component_types, varphi_ub_[l][s - 1u]);
                    }

                    const auto h_last = inst.num_component_types - 1u;
                    rho_ub[l][s][h_last] *= 1 - min_r[l][s];
                }
            }
        }

        /* should we use the upper bounds of varphi */
        std::vector<std::vector<double>> varphi_ub(num_logic_subs);
        for (auto l = 0u; l < num_logic_subs; ++l) {
            varphi_ub[l] = std::vector<double>(num_subs_in_logic[l], 1.0);
        }
        if (general_logx_model_use_varphi_ub) {
            varphi_ub = std::move(varphi_ub_);
        }

        /* model and solve */

        try
        {
            GRBEnv env{};
            GRBModel model{env};

            // variables rho, varphi, x

            auto rho = std::vector<std::vector<std::vector<std::vector<GRBVar>>>>(num_logic_subs);
            for (auto l = 0u; l < num_logic_subs; ++l) {
                rho[l] = std::vector<std::vector<std::vector<GRBVar>>>(num_subs_in_logic[l]);
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    const auto j = subsystems[l][s].first; // j = \delta(s)
                    rho[l][s] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                    for (auto h = 0u; h < inst.num_component_types; ++h) {
                        rho[l][s][h] = std::vector<GRBVar>(num_x[j][h]);
                        for (auto k = 0u; k < num_x[j][h]; ++k) {
                            rho[l][s][h][k] = model.addVar(0.0, rho_ub[l][s][h], 0.0, GRB_CONTINUOUS);
                        }
                    }
                }
            }

            auto varphi = std::vector<std::vector<GRBVar>>(num_logic_subs);
            for (auto l = 0u; l < num_logic_subs; ++l) {
                varphi[l] = std::vector<GRBVar>(num_subs_in_logic[l]);
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    varphi[l][s] = model.addVar(0.0, varphi_ub[l][s], 0.0, GRB_CONTINUOUS);
                }
            }

            auto x = std::vector<std::vector<std::vector<GRBVar>>>(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                x[j] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    x[j][h] = std::vector<GRBVar>(num_x[j][h]);
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        x[j][h][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    }
                }
            }

            /* fixing variable according to domination between different component types */
            if (general_logx_model_use_dominated_components) {
                for (auto [j, h]: inst.dominated_components) {
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        x[j][h][k].set(GRB_DoubleAttr_UB, 0.0);
                    }
                }
            }


            // constraints

            /*
             * boundary condition
             * */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                const auto j = subsystems[l][0u].first;
                model.addConstr(rho[l][0u][0u][0u] == 1 - (1 - compute_jhk(j,0u,num_k[j][0u][0u])) * x[j][0u][0u]);
            }

            /* recursive relation between different subsystems of each logic subsystem */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                const auto j = subsystems[l][0u].first; // j = \delta(s)
                const auto h_last = inst.num_component_types - 1u;
                const auto k_last = num_x[j][h_last] - 1u;

                if (subsystems[l][0u].second) { // R
                    model.addConstr(varphi[l][0u] == 1.0 - rho[l][0u][h_last][k_last]);
                } else { // Q
                    model.addConstr(varphi[l][0u] == rho[l][0u][h_last][k_last]);
                }
            }

            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 1u; s < num_subs_in_logic[l]; ++s) {
                    const auto j = subsystems[l][s].first;
                    const auto h_last = inst.num_component_types - 1u;
                    const auto k_last = num_x[j][h_last] - 1u;

                    if (subsystems[l][s].second) { // R
                        model.addConstr(varphi[l][s] == varphi[l][s - 1u] - rho[l][s][h_last][k_last]);
                    } else { // Q
                        model.addConstr(varphi[l][s] == rho[l][s][h_last][k_last]);
                    }
                }
            }

            /* recursive relation between different subsystems of each logic subsystem */
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 1u; s < num_subs_in_logic[l]; ++s) {
                    const auto j = subsystems[l][s].first; // j = \delta(s)
                    const auto big_m = varphi_ub[l][s - 1u];
                    model.addConstr(rho[l][s][0u][0u] <= varphi[l][s - 1u]);
                    model.addConstr(rho[l][s][0u][0u] <= varphi[l][s - 1u] * compute_jhk(j, 0u, num_k[j][0u][0u]) + big_m * (1 - x[j][0u][0u]));
                    model.addConstr(rho[l][s][0u][0u] >= varphi[l][s - 1u] * compute_jhk(j, 0u, num_k[j][0u][0u]));
                    model.addConstr(rho[l][s][0u][0u] >= varphi[l][s - 1u] - big_m * x[j][0u][0u]);
                }
            }

            // recursive relation within the same subsystem
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    const auto j = subsystems[l][s].first; // j = \delta(s)
                    for (auto h = 1u; h < inst.num_component_types; ++h) {
                        auto k_last = num_x[j][h - 1u] - 1u;
                        auto big_m = 1.0;//compute_comb_big_m_positive(subsystems, l, s, h, k);
                        model.addConstr(rho[l][s][h][0u] <= rho[l][s][h - 1u][k_last]);
                        model.addConstr(rho[l][s][h][0u] <= rho[l][s][h - 1u][k_last] * compute_jhk(j, h, num_k[j][h][0u]) + big_m * (1 - x[j][h][0u]));
                        model.addConstr(rho[l][s][h][0u] >= rho[l][s][h - 1u][k_last] * compute_jhk(j, h, num_k[j][h][0u]));
                        model.addConstr(rho[l][s][h][0u] >= rho[l][s][h - 1u][k_last] - big_m * x[j][h][0u]);
                    }
                }
            }

            // recursive relation within the same subsystem
            for (auto l = 0u; l < num_logic_subs; ++l) {
                for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                    const auto j = subsystems[l][s].first; // j = \delta(s)
                    for (auto h = 0u; h < inst.num_component_types; ++h) {
                        for (auto k = 1u; k < num_x[j][h]; ++k) {
                            auto big_m = 1.0;//compute_comb_big_m_positive(subsystems, l, s, h, k);
                            model.addConstr(rho[l][s][h][k] <= rho[l][s][h][k - 1u]);
                            model.addConstr(rho[l][s][h][k] <= rho[l][s][h][k - 1u] * compute_jhk(j, h, num_k[j][h][k]) + big_m * (1 - x[j][h][k]));
                            model.addConstr(rho[l][s][h][k] >= rho[l][s][h][k - 1u] * compute_jhk(j, h, num_k[j][h][k]));
                            model.addConstr(rho[l][s][h][k] >= rho[l][s][h][k - 1u] - big_m * x[j][h][k]);
                        }
                    }
                }
            }

            // at least one component is used
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                GRBLinExpr expr = 0;

                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        expr += x[j][h][k];
                    }
                }

                model.addConstr(expr >= 1);
            }

            // resource constraints
            add_mixed_log_x_constraints(inst, model, x, num_x, num_k);


            /* valid inequalities */

            if (general_logx_model_use_valid_inequalities) {
                for (auto l1 = 0u; l1 < num_logic_subs; ++l1) {
                    std::set<std::pair<size_t, bool>> a;
                    for (auto s1 = 0u; s1 < num_subs_in_logic[l1]; ++s1) {
                        a.insert(subsystems[l1][s1]);
                        for (auto l2 = 0u; l2 < num_logic_subs; ++l2) {
                            if (l1 == l2) continue;

                            std::set<std::pair<size_t, bool>> b;
                            for (auto s2 = 0u; s2 < num_subs_in_logic[l2]; ++s2) {
                                b.insert(subsystems[l2][s2]);

                                // check if 'a' is a subset of 'b'
                                if (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                                    model.addConstr(varphi[l2][s2] <= varphi[l1][s1]);
                                    break;
                                }
                            }
                        }
                    }
                }
            }

            /* Objective functions */

            {
                GRBLinExpr expr = inst.qr_system_reliability.const_term;

                for (auto l = 0u; l < num_logic_subs; ++l) {
                    const auto last_sub = num_subs_in_logic[l] - 1u;
                    expr += logic_subsystems[l].second * varphi[l][last_sub];
                }

                if (general_logx_model_use_ub_one) {
                    model.addConstr(expr <= 1.0);
                }

                model.setObjective(expr, GRB_MAXIMIZE);
            }

            // solve

            if (has_time_limit) {
                model.set(GRB_DoubleParam_TimeLimit, time_limit);
            }

            auto start_time = std::chrono::high_resolution_clock::now();

            model.optimize();

            auto end_time = std::chrono::high_resolution_clock::now();
            auto cpu_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            const auto status = model.get(GRB_IntAttr_Status);

            if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
            {
                std::cout << "\t-- Model infeasible or unbounded!\n";
                return std::nullopt;
            }

            if (status != GRB_OPTIMAL)
            {
                std::cout << "\t-- Warning: Solution not optimal!\n";
            }

            std::cout << "\t-- OBJ: " << model.get(GRB_DoubleAttr_ObjVal) << "\n";

            // extract solution

            auto s = ProbSolution{inst};
            s.obj = model.get(GRB_DoubleAttr_ObjVal);
            s.cpu_time = cpu_time;

            if (status == GRB_TIME_LIMIT) {
                s.time_limit = 1;
                s.mip_gap = model.get(GRB_DoubleAttr_MIPGap);
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }
            else {
                s.time_limit = 0;
                s.mip_gap = 0.0;
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }

            s.num_components_at_subsystem.resize(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                s.num_components_at_subsystem[j].resize(inst.num_component_types, 0.0);

                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    for (auto k = 0u; k < num_x[j][h]; ++k)
                    {
                        if (x[j][h][k].get(GRB_DoubleAttr_X) > 0.5)
                        {
                            s.num_components_at_subsystem[j][h] += num_k[h][h][k];
                        }
                    }
                }
            }

            s.system_reliability = compute_system_reliability(s.num_components_at_subsystem);
            std::cout << "\t-- OBJ: " << s.system_reliability << "\n";

            s.num_bin_vars = model.get(GRB_IntAttr_NumBinVars);
            s.num_vars = model.get(GRB_IntAttr_NumConstrs);
            s.num_cons = model.get(GRB_IntAttr_NumVars);
            s.num_nzs = model.get(GRB_IntAttr_NumNZs);
            s.num_bb_nodes = model.get(GRB_DoubleAttr_NodeCount);

            return s;
        }
        catch (GRBException &e)
        {
            std::cerr << "GUROBI ERROR:\n" << "\t" << e.getErrorCode() << "\n\t" << e.getMessage() << "\n";
        }
        catch (...)
        {
            std::cout << "OTHER_ERRORS" << "\n";
        }

        return std::nullopt;
    }

    template<typename T>
    std::optional<ProbSolution> ProbSolver::solve_general_tree_model(const std::vector<T>&  systems) const
    {
        std::vector<double> min_r(systems.size());
        for (auto l = 1u; l < systems.size(); ++l) {
            const auto j = systems[l]->m_subsystem_id.value();
            min_r[l] = *std::min_element(inst.component_reliability_at_subsystem[j].begin(),
                                         inst.component_reliability_at_subsystem[j].end());
        }

        const bool general_tree_model_use_varphi_ub = true;
        const bool general_tree_model_use_rho_ub = true;
        const bool general_tree_model_use_valid_inequalities = true;
        const bool general_tree_model_use_aggregated_expression = false;
        const bool general_tree_model_use_dominated_components = true;
        const bool general_tree_model_use_ub_one = false;

        /* calculate the upper bounds of decision variables \varphi_{l} */
        auto varphi_ub_ = std::vector<double>(systems.size(), 1.0);
        for (auto l = 1u; l < systems.size(); ++l) {
            const auto f = systems[l]->m_father->m_id;
            const auto j = systems[l]->m_subsystem_id.value();
            if (systems[l]->m_qr) { // R
                varphi_ub_[l] = varphi_ub_[f];
            } else { // Q
                varphi_ub_[l] = varphi_ub_[f] * (1.0 - min_r[l]);
            }
        }

        /* calculate the upper bounds of decision variables \rho_{lh} */
        auto rho_ub = std::vector<std::vector<double>>(systems.size());
        for (auto l = 0u; l < systems.size(); ++l) {
            rho_ub[l] = std::vector<double>(inst.num_component_types, 1.0);
        }
        if (general_tree_model_use_rho_ub) {
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto f = systems[l]->m_father->m_id;
                rho_ub[l] = std::vector<double>(inst.num_component_types, varphi_ub_[f]);

                const auto h_last = inst.num_component_types - 1u;
                rho_ub[l][h_last] *= 1 - min_r[l];
            }
        }

        /* should we use the upper bounds of varphi */
        std::vector<double> varphi_ub(systems.size(), 1.0);
        if (general_tree_model_use_varphi_ub) {
            varphi_ub = varphi_ub_;
        }

        /* model and solve */

        try
        {
            GRBEnv env;
            GRBModel model{env};

            /* variables */

            auto rho = std::vector<std::vector<GRBVar>>(systems.size());
            for (auto l = 0u; l < systems.size(); ++l) {
                rho[l] = std::vector<GRBVar>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    rho[l][h] = model.addVar(0.0, rho_ub[l][h], 0.0, GRB_CONTINUOUS);
                }
            }

            auto varphi = std::vector<GRBVar>(systems.size());
            for (auto l = 0u; l < systems.size(); ++l) {
                varphi[l] = model.addVar(0.0, varphi_ub[l], 0.0, GRB_CONTINUOUS);
            }

            auto x = std::vector<std::vector<std::vector<GRBVar>>>(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                x[j] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    x[j][h] = std::vector<GRBVar>(inst.component_ub_at_subsystem[j][h] + 1u);
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        x[j][h][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    }
                }
            }

            if (general_tree_model_use_dominated_components) {
                for (auto [j, h]: inst.dominated_components) {
                    x[j][h][0].set(GRB_DoubleAttr_LB, 1.0);
                    for (auto k = 1u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        x[j][h][k].set(GRB_DoubleAttr_UB, 0.0);
                    }
                }
            }

            // constraints

            // root node [l = 0u]
            for (auto h = 0u; h < inst.num_component_types; ++h) {
                rho[0u][h].set(GRB_DoubleAttr_LB, 1.0);
            }
            varphi[0u].set(GRB_DoubleAttr_LB, 1.0);

            /* recursive relation between different subsystems of each logic subsystem */
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto f = systems[l]->m_father->m_id;
                const auto h_last = inst.num_component_types - 1u;

                if (systems[l]->m_qr) { // R node
                    model.addConstr(varphi[l] == varphi[f] - rho[l][h_last]);
                } else { // Q node
                    model.addConstr(varphi[l] == rho[l][h_last]);
                }
            }

            // other nodes [l\ne 0u, h = 0u]
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = systems[l]->m_subsystem_id.value();
                const auto f = systems[l]->m_father->m_id;

                if (general_tree_model_use_aggregated_expression) {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                        GRBLinExpr expr = 0;
                        for (auto k2 = k; k2 <= inst.component_ub_at_subsystem[j][0u]; ++k2) {
                            expr += x[j][0u][k2];
                        }

                        const auto big_m = varphi_ub[f] * (1.0 - compute_jhk(j, 0u, k));
                        model.addConstr(rho[l][0u] <= varphi[f] * compute_jhk(j, 0u, k) + big_m * (1 - expr));
                    }
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                        GRBLinExpr expr = 0;
                        for (auto k2 = 0; k2 <= k; ++k2) {
                            expr += x[j][0u][k2];
                        }

                        const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][0u]);
                        const auto big_m = varphi_ub[f] * (compute_jhk(j, 0u, k) - compute_jhk(j, 0u, K));
                        model.addConstr(rho[l][0u] >= varphi[f] * compute_jhk(j, 0u, k) - big_m * (1 - expr));
                    }
                } else {
                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][0u]; ++k) {
                        {
                            const auto big_m = varphi_ub[f] * (1.0 - compute_jhk(j, 0u, k));
                            model.addConstr(rho[l][0u] <= varphi[f] * compute_jhk(j, 0u, k) + big_m * (1 - x[j][0u][k]));
                        }
                        {
                            const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][0u]);
                            const auto big_m = varphi_ub[f] * (compute_jhk(j, 0u, k) - compute_jhk(j, 0u, K));
                            model.addConstr(rho[l][0u] >= varphi[f] * compute_jhk(j, 0u, k) - big_m * (1 - x[j][0u][k]));
                        }
                    }
                }
            }

            // the propagation within each node
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto f = systems[l]->m_father->m_id;
                const auto j = systems[l]->m_subsystem_id.value();
                for (auto h = 1u; h < inst.num_component_types; ++h) {
                    if (general_tree_model_use_aggregated_expression) {
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                            GRBLinExpr expr = 0;
                            for (auto k2 = k; k2 <= inst.component_ub_at_subsystem[j][h]; ++k2) {
                                expr += x[j][h][k2];
                            }

                            const auto big_m = rho_ub[l][h - 1u] * (1.0 - compute_jhk(j, h, k));
                            model.addConstr(rho[l][h] <= rho[l][h - 1u] * compute_jhk(j, h, k) + big_m * (1 - expr));
                        }
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                            GRBLinExpr expr = 0;
                            for (auto k2 = 0; k2 <= k; ++k2) {
                                expr += x[j][h][k2];
                            }

                            const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][h]);
                            const auto big_m = rho_ub[l][h - 1u] * (compute_jhk(j, h, k) - compute_jhk(j, h, K));
                            model.addConstr(rho[l][h] >= rho[l][h - 1u] * compute_jhk(j, h, k) - big_m * (1 - expr));
                        }
                    } else {
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                            // model.addGenConstrIndicator(x[j][h][k], 1, rho[l][h] == rho[l][h - 1u] * std::pow(1.0 - inst.component_reliability_at_subsystem[j][h], static_cast<double>(k)));
                            {
                                const auto big_m = rho_ub[l][h - 1u] * (1.0 - compute_jhk(j, h, k));
                                model.addConstr(
                                        rho[l][h] <= rho[l][h - 1u] * compute_jhk(j, h, k) + big_m * (1 - x[j][h][k]));
                            }
                            {
                                const auto K = static_cast<double>(inst.component_ub_at_subsystem[j][h]);
                                const auto big_m = rho_ub[l][h - 1u] * (compute_jhk(j, h, k) - compute_jhk(j, h, K));
                                model.addConstr(
                                        rho[l][h] >= rho[l][h - 1u] * compute_jhk(j, h, k) - big_m * (1 - x[j][h][k]));
                            }
                        }
                    }
                }
            }

            // sum_{k=0}^{u_jh} = 1, \forall j, h: exactly a particular number of components of type h is used
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    GRBLinExpr expr = 0;

                    for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += x[j][h][k];
                    }

                    model.addConstr(expr == 1);
                }
            }

            // at least one component is used
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                GRBLinExpr expr = 0;

                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 1u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                        expr += x[j][h][k];
                    }
                }

                model.addConstr(expr >= 1);
            }

            // resource constraints
            if (inst.test_type == TEST_TYPE::MIXED) {
                add_mixed_constraints(inst, model, x);
            }
            else if (inst.test_type == TEST_TYPE::TEST_ONE) {
                add_test_one_constraints(inst, model, x);
            }
            else if (inst.test_type == TEST_TYPE::TEST_TWO) {
                add_test_two_constraints(inst, model, x);
            }


            /* valid inequalities */
            if (general_tree_model_use_valid_inequalities) {
                // upper bound on the reliability of each subsystem after the components are configured
                for (auto l = 1u; l < systems.size(); ++l) {
                    const auto j = systems[l]->m_subsystem_id.value();
                    for (auto h = 0u; h < inst.num_component_types; ++h) {
                        GRBLinExpr expr = 0;
                        for (auto k = 0u; k <= inst.component_ub_at_subsystem[j][h]; ++k) {
                            expr += compute_jhk(j, h, k) * x[j][h][k];
                        }
                        const auto f = systems[l]->m_father->m_id;
                        model.addConstr(rho[l][h] <= varphi_ub[f] * expr);
                    }
                }

                for (auto l1 = 1u; l1 < systems.size(); ++l1) {
                    for (auto l2 = 1u; l2 < systems.size(); ++l2) {
                        if (l1 == l2) continue;

                        // we check starting from leaf nodes
                        if (!systems[l2]->m_children.empty()) continue;

                        // note: the nodes are stored in such a vector that the index of a father node
                        // is always less that of a child node
                        bool on_same_route = false;
                        {
                            size_t lmin = std::min(l1, l2);
                            size_t lmax = std::max(l1, l2);
                            auto f = systems[lmax]->m_father;
                            while (f) {
                                if (f->m_id == lmin) {
                                    on_same_route = true;
                                    break;
                                }
                                f = f->m_father;
                            }
                        }
                        if (on_same_route) continue;

                        std::vector<size_t> aa{l1};
                        std::vector<size_t> bb{l2};

                        auto f1 = systems[l1]->m_father;
                        while (f1 && f1->m_father) {
                            aa.push_back(f1->m_id);
                            f1 = f1->m_father;
                        }

                        auto f2 = systems[l2]->m_father;
                        while (f2 && f2->m_father) {
                            bb.push_back(f2->m_id);
                            f2 = f2->m_father;
                        }

                        std::reverse(bb.begin(), bb.end());

                        std::vector<std::pair<size_t, bool>> a;
                        std::vector<std::pair<size_t, bool>> b;
                        for (auto l: aa) {
                            a.emplace_back(systems[l]->m_subsystem_id.value(), systems[l]->m_qr);
                        }
                        for (auto l: bb) {
                            b.emplace_back(systems[l]->m_subsystem_id.value(), systems[l]->m_qr);
                        }

                        if (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                            auto last = b.back();
                            auto last_id = bb.back();
                            b.pop_back();
                            bb.pop_back();
                            while (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                                last = b.back();
                                last_id = bb.back();
                                b.pop_back();
                                bb.pop_back();
                            }

                            model.addConstr(varphi[last_id] <= varphi[l1]);
                        }
                    }
                }
            }


            /* Objective functions */

            {
                GRBLinExpr expr = 0;
                for (auto l = 0u; l < systems.size(); ++l) {
                    if (systems[l]->m_logic) {
                        expr += systems[l]->m_coeff * varphi[l];
                    }
                }

                if (general_tree_model_use_ub_one) {
                    model.addConstr(expr <= 1.0);
                }
                model.setObjective(expr, GRB_MAXIMIZE);
            }

            // solve

            if (has_time_limit) {
                model.set(GRB_DoubleParam_TimeLimit, time_limit);
            }

            auto start_time = std::chrono::high_resolution_clock::now();

            model.optimize();

            auto end_time = std::chrono::high_resolution_clock::now();
            auto cpu_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            const auto status = model.get(GRB_IntAttr_Status);

            if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
            {
                std::cout << "\t-- Model infeasible or unbounded!\n";
                return std::nullopt;
            }

            if (status != GRB_OPTIMAL)
            {
                std::cout << "\t-- Warning: Solution not optimal!\n";
            }

            std::cout << "\t-- OBJ: " << model.get(GRB_DoubleAttr_ObjVal) << "\n";

            // extract solution

            auto s = ProbSolution{inst};
            s.obj = model.get(GRB_DoubleAttr_ObjVal);
            s.cpu_time = cpu_time;

            if (status == GRB_TIME_LIMIT) {
                s.time_limit = 1;
                s.mip_gap = model.get(GRB_DoubleAttr_MIPGap);
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }
            else {
                s.time_limit = 0;
                s.mip_gap = 0.0;
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }

            s.num_components_at_subsystem.resize(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                s.num_components_at_subsystem[j].resize(inst.num_component_types);

                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    for (auto k = 0u; k < inst.component_ub_at_subsystem[j][h]; ++k)
                    {
                        if (x[j][h][k].get(GRB_DoubleAttr_X) > 0.5)
                        {
                            s.num_components_at_subsystem[j][h] = k;
                        }
                    }
                }
            }

            s.system_reliability = compute_system_reliability(s.num_components_at_subsystem);
            std::cout << "\t-- OBJ: " << s.system_reliability << "\n";

            s.num_bin_vars = model.get(GRB_IntAttr_NumBinVars);
            s.num_vars = model.get(GRB_IntAttr_NumConstrs);
            s.num_cons = model.get(GRB_IntAttr_NumVars);
            s.num_nzs = model.get(GRB_IntAttr_NumNZs);
            s.num_bb_nodes = model.get(GRB_DoubleAttr_NodeCount);

            return s;
        }
        catch (GRBException &e)
        {
            std::cerr << "GUROBI ERROR:" << "\n\t" << e.getErrorCode() << "\n\t" << e.getMessage() << "\n";
        }
        catch (...)
        {
            std::cerr << "OTHER_ERRORS" << "\n";
        }

        return std::nullopt;
    }

    template<typename T>
    std::optional<ProbSolution> ProbSolver::solve_general_tree_logx_model(const std::vector<T>& systems) const
    {
        std::vector<double> min_r(systems.size());
        for (auto l = 1u; l < systems.size(); ++l) {
            const auto j = systems[l]->m_subsystem_id.value();
            min_r[l] = *std::min_element(inst.component_reliability_at_subsystem[j].begin(),
                                         inst.component_reliability_at_subsystem[j].end());
        }

        auto num_x = compute_num_x(inst);
        auto num_k = compute_num_k(inst, num_x);

        const bool general_tree_logx_model_use_varphi_ub = true;
        const bool general_tree_logx_model_use_rho_ub = true;
        const bool general_tree_logx_model_use_valid_inequalities = true;
        const bool general_tree_logx_model_use_dominated_components = true;
        const bool general_tree_logx_model_use_ub_one = false;

        /* calculate the upper bounds of decision variables \varphi_{l} */
        auto varphi_ub_ = std::vector<double>(systems.size(), 1.0);
        for (auto l = 1u; l < systems.size(); ++l) {
            const auto f = systems[l]->m_father->m_id;
            const auto j = systems[l]->m_subsystem_id.value();
            if (systems[l]->m_qr) { // R
                varphi_ub_[l] = varphi_ub_[f];
            } else { // Q
                varphi_ub_[l] = varphi_ub_[f] * (1.0 - min_r[l]);
            }
        }

        /* calculate the upper bounds of decision variables \rho_{lh} */
        auto rho_ub = std::vector<std::vector<double>>(systems.size());
        for (auto l = 0u; l < systems.size(); ++l) {
            rho_ub[l] = std::vector<double>(inst.num_component_types, 1.0);
        }
        if (general_tree_logx_model_use_rho_ub) {
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto f = systems[l]->m_father->m_id;
                rho_ub[l] = std::vector<double>(inst.num_component_types, varphi_ub_[f]);

                const auto h_last = inst.num_component_types - 1u;
                rho_ub[l][h_last] *= 1 - min_r[l];
            }
        }

        /* should we use the upper bounds of varphi */
        std::vector<double> varphi_ub(systems.size(), 1.0);
        if (general_tree_logx_model_use_varphi_ub) {
            varphi_ub = varphi_ub_;
        }

        /* model and solve */

        try
        {
            GRBEnv env;
            GRBModel model{env};

            /* variables */

            auto rho = std::vector<std::vector<std::vector<GRBVar>>>(systems.size());
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = systems[l]->m_subsystem_id.value();
                rho[l] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    rho[l][h] = std::vector<GRBVar>(num_x[j][h]);
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        rho[l][h][k] = model.addVar(0.0, rho_ub[l][h], 0.0, GRB_CONTINUOUS);
                    }
                }
            }

            auto varphi = std::vector<GRBVar>(systems.size());
            for (auto l = 0u; l < systems.size(); ++l) {
                varphi[l] = model.addVar(0.0, varphi_ub[l], 0.0, GRB_CONTINUOUS);
            }

            auto x = std::vector<std::vector<std::vector<GRBVar>>>(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                x[j] = std::vector<std::vector<GRBVar>>(inst.num_component_types);
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    x[j][h] = std::vector<GRBVar>(num_x[j][h]);
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        x[j][h][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
                    }
                }
            }

            if (general_tree_logx_model_use_dominated_components) {
                for (auto [j, h]: inst.dominated_components) {
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        x[j][h][k].set(GRB_DoubleAttr_UB, 0.0);
                    }
                }
            }

            // constraints

            // root node [l = 0u]
            varphi[0u].set(GRB_DoubleAttr_LB, 1.0);

            /* recursive relation between different subsystems of each logic subsystem */
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = systems[l]->m_subsystem_id.value();
                const auto f = systems[l]->m_father->m_id;
                const auto h_last = inst.num_component_types - 1u;
                const auto k_last = num_x[j][h_last] - 1u;

                if (systems[l]->m_qr) { // R node
                    model.addConstr(varphi[l] == varphi[f] - rho[l][h_last][k_last]);
                } else { // Q node
                    model.addConstr(varphi[l] == rho[l][h_last][k_last]);
                }
            }

            // other nodes [l\ne 0u, h = 0u]
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = systems[l]->m_subsystem_id.value();
                const auto f = systems[l]->m_father->m_id;

                const auto big_m = 1.0;//varphi_ub[f];
                model.addConstr(rho[l][0u][0u] <= varphi[f]);
                model.addConstr(rho[l][0u][0u] <= varphi[f] * compute_jhk(j, 0u, num_k[j][0][0u]) + big_m * (1 - x[j][0u][0u]));
                model.addConstr(rho[l][0u][0u] >= varphi[f] * compute_jhk(j, 0u, num_k[j][0][0u]));
                model.addConstr(rho[l][0u][0u] >= varphi[f] - big_m * x[j][0u][0]);
            }

            // the propagation within each node
            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = systems[l]->m_subsystem_id.value();
                for (auto h = 1u; h < inst.num_component_types; ++h) {
                    const auto last_k = num_x[j][h - 1u] - 1u;
                    const auto big_m = 1.0;//rho_ub[l][h - 1u];
                    model.addConstr(rho[l][h][0u] <= rho[l][h - 1u][last_k]);
                    model.addConstr(rho[l][h][0u] <= rho[l][h - 1u][last_k] * compute_jhk(j, h, num_k[j][h][0u]) + big_m * (1 - x[j][h][0u]));
                    model.addConstr(rho[l][h][0u] >= rho[l][h - 1u][last_k] * compute_jhk(j, h, num_k[j][h][0u]));
                    model.addConstr(rho[l][h][0u] >= rho[l][h - 1u][last_k] - big_m * x[j][h][0u]);
                }
            }

            for (auto l = 1u; l < systems.size(); ++l) {
                const auto j = systems[l]->m_subsystem_id.value();
                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 1u; k < num_x[j][h]; ++k) {
                        const auto big_m = 1.0;//rho_ub[l][h];
                        model.addConstr(rho[l][h][k] <= rho[l][h][k - 1u]);
                        model.addConstr(rho[l][h][k] <= rho[l][h][k - 1u] * compute_jhk(j, h, num_k[j][h][k]) + big_m * (1 - x[j][h][k]));
                        model.addConstr(rho[l][h][k] >= rho[l][h][k - 1u] * compute_jhk(j, h, num_k[j][h][k]));
                        model.addConstr(rho[l][h][k] >= rho[l][h][k - 1u] - big_m * x[j][h][k]);
                    }
                }
            }

            // exactly a particular number of components of type h is used
            // - \sum_{k=0}^{u_jh} x_{jh}^{k} >= 1, \forall j\in J, h\in H_j.
            for (auto j = 0u; j < inst.num_subsystems; ++j) {
                GRBLinExpr expr = 0;

                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    for (auto k = 0u; k < num_x[j][h]; ++k) {
                        expr += x[j][h][k];
                    }
                }

                model.addConstr(expr >= 1);
            }

            // resource constraints
            add_mixed_log_x_constraints(inst, model, x, num_x, num_k);

            /* valid inequalities */
            if (general_tree_logx_model_use_valid_inequalities) {
                for (auto l1 = 1u; l1 < systems.size(); ++l1) {
                    for (auto l2 = 1u; l2 < systems.size(); ++l2) {
                        if (l1 == l2) continue;

                        // we check starting from leaf nodes
                        if (!systems[l2]->m_children.empty()) continue;

                        // note: the nodes are stored in such a vector that the index of a father node
                        // is always less that of a child node
                        bool on_same_route = false;
                        {
                            size_t lmin = std::min(l1, l2);
                            size_t lmax = std::max(l1, l2);
                            auto f = systems[lmax]->m_father;
                            while (f) {
                                if (f->m_id == lmin) {
                                    on_same_route = true;
                                    break;
                                }
                                f = f->m_father;
                            }
                        }
                        if (on_same_route) continue;

                        std::vector<size_t> aa{l1};
                        std::vector<size_t> bb{l2};

                        auto f1 = systems[l1]->m_father;
                        while (f1 && f1->m_father) {
                            aa.push_back(f1->m_id);
                            f1 = f1->m_father;
                        }

                        auto f2 = systems[l2]->m_father;
                        while (f2 && f2->m_father) {
                            bb.push_back(f2->m_id);
                            f2 = f2->m_father;
                        }

                        std::reverse(bb.begin(), bb.end());

                        std::vector<std::pair<size_t, bool>> a;
                        std::vector<std::pair<size_t, bool>> b;
                        for (auto l: aa) {
                            a.emplace_back(systems[l]->m_subsystem_id.value(), systems[l]->m_qr);
                        }
                        for (auto l: bb) {
                            b.emplace_back(systems[l]->m_subsystem_id.value(), systems[l]->m_qr);
                        }

                        if (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                            auto last = b.back();
                            auto last_id = bb.back();
                            b.pop_back();
                            bb.pop_back();
                            while (std::includes(b.begin(), b.end(), a.begin(), a.end())) {
                                last = b.back();
                                last_id = bb.back();
                                b.pop_back();
                                bb.pop_back();
                            }

                            model.addConstr(varphi[last_id] <= varphi[l1]);
                        }
                    }
                }
            }


            /* Objective functions */

            {
                GRBLinExpr expr = 0;
                for (auto l = 0u; l < systems.size(); ++l)
                {
                    if (systems[l]->m_logic)
                    {
                        expr += systems[l]->m_coeff * varphi[l];
                    }
                }

                if (general_tree_logx_model_use_ub_one) {
                    model.addConstr(expr <= 1.0);
                }
                model.setObjective(expr, GRB_MAXIMIZE);
            }

            // solve

            if (has_time_limit) {
                model.set(GRB_DoubleParam_TimeLimit, time_limit);
            }

            auto start_time = std::chrono::high_resolution_clock::now();

            model.optimize();

            auto end_time = std::chrono::high_resolution_clock::now();
            auto cpu_time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();

            const auto status = model.get(GRB_IntAttr_Status);

            if (status == GRB_INF_OR_UNBD || status == GRB_INFEASIBLE || status == GRB_UNBOUNDED)
            {
                std::cout << "\t-- Model infeasible or unbounded!\n";
                return std::nullopt;
            }

            if (status != GRB_OPTIMAL)
            {
                std::cout << "\t-- Warning: Solution not optimal!\n";
            }

            std::cout << "\t-- OBJ: " << model.get(GRB_DoubleAttr_ObjVal) << "\n";

            // extract solution

            auto s = ProbSolution{inst};
            s.obj = model.get(GRB_DoubleAttr_ObjVal);
            s.cpu_time = cpu_time;

            if (status == GRB_TIME_LIMIT) {
                s.time_limit = 1;
                s.mip_gap = model.get(GRB_DoubleAttr_MIPGap);
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }
            else {
                s.time_limit = 0;
                s.mip_gap = 0.0;
                s.obj_bound = model.get(GRB_DoubleAttr_ObjBound);
            }

            s.num_components_at_subsystem.resize(inst.num_subsystems);
            for (auto j = 0u; j < inst.num_subsystems; ++j)
            {
                s.num_components_at_subsystem[j].resize(inst.num_component_types, 0.0);

                for (auto h = 0u; h < inst.num_component_types; ++h)
                {
                    for (auto k = 0u; k < num_x[j][h]; ++k)
                    {
                        if (x[j][h][k].get(GRB_DoubleAttr_X) > 0.5)
                        {
                            s.num_components_at_subsystem[j][h] += num_k[j][h][k];
                        }
                    }
                }
            }

            s.system_reliability = compute_system_reliability(s.num_components_at_subsystem);
            std::cout << "\t-- OBJ: " << s.system_reliability << "\n";

            s.num_bin_vars = model.get(GRB_IntAttr_NumBinVars);
            s.num_vars = model.get(GRB_IntAttr_NumConstrs);
            s.num_cons = model.get(GRB_IntAttr_NumVars);
            s.num_nzs = model.get(GRB_IntAttr_NumNZs);
            s.num_bb_nodes = model.get(GRB_DoubleAttr_NodeCount);

            return s;
        }
        catch (GRBException &e)
        {
            std::cerr << "GUROBI ERROR:" << "\n\t" << e.getErrorCode() << "\n\t" << e.getMessage() << "\n";
        }
        catch (...)
        {
            std::cerr << "OTHER_ERRORS" << "\n";
        }

        return std::nullopt;
    }

    /* some useful functions */
    double ProbSolver::compute_system_reliability(const std::vector<std::vector<double>>& x) const {

        const auto num_logic_subs { inst.system_reliability.m_other_terms.size() };

        /* save product terms for the following usage */
        std::vector<std::pair<std::set<int>, int>> logic_subsystems;
        std::vector<size_t> num_subs_in_logic;
        logic_subsystems.reserve(num_logic_subs);
        num_subs_in_logic.reserve(num_logic_subs);
        for (auto& [term, coeff] : inst.system_reliability.m_other_terms){
            logic_subsystems.emplace_back(term, coeff);
            num_subs_in_logic.push_back(term.size());
        }

        /* the subsystems in each logic subsystem */
        std::vector<std::vector<size_t>> subsystems(num_logic_subs);
        for (auto l = 0; l < num_logic_subs; ++l)
        {
            std::copy(logic_subsystems[l].first.begin(), logic_subsystems[l].first.end(),
                      std::back_inserter(subsystems[l]));
        }

        std::vector<double> r(num_logic_subs, 1.0);

        for (auto l = 0u; l < num_logic_subs; ++l) {
            for (auto s = 0u; s < num_subs_in_logic[l]; ++s) {
                auto j = subsystems[l][s]; // j = \delta(x)

                for (auto h = 0u; h < inst.num_component_types; ++h) {
                    r[l] *= std::pow(1.0 - inst.component_reliability_at_subsystem[j][h], x[j][h]);
                }
            }
        }

        double obj = inst.system_reliability.const_term;
        for (auto l = 0u; l < num_logic_subs; ++l) {
            obj += logic_subsystems[l].second * r[l];
        }

        return obj;
    }

    double ProbSolver::compute_jhk(const size_t j, const size_t h, const double k) const {
        return std::pow(1.0 - inst.component_reliability_at_subsystem[j][h], k);
    }

} // namespace