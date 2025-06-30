#pragma once
#include <optional>
#include <vector>
#include "mm_instance.h"
#include "mm_solution.h"

namespace rrap {
    enum struct SOLUTION_METHOD {
        Q_BASIC                   = 110,
        Q_BASIC_LOG_X             = 111,
        Q_TREE                    = 120,
        Q_TREE_LOG_X              = 121,
        Q_OB_TREE                 = 130,
        Q_OB_TREE_LOG_X           = 131,
        GENERAL_BASIC             = 210,
        GENERAL_LOG_X             = 211,
        GENERAL_TREE              = 220,
        GENERAL_TREE_LOG_X        = 221,
        GENERAL_OB_TREE           = 230,
        GENERAL_OB_TREE_LOG_X     = 231,
        ERROR_METHOD              = 999
    };

    SOLUTION_METHOD get_solution_method_code(int method_identifier);
    std::string get_solution_method_name(SOLUTION_METHOD method_code);
    std::string get_solution_method_name(int method_identifier);

    struct ProbSolver {
        explicit ProbSolver(const ProbInstance& instance) : inst{instance} {}

        [[nodiscard]] std::optional<ProbSolution> solve(int method) const;

        void set_time_limit_flag(const bool value) {
            has_time_limit = value;
        }
        void set_time_limit_value(const double value) {
            time_limit = value;
        }

    private:
        const ProbInstance& inst;
        bool has_time_limit = false;
        double time_limit = 0;

        /* Q-form models */
        [[nodiscard]] std::optional<ProbSolution> solve_Q_basic_model() const;
        [[nodiscard]] std::optional<ProbSolution> solve_Q_basic_logx_model() const;
        template<typename T>
        [[nodiscard]] std::optional<ProbSolution> solve_Q_tree_model(const std::vector<T>&  systems) const;
        template<typename T>
        [[nodiscard]] std::optional<ProbSolution> solve_Q_tree_logx_model(const std::vector<T>&  systems) const;

        /* general form models */
        [[nodiscard]] std::optional<ProbSolution> solve_general_basic_model() const;
        [[nodiscard]] std::optional<ProbSolution> solve_general_basic_logx_model() const;
        template<typename T>
        [[nodiscard]] std::optional<ProbSolution> solve_general_tree_model(const std::vector<T>& systems) const;
        template<typename T>
        [[nodiscard]] std::optional<ProbSolution> solve_general_tree_logx_model(const std::vector<T>& systems) const;

        /* helper functions */
        // [[nodiscard]] double compute_lower_bound(size_t l) const;
        [[nodiscard]] double compute_system_reliability(const std::vector<std::vector<double>>& x) const;
        [[nodiscard]] double compute_jhk(size_t j, size_t h, double k ) const;
        };
}
