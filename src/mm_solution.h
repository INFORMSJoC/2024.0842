#pragma once
#include "mm_instance.h"

namespace rrap {
    struct ProbSolution {
        const ProbInstance &instance;

        double obj{0.0};
        double system_reliability{0.0};
        std::vector<std::vector<double>> num_components_at_subsystem;
        double cpu_time{0.0};
        int num_bin_vars{0};
        int num_vars{0};
        int num_cons{0};
        int num_nzs{0};
        double num_bb_nodes{0};  // number of explored branch-and-bound nodes
        int time_limit{0};
        double mip_gap{0.0};
        double obj_bound{0.0};

        explicit ProbSolution(const ProbInstance &instance) :
                instance{instance} {}

        [[nodiscard]] bool check_feasibility() const;
    };

    std::ostream &operator<<(std::ostream &out, const ProbSolution &s);
}
