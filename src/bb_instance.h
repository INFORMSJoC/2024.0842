#pragma once
#include <vector>
#include <string>

enum class TEST_TYPE {
    MIXED       = 0,
    EXAMPLE_ONE = 1,
    EXAMPLE_TWO = 2,
    TEST_ONE    = 3,
    TEST_TWO    = 4
};

struct ProbInstance {
    using SolutionSet = std::vector<int>;

    void read_data();
    ProbInstance(int test_type_, int system_index_, const std::string &inst_file_);

    std::string inst_file;
    int system_index{};
    size_t num_resource_types{};   // m: number of resource types
    size_t num_subsystems{};       // ns: number of m_subsystems_on_route
    size_t num_component_types{};  // nh: number of heterogeneous component types

    TEST_TYPE test_type;

    std::vector<std::vector<int>> component_lb_at_subsystem;
    std::vector<std::vector<int>> component_ub_at_subsystem;
    std::vector<std::vector<double>> component_reliability_at_subsystem;

    std::vector<double> amount_of_each_resources;
    std::vector<std::vector<std::vector<double>>> resource_for_component_at_subsystem;

    std::vector<size_t> num_component_types_at_subsystem;

private:
    void calculate_component_lb_at_subsystem();
    void calculate_component_ub_at_subsystem();
    std::vector<double> calculate_min_resource_usage(int j);
};

std::ostream &operator<<(std::ostream &out, const ProbInstance &inst);