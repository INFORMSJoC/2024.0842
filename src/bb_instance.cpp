#include <fstream>
#include <iostream>
#include <limits>
#include <algorithm>
#include "bb_instance.h"

void ProbInstance::read_data()
{
    std::ifstream ifs{inst_file};
    if (ifs.fail())
    {
        std::cerr << "Error opening ProbInstance file " << inst_file << "\n";
        std::exit(EXIT_FAILURE);
    }

    // 1st line:
    ifs >> num_resource_types;
    ifs >> num_subsystems;
    ifs >> num_component_types;

    // 2nd line: the amount of each resource
    amount_of_each_resources.resize(num_resource_types);
    for (auto i = 0u; i < num_resource_types; ++i)
    {
        ifs >> amount_of_each_resources[i];
    }

    if (test_type == TEST_TYPE::TEST_ONE || test_type == TEST_TYPE::TEST_TWO || test_type == TEST_TYPE::MIXED) {
        // num_subsystems Lines, num_component_types Columns: r_jh
        component_reliability_at_subsystem.resize(num_subsystems);
        for (auto j = 0u; j < num_subsystems; ++j)
        {
            component_reliability_at_subsystem[j].resize(num_component_types);
        }
        for (auto j = 0u; j < num_subsystems; ++j)
        {
            for (auto h = 0u; h < num_component_types; ++h)
            {
                ifs >> component_reliability_at_subsystem[j][h];
            }
        }

        num_component_types_at_subsystem.resize(num_subsystems, 0);
        for (auto j = 0u; j < num_subsystems; ++j) {
            for (auto h = 0u; h < num_component_types; ++h) {
                if (component_reliability_at_subsystem[j][h] < 1e-6) continue;
                ++num_component_types_at_subsystem[j];
            }
        }

        // (num_resource_types * num_subsystems) lines: A matrix
        resource_for_component_at_subsystem.resize(num_resource_types);
        for (auto i = 0u; i < num_resource_types; ++i)
        {
            resource_for_component_at_subsystem[i].resize(num_subsystems);
            for (auto j = 0u; j < num_subsystems; ++j)
            {
                resource_for_component_at_subsystem[i][j].resize(num_component_types);
            }
        }
        for (int i = 0; i < num_resource_types; ++i)
        {
            for (int j = 0; j < num_subsystems; ++j)
            {
                for (int h = 0; h < num_component_types; h++)
                {
                    ifs >> resource_for_component_at_subsystem[i][j][h];
                }
            }
        }
    }

    if (test_type == TEST_TYPE::EXAMPLE_ONE || test_type == TEST_TYPE::EXAMPLE_TWO ||
        test_type == TEST_TYPE::TEST_ONE || test_type == TEST_TYPE::TEST_TWO) {
        num_component_types_at_subsystem.resize(num_subsystems, 0);
        for (auto j = 0u; j < num_subsystems; ++j) {
            num_component_types_at_subsystem[j] = num_component_types;
        }

        component_lb_at_subsystem.resize(num_subsystems);
        for (int j = 0; j < num_subsystems; j++)
        {
            component_lb_at_subsystem[j].resize(num_component_types_at_subsystem[j]);
            for (int h = 0; h < num_component_types_at_subsystem[j]; ++h) {
                ifs >> component_lb_at_subsystem[j][h];
            }
        }

        component_ub_at_subsystem.resize(num_subsystems);
        for (int j = 0; j < num_subsystems; j++)
        {
            component_ub_at_subsystem[j].resize(num_component_types_at_subsystem[j]);
            for (int h = 0; h < num_component_types_at_subsystem[j]; ++h) {
                ifs >> component_ub_at_subsystem[j][h];
            }
        }
    }
    else if (test_type == TEST_TYPE::MIXED) {
        calculate_component_lb_at_subsystem();
        calculate_component_ub_at_subsystem();
    }

    ifs.close();
}

ProbInstance::ProbInstance(const int test_type_, int system_index_, const std::string &inst_file_) {
    if (test_type_ == 1) {
        test_type = TEST_TYPE::EXAMPLE_ONE;
    }
    else if (test_type_ == 2) {
        test_type = TEST_TYPE::EXAMPLE_TWO;
    }
    else if (test_type_ == 3) {
        test_type = TEST_TYPE::TEST_ONE;
    }
    else if (test_type_ == 4) {
        test_type = TEST_TYPE::TEST_TWO;
    }
    else if (test_type_ == 0) {
        test_type = TEST_TYPE::MIXED;
    }
    else {
        std::cerr << "ERROR TYPE" << std::endl;
    }

    inst_file = inst_file_;
    system_index = system_index_;

    read_data();
}

void ProbInstance::calculate_component_lb_at_subsystem() {
    component_lb_at_subsystem.resize(num_subsystems);

    for (int j = 0; j < num_subsystems; j++)
    {
        component_lb_at_subsystem[j].resize(num_component_types_at_subsystem[j], 0);
    }
}


void ProbInstance::calculate_component_ub_at_subsystem()
{
    component_ub_at_subsystem.resize(num_subsystems);

    for (int j = 0; j < num_subsystems; j++)
    {
        const auto minUsage = calculate_min_resource_usage(j);
        component_ub_at_subsystem[j].resize(num_component_types_at_subsystem[j], std::numeric_limits<int>::max());

        for (int h = 0; h < num_component_types_at_subsystem[j]; h++)
        {
            for (int i = 0; i < num_resource_types; i++)
            {
                component_ub_at_subsystem[j][h] = std::min(component_ub_at_subsystem[j][h],
                                                           static_cast<int>(std::floor((amount_of_each_resources[i] - minUsage[i]) / resource_for_component_at_subsystem[i][j][h])));
            }
        }
    }
}

std::vector<double> ProbInstance::calculate_min_resource_usage(const int j)
{
    std::vector<double> min_res_usage(num_resource_types, 0.0);

    for (int i = 0; i < num_resource_types; ++i)
    {
        for (int ja = 0; ja < num_subsystems; ++ja)
        {
            if (ja == j) continue;

            min_res_usage[i] += *std::min_element(resource_for_component_at_subsystem[i][ja].begin(), resource_for_component_at_subsystem[i][ja].end());
        }

        /* original: min_res_usage[i] = Math.Round(min_res_usage[i], 2); */
        min_res_usage[i] = std::round(min_res_usage[i] * 100) / 100.0;
    }

    return min_res_usage;
}

std::ostream &operator<<(std::ostream &out, const ProbInstance &inst)
{
    out << "num_resource_types: " << inst.num_resource_types << "\n";
    out << "num_subsystems: " << inst.num_subsystems << "\n";
    out << "num_component_types: " << inst.num_component_types << "\n";

    out << "amount_of_each_resources: " << "\n";
    for (auto i = 0u; i < inst.num_resource_types; ++i)
    {
        out << "\t resource " << i + 1u << ": " << inst.amount_of_each_resources[i] << "\n";
    }

    if (inst.test_type == TEST_TYPE::TEST_ONE || inst.test_type == TEST_TYPE::TEST_TWO || inst.test_type == TEST_TYPE::MIXED) {
        out << "component_reliability_at_subsystem: " << "\n";
        for (auto j = 0u; j < inst.num_subsystems; ++j)
        {
            for (auto h = 0u; h < inst.num_component_types_at_subsystem[j]; ++h)
            {
                out << "\t component " << h + 1u << ": " << inst.component_reliability_at_subsystem[j][h] << "\t";
            }
            out << "\n";
        }

        out << "resource_for_component_at_subsystem: " << "\n";
        for (int i = 0; i < inst.num_resource_types; ++i)
        {
            out << "\t resource " << i + 1u << "\n";
            for (int j = 0; j < inst.num_subsystems; ++j)
            {
                for (int h = 0; h < inst.num_component_types_at_subsystem[j]; h++)
                {
                    out << "\t[s" << j + 1u << ", c" << h + 1u << "] = " << inst.resource_for_component_at_subsystem[i][j][h] << "\t";
                }
                out << "\n";
            }
        }
    }

    out << "calculate_component_lb_at_subsystem: " << "\n";
    for (int j = 0; j < inst.num_subsystems; ++j)
    {
        for (int h = 0; h < inst.num_component_types_at_subsystem[j]; h++)
        {
            out << "\tL[s" << j << ", c" << h << "] = " << inst.component_lb_at_subsystem[j][h] << "\t";
        }
        out << "\n";
    }

    out << "calculate_component_ub_at_subsystem: " << "\n";
    for (int j = 0; j < inst.num_subsystems; ++j)
    {
        for (int h = 0; h < inst.num_component_types_at_subsystem[j]; h++)
        {
            out << "\tU[s" << j << ", c" << h << "] = " << inst.component_ub_at_subsystem[j][h] << "\t";
        }
        out << "\n";
    }

    return out;
}