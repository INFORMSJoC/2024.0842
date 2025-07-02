#include <cstdlib>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <cxxopts.hpp>
#include "mm_instance.h"
#include "mm_solution.h"
#include "mm_solver.h"

using namespace rrap;

int main(int argc, char *argv[]) {
    /* add program options */

    cxxopts::Options options("ra_mm", "Solving RAP via Exact Mathematical Models.");
    options.add_options()
            ("h,help", "display usage information")
            ("s,system", "system identifier associated with the instance [1,2,...,12]", cxxopts::value<int>())
            ("f,file", "path to the instance data file", cxxopts::value<std::string>())
            ("m,method", "model to use:  110: [Q-BASIC], 120: [Q-TREE], 130: [Q-OB-TREE] 210: [BASIC], 220: [TREE], 230: [OB-TREE]", cxxopts::value<int>())
            ("x,time-limit", "time limit in seconds, 0 for unlimited", cxxopts::value<double>()->default_value("0"))
            ("t,test", "test mode", cxxopts::value<int>()->default_value("0"))
            ("v,version", "print the version", cxxopts::value<std::string>());


    /* parse the options */
    auto result = options.parse(argc, argv);

    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return EXIT_SUCCESS;
    }

    if (result.count("version")) {
        std::cout << "Version: 1.0.0" << std::endl;
        return EXIT_SUCCESS;
    }

    if (!result.count("file")) {
        std::cerr << "ERROR: Data file must be provided." << std::endl;
        std::cerr << options.help() << std::endl;
        return EXIT_FAILURE;
    }

    if (!result.count("system")) {
        std::cerr << "ERROR: system index must be provided." << std::endl;
        std::cerr << options.help() << std::endl;
        return EXIT_FAILURE;
    }

    if (!result.count("method")) {
        std::cerr << "ERROR: Method must be provided." << std::endl;
        std::cerr << options.help() << std::endl;
        return EXIT_FAILURE;
    }

    auto system_identifier = result["system"].as<int>();
    auto data_file         = result["file"].as<std::string>();
    auto method_identifier = result["method"].as<int>();
    auto test_type         = result["test"].as<int>();
    auto method_name  = get_solution_method_name(method_identifier);

    std::cout << "--> system structure identifier: " << system_identifier << std::endl;
    // std::cout << "--> number of subsystems: " << num_subsystems << std::endl;
    std::cout << "--> data file: " << data_file << std::endl;
    std::cout << "--> method name: " << method_name << std::endl;

    const auto time_limit = result["time-limit"].as<double>();

    /* initialize the instance */
    ProbInstance instance(test_type, system_identifier, data_file);
    std::cout << instance << std::endl;

    /* initialize the solver */
    ProbSolver s{instance};

    if (static_cast<int>(time_limit) > 0) {
        s.set_time_limit_flag(true);
        s.set_time_limit_value(time_limit);
    }

    /* solve the instance */
    auto sol = s.solve(method_identifier);

    /* check the status */
    if (sol) {
        std::cout << sol.value() << std::endl;
    } else {
        std::cerr << "ERROR: Failed to solve the instance." << std::endl;
    }

    /* output the solution */
    std::filesystem::path out_dir(data_file);
    auto instance_name = out_dir.stem().string();
    auto out_file_name = "_" + instance_name + "_" + method_name + "_sol.txt";

    if (instance.test_type == TEST_TYPE::MIXED) {
        out_file_name.insert(0, std::string("s") + std::to_string(system_identifier));
        std::cout << "TEST_TYPE::MIXED" << std::endl;
    }
    else if (instance.test_type == TEST_TYPE::EXAMPLE_ONE) {
        out_file_name.insert(0, std::string("e1"));
    }
    else if (instance.test_type == TEST_TYPE::EXAMPLE_TWO) {
        out_file_name.insert(0, std::string("e2"));
    }
    else if (instance.test_type == TEST_TYPE::TEST_ONE) {
        out_file_name.insert(0, std::string("t1"));
        std::cout << "TEST_TYPE::TEST_ONE" << std::endl;
    }
    else if (instance.test_type == TEST_TYPE::TEST_TWO) {
        out_file_name.insert(0, std::string("t2"));
    }
    else {
        std::cerr << "error test type..." << std::endl;
    }

    out_dir.replace_filename(out_file_name);

    std::cout << "--> output directory: " << out_dir.string() << std::endl;

    std::ofstream ofs(out_dir);
    ofs << "instance,"
        << "index_of_system,"
        << "number_of_subsystems,"
        << "number_of_component_types,"
        << "number_of_resource_types,"
        << "number_of_q_terms" << ","
        << "total_length_of_q_terms" << ","
        << "min_length_of_q_terms" << ","
        << "max_length_of_q_terms" << ","
        << "avg_length_of_q_terms"<< ","
        << "number_of_q_nodes,"
        << "num_of_q_tree_branches" << ","
        << "total_length_q_tree" << ","
        << "min_length_q_tree" << ","
        << "max_length_q_tree" << ","
        << "avg_length_q_tree" << ","
        << "number_of_q_sorted_nodes,"
        << "num_of_sorted_q_tree_branches" << ","
        << "total_length_sorted_q_tree" << ","
        << "min_length_sorted_q_tree" << ","
        << "max_length_sorted_q_tree" << ","
        << "avg_length_sorted_q_tree" << ","
        << "number_of_qr_terms,"
        << "total_length_of_qr_terms" << ","
        << "min_length_of_qr_terms" << ","
        << "max_length_of_qr_terms" << ","
        << "avg_length_of_qr_terms" << ","
        << "number_of_qr_nodes,"
        << "num_of_qr_tree_branches" << ","
        << "total_length_qr_tree" << ","
        << "min_length_qr_tree" << ","
        << "max_length_qr_tree" << ","
        << "avg_length_qr_tree" << ","
        << "number_of_qr_sorted_nodes,"
        << "num_of_sorted_qr_tree_branches" << ","
        << "total_length_sorted_qr_tree" << ","
        << "min_length_sorted_qr_tree" << ","
        << "max_length_sorted_qr_tree" << ","
        << "avg_length_sorted_qr_tree" << ","
        << "solution_model,"
        << "number_of_vars,"
        << "number_of_bin_vars,"
        << "number_of_cons,"
        << "number_of_nzs,"
        << "number_of_bb_nodes,"
        << "cpu_time,"
        << "reliability,"
        << "obj,"
        << "time_limit,"
        << "mip_gap,"
        << "obj_bound,"
        << "solution\n";

    // the structure properties of the instance
    auto number_of_q_terms = instance.system_reliability.m_other_terms.size();
    auto total_length_of_q_terms = 0u;
    auto min_length_of_q_terms = std::numeric_limits<size_t>::max();
    auto max_length_of_q_terms = std::numeric_limits<size_t>::min();
    for (const auto&[fst, snd] : instance.system_reliability.m_other_terms) {
        total_length_of_q_terms += fst.size();
        min_length_of_q_terms = std::min(min_length_of_q_terms, fst.size());
        max_length_of_q_terms = std::max(max_length_of_q_terms, fst.size());
    }
    auto avg_length_of_q_terms = static_cast<double>(total_length_of_q_terms) / static_cast<double>(number_of_q_terms);

    auto number_of_qr_terms = instance.qr_system_reliability.m_other_terms.size();
    auto total_length_of_qr_terms = 0u;
    auto min_length_of_qr_terms = std::numeric_limits<size_t>::max();
    auto max_length_of_qr_terms = std::numeric_limits<size_t>::min();
    for (const auto&[fst, snd] : instance.qr_system_reliability.m_other_terms) {
        total_length_of_qr_terms += fst.size();
        min_length_of_qr_terms = std::min(min_length_of_qr_terms, fst.size());
        max_length_of_qr_terms = std::max(max_length_of_qr_terms, fst.size());
    }
    auto avg_length_of_qr_terms = static_cast<double>(total_length_of_qr_terms) / static_cast<double>(number_of_qr_terms);

    unsigned min_length_q_tree;
    unsigned max_length_q_tree;
    unsigned total_length_q_tree;
    unsigned num_of_q_tree_branches;
    double avg_length_q_tree;
    compute_tree_lengths(instance.logic_system_tree.systems, min_length_q_tree, max_length_q_tree,
        avg_length_q_tree, total_length_q_tree, num_of_q_tree_branches);

    unsigned min_length_sorted_q_tree;
    unsigned max_length_sorted_q_tree;
    unsigned total_length_sorted_q_tree;
    unsigned num_of_sorted_q_tree_branches;
    double avg_length_sorted_q_tree;
    compute_tree_lengths(instance.sorted_logic_system_tree.systems, min_length_sorted_q_tree,
        max_length_sorted_q_tree, avg_length_sorted_q_tree, total_length_sorted_q_tree,
        num_of_sorted_q_tree_branches);

    unsigned min_length_qr_tree;
    unsigned max_length_qr_tree;
    unsigned total_length_qr_tree;
    unsigned num_of_qr_tree_branches;
    double avg_length_qr_tree;
    compute_tree_lengths(instance.qr_system_tree.systems, min_length_qr_tree, max_length_qr_tree,
        avg_length_qr_tree, total_length_qr_tree, num_of_qr_tree_branches);

    unsigned min_length_sorted_qr_tree;
    unsigned max_length_sorted_qr_tree;
    unsigned total_length_sorted_qr_tree;
    unsigned num_of_sorted_qr_tree_branches;
    double avg_length_sorted_qr_tree;
    compute_tree_lengths(instance.qr_sorted_system_tree.systems, min_length_sorted_qr_tree, max_length_sorted_qr_tree,
        avg_length_sorted_qr_tree, total_length_sorted_qr_tree, num_of_sorted_qr_tree_branches);

    ofs << instance_name << ","
    << system_identifier << ","
    << instance.num_subsystems << ","
    << instance.num_component_types << ","
    << instance.num_resource_types << ","
    << number_of_q_terms << ","
    << total_length_of_q_terms << ","
    << min_length_of_q_terms << ","
    << max_length_of_q_terms << ","
    << avg_length_of_q_terms << ","
    << instance.logic_system_tree.systems.size() << ","
    << num_of_q_tree_branches << ","
    << total_length_q_tree << ","
    << min_length_q_tree << ","
    << max_length_q_tree << ","
    << avg_length_q_tree << ","
    << instance.sorted_logic_system_tree.systems.size() << ","
    << num_of_sorted_q_tree_branches << ","
    << total_length_sorted_q_tree << ","
    << min_length_sorted_q_tree << ","
    << max_length_sorted_q_tree << ","
    << avg_length_sorted_q_tree << ","
    << instance.qr_system_reliability.m_other_terms.size() << ","
    << total_length_of_qr_terms << ","
    << min_length_of_qr_terms << ","
    << max_length_of_qr_terms << ","
    << avg_length_of_qr_terms << ","
    << instance.qr_system_tree.systems.size() << ","
    << num_of_qr_tree_branches << ","
    << total_length_qr_tree << ","
    << min_length_qr_tree << ","
    << max_length_qr_tree << ","
    << avg_length_qr_tree << ","
    << instance.qr_sorted_system_tree.systems.size() << ","
    << num_of_sorted_qr_tree_branches << ","
    << total_length_sorted_qr_tree << ","
    << min_length_sorted_qr_tree << ","
    << max_length_sorted_qr_tree << ","
    << avg_length_sorted_qr_tree << ","
    << method_name << ","
    << sol->num_vars << ","
    << sol->num_bin_vars << ","
    << sol->num_cons << ","
    << sol->num_nzs << ","
    << sol->num_bb_nodes << ","
    << sol->cpu_time << ","
    << sol->system_reliability << ","
    << sol->obj << ","
    << sol->time_limit << ","
    << sol->mip_gap << ","
    << sol->obj_bound << ",";

    for(auto j = 0u; j < instance.num_subsystems; ++j){
        for(auto h = 0u; h < instance.num_component_types; ++h)
           ofs << sol->num_components_at_subsystem[j][h] << "\t";
    }
    ofs << "\n";

    ofs.close();

    return EXIT_SUCCESS;
}
