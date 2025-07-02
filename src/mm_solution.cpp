#include <cmath>
#include "mm_solution.h"

namespace rrap {

    bool ProbSolution::check_feasibility() const {
        // resource limit
        for(auto i = 0u; i < this->instance.num_resource_types; ++i){
            auto used = 0.0;

            for(auto j = 0u; j < this->instance.num_subsystems; ++j){
                for(auto h = 0u; h < this->instance.num_component_types; ++h){
                   used += this->instance.resource_for_component_at_subsystem[i][j][h] * static_cast<double>(this->num_components_at_subsystem[j][h]);
                }
            }
            if(this->instance.amount_of_each_resources[i] < used){
                std::cerr << "Resource limit is violated!" << std::endl;
                return false;
            }
        }

        // at least one component is used at each subsystem
        for(auto j = 0u; j < this->instance.num_subsystems; ++j){
            auto used = 0ULL;
            for(auto h = 0u; h < this->instance.num_component_types; ++h){
                used += this->num_components_at_subsystem[j][h];
            }
            if(used < 1u){
                std::cout << "At least one component is needed!" << std::endl;
                return false;
            }
        }
        
        std::vector<size_t> v = {1u,1u,0u,0u,0u,3u,0u,0u,0u,0u,1u,0u,1u,0u,0u,0u,0u,1u,0u,0u,1u,0u,0u,0u,0u,0u,1u,0u};
        std::vector<std::vector<size_t>> ncs(this->instance.num_subsystems, std::vector<size_t>(this->instance.num_component_types));

        for(auto j = 0u, c = 0u; j < this->instance.num_subsystems; ++j){
            for(auto h = 0u; h < this->instance.num_component_types; ++h, ++c){
                ncs[j][h] = v[c];
                std::cout << "===> c = " << c << std::endl;
            }
        }

        for(auto j = 0u, c = 0u; j < this->instance.num_subsystems; ++j){
            for(auto h = 0u; h < this->instance.num_component_types; ++h, ++c){
                ncs[j][h] = v[c];
                std::cout << "===> JH = " << ncs[j][h] << std::endl;
            }
        }

        
        // calculate the system reliability
        std::vector<double> reliabilities(this->instance.num_subsystems, 1.0);
        std::vector<double> unreliability(this->instance.num_subsystems, 1.0);
        for(auto j = 0u; j < this->instance.num_subsystems; ++j){
            for(auto h = 0u; h < this->instance.num_component_types; ++h){
                unreliability[j] *= std::pow(1.0 - this->instance.component_reliability_at_subsystem[j][h], ncs[j][h]);//this->num_components_at_subsystem[j][h]);
            }
            reliabilities[j] -= unreliability[j];
            std::cout << "==> " <<  reliabilities[j] << "\t" << unreliability[j] << std::endl;
        }

        
        auto Q345 = unreliability[2] + reliabilities[2] * unreliability[3] * unreliability[4];
        auto R345 = 1.0 - Q345;
        auto S4 = unreliability[6] * (1.0 - (1.0 - reliabilities[0] * reliabilities[1]) * (1.0 - R345 * reliabilities[5])) 
            + reliabilities[6] * (1.0 - unreliability[0] * Q345) * (1.0 - unreliability[1] * unreliability[5]);
        std::cout << "Objective is " << S4 << std::endl;


        if (std::abs(S4 - this->system_reliability) > 1e-6){
            std::cout << "Objective is incorrect!" << std::endl;
            return false;
        }

        return true;
    }


    std::ostream& operator<<(std::ostream& out, const ProbSolution& s) {
        out << s.cpu_time << ",";
        out << s.num_vars << ",";
        out << s.num_bin_vars << ",";
        out << s.num_cons << ",";

        out << s.system_reliability <<  ",";
        for(auto j = 0u; j < s.instance.num_subsystems; ++j){
            for(auto h = 0u; h < s.instance.num_component_types; ++h)
            out << s.num_components_at_subsystem[j][h] << "\t";
        }

        return out;
    }
}