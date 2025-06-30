#include <algorithm>
#include <cassert>
#include "mm_reliability.h"

namespace rrap {
        
    // constructors
    QSystemReliability::QSystemReliability(const size_t num_subsystems, const int constant):
            num_subsystems{num_subsystems}, const_term{constant} {}

    QSystemReliability::QSystemReliability(std::initializer_list<int> subsystems, int constant): const_term(constant){
        num_subsystems = subsystems.size();

        std::set<int> term;

        for(auto& index : subsystems){
            term.insert(index);
        }

        m_other_terms[term] += 1;
    }

    void clear_zero_terms(QSystemReliability& s){
        for(auto iter = s.m_other_terms.begin(); iter != s.m_other_terms.end();
            iter = iter->second == 0 ? s.m_other_terms.erase(iter) : std::next(iter)) {
        }
    }

    // operators
    QSystemReliability operator+(const QSystemReliability& left, const QSystemReliability& right)
    {
        assert(left.num_subsystems == right.num_subsystems);

        auto result = left;
        
        result.const_term += right.const_term;

        for(const auto& term : right.m_other_terms){
            result.m_other_terms[term.first] += term.second;
        }

        clear_zero_terms(result);

        return result;
    }

    QSystemReliability operator-(const QSystemReliability& left, const QSystemReliability& right){
        assert(left.num_subsystems == right.num_subsystems);

        auto result = left;
        
        result.const_term -= right.const_term;

        for(const auto& term : right.m_other_terms){
            result.m_other_terms[term.first] -= term.second;
        }

        clear_zero_terms(result);

        return result;
    }

    bool has_common_element(const QSystemReliability& left, const QSystemReliability& right)
    {
        auto first1 = left.m_other_terms.begin();
        auto last1 =left.m_other_terms.end();
        auto first2 = right.m_other_terms.begin();
        auto last2 =right.m_other_terms.end();

        while (first1 != last1 && first2 != last2)
        {
            if (*first1 < *first2) { ++first1; }
            else if (*first2 < *first1) { ++first2; }
            else { return true; }
        }

        return false;
    }

    QSystemReliability operator*(const QSystemReliability& left, const QSystemReliability& right){
        assert(left.num_subsystems == right.num_subsystems);
        assert(!has_common_element(left, right));

        auto result = QSystemReliability(left.num_subsystems);

        result.const_term = left.const_term * right.const_term;

        if (left.const_term != 0) {
            for(auto& term : right.m_other_terms) {
                result.m_other_terms[term.first] += left.const_term * term.second;
            }
        }

        if (right.const_term != 0) {
            for(auto& term : left.m_other_terms) {
                result.m_other_terms[term.first] += right.const_term * term.second;
            }
        }

        for(auto& term1 : right.m_other_terms){
            for(auto& term2 : left.m_other_terms){
                std::set<int> term;
                std::set_union(term1.first.begin(), term1.first.end(),term2.first.begin(), term2.first.end(),
                               std::inserter(term, term.end()));

                result.m_other_terms[term] += term1.second * term2.second;
            }
        }

        clear_zero_terms(result);

        return result;
    } 


    QSystemReliability operator*(int num, const QSystemReliability& other){

        if(num == 0) return QSystemReliability(other.num_subsystems);

        auto result = other;
        result.const_term *= num;

        for (auto& term : result.m_other_terms){
            term.second *= num;
        }

        return result;
    }

    QSystemReliability operator*(const QSystemReliability& other, int num) {
        return num * other;
    }

    QSystemReliability operator+(int num, const QSystemReliability& s){
        auto result = s;
        result.const_term += num;
        return result;
    }


    QSystemReliability operator+(const QSystemReliability& other, int num){
        return (num + other);
    }


    QSystemReliability operator-(int num, const QSystemReliability& s){
        auto result = s;
        result.const_term = num - result.const_term;

        for (auto& term : result.m_other_terms){
            term.second = -term.second;
        }

        return result;
    }


    QSystemReliability operator-(const QSystemReliability& s, int num){
        auto result = s;
        result.const_term -= num;
        return result;
    }

   std::ostream& operator<<(std::ostream& out, const QSystemReliability& s) {
        
        // out << s.other_terms.size() << "TERMS" << "\n";
        out << s.const_term;

       for (auto& term : s.m_other_terms){
           if (term.second > 0){
               out << " + " << term.second << "Q_";
           } else if (term.second < 0){
               out << " - " << -term.second << "Q_";
           }
           else{
               out << "ERROR_TERM: " << term.second << "\n";
           }

           for(const auto& t : term.first){
               out << t;
           }
       }

        return out;
    }

    QSystemReliability& QSystemReliability::operator+=(const QSystemReliability& rhs){
        assert(this->num_subsystems == rhs.num_subsystems);

        this->const_term += rhs.const_term;

        for(const auto& term : rhs.m_other_terms){
            this->m_other_terms[term.first] += term.second;
        }

        clear_zero_terms(*this);

        return (*this);
    }

    QSystemReliability& QSystemReliability::operator-=(const QSystemReliability& rhs){
        assert(this->num_subsystems == rhs.num_subsystems);

        this->const_term -= rhs.const_term;

        for(const auto& term : rhs.m_other_terms){
            this->m_other_terms[term.first] -= term.second;
        }

        clear_zero_terms(*this);

        return (*this);
    }

    QSystemReliability& QSystemReliability::operator*=(int multi){
        if(multi == 0) {
            this->const_term = 0;
            this->m_other_terms.clear();
        }
        else{
            this->const_term *= multi;

            for (auto& term : this->m_other_terms){
                term.second *= multi;
            }
        }

        return (*this);
    }


    QSystemReliability::QSystemReliability(const QSystemReliability &from) {
        num_subsystems = from.num_subsystems;
        const_term = from.const_term;
        m_other_terms = from.m_other_terms;
    }

    QSystemReliability& QSystemReliability::operator=(const QSystemReliability &from) {
        if(this != &from){
            num_subsystems = from.num_subsystems;
            const_term = from.const_term;
            m_other_terms = from.m_other_terms;
        }

        return *this;
    }


    // =================================================================

    QRSystemReliability::QRSystemReliability(const size_t num_subsystems, const int constant) :
            num_subsystems{num_subsystems}, const_term{constant} {}


    QRSystemReliability::QRSystemReliability(size_t num_subsystems, int s, bool qr, int coeff, int constant) :
            num_subsystems {num_subsystems}, const_term {constant} {

        std::set<std::pair<int, bool>> term { std::make_pair(s, qr) };
        m_other_terms[term] = coeff;
    }

    QRSystemReliability::QRSystemReliability(size_t num_subsystems, const std::vector<int> &ss,
                                             const std::vector<bool> &qrs, const std::vector<int> &coeffs, int constant)
            :
            num_subsystems {num_subsystems}, const_term {constant} {

        assert(ss.size() == qrs.size() && ss.size() == coeffs.size());

        for (auto s = 0; s < ss.size(); ++s) {
            std::set<std::pair<int, bool>> term{ std::make_pair(ss[s], qrs[s]) };
            m_other_terms[term] = coeffs[s];
        }
    }

    QRSystemReliability::QRSystemReliability(std::initializer_list<std::pair<int, bool>> subsystems, int constant) : const_term(constant){
        num_subsystems = subsystems.size();

        std::set<std::pair<int,bool>> term;

        for(auto& index : subsystems){
            term.insert(index);
        }

        m_other_terms[term] += 1;
    }

    QRSystemReliability::QRSystemReliability(const QRSystemReliability &from) {
        num_subsystems = from.num_subsystems;
        const_term = from.const_term;
        m_other_terms = from.m_other_terms;
    }

    QRSystemReliability &QRSystemReliability::operator=(const QRSystemReliability &from) {
        if(this != &from){
            num_subsystems = from.num_subsystems;
            const_term = from.const_term;
            m_other_terms = from.m_other_terms;
        }

        return *this;
    }

    void clear_zero_terms(QRSystemReliability& s){
        for(auto iter = s.m_other_terms.begin(); iter != s.m_other_terms.end();
            iter = iter->second == 0 ? s.m_other_terms.erase(iter) : std::next(iter)) {
        }
    }



    // QR

    QRSystemReliability operator+(const QRSystemReliability& left, const QRSystemReliability& right)
    {
        assert(left.num_subsystems == right.num_subsystems);

        auto result = left;

        result.const_term += right.const_term;

        for(const auto& term : right.m_other_terms){
            result.m_other_terms[term.first] += term.second;
        }

        clear_zero_terms(result);

        return result;
    }

    QRSystemReliability operator-(const QRSystemReliability& left, const QRSystemReliability& right){
        assert(left.num_subsystems == right.num_subsystems);

        auto result = left;

        result.const_term -= right.const_term;

        for(const auto& term : right.m_other_terms){
            result.m_other_terms[term.first] -= term.second;
        }

        clear_zero_terms(result);

        return result;
    }

    bool has_common_element(const QRSystemReliability& left, const QRSystemReliability& right)
    {
        auto first1 = left.m_other_terms.begin();
        auto last1 =left.m_other_terms.end();
        auto first2 = right.m_other_terms.begin();
        auto last2 =right.m_other_terms.end();

        while (first1 != last1 && first2 != last2)
        {
            if (*first1 < *first2) { ++first1; }
            else if (*first2 < *first1) { ++first2; }
            else { return true; }
        }

        return false;
    }

    QRSystemReliability operator*(const QRSystemReliability& left, const QRSystemReliability& right){
        assert(left.num_subsystems == right.num_subsystems);
        assert(!has_common_element(left, right));

        auto result = QRSystemReliability(left.num_subsystems);

        result.const_term = left.const_term * right.const_term;

        if (left.const_term != 0) {
            for(auto& term : right.m_other_terms) {
                result.m_other_terms[term.first] += left.const_term * term.second;
            }
        }

        if (right.const_term != 0) {
            for(auto& term : left.m_other_terms) {
                result.m_other_terms[term.first] += right.const_term * term.second;
            }
        }

        for(auto& term1 : right.m_other_terms){
            for(auto& term2 : left.m_other_terms){
                std::set<std::pair<int,bool>> term;
                std::set_union(term1.first.begin(), term1.first.end(),term2.first.begin(), term2.first.end(),
                               std::inserter(term, term.end()));

                result.m_other_terms[term] += term1.second * term2.second;
            }
        }

        clear_zero_terms(result);

        return result;
    }


    QRSystemReliability operator*(int num, const QRSystemReliability& other){

        if(num == 0) return QRSystemReliability(other.num_subsystems);

        auto result = other;
        result.const_term *= num;

        for (auto& term : result.m_other_terms){
            term.second *= num;
        }

        return result;
    }

    QRSystemReliability operator*(const QRSystemReliability& other, int num) {
        return num * other;
    }

    QRSystemReliability operator+(int num, const QRSystemReliability& s){
        auto result = s;
        result.const_term += num;
        return result;
    }


    QRSystemReliability operator+(const QRSystemReliability& other, int num){
        return (num + other);
    }


    QRSystemReliability operator-(int num, const QRSystemReliability& s){
        auto result = s;
        result.const_term = num - result.const_term;

        for (auto& term : result.m_other_terms){
            term.second = -term.second;
        }

        return result;
    }



    QRSystemReliability operator-(const QRSystemReliability& s, int num){
        auto result = s;
        result.const_term -= num;
        return result;
    }

    QRSystemReliability &QRSystemReliability::operator+=(const QRSystemReliability &rhs) {
        assert(this->num_subsystems == rhs.num_subsystems);

        this->const_term += rhs.const_term;

        for(const auto& term : rhs.m_other_terms){
            this->m_other_terms[term.first] += term.second;
        }

        clear_zero_terms(*this);

        return (*this);
    }

    QRSystemReliability &QRSystemReliability::operator-=(const QRSystemReliability &rhs) {
        assert(this->num_subsystems == rhs.num_subsystems);

        this->const_term -= rhs.const_term;

        for(const auto& term : rhs.m_other_terms){
            this->m_other_terms[term.first] -= term.second;
        }

        clear_zero_terms(*this);

        return (*this);
    }

    QRSystemReliability &QRSystemReliability::operator*=(int multi) {
        if(multi == 0) {
            this->const_term = 0;
            this->m_other_terms.clear();
        }
        else{
            this->const_term *= multi;

            for (auto& term : this->m_other_terms){
                term.second *= multi;
            }
        }

        return (*this);
    }

    std::ostream& operator<<(std::ostream& out, const QRSystemReliability& s) {

        // out << s.other_terms.size() << "TERMS" << "\n";
        out << s.const_term;

        for (auto& term : s.m_other_terms){
            if (term.second > 0){
                out << " + " << term.second;
            } else if (term.second < 0){
                out << " - " << -term.second;
            }
            else{
                out << "ERROR_TERM: " << term.second << "\n";
            }

            for(const auto& t : term.first){
                if (t.second) {
                    out << "R" << t.first;
                }
                else{
                    out << "Q" << t.first;
                }
            }
        }

        return out;
    }
}