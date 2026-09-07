#include "r410a.hpp"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main() {
    // CoolProp 7.2.0 HEOS fixed-composition R410A, obtained independently through
    // PropsSI(D,T,t,P,p,R410A). This verifies implementation, not oil-mixture data.
    constexpr std::array<std::array<double,3>,15> references{{
        {293.15,1e5,3.019666669276735}, {293.15,1e6,35.23441309509108},
        {313.15,1e5,2.818234879413699}, {313.15,1e6,31.549508396000213}, {313.15,2e6,75.71094588895434},
        {333.15,1e5,2.643202061694997}, {333.15,1e6,28.81478249367936}, {333.15,2e6,65.09601484848768},
        {353.15,1e5,2.4893850349057947}, {353.15,1e6,26.642335228450023}, {353.15,2e6,58.24076465999728},
        {313.15,2.2e6,87.68721787568509}, {333.15,3e6,116.54814841145676},
        {353.15,5e6,241.2714995900332}, {353.15,1e7,840.191953522374}
    }};
    double worst=0;
    constexpr std::array<std::array<double,4>,4> caloric{{
        {293.15,1e5,422667.20917943725,692.9550441891065},
        {313.15,2e6,411317.29549059365,900.2977111820805},
        {333.15,1e6,443288.19531610294,793.9692552712804},
        {353.15,2e6,451126.7486514086,850.2591188349703}
    }};
    for (const auto& row : caloric) {
        const double rho=R410A::gas_density(row[1],row[0]);
        if (std::abs(R410A::internal_energy(rho,row[0])/row[2]-1)>1e-10 ||
            std::abs(R410A::heat_capacity_cv(rho,row[0])/row[3]-1)>1e-10) return 4;
        const double cv_fd=(R410A::internal_energy(rho,row[0]+.001)-R410A::internal_energy(rho,row[0]-.001))/.002;
        if (std::abs(cv_fd/row[3]-1)>1e-8) return 5;
    }
    for (const auto& row : references) {
        const double rho=R410A::gas_density(row[1],row[0]);
        worst=std::max(worst,std::abs(rho/row[2]-1));
        const double eps=rho*1e-5;
        const double fd=(R410A::pressure(rho+eps,row[0])-R410A::pressure(rho-eps,row[0]))/(2*eps);
        if (std::abs(fd/R410A::pressure_derivative(rho,row[0])-1)>1e-8) return 2;
    }
    for (auto invalid : {std::array<double,2>{2e6,293.15}, {-1,313.15}, {1e5,199}, {1e5,501}}) {
        bool rejected=false;
        try { R410A::gas_density(invalid[0],invalid[1]); } catch (const std::domain_error&) { rejected=true; }
        if (!rejected) return 3;
    }
    std::cout << "Maximum R410A density relative difference from CoolProp: " << worst << '\n';
    return worst<1e-10 ? 0 : 1;
}
