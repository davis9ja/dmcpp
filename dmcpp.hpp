#ifndef DMC_HPP_
#define DMC_HPP_

#include <omp.h>
#include <math.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <algorithm>
#include <fstream>
#include <iostream>

void creator(boost::numeric::ublas::vector<int> &phi, int p);
void annihilator(boost::numeric::ublas::vector<int> &phi, int p);
int inner_product(boost::numeric::ublas::vector<int> phi_bra, 
                       boost::numeric::ublas::vector<int> phi_ket);
double trace_1b(boost::numeric::ublas::vector<double> data);
double trace_2b(boost::numeric::ublas::vector<double> data);

std::string state2occupation(std::string state);
boost::numeric::ublas::matrix<int> gen_basis(int nholes, int nparticles);
boost::numeric::ublas::matrix<int> readBasisFromFile(std::string file_path);

boost::numeric::ublas::vector<double> density_1b(int nholes, int nparticles, boost::numeric::ublas::vector<double> weights, std::string basis_path);
boost::numeric::ublas::vector<double> density_2b(int nholes, int nparticles, boost::numeric::ublas::vector<double> weights, std::string basis_path);
boost::numeric::ublas::vector<double> density_3b(int nholes, int nparticles, boost::numeric::ublas::vector<double> weights, std::string basis_path);

inline bool fexists (const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    } 
}

#endif
