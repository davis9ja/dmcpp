//////////////////////////////////////////////////////////////////////////
// The purpose of this program is to compute 1b/2b/3b density matrices, //
// where an eigenstate is expanded in a basis of Slater determinant     //
// configurations in an occupation number representation.               //
//                                                                      //
//     Author: Jacob Davison                                            //
//     Date: 08/05/2021                                                 //
//////////////////////////////////////////////////////////////////////////
#include "dmcpp.hpp"

using namespace boost::numeric::ublas;
//using namespace std;

void creator(vector<int> &phi, int p) {

    if ( phi[0] != -2 ) {
        

        int idx = p+1;
        int phase_exp = 0;
        int phase;
    
        for (int i = 1; i < idx; i++)
            phase_exp += phi[i];

        phase = pow(-1, phase_exp);
        
        if (phi[idx] == 0) {
            phi[idx] = 1;
            phi[0] = phi[0]*phase;
        } else {
            phi.resize(1, false);
            phi[0] = -2;
        }
    }
}

void annihilator(vector<int> &phi, int p) {

    if ( phi[0] != -2 ) {
     
        int idx = p+1;
        int phase_exp = 0;
        int phase;
    
        for (int i = 1; i < idx; i++)
            phase_exp += phi[i];

        phase = pow(-1, phase_exp);

        if (phi[idx] == 1) {
            phi[idx] = 0;
            phi[0] = phi[0]*phase;
        } else {
            phi.resize(1, false);
            phi[0] = -2;
        }
    }
    
}

int inner_product(vector<int> phi_bra, vector<int> phi_ket) {

    if ( (phi_bra[0] == -2) || (phi_ket[0] == -2) )
        return 0;

    vector_range<vector<int>> bra (phi_bra, range(1,phi_bra.size()));
    vector_range<vector<int>> ket (phi_ket, range(1,phi_ket.size()));
    
    bool match = std::equal(bra.begin(), bra.end(), ket.begin());
    
    if (match)
        return 1*phi_bra[0]*phi_ket[0];
    else
        return 0;
}

/*
  Custom sort function.
 */
struct {
    bool operator()(std::string s1, std::string s2) const { 
        int half = (int)s1.length()/2;
        int pairs1 = 0;
        int pairs2 = 0;

        std::string up1 = s1.substr(0,half);
        std::string dn1 = s1.substr(half,s1.length());
        for (int i = 0; i < up1.length(); i++) 
            if ( up1[i] == dn1[i] == 1 )
                pairs1 += 1;

        std::string up2 = s2.substr(0,half);
        std::string dn2 = s2.substr(half,s1.length());
        for (int i = 0; i < up2.length(); i++) 
            if ( up2[i] == dn2[i] == 1 )
                pairs2 += 1;
        
        return pairs1 < pairs2;
    }
} pair_comp;

std::string state2occupation(std::string state) {

    int half = (int)state.length()/2;
    std::string occupation = "";
    std::string up1 = state.substr(0,half);
    std::string dn1 = state.substr(half,state.length());

    for (int i = 0; i < up1.length(); i++)  {
        occupation += up1[i];
        occupation += dn1[i];
    }

    return occupation;
}


/*
  Generate only Sz = 0 block. Run through permuations of
  occupation string. Ordering ends up not intuitive; not sure
  how to fix this yet.
 */
matrix<int> gen_basis(int nholes, int nparticles) {

    int numStates = nholes+nparticles;
    std::string state;
    int total = 0;
    std::vector<std::string> permutations;
    std::vector<std::string> states_str_vec;

    for (int i = 0; i < (int)nholes/2; i++)
        state += "1";
    for (int i = 0; i < (int)nparticles/2; i++)
        state += "0";

    
    std::sort(state.begin(), state.end());
    do {
        total += 1;
        permutations.push_back(state);
        
    } while(std::next_permutation(state.begin(), state.end()));

    for (int i = 0; i < permutations.size(); i++)
        for (int j = 0; j < permutations.size(); j++)
            states_str_vec.push_back(permutations[i]+permutations[j]);
    std::sort(states_str_vec.begin(), states_str_vec.end(), pair_comp);

    std::string occupation_str;
    matrix<int> basis(states_str_vec.size(), numStates+1);
    for (int i = 0; i < states_str_vec.size(); i++) {
        occupation_str = state2occupation(states_str_vec[states_str_vec.size()-1-i]);
        occupation_str = "1" + occupation_str;

        for (int j = 0; j < occupation_str.length(); j++) {
            const char& s = occupation_str[j];
            basis(i,j) = (s-'0');
        }
    }

    return basis;
}

/*
  Read SD basis from a file. The first line is assumed to be
      <# s.p. states> <# states>
  which gives the columns and rows for the matrix<int> that will
  be returned as the basis.
 */
matrix<int> readBasisFromFile(std::string file_path) {
    std::string line;
    std::ifstream myfile (file_path);
    
    int num_sp, num_states, count;

    if (myfile.is_open()) {
        getline(myfile,line);
        size_t pos = 0;
        std::string delim = " ";
        std::string token;
        

        pos = line.find(delim);
        token = line.substr(0, pos);
        line.erase(0, pos + delim.length());
        
        num_sp = std::stoi(token);
        num_states = std::stoi(line);
    }

    std::cout << num_sp << " " << num_states << std::endl;

    matrix<int> M(num_states, num_sp+1);
    count = 0;
    if (myfile.is_open()) {            
            while ( getline (myfile,line) )
                {
                    for (int i = 0; i < num_sp; i++) {
                        const char& s = line[i];
                        M(count,i) = (s-'0');
                    }
                    count++;
                        
                }
            myfile.close();
        }

    else std::cout << "Unable to open file"; 

    return M;
}


/*
  Build density matrix from gen_basis() function. Ordering is not logical.
 */
vector<double> density_1b(int nholes, int nparticles, vector<double> weights, std::string basis_path = "") {
    matrix<int> basis;
    std::string def = "";
    int numSP = nholes+nparticles;
    
    if (basis_path == def) {
        basis = gen_basis(nholes, nparticles);
    }
    else if (fexists(basis_path)) {
        basis = readBasisFromFile(basis_path);
    }
    else {
        std::cout << "ARGS MALFORMED OR FILE DOES NOT EXIST\n";
        exit(1);
    }

    std::cout << basis << std::endl;

    int numBasisStates = basis.size1();

    vector<double> rho1b(numSP*numSP);
    for (int i = 0; i < rho1b.size(); i++)
        rho1b(i) = 0.0;

    vector<int> state_bra, state_ket;
    double coeff_bra, coeff_ket;

    int p,q,i,j;

    for (p = 0; p < numSP; p++) {
        for (q = 0; q < numSP; q++) {
            for (i = 0; i < numBasisStates; i++) {
                state_ket = matrix_row<matrix<int>>(basis, i);
                coeff_ket = weights[i];
                annihilator(state_ket, q);
                creator(state_ket, p);
                //std::cout << state_ket << std::endl;
                
                for (j = 0; j < numBasisStates; j++) {
                    state_bra = matrix_row<matrix<int>>(basis, j);
                    coeff_bra = weights[j];

                    rho1b[p*numSP + q] += coeff_bra*coeff_ket*inner_product(state_bra, state_ket);

                    // if ( !(std::abs(rho1b[p*numSP + q]) < 1e-16) )
                    //     std::cout << rho1b[p*numSP + q] << state_bra << " " << state_ket << std::endl;
                }
            }
        }
    }
    
    return rho1b;

}


/*
  Build density matrix from gen_basis(). Ordering is not logical.
 */
vector<double> density_2b(int nholes, int nparticles, vector<double> weights, std::string basis_path = "") {

    matrix<int> basis;
    std::string def = "";
    int numSP = nholes+nparticles;
    
    if (basis_path == def) {
        basis = gen_basis(nholes, nparticles);
    }
    else if (fexists(basis_path)) {
        basis = readBasisFromFile(basis_path);
    }
    else {
        std::cout << "ARGS MALFORMED OR FILE DOES NOT EXIST\n";
        exit(1);
    }

    int numBasisStates = basis.size1();

    vector<double> rho2b(numSP*numSP*numSP*numSP);
    for (int i = 0; i < rho2b.size(); i++)
        rho2b[i] = 0.0;

    vector<int> state_bra, state_ket;
    double coeff_bra, coeff_ket, result;
    
    int p,q,r,s,i,j;

    for (q = 0; q < numSP; q++) {
        for (p = 0; p < q; p++) {
            for (s = 0; s < numSP; s++) {
                for (r = 0; r < s; r++) {
                    for (i = 0; i < numBasisStates; i++) {
                        state_ket = matrix_row<matrix<int>>(basis, i);
                        coeff_ket = weights[i];
                        
                        annihilator(state_ket, r);
                        annihilator(state_ket, s);
                        creator(state_ket, q);
                        creator(state_ket, p);

                
                        for (j = 0; j < numBasisStates; j++) {
                            state_bra = matrix_row<matrix<int>>(basis, j);
                            coeff_bra = weights[j];
                            result = coeff_bra*coeff_ket*inner_product(state_bra, state_ket);

                            rho2b[p*numSP*numSP*numSP + q*numSP*numSP + r*numSP + s] += result;
                            rho2b[q*numSP*numSP*numSP + p*numSP*numSP + r*numSP + s] = -result;
                            rho2b[p*numSP*numSP*numSP + q*numSP*numSP + s*numSP + r] = -result;
                            rho2b[q*numSP*numSP*numSP + p*numSP*numSP + s*numSP + r] = result;

                            // if ( !(std::abs(result) < 1e-16) )
                            //     std::cout << result << state_bra << " " << state_ket << std::endl;
                        }
                    }
                }
            }
        }
    }
    
    return rho2b;

}


double trace_1b(vector<double> data) {
    double sum = 0.;
    
    for (int i = 0; i < 8; i++)
            sum += data[i*8+i];
    return sum;
}

double trace_2b(vector<double> data) {
    double sum = 0.;

    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 8; j++)
            sum += data[i*8*8*8 + j*8*8 + i*8 + j];
}

// int main(int argc, char** argv) {
//     int nholes = std::atoi(argv[1]);
//     int nparticles = std::atoi(argv[2]);
//     std::string basis_path = "";

//     if (argc > 3)
//         basis_path = argv[3];
                  
//     vector<double> weights(36);
//     for(int i = 0; i < 36; i++) {
//         weights[i] = 0.0;
//         if (i == 0)
//             weights[i] = sqrt(0.8);
//         if (i == 1)
//             weights[i] = sqrt(0.2);
//     }
        
//     vector<double> rho1b = density_1b(nholes, nparticles, weights, basis_path);
//     vector<double> rho2b = density_2b(nholes, nparticles, weights, basis_path);
//     std::cout << rho1b << std::endl;
//     int idx;
//     for (int i = 0; i < 8; i++) {
//         for (int j = 0; j < 8; j++) {
//             idx = i*8 + j;
//             printf("%0.2f  ", rho1b[idx]);
//         }
//         std::cout << "\n";
//     }

        
    
//     double trace1b = trace_1b(rho1b);
//     double trace2b = trace_2b(rho2b);
//     std::cout << trace2b << std::endl;
// }
