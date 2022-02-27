/*
 * LU_Decomp.hpp
 *
 *  Created on: Mar 5, 2018
 *      Author: lindorfer
 */

#ifndef LU_DECOMP_HPP_
#define LU_DECOMP_HPP_

#include "Utilities.hpp"
#include <iostream>
#include <initializer_list>
#include <list>

using namespace std;

template <typename matrix_type, typename vector_type> vector_type solve_lineq(matrix_type& A, vector_type& b){

	//----Solve system of linear equations A.x = b using boost::ublas-----

	matrix_type A_factorized = A;
	ublas::permutation_matrix<size_t> pm(b.size());

	int singular = lu_factorize(A_factorized, pm);
	if (singular) throw std::runtime_error("[solve_lineq()] A is singular.");

	vector_type result(b);
	lu_substitute(A_factorized, pm, result);

	return result;
}


template<typename T, typename Vec_Type, typename limit_type>
double Gauss_Chebychev_Quadrature(T& function, Vec_Type& nodes, Vec_Type& weights, const limit_type& lower, const limit_type& upper, const bool& init_grid = false){

	//-----Inititalize a Sample Grid of 1e6 Points if Grid is NOT Pre-Defined, i.e. init_grid != true-----
	if(init_grid == true){
		double n = 100000;
		double pi = 3.1415926535897932384626434;

		nodes.resize(n);
		weights.resize(n);

		for (long double i = 1; i <= n; i++) {
			nodes(i - 1) = cos((2.0 * i - 1.0) / (2.0 * n) * pi);
			weights(i - 1) = pi / n;
		}
	}

	//-----Do the Gauss-Chebychev Integration-----
	double sum = 0.0;
	int n = nodes.size();

	//-----Transform the Upper / Lower Limits to [-1, 1]-----
	auto f = [&function](const double& xpara, const double& a, const double& b) {
		double x = 0;
		x = (b - a) / 2.0 * xpara + (b + a) / 2.0;
		return ((b - a) / 2.0 * function(x));
	};

	//-----Perform the Gauss-Chebychev-Quadrature Summation-----
	for (int i = 0; i < n; i++) {
		sum += weights(i) * f(nodes(i), lower, upper) * sqrt(1 - nodes(i) * nodes(i));
	}

	return sum;
}



//void foo(const std::list<std::string> & myArguments) {
//
//	for(auto it = myArguments.begin(); it != myArguments.end(); it++){
//		cout << *it;
//	}
//	cout << endl << "Number of Elements in List: " << myArguments.size();
//}

//template<list list_type> void FillVec(list<list_type> mylist){
//
//	for (auto it = mylist.begin(); it != mylist.end(); it++) {
//		cout << *it;
//	}
//	cout << endl << "Number of Elements in List: " << mylist.size();
//}


#endif /* LU_DECOMP_HPP_ */
