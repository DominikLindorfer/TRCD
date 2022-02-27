//============================================================================
// Name        : Utilities.hpp
// Author      : lindorfer
// Version     :
// Copyright   : Your copyright notice
// Description : Simplistic C++ Utilities to make life easier
//============================================================================

#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include <bits/stdc++.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

//-----This stuff is for LU Decomposition-----
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/assign.hpp>
#include <boost/numeric/ublas/operations.hpp>
#include <boost/numeric/ublas/assignment.hpp>

template <typename T>
using Vektor = boost::numeric::ublas::vector<T>;

template <typename T>
using Matrix = boost::numeric::ublas::matrix<T>;

namespace ublas = boost::numeric::ublas;

//-----Push_Back Routine for UBLAS Vektor-----
template <typename T, typename A>
void push_back(Vektor<T>& vec, A item) {

	vec.resize((int)vec.size() + 1);
	vec(vec.size() - 1) = item;
}

//-----Initialize Ublas-Vektor by initializer list-----
template<typename T, typename A>
void vec_append(Vektor<A>& v, std::initializer_list<T> l) {

	for(auto it = l.begin(); it != l.end(); ++it) {
	    push_back(v, *it);
	}
}

//-----Initialize Ublas-Matrix by initializer lists-----
template<typename T, typename A>
void mat_append(Matrix<A>& mat, std::initializer_list< std::initializer_list<T> > st) {

	//-----Resize Matrix to initializer_list sizes & fill matrix after-----
	mat.resize((int)st.size(), (int)st.begin()->size());


	int i = 0;
	int j = 0;

	for (const auto& l : st) {
		j = 0;
		for (const auto& v : l) {

			mat(i, j) = v;
			j++;
		}
		i++;
	}
}
using namespace std;

//-----Range-----
//-----The data-type of the step is used for the data-type of the vector-----
template <typename start_type, typename stop_type, typename step_type>
vector<step_type> range(start_type start, stop_type stop, step_type step){
  if (step == step_type(0)){
    throw std::invalid_argument("step for range must be non-zero");
  }

  if(typeid(step) != typeid(start)){

  }

  vector<step_type> result;
  step_type i = start;

  while ((step > 0) ? (i < stop) : (i > stop)){
    result.push_back(i);
    i += step;
  }

  return result;
}

template <typename T>
vector<T> range(T start, T stop){
  return range(start, stop, T(1));
}

template <typename T>
vector<T> range(T stop){
  return range(T(0), stop, T(1));
}

//-----std::vector cout-----
template<typename T>
ostream& operator<< (ostream& out, const vector<T>& v) {
    out << "{";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
        out << v[i];
        if (i != last)
            out << ", ";
    }
    out << "}";
    return out;
}

//-----std::tuple cout-----
template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I == sizeof...(Tp), string>::type
stringval(const std::tuple<Tp...> & t)
{
  std::stringstream buffer;
  buffer << "}";
  return buffer.str();
}

template<std::size_t I = 0, typename... Tp>
inline typename std::enable_if<I < sizeof...(Tp), string>::type
stringval(const std::tuple<Tp...> & t)
{
  std::stringstream buffer;
  size_t len = sizeof...(Tp);
  if(I==0)
      buffer << "{";
  buffer << std::get<I>(t);
  if(I < len - 1)
    buffer << ", ";
  buffer << stringval<I + 1, Tp...>(t);
  return buffer.str();
}

template<typename... Tp> ostream& operator <<(ostream& out, const std::tuple<Tp...> & t)
{
  out << stringval(t);
  return out;
}

//template <typename t, template<typename> typename vec> void move_element(vec<t>& v, size_t oldIndex, size_t newIndex){
template <typename vec> void move_element(vec& v, size_t oldIndex, size_t newIndex){
    if (oldIndex > newIndex)
        std::rotate(v.rend() - oldIndex - 1, v.rend() - oldIndex, v.rend() - newIndex);
    else
        std::rotate(v.begin() + oldIndex, v.begin() + oldIndex + 1, v.begin() + newIndex + 1);
}

template <typename vec> void erase_element(vec& v, size_t idx){

	assert(idx < v.size());
	for (uint i = idx; i < v.size() - 1; i++) {
	        v(i) = v(i + 1);
	    }
	    v.resize(v.size() - 1);
}


#endif /* UTILITIES_HPP_ */
