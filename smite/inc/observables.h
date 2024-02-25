#ifndef OBSERV_H_
#define OBSERV_H_

#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime> 
#include <cstdlib> // Need it for rand() and etc
#include <functional> // Need it for function as argument
#include <numeric>

#include "cross_sections.h"
#include "integrator.h"

std::vector<long double> S_to_LAB(const std::vector<long double> &x, const std::vector<long double> &params);
std::vector<long double> S_to_LAB_params(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x);

std::vector<std::vector<long double>> P_Sframe(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x);
std::vector<std::vector<long double>> P_LAB(const long double mjarg, const long double mkarg, const long double Mnarg, const std::vector<long double> &x);

#endif