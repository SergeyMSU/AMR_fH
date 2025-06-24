#pragma once

// Boost библиотека (надо подключать к компилятору отдельно)
#include "boost/multi_array.hpp"
//#include <boost/parameter.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <array>


#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <list>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <map>
#include <math.h>
#include <cmath>
#include <limits>
#include <iterator>
#include <cstdlib>
#include <optional>
#include <omp.h>
#include <chrono>


const double const_pi = 3.14159265358979323846;
const double cpi4 = 4.0 * const_pi;
const double cpi8 = 8.0 * const_pi;
const double spi4 = sqrt(cpi4);
const double eps = 1E-12;
const double epsb = 1E-4;
const double eps_p = 1E-6;
const double eps_d = 1E-3;

#define kv(x) ((x) * (x))
#define kvg(x) (pow(x, 2.0 * this->phys_param->gamma))
#define kyb(x) ((x) * (x) * (x))
#define kvv(x, y, z) ((x) * (x) + (y) * (y) + (z) * (z))
#define norm2(x, y, z) (sqrt(kvv(x, y, z)))
#define whach(x) cout << #x <<": " << (x) << endl

class AMR_f;
class AMR_cell;

#include "AMR_cell.h"
#include "AMR_f.h"

using namespace std;


void Print(std::vector<std::array<unsigned int, 3>>& numbers);