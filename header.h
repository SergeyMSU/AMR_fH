#pragma once

// Boost библиотека (надо подключать к компилятору отдельно)
#include "boost/multi_array.hpp"
#include <boost/parameter.hpp>
#include <iostream>
#include <vector>
#include <array>

#include <string>
#include <sstream>
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

#define whach(x) cout << #x <<": " << (x) << endl

class AMR_f;
class AMR_cell;

#include "AMR_cell.h"
#include "AMR_f.h"

using namespace std;


void Print(std::vector<std::array<unsigned int, 3>>& numbers);