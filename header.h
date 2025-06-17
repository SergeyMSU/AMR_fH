#pragma once

// Boost библиотека (надо подключать к компилятору отдельно)
#include "boost/multi_array.hpp"
#include <boost/parameter.hpp>
#include <iostream>
#include <vector>
#include <array>

class AMR_f;
class AMR_cell;

#include "AMR_cell.h"
#include "AMR_f.h"

using namespace std;


void Print(std::vector<std::array<unsigned int, 3>>& numbers);