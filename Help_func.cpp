#include "header.h"


void Print(std::vector<std::array<unsigned int, 3>>& numbers)
{
	cout << "--------------" << endl;
	for (auto& i : numbers)
	{
		cout << i[0] << " " << i[1] << " " << i[2] << endl;
	}
	cout << "--------------" << endl;
}