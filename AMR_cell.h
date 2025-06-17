#pragma once
#include "Header.h"

class AMR_cell
{
public:
	double f;
	AMR_cell* I_self;            // ”казатель на себ€


	unsigned short int level = 0;
	AMR_cell* parent = nullptr;    // ячейка - родитель
	unsigned int nx = 0;           // Ќомер данной €чейки в €чейке-родителе
	unsigned int ny = 0;
	unsigned int nz = 0;

	bool is_divided = false;
	boost::multi_array<AMR_cell*, 3> cells;  // ячейки - дети

	AMR_cell();

	void divide(unsigned short int n1, unsigned short int n2, unsigned short int n3); // –азделить €чейку

	AMR_cell* find_cell(const double& x, const double& y, const double& z, const double& xL, 
		const double& xR, const double& yL, const double& yR, const double& zL, const double& zR);
	// »щет €чейку по еЄ сосед€м

	AMR_cell* get_sosed(AMR_f* AMR, short int nn);
	// nn = 0, 1, 2, 3, 4, 5
	//     по х вперЄд - назад, по y ...

	void Print_info(void);

	void Get_index(std::vector<std::array<unsigned int, 3>>& numbers);
	// ѕолучить индекс €чейки

	void Get_Center(AMR_f* AMR, std::array<double, 3>& center);

};

