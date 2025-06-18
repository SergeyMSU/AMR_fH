#pragma once
#include "Header.h"

class AMR_f
{
public:
	double xL;
	double xR;
	double yL;
	double yR;
	double zL;
	double zR;
	unsigned int xn;
	unsigned int yn;
	unsigned int zn;

	boost::multi_array<AMR_cell*, 3> cells;

	AMR_f(const double& xL, const double& xR, const double& yL, const double& yR, const double& zL,
		const double& zR, unsigned int xn, unsigned int yn, unsigned int zn);


	AMR_cell* find_cell(const double& x, const double& y, const double& z);
	// »щет €чейку (указатель на неЄ) по координате

	void Get_all_cells(vector< AMR_cell*>& cells); // ѕолучить список действительных €чеек (неразделЄнных)


	void Print_info(void);

	void Print_all_center_Tecplot(AMR_f* AMR);
	void Print_all_sosed_Tecplot(AMR_f* AMR);
};

