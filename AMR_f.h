#pragma once
#include "Header.h"
#include <vector>

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
	// ���� ������ (��������� �� ��) �� ����������

	void Get_all_cells(std::vector<AMR_cell*>& cells); // �������� ������ �������������� ����� (������������)


	void Print_info(void);

	void Print_all_center_Tecplot(AMR_f* AMR);
	void Print_slice_Tecplot(AMR_f* AMR, const double& a, const double& b, const double& c, const double& d);
	// ���������  a x + b y + c z + d = 0

	void Print_all_sosed_Tecplot(AMR_f* AMR);
};

