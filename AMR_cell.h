#pragma once
#include "Header.h"

class AMR_cell
{
public:
	double f;
	AMR_cell* I_self;            // ��������� �� ����


	unsigned short int level = 0;
	AMR_cell* parent = nullptr;    // ������ - ��������
	unsigned int nx = 0;           // ����� ������ ������ � ������-��������
	unsigned int ny = 0;
	unsigned int nz = 0;

	bool is_divided = false;
	boost::multi_array<AMR_cell*, 3> cells;  // ������ - ����

	AMR_cell();

	void divide(unsigned short int n1, unsigned short int n2, unsigned short int n3); // ��������� ������

	AMR_cell* find_cell(const double& x, const double& y, const double& z, const double& xL, 
		const double& xR, const double& yL, const double& yR, const double& zL, const double& zR);
	// ���� ������ �� � �������

	AMR_cell* get_sosed(AMR_f* AMR, short int nn);
	// nn = 0, 1, 2, 3, 4, 5
	//     �� � ����� - �����, �� y ...

	void Print_info(void);

	void Get_index(std::vector<std::array<unsigned int, 3>>& numbers);
	// �������� ������ ������

	void Get_Center(AMR_f* AMR, std::array<double, 3>& center); // �������� ����� ������ (���� ���� ��� �������)
	void Get_Center(AMR_f* AMR, std::array<double, 3>& center, std::array<double, 3>& razmer); // �������� ����� ������ (���� ���� ��� �������)

	void Get_Centers(AMR_f* AMR, std::vector<std::array<double, 3>>& centers); // �������� ������ ������ (������� ������ ���������)
	void Get_all_cells(vector< AMR_cell*>& cells); // �������� ������ �������������� ����� (������������)
};

