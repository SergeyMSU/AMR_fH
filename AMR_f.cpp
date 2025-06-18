#include "AMR_f.h"
using namespace std;

AMR_f::AMR_f(const double& xL, const double& xR, const double& yL, const double& yR, const double& zL,
	const double& zR, unsigned int xn, unsigned int yn, unsigned int zn)
{
	this->xL = xL;
	this->xR = xR;

	this->yL = yL;
	this->yR = yR;

	this->zL = zL;
	this->zR = zR;

	this->xn = xn;
	this->yn = yn;
	this->zn = zn;


	this->cells.resize(boost::extents[xn][yn][zn]);

	for (int i = 0; i < xn; ++i) {
		for (int j = 0; j < yn; ++j) {
			for (int k = 0; k < zn; ++k) {
				auto A = new AMR_cell();
				A->nx = i;
				A->ny = j;
				A->nz = k;
				A->parent = nullptr;
				A->I_self = A;
				A->level = 0;
				this->cells[i][j][k] = A;
			}
		}
	}
}

AMR_cell* AMR_f::find_cell(const double& x, const double& y, const double& z)
{
	if(x < this->xL || x > this->xR) return nullptr;
	if(y < this->yL || y > this->yR) return nullptr;
	if(z < this->zL || z > this->zR) return nullptr;

	double dx = (this->xR - this->xL) / this->xn;
	int index1 = static_cast<int>((x - this->xL) / dx);
	if (index1 == this->xn) index1 = this->xn - 1;

	double dy = (this->yR - this->yL) / this->yn;
	int index2 = static_cast<int>((y - this->yL) / dy);
	if (index2 == this->yn) index2 = this->yn - 1;

	double dz = (this->zR - this->zL) / this->zn;
	int index3 = static_cast<int>((z - this->zL) / dz);
	if (index3 == this->zn) index3 = this->zn - 1;

	auto A = this->cells[index1][index2][index3];

	if (A->is_divided == false)
	{
		return A;
	}
	else
	{
		return A->find_cell(x, y, z,
			xL + index1 * dx, xL + (index1 + 1) * dx,
			yL + index2 * dy, yL + (index2 + 1) * dy,
			zL + index3 * dz, zL + (index3 + 1) * dz);
	}
}

void AMR_f::Get_all_cells(vector<AMR_cell*>& cells)
{
	const auto& shape = this->cells.shape();
	const size_t nx = shape[0];
	const size_t ny = shape[1];
	const size_t nz = shape[2];
	for (size_t i = 0; i < nx; ++i)
	{
		for (size_t j = 0; j < ny; ++j)
		{
			for (size_t k = 0; k < nz; ++k)
			{
				AMR_cell* cell = this->cells[i][j][k];
				if (cell->is_divided == false) {
					cells.push_back(cell);
				}
				else
				{
					cell->Get_all_cells(cells);
				}
			}
		}
	}
}

void AMR_f::Print_info(void)
{
	const auto& shape = this->cells.shape();
	const size_t nx = shape[0];
	const size_t ny = shape[1];
	const size_t nz = shape[2];
	for (size_t i = 0; i < nx; ++i) 
	{
		for (size_t j = 0; j < ny; ++j) 
		{
			for (size_t k = 0; k < nz; ++k) 
			{
				AMR_cell* cell = cells[i][j][k];
				if (cell != nullptr) {
					cell->Print_info();
				}
			}
		}
	}
}

void AMR_f::Print_all_center_Tecplot(AMR_f* AMR)
{
	ofstream fout;
	string name_f = "Tecplot_print_cell_center_3D.txt";
	std::vector<std::array<double, 3>> centers;

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	const size_t dim1 = this->cells.shape()[0];
	const size_t dim2 = this->cells.shape()[1];
	const size_t dim3 = this->cells.shape()[2];

	for (size_t i = 0; i < dim1; ++i) {
		for (size_t j = 0; j < dim2; ++j) {
			for (size_t k = 0; k < dim3; ++k) {
				AMR_cell* cell = cells[i][j][k];
				cell->Get_Centers(AMR, centers);
			}
		}
	}

	for (auto& i : centers)
	{
		fout << i[0] << " " << i[1] << " " <<  i[2] << endl;
	}

	fout.close();
}

void AMR_f::Print_slice_Tecplot(AMR_f* AMR, const double& a, const double& b, const double& c, const double& d)
{
	std::vector<AMR_cell*> all_cells;
	this->Get_all_cells(all_cells);

	for (auto& i : all_cells)
	{

	}
}

void AMR_f::Print_all_sosed_Tecplot(AMR_f* AMR)
{
	ofstream fout;
	string name_f = "Tecplot_print_sosed_3D.txt";

	fout.open(name_f);
	fout << "TITLE = HP  VARIABLES = X, Y, Z" << endl;

	vector< AMR_cell*> all_cells;
	this->Get_all_cells(all_cells);

	fout << "ZONE T=HP, N = " << all_cells.size() * 6 * 2 << ", E = " << all_cells.size() * 6 << ", F=FEPOINT, ET=LINESEG" << endl;

	for (auto& cc : all_cells)
	{
		std::array<double, 3> center;
		cc->Get_Center(AMR, center);
		for (size_t ii = 0; ii < 6; ii++)
		{
			auto ss = cc->get_sosed(AMR, ii);
			if (ss == nullptr) ss = cc;

			std::array<double, 3> center2;
			ss->Get_Center(AMR, center2);

			fout << center[0] << " " << center[1] << " " << center[2] << endl;
			fout << (center[0] + center2[0]) / 2.0 << " " << (center[1] + center2[1]) / 2.0
				<< " " << (center[2] + center2[2]) / 2.0 << endl;
		}
	}

	for (int m = 0; m < all_cells.size() * 6; m++)
	{
		fout << 2 * m + 1 << " " << 2 * m + 2 << endl;
	}


	fout.close();
}

