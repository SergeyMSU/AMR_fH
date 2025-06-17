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

