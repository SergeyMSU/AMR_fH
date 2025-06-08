#include "AMR_cell.h"

AMR_cell::AMR_cell()
{
	this->f = 0.0;
}

void AMR_cell::divide(unsigned short int n1, unsigned short int n2, unsigned short int n3)
{
	this->is_divided = true;

	this->cells.resize(boost::extents[n1][n2][n3]);
	for (int i = 0; i < n1; ++i) {
		for (int j = 0; j < n2; ++j) {
			for (int k = 0; k < n3; ++k)
			{
				auto A = new AMR_cell();
				A->nx = i;
				A->ny = j;
				A->nz = k;
				A->parent = this->I_self;
				A->level = this->level + 1;
				this->cells[i][j][k] = A;
				A->I_self = A;

				A->f = this->f;                        //  Просто сносим значение
			}
		}
	}

}

AMR_cell* AMR_cell::find_cell(const double& x, const double& y, const double& z, const double& xL,
	const double& xR, const double& yL, const double& yR, const double& zL, const double& zR)
{
	unsigned int xn = this->cells.shape()[0];
	unsigned int yn = this->cells.shape()[1];
	unsigned int zn = this->cells.shape()[2];

	double dx = (xR - xL) / xn;
	int index1 = static_cast<int>((x - xL) / dx);
	if (index1 == xn) index1 = xn - 1;

	double dy = (yR - yL) / yn;
	int index2 = static_cast<int>((y - yL) / dy);
	if (index2 == yn) index2 = yn - 1;

	double dz = (zR - zL) / zn;
	int index3 = static_cast<int>((z - zL) / dz);
	if (index3 == zn) index3 = zn - 1;

	auto A = this->cells[index1][index2][index3];

	if (A->is_divided == false)
	{
		return A;
	}
	else
	{
		return A->find_cell(x, y, z,
			xL + index1 * dx, xR + (index1 + 1) * dx,
			yL + index2 * dy, yR + (index2 + 1) * dy,
			zL + index3 * dz, zR + (index3 + 1) * dz);
	}
}

AMR_cell* AMR_cell::get_sosed(AMR_f* AMR, short int nn)
{
	// Нужно проверить, вдруг это основная ячейка, у которой нет родителя
	
	if (this->level == 0)
	{
		if (this->nx > 0)
		{
			return AMR->cells[this->nx - 1][this->ny][this->nz];
		}
		else
		{
			return nullptr;
		}
	}


	// Ищем соседа слева
	if (nn == 0)
	{

	}
	else if (nn == 1)
	{
		if (this->nx > 0)
		{
			return this->parent->cells[this->nx - 1][this->ny][this->nz];
		}
		else
		{
			return this->parent->get_sosed(AMR, nn);
		}
	}


	return nullptr;
}

void AMR_cell::Print_info(void)
{
	cout << "level: " << this->level << endl;
	cout << "nx: " << this->nx << endl;
	cout << "ny: " << this->ny << endl;
	cout << "nz: " << this->nz << endl;
	cout << "_________________________________" <<  endl;
}
