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
	std::vector<std::array<unsigned int, 3>> numbers;
	numbers.resize(this->level + 1);
	this->Get_index(numbers);

	std::array<double, 3> center;
	std::array<double, 3> razmer;
	this->Get_Center(AMR, center, razmer);


	auto nn1 = this->cells.shape()[0];
	auto nn2 = this->cells.shape()[1];
	auto nn3 = this->cells.shape()[2];
	
	if (this->level == 0)
	{
		switch (nn) {
		case 0:
			if (this->nx < nn1 - 1)
			{
				return AMR->cells[this->nx + 1][this->ny][this->nz];
			}
			else { return nullptr; }
			break;
		case 1:
			if (this->nx > 0)
			{
				return AMR->cells[this->nx - 1][this->ny][this->nz];
			}
			else { return nullptr; }
			break;
		case 2:
			if (this->ny < nn2 - 1)
			{
				return AMR->cells[this->nx][this->ny + 1][this->nz];
			}
			else { return nullptr; }
			break;
		case 3:
			if (this->ny > 0)
			{
				return AMR->cells[this->nx][this->ny - 1][this->nz];
			}
			else { return nullptr; }
			break;
		case 4:
			if (this->nz < nn3 - 1)
			{
				return AMR->cells[this->nx][this->ny][this->nz + 1];
			}
			else { return nullptr; }
			break;
		case 5:
			if (this->nz > 0)
			{
				return AMR->cells[this->nx][this->ny][this->nz - 1];
			}
			else { return nullptr; }
			break;
		default:
			cout << "ERROR 874658767843659837459" << endl;
			exit(-1);
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
			return this->parent->get_sosed(AMR, nn); // Это не будет работать
		}
	}


	return nullptr;
}

void AMR_cell::Print_info(void)
{
	if (this->is_divided == false)
	{
		cout << "_________________________________" << endl;
		cout << "level: " << this->level << endl;
		cout << "is devided?: " << this->is_divided << endl;
		cout << "nx: " << this->nx << endl;
		cout << "ny: " << this->ny << endl;
		cout << "nz: " << this->nz << endl;
		cout << "_________________________________" << endl;
	}
	else
	{
		cout << "_________________________________" << endl;
		cout << "level: " << this->level << endl;
		cout << "is devided?: " << this->is_divided << endl;
		cout << "nx: " << this->nx << endl;
		cout << "ny: " << this->ny << endl;
		cout << "nz: " << this->nz << endl;
		cout << "_______Spusk_______" << endl;

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
}

void AMR_cell::Get_index(std::vector<std::array<unsigned int, 3>>& numbers)
{
	//numbers.resize(this->level + 1);
	numbers[this->level] = std::array<unsigned int, 3>{this->nx, this->ny, this->nz};
	if (this->level == 0) return;

	this->parent->Get_index(numbers);
}

void AMR_cell::Get_Center(AMR_f* AMR, std::array<double, 3>& center)
{
	std::vector<std::array<unsigned int, 3>> numbers;
	numbers.resize(this->level + 1);
	this->Get_index(numbers);

	center[2] = center[1] = center[0] = 0.0;

	unsigned int xn = AMR->xn;
	unsigned int yn = AMR->yn;
	unsigned int zn = AMR->zn;

	double xL = AMR->xL;
	double xR = AMR->xR;

	double yL = AMR->yL;
	double yR = AMR->yR;

	double zL = AMR->zL;
	double zR = AMR->zR;

	AMR_cell* cell = nullptr;
	short int stk = 0;

	for (auto& i : numbers)
	{
		stk++;
		double dx = (xR - xL) / xn;
		double dy = (yR - yL) / yn;
		double dz = (zR - zL) / zn;
		xR = xL + (i[0] + 1) * dx;
		xL = xL + i[0] * dx;
		yR = yL + (i[1] + 1) * dy;
		yL = yL + i[1] * dy;
		zR = zL + (i[2] + 1) * dz;
		zL = zL + i[2] * dz;
		

		if (stk == 1)
		{
			cell = AMR->cells[i[0]][i[1]][i[2]];
		}
		else
		{
			/*cout << "A " << endl;
			whach(i[0]);
			whach(i[1]);
			whach(i[2]);
			whach(xL);
			whach(xR);
			whach(dx);
			whach(xn);*/
			cell = cell->cells[i[0]][i[1]][i[2]];
		}

		if (cell->is_divided == true)
		{
			xn = cell->cells.shape()[0];
			yn = cell->cells.shape()[1];
			zn = cell->cells.shape()[2];
		}
	}

	center[0] = (xR + xL) / 2.0;
	center[1] = (yR + yL) / 2.0;
	center[2] = (zR + zL) / 2.0;

	return;
}

void AMR_cell::Get_Center(AMR_f* AMR, std::array<double, 3>& center, std::array<double, 3>& razmer)
{
	std::vector<std::array<unsigned int, 3>> numbers;
	numbers.resize(this->level + 1);
	this->Get_index(numbers);

	center[2] = center[1] = center[0] = 0.0;

	unsigned int xn = AMR->xn;
	unsigned int yn = AMR->yn;
	unsigned int zn = AMR->zn;

	double xL = AMR->xL;
	double xR = AMR->xR;

	double yL = AMR->yL;
	double yR = AMR->yR;

	double zL = AMR->zL;
	double zR = AMR->zR;

	AMR_cell* cell = nullptr;
	short int stk = 0;

	for (auto& i : numbers)
	{
		stk++;
		double dx = (xR - xL) / xn;
		double dy = (yR - yL) / yn;
		double dz = (zR - zL) / zn;
		xR = xL + (i[0] + 1) * dx;
		xL = xL + i[0] * dx;
		yR = yL + (i[1] + 1) * dy;
		yL = yL + i[1] * dy;
		zR = zL + (i[2] + 1) * dz;
		zL = zL + i[2] * dz;


		if (stk == 1)
		{
			cell = AMR->cells[i[0]][i[1]][i[2]];
		}
		else
		{
			cell = cell->cells[i[0]][i[1]][i[2]];
		}

		if (cell->is_divided == true)
		{
			xn = cell->cells.shape()[0];
			yn = cell->cells.shape()[1];
			zn = cell->cells.shape()[2];
		}
	}

	center[0] = (xR + xL) / 2.0;
	center[1] = (yR + yL) / 2.0;
	center[2] = (zR + zL) / 2.0;

	razmer[0] = xR - xL;
	razmer[1] = yR - yL;
	razmer[2] = zR - zL;

	return;
}

void AMR_cell::Get_Centers(AMR_f* AMR, std::vector<std::array<double, 3>>& centers)
{
	if (this->is_divided == false)
	{
		std::array<double, 3> center;
		this->Get_Center(AMR, center);
		centers.push_back(center);
	}
	else
	{
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
	}
}

void AMR_cell::Get_all_cells(vector<AMR_cell*>& cells)
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

