#include <iostream>
#include "header.h"

int main()
{
    AMR_f* S = new AMR_f(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 5, 5, 1);

    auto A = S->find_cell(0.0, 0.0, 0.0);
    //S->Print_info();
    A->divide(3, 3, 1);
    //S->Print_info();
    A = S->find_cell(0.0, 0.0, 0.0);
    A->divide(5, 5, 1);

    std::vector<std::array<unsigned int, 3>> numbers;
    numbers.resize(A->level + 1);
    A->Get_index(numbers);

    std::array<double, 3> center;
    A->Get_Center(S, center);
    cout << "center = " << center[0] << " " << center[1] << " " << center[2] << endl;

    Print(numbers);

    S->Print_all_center_Tecplot(S);
}
