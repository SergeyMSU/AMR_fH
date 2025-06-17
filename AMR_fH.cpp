#include <iostream>
#include "header.h"

int main()
{
    AMR_f* S = new AMR_f(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 5, 5, 1);

    auto A = S->find_cell(0.1, 0.0, 0.0);
    //S->Print_info();
    A->divide(3, 3, 1);
    //S->Print_info();
    A = S->find_cell(0.18, 0.01, 0.01);
    auto B = A->get_sosed(S, 5);

    std::vector<std::array<unsigned int, 3>> numbers;
    numbers.resize(A->level + 1);
    A->Get_index(numbers);

    Print(numbers);

}
