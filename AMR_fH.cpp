#include <iostream>
#include "header.h"

int main()
{
    AMR_f* S = new AMR_f(-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 5, 5, 5);

    auto A = S->find_cell(0.0, 0.0, 0.0);
    A->divide(3, 3, 3);
    A = S->find_cell(0.18, 0.01, 0.01);
    auto B = A->get_sosed(S, 1);


    A->Print_info();
    B->Print_info();
}
