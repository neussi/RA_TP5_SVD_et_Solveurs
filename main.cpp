#include <iostream>
#include <cstdio>
#include "SVD.h"

int main() {
    using namespace NkMath;
    printf("=== TP5: SVD (Singular Value Decomposition) ===\n\n");

    MatXd A(3, 2);
    A(0, 0) = 1; A(0, 1) = 2;
    A(1, 0) = 3; A(1, 1) = 4;
    A(2, 0) = 5; A(2, 1) = 6;

    printf("Matrice A (3x2) :\n  1 2\n  3 4\n  5 6\n\n");

    SVDResult res = SVD(A);

    printf("Valeurs singulieres (sigma) :\n");
    for (double s : res.sigma) printf("  %.4f", s);
    printf("\n\n");

    printf("  SVD calculee avec succes via Jacobi one-sided.\n");
    return 0;
}
