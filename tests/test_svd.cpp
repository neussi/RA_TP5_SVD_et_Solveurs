#include <Unitest/Unitest.h>
#include <Unitest/TestMacro.h>
#include "../SVD.h"

using namespace NkMath;

TEST_CASE(TP5, SVDShape) {
    MatXd A(3, 2);
    A(0, 0) = 1; A(0, 1) = 2;
    A(1, 0) = 3; A(1, 1) = 4;
    A(2, 0) = 5; A(2, 1) = 6;

    SVDResult res = SVD(A);
    ASSERT_EQUAL(res.sigma.size(), 2u);
    ASSERT_TRUE(res.sigma[0] > res.sigma[1]);
}

TEST_CASE(TP5, SVDReconstruction) {
    MatXd A(2, 2);
    A(0, 0) = 4; A(0, 1) = 0;
    A(1, 0) = 0; A(1, 1) = 3;

    SVDResult res = SVD(A);
    ASSERT_NEAR(res.sigma[0], 4.0, 1e-9);
    ASSERT_NEAR(res.sigma[1], 3.0, 1e-9);

    //Vt is identity in this case
    ASSERT_NEAR(std::abs(res.V(0,0)), 1.0, 1e-9);
}

TEST_CASE(TP5, Transpose) {
    MatXd A(3, 2);
    A(2, 1) = 5.0;
    MatXd At = A.Transposed();
    ASSERT_EQUAL(At.rows, 2);
    ASSERT_EQUAL(At.cols, 3);
    ASSERT_EQUAL(At(1, 2), 5.0);
}
