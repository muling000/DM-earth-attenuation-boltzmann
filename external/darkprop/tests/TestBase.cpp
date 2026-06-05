#include "Test.hpp"

TEST_CASE("rotation_axis_angle") {
    SUBCASE("no rotation") {
        Vector3 n {0, 0, 1};
        Vector3 v {0, 1, 0};
        Vector3 vnew = v;
        rotation_axis_angle(vnew, n, 0);
        CHECK_EQ(v, vnew);
    }

    SUBCASE("rotation along axis") {
        Vector3 n {1.0 / sqrt2, 1.0 / sqrt2, 0};
        Vector3 v {1, 1, 0};
        Vector3 vnew = v;
        rotation_axis_angle(vnew, n, 2*M_PI);
        CHECK(vnew.isApprox(v, PREC));
    }

    SUBCASE("rotate back") {
        Vector3 n {1, 1, 1};
        n.normalize();
        Vector3 v {1, 0, 0};
        Vector3 vnew = v;
        rotation_axis_angle(vnew, n, M_PI);
        rotation_axis_angle(vnew, n, -M_PI);
        CHECK(vnew.isApprox(v, PREC));
    }

    SUBCASE("roate 2pi") {
        Vector3 n {1, 2, 3};
        n.normalize();
        Vector3 v {1, 0, 0};
        Vector3 vnew = v;
        rotation_axis_angle(vnew, n, 2*M_PI);
        CHECK(vnew.isApprox(v, PREC));
    }

    SUBCASE("roate pi/2") {
        Vector3 ex {1, 0, 0};
        Vector3 ey {0, 1, 0};
        Vector3 ez {0, 0, 1};
        Vector3 xnew = ex;
        Vector3 ynew = ey;
        Vector3 znew = ez;
        rotation_axis_angle(xnew, ez, M_PI/2);
        CHECK(xnew.isApprox(ey, PREC));
        rotation_axis_angle(ynew, ex, M_PI/2);
        CHECK(ynew.isApprox(ez, PREC));
        rotation_axis_angle(znew, ey, M_PI/2);
        CHECK(znew.isApprox(ex, PREC));
    }
}
