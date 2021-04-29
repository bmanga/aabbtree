#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <abt/aabb_tree.hpp>

using namespace abt;
TEST_CASE("point")
{
  auto pt2d = point<2, double>{3, 5};
  auto pt2i = point<2, int>{3, 5};
  auto pt3d = point<3, double>{3, 5};
  auto pt3i = point<3, int>{3, 5};
  // Ensure remaining points are zero-initialized.
  REQUIRE(pt3d[2] == 0.0);
  REQUIRE(pt3i[2] == 0);
}

TEST_CASE_TEMPLATE("aabb 2d", T, double, float, int, short)
{
  using aabb = aabb<2, T>;
  using point = aabb::point;
  auto bb2d = aabb();
  REQUIRE(bb2d.lowerBound.x() == 0);
  REQUIRE(bb2d.upperBound.x() == 0);

  bb2d = aabb({0, 0}, {4, 4});
  REQUIRE(bb2d.lowerBound == point{0, 0});
  REQUIRE(bb2d.upperBound == point{4, 4});
  // Not sure.
  REQUIRE(bb2d.surfaceArea == 16);
  REQUIRE(bb2d.centre == point{2, 2});

  auto bb2d1 = aabb({2, 2}, {4, 4});
  auto bb2d2 = aabb({3, 3}, {5, 5});
  SUBCASE("merge")
  {
    bb2d.merge(bb2d1, bb2d2);
    REQUIRE(bb2d.lowerBound == point{2, 2});
    REQUIRE(bb2d.upperBound == point{5, 5});
    // Not sure.
    REQUIRE(bb2d.surfaceArea == 12);
    REQUIRE(bb2d.centre == point{3.5, 3.5});
  }
  SUBCASE("contains")
  {
    REQUIRE(bb2d.contains(bb2d1));
    REQUIRE(!bb2d.contains(bb2d2));
  }
  SUBCASE("overlaps")
  {
    REQUIRE(bb2d.overlaps(bb2d1, true));
    REQUIRE(bb2d.overlaps(bb2d2, true));
  }
}

TEST_CASE_TEMPLATE("tree 2d", T, double, float, int, short)
{
  using tree = tree<2, T>;
  using aabb = tree::aabb;
  using point = tree::point;

  tree t;
  t.insertParticle(1, {0, 0}, {2, 2});
  t.insertParticle(2, {1, 1}, {3, 3});
  t.insertParticle(3, {2, 2}, {4, 4});
  t.insertParticle(4, {5, 5}, {7, 7});
  auto intersections = t.get_overlaps(aabb{{1, 1}, {2, 2}});
  REQUIRE(intersections.size() == 3);
  REQUIRE(t.any_overlap(aabb{{1, 1}, {2, 2}}));
  intersections = t.get_overlaps(aabb{{1, 1}, {2, 2}}, false);
  REQUIRE(intersections.size() == 2);
  intersections = t.get_overlaps(point{1, 1});
  REQUIRE(intersections.size() == 2);
  REQUIRE(t.any_overlap(point{1, 1}));
  intersections = t.get_overlaps(point{1, 1}, false);
  REQUIRE(intersections.size() == 1);
}
