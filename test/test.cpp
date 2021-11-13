#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include <abt/aabb_tree.hpp>

using namespace abt;
TEST_CASE("point")
{
  auto pt2d = point<2, double>{3, 5};
  auto pt2i = point<2, int>{3, 5};
  auto pt3d = point<3, double>{3, 5, 7};
  auto pt3i = point<3, int>{3, 5, 7};
  REQUIRE(pt3d[2] == 7);
  REQUIRE(pt3i[2] == 7);
}

TEST_CASE_TEMPLATE("aabb 2d", T, double, float, int, short)
{
  using aabb = aabb<2, T>;
  using point = aabb::point;
  auto bb2d = aabb();
  REQUIRE(bb2d.lowerBound[0] == 0);
  REQUIRE(bb2d.upperBound[0] == 0);
  REQUIRE(bb2d == aabb({0, 0}, {0, 0}));
  REQUIRE(bb2d != aabb({1, 0}, {0, 0}));

  bb2d = aabb({0, 0}, {4, 4});
  REQUIRE(all_of(bb2d.lowerBound == point{0, 0}));
  REQUIRE(all_of(bb2d.upperBound == point{4, 4}));
  // Not sure.
  REQUIRE(bb2d.surfaceArea == 32);
  REQUIRE(all_of(abt::compute_center(bb2d) == point{2, 2}));

  auto bb2d1 = aabb({2, 2}, {4, 4});
  auto bb2d2 = aabb({3, 3}, {5, 5});
  SUBCASE("merge")
  {
    bb2d.merge(bb2d1, bb2d2);
    REQUIRE(all_of(bb2d.lowerBound == point{2, 2}));
    REQUIRE(all_of(bb2d.upperBound == point{5, 5}));
    // Not sure.
    REQUIRE(bb2d.surfaceArea == 18);
    REQUIRE(all_of(abt::compute_center(bb2d) == point{(T)3.5, (T)3.5}));
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

TEST_CASE_TEMPLATE("aabb 3d", T, double, float, int, short)
{
  using aabb = aabb<3, T>;
  using point = aabb::point;
  auto bb3d = aabb();
  REQUIRE(bb3d.lowerBound[0] == 0);
  REQUIRE(bb3d.upperBound[0] == 0);
  REQUIRE(bb3d == aabb({0, 0, 0}, {0, 0, 0}));
  REQUIRE(bb3d != aabb({1, 0, 0}, {0, 0, 0}));

  bb3d = aabb({0, 0, 0}, {4, 4, 4});
  REQUIRE(all_of(bb3d.lowerBound == point{0, 0, 0}));
  REQUIRE(all_of(bb3d.upperBound == point{4, 4, 4}));
  // Not sure.
  REQUIRE(bb3d.surfaceArea == 48);
  REQUIRE(all_of(abt::compute_center(bb3d) == point{2, 2, 2}));

  /*
  auto bb2d1 = aabb({2, 2}, {4, 4});
  auto bb2d2 = aabb({3, 3}, {5, 5});
  SUBCASE("merge")
  {
    bb2d.merge(bb2d1, bb2d2);
    REQUIRE(all_of(bb2d.lowerBound == point{2, 2}));
    REQUIRE(all_of(bb2d.upperBound == point{5, 5}));
    // Not sure.
    REQUIRE(bb2d.surfaceArea == 9);
    REQUIRE(all_of(abt::compute_center(bb2d) == point{(T)3.5, (T)3.5}));
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
  */
}

TEST_CASE_TEMPLATE("tree 2d", T, double, float, int32_t, uint32_t, int16_t, uint16_t)
{
  using tree = tree<2, T>;
  using aabb = tree::aabb;
  using point = tree::point;

  tree t;
  t.insert({{0, 0}, {2, 2}});
  t.insert({{1, 1}, {3, 3}});
  t.insert({{2, 2}, {4, 4}});
  t.insert({{5, 5}, {7, 7}});
  REQUIRE(t.size() == 4);
  REQUIRE(t.get_height() == 2);
  t.rebuild();
  REQUIRE(t.get_height() == 3);
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

  REQUIRE(t.any_overlap(point{1, 1}, [] { return true; }));
  REQUIRE(t.any_overlap(point{1, 1}, [] (tree::node_id) { return true; }));
  REQUIRE(t.any_overlap(point{1, 1}, [] (aabb) { return true; }));
  REQUIRE(t.any_overlap(point{1, 1}, [](tree::node_id, aabb) { return true; }));

  t.clear();
  REQUIRE(t.size() == 0);
  t.insert({{2, 2}, {4, 4}});
  REQUIRE(t.size() == 1);
}

TEST_CASE_TEMPLATE("periodic tree 2d",
                   T,
                   double,
                   float,
                   int32_t,
                   int16_t)
{
  using tree = tree<2, T>;
  using aabb = tree::aabb;
  using point = tree::point;

  Vc::fixed_size_simd<T, 2> bounds = {{10, 10}};

  tree t;
  auto n1 = t.insert({{0, 0}, {2, 2}});
  auto n2 = t.insert({{1, 1}, {3, 3}});
  auto n3 = t.insert({{2, 2}, {4, 4}});
  auto n4 = t.insert({{5, 5}, {7, 7}});

  REQUIRE(t.get_aabb(n1) == aabb{{0, 0}, {2, 2}});
  REQUIRE(t.get_aabb(n2) == aabb{{1, 1}, {3, 3}});
  REQUIRE(t.get_aabb(n3) == aabb{{2, 2}, {4, 4}});
  REQUIRE(t.get_aabb(n4) == aabb{{5, 5}, {7, 7}});

  auto intersections = t.get_overlaps(point{10, 1}, true, bounds);
  REQUIRE(intersections.size() == 1);

  intersections = t.get_overlaps(point{10, 11}, true, bounds);
  REQUIRE(intersections.size() == 1);

  intersections = t.get_overlaps(point{11, 1}, true, bounds);
  REQUIRE(intersections.size() == 2);

  intersections = t.get_overlaps(point{11, 11}, true, bounds);
  REQUIRE(intersections.size() == 2);
}

TEST_CASE_TEMPLATE("optimal tree 2d", T, double, float, int, short)
{
  using tree = tree<2, T>;
  using aabb = tree::aabb;
  using point = tree::point;

  std::array<aabb, 4> bbs = {
      {{{0, 0}, {2, 2}}, {{1, 1}, {3, 3}}, {{2, 2}, {4, 4}}, {{5, 5}, {7, 7}}}};

  tree t(bbs);
  REQUIRE(t.size() == 4);
  REQUIRE(t.get_height() == 3);
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

  REQUIRE(t.any_overlap(point{1, 1}, [] { return true; }));
  REQUIRE(t.any_overlap(point{1, 1}, [](tree::node_id) { return true; }));
  REQUIRE(t.any_overlap(point{1, 1}, [](aabb) { return true; }));
  REQUIRE(t.any_overlap(point{1, 1}, [](tree::node_id, aabb) { return true; }));

  t.clear();
  REQUIRE(t.size() == 0);
  t.insert({{2, 2}, {4, 4}});
  REQUIRE(t.size() == 1);
}
