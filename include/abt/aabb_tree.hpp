/*
  Copyright (c) 2009 Erin Catto http://www.box2d.org
  Copyright (c) 2016-2018 Lester Hedges <lester.hedges+aabbcc@gmail.com>

  This software is provided 'as-is', without any express or implied
  warranty. In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.

  3. This notice may not be removed or altered from any source distribution.

  This code was adapted from parts of the Box2D Physics Engine,
  http://www.box2d.org
*/

#ifndef _AABB_H
#define _AABB_H

#include <algorithm>
#include <array>
#include <cassert>
#include <stdexcept>
#include <vector>

#include <limits>

#include <unordered_map>

//#define vector moo

/// Null node flag.
const unsigned int NULL_NODE = 0xffffffff;

namespace abt {
template <unsigned Dim, typename ValTy = double>
struct point {
  using value_type = ValTy;

  template <typename... Values>
  point(Values... values) : values{value_type(values)...}
  {
  }

  bool operator==(const point &other) const { return values == other.values; }

  value_type x() { return values[0]; }
  value_type y() { return values[1]; }
  value_type z() { return values[2]; }

  constexpr size_t size() const { return Dim; }

  value_type &operator[](size_t idx) { return values[idx]; }
  const value_type &operator[](size_t idx) const { return values[idx]; }

  std::array<value_type, Dim> values = {};
};
/*! \brief The axis-aligned bounding box object.

    Axis-aligned bounding boxes (AABBs) store information for the minimum
    orthorhombic bounding-box for an object. Support is provided for
    dimensions >= 2. (In 2D the bounding box is either a rectangle,
    in 3D it is a rectangular prism.)

    Class member functions provide functionality for merging AABB objects
    and testing overlap with other AABBs.
 */
template <unsigned Dim, typename ValTy = double>
class aabb {
 public:
  using point = abt::point<Dim, ValTy>;
  using value_type = ValTy;

  /// Constructor.
  aabb() = default;

  //! Constructor.
  /*! \param lowerBound_
          The lower bound in each dimension.

      \param upperBound_
          The upper bound in each dimension.
   */
  aabb(const point &lower_bound, const point &upper_bound)
      : lowerBound(lower_bound), upperBound(upper_bound)
  {
    surfaceArea = compute_surface_area();
    centre = compute_center();
  }

  static aabb of_sphere(const point &center, value_type radius)
  {
    point lb, ub;
    for (unsigned int i = 0; i < Dim; i++) {
      lb[i] = center[i] - radius;
      ub[i] = center[i] + radius;
    }
    return {lb, ub};
  }

  /// Compute the surface area of the box.
  value_type compute_surface_area() const
  {
    // Sum of "area" of all the sides.
    value_type sum = 0;

    // General formula for one side: hold one dimension constant
    // and multiply by all the other ones.
    for (unsigned int d1 = 0; d1 < Dim; d1++) {
      // "Area" of current side.
      value_type product = 1;

      for (unsigned int d2 = 0; d2 < Dim; d2++) {
        if (d1 == d2)
          continue;

        value_type dx = upperBound[d2] - lowerBound[d2];
        product *= dx;
      }

      // Update the sum.
      sum += product;
    }

    return 2 * sum;
  }

  /// Get the surface area of the box.
  value_type get_surface_area() const { return surfaceArea; }

  //! Merge two AABBs into this one.
  /*! \param aabb1
          A reference to the first AABB.

      \param aabb2
          A reference to the second AABB.
   */
  void merge(const aabb &aabb1, const aabb &aabb2)
  {
    for (unsigned int i = 0; i < lowerBound.size(); i++) {
      lowerBound[i] = std::min(aabb1.lowerBound[i], aabb2.lowerBound[i]);
      upperBound[i] = std::max(aabb1.upperBound[i], aabb2.upperBound[i]);
    }

    surfaceArea = compute_surface_area();
    centre = compute_center();
  }

  //! Test whether the AABB is contained within this one.
  /*! \param aabb
          A reference to the AABB.

      \return
          Whether the AABB is fully contained.
   */
  bool contains(const aabb &aabb) const
  {
    for (unsigned int i = 0; i < lowerBound.size(); i++) {
      if (aabb.lowerBound[i] < lowerBound[i])
        return false;
      if (aabb.upperBound[i] > upperBound[i])
        return false;
    }

    return true;
  }

  //! Test whether the AABB overlaps this one.
  /*! \param aabb
          A reference to the AABB.

      \param touchIsOverlap
          Does touching constitute an overlap?

      \return
          Whether the AABB overlaps.
   */
  bool overlaps(const aabb &aabb, bool touchIsOverlap) const
  {
    if (touchIsOverlap) {
      for (unsigned int i = 0; i < lowerBound.size(); ++i) {
        if (aabb.upperBound[i] < lowerBound[i] ||
            aabb.lowerBound[i] > upperBound[i]) {
          return false;
        }
      }
    }
    else {
      for (unsigned int i = 0; i < lowerBound.size(); ++i) {
        if (aabb.upperBound[i] <= lowerBound[i] ||
            aabb.lowerBound[i] >= upperBound[i]) {
          return false;
        }
      }
    }

    return true;
  }

  //! Test whether the point overlaps this one.
  /*! \param point
          A reference to the point.

      \param touchIsOverlap
          Does touching constitute an overlap?

      \return
          Whether the AABB overlaps.
   */
  bool overlaps(const point &pt, bool touchIsOverlap) const
  {
    if (touchIsOverlap) {
      for (unsigned int i = 0; i < lowerBound.size(); ++i) {
        if (pt[i] < lowerBound[i] || pt[i] > upperBound[i]) {
          return false;
        }
      }
    }
    else {
      for (unsigned int i = 0; i < lowerBound.size(); ++i) {
        if (pt[i] <= lowerBound[i] || pt[i] >= upperBound[i]) {
          return false;
        }
      }
    }

    return true;
  }

  //! Compute the centre of the AABB.
  /*! \returns
          The position vector of the AABB centre.
   */
  point compute_center()
  {
    point center;

    for (unsigned int i = 0; i < Dim; i++)
      center[i] = 0.5 * (lowerBound[i] + upperBound[i]);

    return center;
  }

  /// Lower bound of AABB in each dimension.
  point lowerBound;

  /// Upper bound of AABB in each dimension.
  point upperBound;

  /// The position of the AABB centre.
  point centre;

  /// The AABB's surface area.
  value_type surfaceArea = 0;
};

enum visit_action : char { visit_stop, visit_continue };

/*! \brief The dynamic AABB tree.

    The dynamic AABB tree is a hierarchical data structure that can be used
    to efficiently query overlaps between objects of arbitrary shape and
    size that lie inside of a simulation box. Support is provided for
    periodic and non-periodic boxes, as well as boxes with partial
    periodicity, e.g. periodic along specific axes.
 */
template <unsigned Dim, typename ValTy = double>
class tree {
 public:
  using value_type = ValTy;
  using aabb = abt::aabb<Dim, value_type>;
  using point = abt::point<Dim, value_type>;
  template <typename Ty>
  using vec = std::array<Ty, Dim>;

 private:
  /*! \brief A node of the AABB tree.

   Each node of the tree contains an AABB object which corresponds to a
   entry, or a group of entrys, in the simulation box. The AABB
   objects of individual entrys are "fattened" before they are stored
   to avoid having to continually update and rebalance the tree when
   displacements are small.

   Nodes are aware of their position within in the tree. The isLeaf member
   function allows the tree to query whether the node is a leaf, i.e. to
   determine whether it holds a single entry.
  */
  struct node {
    /// Constructor.
    node() = default;

    /// The fattened axis-aligned bounding box.
    aabb bb;

    /// Index of the parent node.
    unsigned int parent = 0;

    /// Index of the next node.
    unsigned int next = 0;

    /// Index of the left-hand child.
    unsigned int left = 0;

    /// Index of the right-hand child.
    unsigned int right = 0;

    /// Height of the node. This is 0 for a leaf and -1 for a free node.
    int height = 0;

    /// The id of the entry that the node contains (leaf nodes only).
    unsigned int id = 0;

    //! Test whether the node is a leaf.
    /*! \return
            Whether the node is a leaf node.
     */
    bool isLeaf() const { return (left == NULL_NODE); }
  };

 public:
  //! Constructor (non-periodic).
  /*! \param skin_thickness
          The skin thickness for fattened AABBs, as a fraction
          of the AABB base length.

      \param nentrys
          The number of entrys (for fixed entry number systems).

      \param touchIsOverlap
          Does touching count as overlapping in query operations?
   */
  tree(value_type skin_thickness = 0, unsigned int initial_size = 16)
      : m_is_periodic(false), m_skin_thickness(skin_thickness)
  {
    // Initialise the periodicity vector.
    std::fill(m_periodicity.begin(), m_periodicity.end(), false);

    // Initialise the tree.
    m_root = NULL_NODE;
    m_node_count = 0;
    m_node_capacity = initial_size;
    m_nodes.resize(m_node_capacity);

    // Build a linked list for the list of free nodes.
    for (unsigned int i = 0; i < m_node_capacity - 1; i++) {
      m_nodes[i].next = i + 1;
      m_nodes[i].height = -1;
    }
    m_nodes[m_node_capacity - 1].next = NULL_NODE;
    m_nodes[m_node_capacity - 1].height = -1;

    // Assign the index of the first free node.
    m_free_list = 0;
  }

  //! Constructor (custom periodicity).
  /*! \param skinThickness
          The skin thickness for fattened AABBs, as a fraction
          of the AABB base length.

      \param periodicity
          Whether the system is periodic in each dimension.

      \param boxSize
          The size of the simulation box in each dimension.

      \param nentrys
          The number of entrys (for fixed entry number systems).

      \param touchIsOverlap
          Does touching count as overlapping in query operations?
   */
  tree(double skinThickness,
       const vec<bool> &periodicity,
       const vec<double> &boxSize,
       unsigned int initial_size = 16)
      : m_skin_thickness(skinThickness),
        m_periodicity(periodicity),
        m_box_size(boxSize)
  {
    // Initialise the tree.
    m_root = NULL_NODE;
    m_node_count = 0;
    m_node_capacity = initial_size;
    m_nodes.resize(m_node_capacity);

    // Build a linked list for the list of free nodes.
    for (unsigned int i = 0; i < m_node_capacity - 1; i++) {
      m_nodes[i].next = i + 1;
      m_nodes[i].height = -1;
    }
    m_nodes[m_node_capacity - 1].next = NULL_NODE;
    m_nodes[m_node_capacity - 1].height = -1;

    // Assign the index of the first free node.
    m_free_list = 0;

    // Check periodicity.
    m_is_periodic = false;

    for (unsigned int i = 0; i < Dim; i++) {
      m_pos_min_image[i] = 0.5 * boxSize[i];
      m_neg_min_image[i] = -0.5 * boxSize[i];

      if (periodicity[i])
        m_is_periodic = true;
    }
  }

  //! Set the periodicity of the simulation box.
  /*! \param periodicity_
          Whether the system is periodic in each dimension.
   */
  void set_periodicity(const vec<bool> &periodicity)
  {
    m_periodicity = periodicity;
  }

  //! Set the size of the simulation box.
  /*! \param box_size
          The size of the simulation box in each dimension.
   */
  void set_box_size(const vec<double> &box_size) { m_box_size = box_size; }

  //! Insert a entry into the tree (arbitrary shape with bounding box).
  /*! \param index
          The index of the entry.

      \param lowerBound
          The lower bound in each dimension.

      \param upperBound
          The upper bound in each dimension.
   */
  void insert(unsigned int id, const aabb &bb)
  {
    // Make sure the entry doesn't already exist.
    if (m_id_map.count(id) != 0) {
      throw std::invalid_argument("[ERROR]: entry already exists in tree!");
    }

    // Allocate a new node for the entry.
    unsigned int node_idx = allocate_node();
    auto &node = m_nodes[node_idx];
    node.id = id;
    node.bb = bb;

    // Fatten the AABB.
    for (unsigned int i = 0; i < Dim; i++) {
      auto sz = bb.upperBound[i] - bb.lowerBound[i];
      node.bb.lowerBound[i] -= m_skin_thickness * sz;
      node.bb.upperBound[i] += m_skin_thickness * sz;
    }
    node.bb.surfaceArea = node.bb.compute_surface_area();
    node.bb.centre = node.bb.compute_center();

    // Zero the height.
    node.height = 0;

    // Insert a new leaf into the tree.
    insert_leaf(node_idx);

    // Add the new entry to the map.
    m_id_map.insert(std::unordered_map<unsigned int, unsigned int>::value_type(
        id, node_idx));
  }

  /// Return the number of entrys in the tree.
  unsigned int size() const { return m_id_map.size(); }

  //! Remove a entry from the tree.
  /*! \param entry
          The entry index (entryMap will be used to map the node).
   */
  void remove(unsigned int id)
  {
    // Map iterator.
    std::unordered_map<unsigned int, unsigned int>::iterator it;

    // Find the entry.
    it = m_id_map.find(id);

    // The entry doesn't exist.
    if (it == m_id_map.end()) {
      throw std::invalid_argument("[ERROR]: Invalid entry index!");
    }

    // Extract the node index.
    unsigned int node = it->second;

    // Erase the entry from the map.
    m_id_map.erase(it);

    assert(node < m_node_capacity);
    assert(m_nodes[node].isLeaf());

    remove_leaf(node);
    free_node(node);
  }

  /// Remove all entrys from the tree.
  void clear()
  {
    // Iterator pointing to the start of the entry map.
    std::unordered_map<unsigned int, unsigned int>::iterator it =
        m_id_map.begin();

    // Iterate over the map.
    while (it != m_id_map.end()) {
      // Extract the node index.
      unsigned int node = it->second;

      assert(node < m_node_capacity);
      assert(m_nodes[node].isLeaf());

      remove_leaf(node);
      free_node(node);

      it++;
    }

    // Clear the entry map.
    m_id_map.clear();
  }

  //! Update the tree if a entry moves outside its fattened AABB.
  /*! \param entry
          The entry index (entryMap will be used to map the node).

      \param lowerBound
          The lower bound in each dimension.

      \param upperBound
          The upper bound in each dimension.

      \param alwaysReinsert
          Always reinsert the entry, even if it's within its old AABB
     (default: false)
   */
  bool update(unsigned int id, aabb bb, bool always_reinsert = false)
  {
    // Map iterator.
    std::unordered_map<unsigned int, unsigned int>::iterator it;

    // Find the entry.
    it = m_id_map.find(id);

    // The entry doesn't exist.
    if (it == m_id_map.end()) {
      throw std::invalid_argument("[ERROR]: Invalid entry index!");
    }

    // Extract the node index.
    unsigned int node_idx = it->second;
    auto &node = m_nodes[node_idx];

    assert(node_idx < m_node_capacity);
    assert(node.isLeaf());

    // No need to update if the entry is still within its fattened AABB.
    if (!always_reinsert && node.bb.contains(bb))
      return false;

    // Remove the current leaf.
    remove_leaf(node_idx);

    // Fatten the new AABB.
    for (unsigned int i = 0; i < Dim; i++) {
      auto sz = bb.upperBound[i] - bb.lowerBound[i];
      bb.lowerBound[i] -= m_skin_thickness * sz;
      bb.upperBound[i] += m_skin_thickness * sz;
    }

    // Assign the new AABB.
    node.bb = bb;

    // Update the surface area and centroid.
    node.bb.surfaceArea = node.bb.compute_surface_area();
    node.bb.centre = node.bb.compute_center();

    // Insert a new leaf node.
    insert_leaf(node_idx);

    return true;
  }

  //! Query the tree to find candidate interactions for a entry.
  /*! \param entry
          The entry index.

      \return entrys
          A vector of entry indices.
   */
  std::vector<unsigned int> get_overlaps(unsigned int id)
  {
    // Make sure that this is a valid entry.
    if (m_id_map.count(id) == 0) {
      throw std::invalid_argument("[ERROR]: Invalid entry index!");
    }

    // Test overlap of entry AABB against all other entrys.
    return query(m_nodes[m_id_map.find(id)->second].bb);
  }

  //! Query the tree to find candidate interactions for an AABB.
  /*! \param aabb
          The AABB.
      \param out
          An output iterator

      \return entrys
          A vector of entry indices.
   */
  template <class Query, class OutputIterator>
  void get_overlaps(const Query &query,
                    OutputIterator out,
                    bool include_touch = true) const
  {
    visit_overlaps(
        query, [&](unsigned int id) { *out++ = id; }, include_touch);
  }

  //! Query the tree to find candidate interactions for an AABB.
  /*! \param aabb
          The AABB.

      \return entrys
          A vector of entry indices.
   */
  template <class Query>
  std::vector<unsigned int> get_overlaps(const Query &query,
                                         bool include_touch = true) const
  {
    std::vector<unsigned int> overlaps;
    visit_overlaps(
        query, [&](unsigned int id) { overlaps.push_back(id); }, include_touch);
    return overlaps;
  }

  template <class Query>
  bool any_overlap(const Query &query, bool include_touch = true) const
  {
    return any_overlap(
        query, [](unsigned) { return true; }, include_touch);
  }

  template <class Query, class Fn>
  bool any_overlap(const Query &query, Fn &&fn, bool include_touch = true) const
  {
    constexpr bool fn_with_bb = std::is_invocable_v<Fn, unsigned, aabb>;
    static_assert(std::is_invocable_v<Fn, unsigned int> || fn_with_bb,
                  "Wrong function signature");
    bool overlap = false;
    if constexpr (fn_with_bb) {
      auto wrap_fn = [&overlap, &fn](unsigned int id, const aabb &bb) {
        bool success = std::forward<Fn>(fn)(id, bb);
        overlap |= success;
        return success ? visit_stop : visit_continue;
      };
      visit_overlaps(query, wrap_fn, include_touch);
    }
    else {
      auto wrap_fn = [&overlap, &fn](unsigned int id) {
        bool success = std::forward<Fn>(fn)(id);
        overlap |= success;
        return success ? visit_stop : visit_continue;
      };
      visit_overlaps(query, wrap_fn, include_touch);
    }

    return overlap;
  }

  template <class Query, class Fn>
  void visit_overlaps(const Query &query,
                      Fn &&fn,
                      bool include_touch = true) const
  {
    constexpr bool query_is_point = std::is_same_v<Query, point>;
    constexpr bool query_is_aabb = std::is_same_v<Query, aabb>;
    static_assert(query_is_point || query_is_aabb,
                  "Only point or aabb queries are supported");
    constexpr bool fn_with_bb = std::is_invocable_v<Fn, unsigned, aabb>;

    using rt = decltype(get_return_type<fn_with_bb>(std::forward<Fn>(fn)));
    constexpr bool fn_returns_action = std::is_convertible_v<rt, visit_action>;
    static_assert(fn_returns_action || std::is_same_v<rt, void>,
                  "Only void or visit_action return types are allowed");

    // Make sure the tree isn't empty.
    if (m_id_map.size() == 0) {
      return;
    }

    static thread_local std::vector<unsigned int> stack(64);
    stack.clear();

    stack.push_back(m_root);

    while (stack.size() > 0) {
      unsigned int node = stack.back();
      stack.pop_back();

      // Copy the AABB.
      auto nodeAABB = m_nodes[node].bb;

      if (node == NULL_NODE)
        continue;

      if (m_is_periodic) {
        vec<double> separation;
        vec<double> shift;
        point center;
        if constexpr (query_is_point) {
          center = query;
        }
        else {
          center = query.centre;
        }
        for (unsigned int i = 0; i < Dim; i++)
          separation[i] = nodeAABB.centre[i] - center[i];

        bool isShifted = minimum_image(separation, shift);

        // Shift the AABB.
        if (isShifted) {
          for (unsigned int i = 0; i < Dim; i++) {
            nodeAABB.lowerBound[i] += shift[i];
            nodeAABB.upperBound[i] += shift[i];
          }
        }
      }

      // Test for overlap between the AABBs.
      if (nodeAABB.overlaps(query, include_touch)) {
        // Check that we're at a leaf node.
        if (m_nodes[node].isLeaf()) {
          const auto &n = m_nodes[node];
          if constexpr (fn_returns_action) {
            visit_action visit_act;
            if constexpr (fn_with_bb) {
              visit_act = std::forward<Fn>(fn)(n.id, n.bb);
            }
            else {
              visit_act = std::forward<Fn>(fn)(n.id);
            }
            if (visit_act == visit_stop) {
              return;
            }
          }
          else {
            if constexpr (fn_with_bb) {
              std::forward<Fn>(fn)(n.id, n.bb);
            }
            else {
              std::forward<Fn>(fn)(n.id);
            }
          }
        }
        else {
          stack.push_back(m_nodes[node].left);
          stack.push_back(m_nodes[node].right);
        }
      }
    }
  }

  //! Get a entry AABB.
  /*! \param entry
          The entry index.
   */
  const aabb &get_aabb(unsigned int id) const
  {
    return m_nodes[m_id_map.at(id)].bb;
  }

  //! Get the height of the tree.
  /*! \return
          The height of the binary tree.
   */
  unsigned int get_height() const
  {
    if (m_root == NULL_NODE)
      return 0;
    return m_nodes[m_root].height;
  }

  /// Validate the tree.
  void validate() const
  {
#ifndef NDEBUG
    validate_structure(m_root);
    validate_metrics(m_root);

    unsigned int freeCount = 0;
    unsigned int freeIndex = m_free_list;

    while (freeIndex != NULL_NODE) {
      assert(freeIndex < m_node_capacity);
      freeIndex = m_nodes[freeIndex].next;
      freeCount++;
    }

    assert(get_height() == compute_height());
    assert((m_node_count + freeCount) == m_node_capacity);
#endif
  }

  /// Rebuild an optimal tree.
  void rebuild()
  {
    std::vector<unsigned int> nodeIndices(m_node_count);
    unsigned int count = 0;

    for (unsigned int i = 0; i < m_node_capacity; i++) {
      // Free node.
      if (m_nodes[i].height < 0)
        continue;

      if (m_nodes[i].isLeaf()) {
        m_nodes[i].parent = NULL_NODE;
        nodeIndices[count] = i;
        count++;
      }
      else
        free_node(i);
    }

    while (count > 1) {
      double minCost = std::numeric_limits<double>::max();
      int iMin = -1, jMin = -1;

      for (unsigned int i = 0; i < count; i++) {
        aabb aabbi = m_nodes[nodeIndices[i]].bb;

        for (unsigned int j = i + 1; j < count; j++) {
          aabb aabbj = m_nodes[nodeIndices[j]].bb;
          aabb aabb;
          aabb.merge(aabbi, aabbj);
          double cost = aabb.get_surface_area();

          if (cost < minCost) {
            iMin = i;
            jMin = j;
            minCost = cost;
          }
        }
      }

      unsigned int index1 = nodeIndices[iMin];
      unsigned int index2 = nodeIndices[jMin];

      unsigned int parent = allocate_node();
      m_nodes[parent].left = index1;
      m_nodes[parent].right = index2;
      m_nodes[parent].height =
          1 + std::max(m_nodes[index1].height, m_nodes[index2].height);
      m_nodes[parent].bb.merge(m_nodes[index1].bb, m_nodes[index2].bb);
      m_nodes[parent].parent = NULL_NODE;

      m_nodes[index1].parent = parent;
      m_nodes[index2].parent = parent;

      nodeIndices[jMin] = nodeIndices[count - 1];
      nodeIndices[iMin] = parent;
      count--;
    }

    m_root = nodeIndices[0];

    validate();
  }

 private:
  /// The index of the root node.
  unsigned int m_root;

  /// The dynamic tree.
  std::vector<node> m_nodes;

  /// The current number of nodes in the tree.
  unsigned int m_node_count;

  /// The current node capacity.
  unsigned int m_node_capacity;

  /// The position of node at the top of the free list.
  unsigned int m_free_list;

  /// Whether the system is periodic along at least one axis.
  bool m_is_periodic;

  /// The skin thickness of the fattened AABBs, as a fraction of the AABB base
  /// length.
  double m_skin_thickness;

  /// Whether the system is periodic along each axis.
  vec<bool> m_periodicity;

  /// The size of the system in each dimension.
  vec<double> m_box_size;

  /// The position of the negative minimum image.
  point m_neg_min_image;

  /// The position of the positive minimum image.
  point m_pos_min_image;

  /// A map between entry and node indices.
  std::unordered_map<unsigned int, unsigned int> m_id_map;

 private:
  template <bool with_aabb, class Fn>
  static decltype(auto) get_return_type(Fn &&fn)
  {
    if constexpr (with_aabb) {
      return std::forward<Fn>(fn)(0, aabb());
    }
    else {
      return std::forward<Fn>(fn)(0);
    }
  }
  //! Allocate a new node.
  /*! \return
          The index of the allocated node.
   */
  unsigned int allocate_node()
  {
    // Exand the node pool as needed.
    if (m_free_list == NULL_NODE) {
      assert(m_node_count == m_node_capacity);

      // The free list is empty. Rebuild a bigger pool.
      m_node_capacity *= 2;
      m_nodes.resize(m_node_capacity);

      // Build a linked list for the list of free nodes.
      for (unsigned int i = m_node_count; i < m_node_capacity - 1; i++) {
        m_nodes[i].next = i + 1;
        m_nodes[i].height = -1;
      }
      m_nodes[m_node_capacity - 1].next = NULL_NODE;
      m_nodes[m_node_capacity - 1].height = -1;

      // Assign the index of the first free node.
      m_free_list = m_node_count;
    }

    // Peel a node off the free list.
    unsigned int node = m_free_list;
    m_free_list = m_nodes[node].next;
    m_nodes[node].parent = NULL_NODE;
    m_nodes[node].left = NULL_NODE;
    m_nodes[node].right = NULL_NODE;
    m_nodes[node].height = 0;
    m_node_count++;

    return node;
  }

  //! Free an existing node.
  /*! \param node
          The index of the node to be freed.
   */
  void free_node(unsigned int node)
  {
    assert(node < m_node_capacity);
    assert(0 < m_node_count);

    m_nodes[node].next = m_free_list;
    m_nodes[node].height = -1;
    m_free_list = node;
    m_node_count--;
  }

  //! Insert a leaf into the tree.
  /*! \param leaf
          The index of the leaf node.
   */
  void insert_leaf(unsigned int leaf)
  {
    if (m_root == NULL_NODE) {
      m_root = leaf;
      m_nodes[m_root].parent = NULL_NODE;
      return;
    }

    // Find the best sibling for the node.

    aabb leafAABB = m_nodes[leaf].bb;
    unsigned int index = m_root;

    while (!m_nodes[index].isLeaf()) {
      // Extract the children of the node.
      unsigned int left = m_nodes[index].left;
      unsigned int right = m_nodes[index].right;

      double surfaceArea = m_nodes[index].bb.get_surface_area();

      aabb combinedAABB;
      combinedAABB.merge(m_nodes[index].bb, leafAABB);
      double combinedSurfaceArea = combinedAABB.get_surface_area();

      // Cost of creating a new parent for this node and the new leaf.
      double cost = 2.0 * combinedSurfaceArea;

      // Minimum cost of pushing the leaf further down the tree.
      double inheritanceCost = 2.0 * (combinedSurfaceArea - surfaceArea);

      // Cost of descending to the left.
      double costLeft;
      if (m_nodes[left].isLeaf()) {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[left].bb);
        costLeft = aabb.get_surface_area() + inheritanceCost;
      }
      else {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[left].bb);
        double oldArea = m_nodes[left].bb.get_surface_area();
        double newArea = aabb.get_surface_area();
        costLeft = (newArea - oldArea) + inheritanceCost;
      }

      // Cost of descending to the right.
      double costRight;
      if (m_nodes[right].isLeaf()) {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[right].bb);
        costRight = aabb.get_surface_area() + inheritanceCost;
      }
      else {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[right].bb);
        double oldArea = m_nodes[right].bb.get_surface_area();
        double newArea = aabb.get_surface_area();
        costRight = (newArea - oldArea) + inheritanceCost;
      }

      // Descend according to the minimum cost.
      if ((cost < costLeft) && (cost < costRight))
        break;

      // Descend.
      if (costLeft < costRight)
        index = left;
      else
        index = right;
    }

    unsigned int sibling = index;

    // Create a new parent.
    unsigned int oldParent = m_nodes[sibling].parent;
    unsigned int newParent = allocate_node();
    m_nodes[newParent].parent = oldParent;
    m_nodes[newParent].bb.merge(leafAABB, m_nodes[sibling].bb);
    m_nodes[newParent].height = m_nodes[sibling].height + 1;

    // The sibling was not the root.
    if (oldParent != NULL_NODE) {
      if (m_nodes[oldParent].left == sibling)
        m_nodes[oldParent].left = newParent;
      else
        m_nodes[oldParent].right = newParent;

      m_nodes[newParent].left = sibling;
      m_nodes[newParent].right = leaf;
      m_nodes[sibling].parent = newParent;
      m_nodes[leaf].parent = newParent;
    }
    // The sibling was the root.
    else {
      m_nodes[newParent].left = sibling;
      m_nodes[newParent].right = leaf;
      m_nodes[sibling].parent = newParent;
      m_nodes[leaf].parent = newParent;
      m_root = newParent;
    }

    // Walk back up the tree fixing heights and AABBs.
    index = m_nodes[leaf].parent;
    while (index != NULL_NODE) {
      index = balance(index);

      unsigned int left = m_nodes[index].left;
      unsigned int right = m_nodes[index].right;

      assert(left != NULL_NODE);
      assert(right != NULL_NODE);

      m_nodes[index].height =
          1 + std::max(m_nodes[left].height, m_nodes[right].height);
      m_nodes[index].bb.merge(m_nodes[left].bb, m_nodes[right].bb);

      index = m_nodes[index].parent;
    }
  }

  //! Remove a leaf from the tree.
  /*! \param leaf
          The index of the leaf node.
   */
  void remove_leaf(unsigned int leaf)
  {
    if (leaf == m_root) {
      m_root = NULL_NODE;
      return;
    }

    unsigned int parent = m_nodes[leaf].parent;
    unsigned int grandParent = m_nodes[parent].parent;
    unsigned int sibling;

    if (m_nodes[parent].left == leaf)
      sibling = m_nodes[parent].right;
    else
      sibling = m_nodes[parent].left;

    // Destroy the parent and connect the sibling to the grandparent.
    if (grandParent != NULL_NODE) {
      if (m_nodes[grandParent].left == parent)
        m_nodes[grandParent].left = sibling;
      else
        m_nodes[grandParent].right = sibling;

      m_nodes[sibling].parent = grandParent;
      free_node(parent);

      // Adjust ancestor bounds.
      unsigned int index = grandParent;
      while (index != NULL_NODE) {
        index = balance(index);

        unsigned int left = m_nodes[index].left;
        unsigned int right = m_nodes[index].right;

        m_nodes[index].bb.merge(m_nodes[left].bb, m_nodes[right].bb);
        m_nodes[index].height =
            1 + std::max(m_nodes[left].height, m_nodes[right].height);

        index = m_nodes[index].parent;
      }
    }
    else {
      m_root = sibling;
      m_nodes[sibling].parent = NULL_NODE;
      free_node(parent);
    }
  }

  //! Balance the tree.
  /*! \param leaf
          The index of the node.
   */
  unsigned int balance(unsigned int node)
  {
    assert(node != NULL_NODE);

    if (m_nodes[node].isLeaf() || (m_nodes[node].height < 2))
      return node;

    unsigned int left = m_nodes[node].left;
    unsigned int right = m_nodes[node].right;

    assert(left < m_node_capacity);
    assert(right < m_node_capacity);

    int currentBalance = m_nodes[right].height - m_nodes[left].height;

    // Rotate right branch up.
    if (currentBalance > 1) {
      unsigned int rightLeft = m_nodes[right].left;
      unsigned int rightRight = m_nodes[right].right;

      assert(rightLeft < m_node_capacity);
      assert(rightRight < m_node_capacity);

      // Swap node and its right-hand child.
      m_nodes[right].left = node;
      m_nodes[right].parent = m_nodes[node].parent;
      m_nodes[node].parent = right;

      // The node's old parent should now point to its right-hand child.
      if (m_nodes[right].parent != NULL_NODE) {
        if (m_nodes[m_nodes[right].parent].left == node)
          m_nodes[m_nodes[right].parent].left = right;
        else {
          assert(m_nodes[m_nodes[right].parent].right == node);
          m_nodes[m_nodes[right].parent].right = right;
        }
      }
      else
        m_root = right;

      // Rotate.
      if (m_nodes[rightLeft].height > m_nodes[rightRight].height) {
        m_nodes[right].right = rightLeft;
        m_nodes[node].right = rightRight;
        m_nodes[rightRight].parent = node;
        m_nodes[node].bb.merge(m_nodes[left].bb, m_nodes[rightRight].bb);
        m_nodes[right].bb.merge(m_nodes[node].bb, m_nodes[rightLeft].bb);

        m_nodes[node].height =
            1 + std::max(m_nodes[left].height, m_nodes[rightRight].height);
        m_nodes[right].height =
            1 + std::max(m_nodes[node].height, m_nodes[rightLeft].height);
      }
      else {
        m_nodes[right].right = rightRight;
        m_nodes[node].right = rightLeft;
        m_nodes[rightLeft].parent = node;
        m_nodes[node].bb.merge(m_nodes[left].bb, m_nodes[rightLeft].bb);
        m_nodes[right].bb.merge(m_nodes[node].bb, m_nodes[rightRight].bb);

        m_nodes[node].height =
            1 + std::max(m_nodes[left].height, m_nodes[rightLeft].height);
        m_nodes[right].height =
            1 + std::max(m_nodes[node].height, m_nodes[rightRight].height);
      }

      return right;
    }

    // Rotate left branch up.
    if (currentBalance < -1) {
      unsigned int leftLeft = m_nodes[left].left;
      unsigned int leftRight = m_nodes[left].right;

      assert(leftLeft < m_node_capacity);
      assert(leftRight < m_node_capacity);

      // Swap node and its left-hand child.
      m_nodes[left].left = node;
      m_nodes[left].parent = m_nodes[node].parent;
      m_nodes[node].parent = left;

      // The node's old parent should now point to its left-hand child.
      if (m_nodes[left].parent != NULL_NODE) {
        if (m_nodes[m_nodes[left].parent].left == node)
          m_nodes[m_nodes[left].parent].left = left;
        else {
          assert(m_nodes[m_nodes[left].parent].right == node);
          m_nodes[m_nodes[left].parent].right = left;
        }
      }
      else
        m_root = left;

      // Rotate.
      if (m_nodes[leftLeft].height > m_nodes[leftRight].height) {
        m_nodes[left].right = leftLeft;
        m_nodes[node].left = leftRight;
        m_nodes[leftRight].parent = node;
        m_nodes[node].bb.merge(m_nodes[right].bb, m_nodes[leftRight].bb);
        m_nodes[left].bb.merge(m_nodes[node].bb, m_nodes[leftLeft].bb);

        m_nodes[node].height =
            1 + std::max(m_nodes[right].height, m_nodes[leftRight].height);
        m_nodes[left].height =
            1 + std::max(m_nodes[node].height, m_nodes[leftLeft].height);
      }
      else {
        m_nodes[left].right = leftRight;
        m_nodes[node].left = leftLeft;
        m_nodes[leftLeft].parent = node;
        m_nodes[node].bb.merge(m_nodes[right].bb, m_nodes[leftLeft].bb);
        m_nodes[left].bb.merge(m_nodes[node].bb, m_nodes[leftRight].bb);

        m_nodes[node].height =
            1 + std::max(m_nodes[right].height, m_nodes[leftLeft].height);
        m_nodes[left].height =
            1 + std::max(m_nodes[node].height, m_nodes[leftRight].height);
      }

      return left;
    }

    return node;
  }

  //! Compute the height of the tree.
  /*! \return
          The height of the entire tree.
   */
  unsigned int compute_height() const { return compute_height(m_root); }

  //! Compute the height of a sub-tree.
  /*! \param node
          The index of the root node.

      \return
          The height of the sub-tree.
   */
  unsigned int compute_height(unsigned int node) const
  {
    assert(node < m_node_capacity);

    if (m_nodes[node].isLeaf())
      return 0;

    unsigned int height1 = computeHeight(m_nodes[node].left);
    unsigned int height2 = computeHeight(m_nodes[node].right);

    return 1 + std::max(height1, height2);
  }

  //! Assert that the sub-tree has a valid structure.
  /*! \param node
          The index of the root node.
   */
  void validate_structure(unsigned int node) const
  {
    if (node == NULL_NODE)
      return;

    if (node == m_root)
      assert(m_nodes[node].parent == NULL_NODE);

    unsigned int left = m_nodes[node].left;
    unsigned int right = m_nodes[node].right;

    if (m_nodes[node].isLeaf()) {
      assert(left == NULL_NODE);
      assert(right == NULL_NODE);
      assert(m_nodes[node].height == 0);
      return;
    }

    assert(left < m_node_capacity);
    assert(right < m_node_capacity);

    assert(m_nodes[left].parent == node);
    assert(m_nodes[right].parent == node);

    validate_structure(left);
    validate_structure(right);
  }

  //! Assert that the sub-tree has valid metrics.
  /*! \param node
          The index of the root node.
   */
  void validate_metrics(unsigned int node) const
  {
    if (node == NULL_NODE)
      return;

    unsigned int left = m_nodes[node].left;
    unsigned int right = m_nodes[node].right;

    if (m_nodes[node].isLeaf()) {
      assert(left == NULL_NODE);
      assert(right == NULL_NODE);
      assert(m_nodes[node].height == 0);
      return;
    }

    assert(left < m_node_capacity);
    assert(right < m_node_capacity);

    int height1 = m_nodes[left].height;
    int height2 = m_nodes[right].height;
    int height = 1 + std::max(height1, height2);
    (void)height;  // Unused variable in Release build
    assert(m_nodes[node].height == height);

    aabb aabb;
    aabb.merge(m_nodes[left].bb, m_nodes[right].bb);

    for (unsigned int i = 0; i < Dim; i++) {
      assert(aabb.lowerBound[i] == m_nodes[node].bb.lowerBound[i]);
      assert(aabb.upperBound[i] == m_nodes[node].bb.upperBound[i]);
    }

    validate_metrics(left);
    validate_metrics(right);
  }

  //! Apply periodic boundary conditions.
  /* \param position
          The position vector.
   */
  void periodic_boundaries(vec<double> &position)
  {
    for (unsigned int i = 0; i < Dim; i++) {
      if (position[i] < 0) {
        position[i] += m_box_size[i];
      }
      else {
        if (position[i] >= m_box_size[i]) {
          position[i] -= m_box_size[i];
        }
      }
    }
  }

  //! Compute minimum image separation.
  /*! \param separation
          The separation vector.

      \param shift
          The shift vector.

      \return
          Whether a periodic shift has been applied.
   */
  bool minimum_image(vec<double> &separation, vec<double> &shift) const
  {
    bool isShifted = false;

    for (unsigned int i = 0; i < Dim; i++) {
      if (separation[i] < m_neg_min_image[i]) {
        separation[i] += m_periodicity[i] * m_box_size[i];
        shift[i] = m_periodicity[i] * m_box_size[i];
        isShifted = true;
      }
      else {
        if (separation[i] >= m_pos_min_image[i]) {
          separation[i] -= m_periodicity[i] * m_box_size[i];
          shift[i] = -(m_periodicity[i] * m_box_size[i]);
          isShifted = true;
        }
      }
    }

    return isShifted;
  }
};

// Utility typedefs
#define TYPEDEFS(suffix, dim, type)       \
  using point##suffix = point<dim, type>; \
  using aabb##suffix = aabb<dim, type>;   \
  using tree##suffix = tree<dim, type>

TYPEDEFS(2d, 2, double);
TYPEDEFS(2f, 2, float);
TYPEDEFS(2i, 2, int);
TYPEDEFS(3d, 3, double);
TYPEDEFS(3f, 3, float);
TYPEDEFS(3i, 3, int);

#undef TYPEDEFS

}  // namespace abt

#endif /* _AABB_H */
