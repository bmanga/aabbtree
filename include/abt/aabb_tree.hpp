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

  std::array<value_type, Dim> values;
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
    surfaceArea = computeSurfaceArea();
    centre = computeCentre();
  }

  /// Compute the surface area of the box.
  value_type computeSurfaceArea() const
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
  value_type getSurfaceArea() const { return surfaceArea; }

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

    surfaceArea = computeSurfaceArea();
    centre = computeCentre();
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
  point computeCentre()
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
  value_type surfaceArea;
};

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
   particle, or a group of particles, in the simulation box. The AABB
   objects of individual particles are "fattened" before they are stored
   to avoid having to continually update and rebalance the tree when
   displacements are small.

   Nodes are aware of their position within in the tree. The isLeaf member
   function allows the tree to query whether the node is a leaf, i.e. to
   determine whether it holds a single particle.
  */
  struct node {
    /// Constructor.
    node() = default;

    /// The fattened axis-aligned bounding box.
    aabb bb;

    /// Index of the parent node.
    unsigned int parent;

    /// Index of the next node.
    unsigned int next;

    /// Index of the left-hand child.
    unsigned int left;

    /// Index of the right-hand child.
    unsigned int right;

    /// Height of the node. This is 0 for a leaf and -1 for a free node.
    int height;

    /// The index of the particle that the node contains (leaf nodes only).
    unsigned int particle;

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

      \param nParticles
          The number of particles (for fixed particle number systems).

      \param touchIsOverlap
          Does touching count as overlapping in query operations?
   */
  tree(value_type skin_thickness = 0,
       unsigned int nParticles = 16,
       bool touchIsOverlap = true)
      : m_is_periodic(false),
        m_skin_thickness(skin_thickness),
        m_touch_is_overlap(touchIsOverlap)
  {
    // Initialise the periodicity vector.
    std::fill(m_periodicity.begin(), m_periodicity.end(), false);

    // Initialise the tree.
    m_root = NULL_NODE;
    m_node_count = 0;
    m_node_capacity = nParticles;
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

      \param nParticles
          The number of particles (for fixed particle number systems).

      \param touchIsOverlap
          Does touching count as overlapping in query operations?
   */
  tree(double skinThickness,
       const vec<bool> &periodicity,
       const vec<double> &boxSize,
       unsigned int nParticles = 16,
       bool touchIsOverlap = true)
      : m_skin_thickness(skinThickness),
        m_periodicity(periodicity),
        m_box_size(boxSize),
        m_touch_is_overlap(touchIsOverlap)
  {
    // Initialise the tree.
    m_root = NULL_NODE;
    m_node_count = 0;
    m_node_capacity = nParticles;
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
  void setPeriodicity(const vec<bool> &periodicity)
  {
    m_periodicity = periodicity;
  }

  //! Set the size of the simulation box.
  /*! \param box_size
          The size of the simulation box in each dimension.
   */
  void setBoxSize(const vec<double> &box_size) { m_box_size = box_size; }

  //! Insert a particle into the tree (point particle).
  /*! \param index
          The index of the particle.

      \param position
          The position vector of the particle.

      \param radius
          The radius of the particle.
   */
  void insertParticle(unsigned int particle,
                      const point &position,
                      double radius)
  {
    // Make sure the particle doesn't already exist.
    if (m_particle_map.count(particle) != 0) {
      throw std::invalid_argument("[ERROR]: Particle already exists in tree!");
    }

    // Allocate a new node for the particle.
    unsigned int node = allocateNode();

    // AABB size in each dimension.
    vec<double> size;

    // Compute the AABB limits.
    for (unsigned int i = 0; i < Dim; i++) {
      m_nodes[node].bb.lowerBound[i] = position[i] - radius;
      m_nodes[node].bb.upperBound[i] = position[i] + radius;
      size[i] = m_nodes[node].bb.upperBound[i] - m_nodes[node].bb.lowerBound[i];
    }

    // Fatten the AABB.
    for (unsigned int i = 0; i < Dim; i++) {
      m_nodes[node].bb.lowerBound[i] -= m_skin_thickness * size[i];
      m_nodes[node].bb.upperBound[i] += m_skin_thickness * size[i];
    }
    m_nodes[node].bb.surfaceArea = m_nodes[node].bb.computeSurfaceArea();
    m_nodes[node].bb.centre = m_nodes[node].bb.computeCentre();

    // Zero the height.
    m_nodes[node].height = 0;

    // Insert a new leaf into the tree.
    insertLeaf(node);

    // Add the new particle to the map.
    m_particle_map.insert(
        std::unordered_map<unsigned int, unsigned int>::value_type(particle,
                                                                   node));

    // Store the particle index.
    m_nodes[node].particle = particle;
  }

  //! Insert a particle into the tree (arbitrary shape with bounding box).
  /*! \param index
          The index of the particle.

      \param lowerBound
          The lower bound in each dimension.

      \param upperBound
          The upper bound in each dimension.
   */
  void insertParticle(unsigned int particle,
                      const point &lowerBound,
                      const point &upperBound)
  {
    // Make sure the particle doesn't already exist.
    if (m_particle_map.count(particle) != 0) {
      throw std::invalid_argument("[ERROR]: Particle already exists in tree!");
    }

    // Allocate a new node for the particle.
    unsigned int node = allocateNode();

    // AABB size in each dimension.
    vec<double> size;

    // Compute the AABB limits.
    for (unsigned int i = 0; i < Dim; i++) {
      // Validate the bound.
      if (lowerBound[i] > upperBound[i]) {
        throw std::invalid_argument(
            "[ERROR]: AABB lower bound is greater than the upper bound!");
      }

      m_nodes[node].bb.lowerBound[i] = lowerBound[i];
      m_nodes[node].bb.upperBound[i] = upperBound[i];
      size[i] = upperBound[i] - lowerBound[i];
    }

    // Fatten the AABB.
    for (unsigned int i = 0; i < Dim; i++) {
      m_nodes[node].bb.lowerBound[i] -= m_skin_thickness * size[i];
      m_nodes[node].bb.upperBound[i] += m_skin_thickness * size[i];
    }
    m_nodes[node].bb.surfaceArea = m_nodes[node].bb.computeSurfaceArea();
    m_nodes[node].bb.centre = m_nodes[node].bb.computeCentre();

    // Zero the height.
    m_nodes[node].height = 0;

    // Insert a new leaf into the tree.
    insertLeaf(node);

    // Add the new particle to the map.
    m_particle_map.insert(
        std::unordered_map<unsigned int, unsigned int>::value_type(particle,
                                                                   node));

    // Store the particle index.
    m_nodes[node].particle = particle;
  }

  /// Return the number of particles in the tree.
  unsigned int nParticles() { return m_particle_map.size(); }

  //! Remove a particle from the tree.
  /*! \param particle
          The particle index (particleMap will be used to map the node).
   */
  void removeParticle(unsigned int particle)
  {
    // Map iterator.
    std::unordered_map<unsigned int, unsigned int>::iterator it;

    // Find the particle.
    it = m_particle_map.find(particle);

    // The particle doesn't exist.
    if (it == m_particle_map.end()) {
      throw std::invalid_argument("[ERROR]: Invalid particle index!");
    }

    // Extract the node index.
    unsigned int node = it->second;

    // Erase the particle from the map.
    m_particle_map.erase(it);

    assert(node < m_node_capacity);
    assert(m_nodes[node].isLeaf());

    removeLeaf(node);
    freeNode(node);
  }

  /// Remove all particles from the tree.
  void removeAll()
  {
    // Iterator pointing to the start of the particle map.
    std::unordered_map<unsigned int, unsigned int>::iterator it =
        m_particle_map.begin();

    // Iterate over the map.
    while (it != m_particle_map.end()) {
      // Extract the node index.
      unsigned int node = it->second;

      assert(node < m_node_capacity);
      assert(m_nodes[node].isLeaf());

      removeLeaf(node);
      freeNode(node);

      it++;
    }

    // Clear the particle map.
    m_particle_map.clear();
  }

  //! Update the tree if a particle moves outside its fattened AABB.
  /*! \param particle
          The particle index (particleMap will be used to map the node).

      \param position
          The position vector of the particle.

      \param radius
          The radius of the particle.

      \param alwaysReinsert
          Always reinsert the particle, even if it's within its old AABB
     (default:false)

      \return
          Whether the particle was reinserted.
   */
  bool updateParticle(unsigned int particle,
                      const point &position,
                      double radius,
                      bool alwaysReinsert = false)
  {
    // AABB bounds vectors.
    point lowerBound;
    point upperBound;

    // Compute the AABB limits.
    for (unsigned int i = 0; i < Dim; i++) {
      lowerBound[i] = position[i] - radius;
      upperBound[i] = position[i] + radius;
    }

    // Update the particle.
    return updateParticle(particle, lowerBound, upperBound, alwaysReinsert);
  }

  //! Update the tree if a particle moves outside its fattened AABB.
  /*! \param particle
          The particle index (particleMap will be used to map the node).

      \param lowerBound
          The lower bound in each dimension.

      \param upperBound
          The upper bound in each dimension.

      \param alwaysReinsert
          Always reinsert the particle, even if it's within its old AABB
     (default: false)
   */
  bool updateParticle(unsigned int particle,
                      const point &lowerBound,
                      const point &upperBound,
                      bool alwaysReinsert = false)
  {
    // Map iterator.
    std::unordered_map<unsigned int, unsigned int>::iterator it;

    // Find the particle.
    it = m_particle_map.find(particle);

    // The particle doesn't exist.
    if (it == m_particle_map.end()) {
      throw std::invalid_argument("[ERROR]: Invalid particle index!");
    }

    // Extract the node index.
    unsigned int node = it->second;

    assert(node < m_node_capacity);
    assert(m_nodes[node].isLeaf());

    // AABB size in each dimension.
    vec<double> size;

    // Compute the AABB limits.
    for (unsigned int i = 0; i < Dim; i++) {
      // Validate the bound.
      if (lowerBound[i] > upperBound[i]) {
        throw std::invalid_argument(
            "[ERROR]: AABB lower bound is greater than the upper bound!");
      }

      size[i] = upperBound[i] - lowerBound[i];
    }

    // Create the new AABB.
    aabb aabb(lowerBound, upperBound);

    // No need to update if the particle is still within its fattened AABB.
    if (!alwaysReinsert && m_nodes[node].bb.contains(aabb))
      return false;

    // Remove the current leaf.
    removeLeaf(node);

    // Fatten the new AABB.
    for (unsigned int i = 0; i < Dim; i++) {
      aabb.lowerBound[i] -= m_skin_thickness * size[i];
      aabb.upperBound[i] += m_skin_thickness * size[i];
    }

    // Assign the new AABB.
    m_nodes[node].bb = aabb;

    // Update the surface area and centroid.
    m_nodes[node].bb.surfaceArea = m_nodes[node].bb.computeSurfaceArea();
    m_nodes[node].bb.centre = m_nodes[node].bb.computeCentre();

    // Insert a new leaf node.
    insertLeaf(node);

    return true;
  }

  //! Query the tree to find candidate interactions for a particle.
  /*! \param particle
          The particle index.

      \return particles
          A vector of particle indices.
   */
  std::vector<unsigned int> query(unsigned int particle)
  {
    // Make sure that this is a valid particle.
    if (m_particle_map.count(particle) == 0) {
      throw std::invalid_argument("[ERROR]: Invalid particle index!");
    }

    // Test overlap of particle AABB against all other particles.
    return query(m_nodes[m_particle_map.find(particle)->second].bb);
  }

  //! Query the tree to find candidate interactions for an AABB.
  /*! \param aabb
          The AABB.

      \return particles
          A vector of particle indices.
   */

  template <class Query>
  std::vector<unsigned int> query(const Query &query) const
  {
    constexpr bool query_is_point = std::is_same_v<Query, point>;
    constexpr bool query_is_aabb = std::is_same_v<Query, aabb>;
    static_assert(query_is_point || query_is_aabb,
                  "Only point or aabb queries are supported");
    // Make sure the tree isn't empty.
    if (m_particle_map.size() == 0) {
      return std::vector<unsigned int>();
    }

    std::vector<unsigned int> stack;
    stack.reserve(256);
    stack.push_back(m_root);

    std::vector<unsigned int> particles;

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

        bool isShifted = minimumImage(separation, shift);

        // Shift the AABB.
        if (isShifted) {
          for (unsigned int i = 0; i < Dim; i++) {
            nodeAABB.lowerBound[i] += shift[i];
            nodeAABB.upperBound[i] += shift[i];
          }
        }
      }

      // Test for overlap between the AABBs.
      if (nodeAABB.overlaps(query, m_touch_is_overlap)) {
        // Check that we're at a leaf node.
        if (m_nodes[node].isLeaf()) {
          particles.push_back(m_nodes[node].particle);
        }
        else {
          stack.push_back(m_nodes[node].left);
          stack.push_back(m_nodes[node].right);
        }
      }
    }

    return particles;
  }

  //! Get a particle AABB.
  /*! \param particle
          The particle index.
   */
  const aabb &getAABB(unsigned int particle) const
  {
    return m_nodes[m_particle_map.at(particle)].bb;
  }

  //! Get the height of the tree.
  /*! \return
          The height of the binary tree.
   */
  unsigned int getHeight() const
  {
    if (m_root == NULL_NODE)
      return 0;
    return m_nodes[m_root].height;
  }

  //! Get the number of nodes in the tree.
  /*! \return
          The number of nodes in the tree.
   */
  unsigned int getNodeCount() const { return m_node_count; }

  //! Compute the maximum balancance of the tree.
  /*! \return
          The maximum difference between the height of two
          children of a node.
   */
  unsigned int computeMaximumBalance() const
  {
    unsigned int maxBalance = 0;
    for (unsigned int i = 0; i < m_node_capacity; i++) {
      if (m_nodes[i].height <= 1)
        continue;

      assert(m_nodes[i].isLeaf() == false);

      unsigned int balance = std::abs(m_nodes[m_nodes[i].left].height -
                                      m_nodes[m_nodes[i].right].height);
      maxBalance = std::max(maxBalance, balance);
    }

    return maxBalance;
  }

  void set_touch_is_overlap(bool is_overlap)
  {
    m_touch_is_overlap = is_overlap;
  }

  bool is_touch_overlap() const { return m_touch_is_overlap; }

  //! Compute the surface area ratio of the tree.
  /*! \return
          The ratio of the sum of the node surface area to the surface
          area of the root node.
   */
  double computeSurfaceAreaRatio() const
  {
    if (m_root == NULL_NODE)
      return 0.0;

    double rootArea = m_nodes[m_root].bb.computeSurfaceArea();
    double totalArea = 0.0;

    for (unsigned int i = 0; i < m_node_capacity; i++) {
      if (m_nodes[i].height < 0)
        continue;

      totalArea += m_nodes[i].bb.computeSurfaceArea();
    }

    return totalArea / rootArea;
  }

  /// Validate the tree.
  void validate() const
  {
#ifndef NDEBUG
    validateStructure(m_root);
    validateMetrics(m_root);

    unsigned int freeCount = 0;
    unsigned int freeIndex = m_free_list;

    while (freeIndex != NULL_NODE) {
      assert(freeIndex < m_node_capacity);
      freeIndex = m_nodes[freeIndex].next;
      freeCount++;
    }

    assert(getHeight() == computeHeight());
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
        freeNode(i);
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
          double cost = aabb.getSurfaceArea();

          if (cost < minCost) {
            iMin = i;
            jMin = j;
            minCost = cost;
          }
        }
      }

      unsigned int index1 = nodeIndices[iMin];
      unsigned int index2 = nodeIndices[jMin];

      unsigned int parent = allocateNode();
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

  /// A map between particle and node indices.
  std::unordered_map<unsigned int, unsigned int> m_particle_map;

  /// Does touching count as overlapping in tree queries?
  bool m_touch_is_overlap;

  //! Allocate a new node.
  /*! \return
          The index of the allocated node.
   */
  unsigned int allocateNode()
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
  void freeNode(unsigned int node)
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
  void insertLeaf(unsigned int leaf)
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

      double surfaceArea = m_nodes[index].bb.getSurfaceArea();

      aabb combinedAABB;
      combinedAABB.merge(m_nodes[index].bb, leafAABB);
      double combinedSurfaceArea = combinedAABB.getSurfaceArea();

      // Cost of creating a new parent for this node and the new leaf.
      double cost = 2.0 * combinedSurfaceArea;

      // Minimum cost of pushing the leaf further down the tree.
      double inheritanceCost = 2.0 * (combinedSurfaceArea - surfaceArea);

      // Cost of descending to the left.
      double costLeft;
      if (m_nodes[left].isLeaf()) {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[left].bb);
        costLeft = aabb.getSurfaceArea() + inheritanceCost;
      }
      else {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[left].bb);
        double oldArea = m_nodes[left].bb.getSurfaceArea();
        double newArea = aabb.getSurfaceArea();
        costLeft = (newArea - oldArea) + inheritanceCost;
      }

      // Cost of descending to the right.
      double costRight;
      if (m_nodes[right].isLeaf()) {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[right].bb);
        costRight = aabb.getSurfaceArea() + inheritanceCost;
      }
      else {
        aabb aabb;
        aabb.merge(leafAABB, m_nodes[right].bb);
        double oldArea = m_nodes[right].bb.getSurfaceArea();
        double newArea = aabb.getSurfaceArea();
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
    unsigned int newParent = allocateNode();
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
  void removeLeaf(unsigned int leaf)
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
      freeNode(parent);

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
      freeNode(parent);
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
  unsigned int computeHeight() const { return computeHeight(m_root); }

  //! Compute the height of a sub-tree.
  /*! \param node
          The index of the root node.

      \return
          The height of the sub-tree.
   */
  unsigned int computeHeight(unsigned int node) const
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
  void validateStructure(unsigned int node) const
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

    validateStructure(left);
    validateStructure(right);
  }

  //! Assert that the sub-tree has valid metrics.
  /*! \param node
          The index of the root node.
   */
  void validateMetrics(unsigned int node) const
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

    validateMetrics(left);
    validateMetrics(right);
  }

  //! Apply periodic boundary conditions.
  /* \param position
          The position vector.
   */
  void periodicBoundaries(vec<double> &position)
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
  bool minimumImage(vec<double> &separation, vec<double> &shift) const
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
          shift[i] = -m_periodicity[i] * m_box_size[i];
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
