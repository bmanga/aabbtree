/*
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
*/

#include <fstream>
#include <iostream>

#include "abt/aabb_tree.hpp"
#include "MersenneTwister.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

using namespace abt;
template <class T>
using vec = abt::tree2d::vec<T>;

/*! \file hard_disc.cpp

  An example showing the use of AABB trees for simulating the dynamics
  of a binary hard disc system where there is a large size asymmetry (10:1)
  between the particle species.
*/

// FUNCTION DEFINITIONS
void minimumImage(vec<double> &separation,
                  const vec<bool> &periodicity,
                  const vec<double> &boxSize)
{
  for (unsigned int i = 0; i < 2; i++) {
    if (separation[i] < -0.5 * boxSize[i]) {
      separation[i] += periodicity[i] * boxSize[i];
    }
    else {
      if (separation[i] >= 0.5 * boxSize[i]) {
        separation[i] -= periodicity[i] * boxSize[i];
      }
    }
  }
}

bool overlaps(const point2d &position1,
              const point2d &position2,
              const vec<bool> &periodicity,
              const vec<double> &boxSize,
              double cutOff)
{
  // Calculate particle separation.
  vec<double> separation;
  separation[0] = (position1[0] - position2[0]);
  separation[1] = (position1[1] - position2[1]);

  // Calculate minimum image separation.
  minimumImage(separation, periodicity, boxSize);

  double rSqd = separation[0] * separation[0] + separation[1] * separation[1];

  if (rSqd < cutOff)
    return true;
  else
    return false;
}

void periodicBoundaries(point2d &position,
                        const vec<bool> &periodicity,
                        const vec<double> &boxSize)
{
  for (unsigned int i = 0; i < 2; i++) {
    if (position[i] < 0) {
      position[i] += periodicity[i] * boxSize[i];
    }
    else {
      if (position[i] >= boxSize[i]) {
        position[i] -= periodicity[i] * boxSize[i];
      }
    }
  }
}

void printVMD(const std::string &fileName,
              const std::vector<point2d> &positionsSmall,
              const std::vector<point2d> &positionsLarge)
{
  FILE *pFile;
  pFile = fopen(fileName.c_str(), "a");

  fprintf(pFile, "%lu\n\n", positionsSmall.size() + positionsLarge.size());
  for (unsigned int i = 0; i < positionsSmall.size(); i++)
    fprintf(pFile, "0 %lf %lf 0\n", positionsSmall[i][0], positionsSmall[i][1]);
  for (unsigned int i = 0; i < positionsLarge.size(); i++)
    fprintf(pFile, "1 %lf %lf 0\n", positionsLarge[i][0], positionsLarge[i][1]);

  fclose(pFile);
}


// MAIN FUNCTION

int main(int argc, char **argv)
{
  // Print git commit info, if present.
#ifdef COMMIT
  std::cout << "Git commit: " << COMMIT << "\n";
#endif

  // Print git branch info, if present.
#ifdef BRANCH
  std::cout << "Git branch: " << BRANCH << "\n";
#endif

  /*****************************************************************/
  /*      Set parameters, initialise variables and objects.        */
  /*****************************************************************/

  unsigned int nSweeps = 100000;      // The number of Monte Carlo sweeps.
  unsigned int sampleInterval = 100;  // The number of sweeps per sample.
  unsigned int nSmall = 1000;         // The number of small particles.
  unsigned int nLarge = 100;          // The number of large particles.
  double diameterSmall = 1;           // The diameter of the small particles.
  double diameterLarge = 10;          // The diameter of the large particles.
  double density = 0.1;               // The system density.
  double maxDisp = 0.1;  // Maximum trial displacement (in units of diameter).

  // Total particles.
  unsigned int nParticles = nSmall + nLarge;

  // Number of samples.
  unsigned int nSamples = nSweeps / sampleInterval;

  // Particle radii.
  double radiusSmall = 0.5 * diameterSmall;
  double radiusLarge = 0.5 * diameterLarge;

  // Output formatting flag.
  unsigned int format = std::floor(std::log10(nSamples));

  // Work out base length of simulation box.
  double baseLength =
      std::pow((M_PI * (nSmall * diameterSmall + nLarge * diameterLarge)) /
                   (4.0 * density),
               1.0 / 2.0);
  vec<double> boxSize({baseLength, baseLength});
  vec<bool> periodicity{true, true};

  // Initialise the random number generator.
  MersenneTwister rng;

  // Initialise the AABB trees.
  abt::tree<2, double> treeSmall(nSmall);
  abt::tree<2, double> treeLarge(nLarge);


    // Initialise particle position vectors.
  std::unordered_map<unsigned, point2d> positionsSmall(nSmall);
  std::unordered_map<unsigned, point2d> positionsLarge(nLarge);
  /*****************************************************************/
  /*             Generate the initial AABB trees.                  */
  /*****************************************************************/

  // First the large particles.

  std::cout << "\nInserting large particles into AABB tree ...\n";
  for (unsigned int i = 0; i < nLarge; i++) {
    // Initialise the particle position vector.
    point2d position;
    aabb2d aabb;

    // Insert the first particle directly.
    if (i == 0) {
      // Generate a random particle position.
      position[0] = boxSize[0] * rng();
      position[1] = boxSize[1] * rng();
      aabb = fattened(aabb2d::of_sphere(position, radiusLarge), maxDisp);
    }

    // Check for overlaps.
    else {
      // Keep trying until there is no overlap.
      do {
        position[0] = boxSize[0] * rng();
        position[1] = boxSize[1] * rng();

        aabb = fattened(abt::aabb2d::of_sphere(position, radiusLarge), maxDisp);
      } while (treeLarge.any_overlap(
        aabb, [&](unsigned id, const aabb2d& overlap) {
          // Cut-off distance.
          double cutOff = 2.0 * radiusLarge;
          cutOff *= cutOff;

          // Particles overlap.
          return overlaps(position, overlap.centre, periodicity, boxSize, cutOff);
        }, true, boxSize));
    }

    // Insert the particle into the tree.
    treeLarge.insert(aabb);
  }
  std::cout << "Tree generated!\n";

  // Now fill the gaps with small particles.

  std::cout << "\nInserting small particles into AABB tree ...\n";
  for (unsigned int i = 0; i < nSmall; i++) {
    // Initialise the particle position vector.
    point2d position;
    aabb2d aabb;

    do {
      position[0] = boxSize[0] * rng();
      position[1] = boxSize[1] * rng();
      aabb = fattened(aabb2d::of_sphere(position, radiusSmall), maxDisp);
    } while (
        treeLarge.any_overlap(aabb,
                              [&](unsigned id, const aabb2d &overlap) {
                                // Cut-off distance.
                                double cutOff = radiusSmall + radiusLarge;
                                cutOff *= cutOff;

                                // Particles overlap.
                                return overlaps(position, overlap.centre,
                                                periodicity, boxSize, cutOff);
            },
            true, boxSize) ||
        treeSmall.any_overlap(aabb,
                              [&](unsigned id, const aabb2d &overlap) {
                                // Cut-off distance.
                                double cutOff = 2.0 * radiusSmall;
                                cutOff *= cutOff;

                                // Particles overlap.
                                return overlaps(position, overlap.centre,
                                                periodicity, boxSize, cutOff);
                 },
                 true, boxSize)

    );
    // Insert particle into tree.
    treeSmall.insert(aabb);
  }
  std::cout << "Tree generated!\n";

  /*****************************************************************/
  /*      Perform the dynamics, updating the tree as we go.        */
  /*****************************************************************/

  // Clear the trajectory file.
  FILE *pFile;
  pFile = fopen("trajectory.xyz", "w");
  fclose(pFile);

  unsigned int sampleFlag = 0;
  unsigned int nSampled = 0;

  std::vector<unsigned> largeParticle_ids;
  std::vector<unsigned> smallParticle_ids;

  treeLarge.for_each([&](unsigned id, const auto &) { largeParticle_ids.push_back(id); });
  treeSmall.for_each([&](unsigned id, const auto &) { smallParticle_ids.push_back(id); });

  std::cout << "\nRunning dynamics ...\n";
  for (unsigned int i = 0; i < nSweeps; i++) {
    for (unsigned int j = 0; j < nParticles; j++) {
      // Choose a random particle.
      unsigned int particle = rng.integer(0, nParticles - 1);

      // Determine the particle type.
      unsigned int particleType = (particle < nSmall) ? 0 : 1;

      // Determine the radius of the particle.
      double radius = (particleType == 0) ? radiusSmall : radiusLarge;

      if (particleType == 1)
        particle -= nSmall;

      // Initialise vectors.
      vec<double> displacement;
      point2d position;

      // Calculate the new particle position and displacement
      if (particleType == 0) {
        particle = smallParticle_ids[particle];
        displacement[0] = maxDisp * diameterSmall * (2.0 * rng() - 1.0);
        displacement[1] = maxDisp * diameterSmall * (2.0 * rng() - 1.0);
        position = treeSmall.get_aabb(particle).centre + displacement;
      }
      else {
        particle = largeParticle_ids[particle];
        displacement[0] = maxDisp * diameterLarge * (2.0 * rng() - 1.0);
        displacement[1] = maxDisp * diameterLarge * (2.0 * rng() - 1.0);
        position = treeLarge.get_aabb(particle).centre + displacement;
      }

      // Apply periodic boundary conditions.
      periodicBoundaries(position, periodicity, boxSize);

      // Generate the AABB.
      auto aabb = fattened(aabb2d::of_sphere(position, radius), maxDisp);

      if (!treeLarge.any_overlap(aabb,
                                [&](unsigned id, const aabb2d &overlap) {
                                  // Cut-off distance.
                                  double cutOff = radius + radiusLarge;
                                  cutOff *= cutOff;

                                  if (id == particle) {
                                    // Self overlap, ignore.
                                    return false;
                                  }

                                  // Particles overlap.
                                  return overlaps(position, overlap.centre,
                                                  periodicity, boxSize, cutOff);
              },
              true, boxSize) &&
          !treeSmall.any_overlap(aabb,
                                [&](unsigned id, const aabb2d &overlap) {
                                  // Cut-off distance.
                                  double cutOff = radius + radiusSmall;
                                  cutOff *= cutOff;

                                  if (id == particle) {
                                    // Self overlap, ignore.
                                    return false;
                                  }

                                  // Particles overlap.
                                  return overlaps(position, overlap.centre,
                                                  periodicity, boxSize, cutOff);
              },
              true, boxSize)

      ) {
        // Accept the move.
        if (particleType == 0) {
          positionsSmall[particle] = position;
          treeSmall.update(particle, aabb);
        }
        else {
          positionsLarge[particle] = position;
          treeLarge.update(particle, aabb);
        }
      }
    }

    sampleFlag++;

    if (sampleFlag == sampleInterval) {
      sampleFlag = 0;
      nSampled++;

      std::vector<point2d> smallVec;
      for (unsigned x = 0; x < positionsSmall.size(); ++x) {
        smallVec.push_back(positionsSmall[x]);
      }

      std::vector<point2d> largeVec;
      for (unsigned x = 0; x < positionsLarge.size(); ++x) {
        largeVec.push_back(positionsLarge[x]);
      }

      printVMD("trajectory.xyz", smallVec, largeVec);

      if (format == 1)
        printf("Saved configuration %2d of %2d\n", nSampled, nSamples);
      else if (format == 2)
        printf("Saved configuration %3d of %3d\n", nSampled, nSamples);
      else if (format == 3)
        printf("Saved configuration %4d of %4d\n", nSampled, nSamples);
      else if (format == 4)
        printf("Saved configuration %5d of %5d\n", nSampled, nSamples);
      else if (format == 5)
        printf("Saved configuration %6d of %6d\n", nSampled, nSamples);
    }
  }

  std::cout << "Done!\n";

  return (EXIT_SUCCESS);
}