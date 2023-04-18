#ifndef DENDRO_KT_FILTER_FUNCTION_H
#define DENDRO_KT_FILTER_FUNCTION_H

#include <functional>

namespace ibm
{

/**
 * In collaboration with Baskar's group, 2020-07-23.
 *
 * Convention: Want to carve out and discard a coarse representation of "In,"
 *             while retaining the coarse representation of "Out,"
 *             and be able to identify boundary elements and nodes.
 *
 * All of space is partitioned into two complementary and disjoint sets, called "In" and "Out".
 *
 * So that all boundary nodes lie in "In," the "In" set must be closed and the "Out" set must be open.
 *
 * The filter function f(const double * p, double sz) takes in a cube "c" having zero or positive side length.
 * f returns In if the closure of "c" is contained in "In".
 * f returns Out if the closure of "c" is contained in "Out".
 * f returns Intercepted otherwise.
 *
 *
 * Properties needed for algorithm:
 * Tree:
 *   - The union (finite) of the closures of all In octants is contained in the "In" set.
 *
 *   - Hence the interior of the tree formed from all Out or Intercepted octants contains the "Out" set.
 *
 *   - The tree is minimal, for if any octant could be removed, that octant would be In .
 *
 * Nodes:
 *   - Since the interior of the tree contains the "Out" set,
 *     all nodes on the boundary of the tree are contained
 *     in the "In" set (->no false negatives).
 *
 *   - All nodes contained in the "In" set cannot be on the boundary of an Out octant
 *     (->the Intercepted octants are the boundary octants).
 *
 *
 *
 * Then @maksbh's algorithm can be proved correct.
 *
 * Algorithm:
 * 1. First pass mark the elements as In, Out and Intercepted.
 *    Instead of returning true/false, we will return the flags: In, Out or Intercepted.
 *    The Intercepted elements are marked as boundary elements.
 *    The In elements are thrown away.
 *
 * 2. Within the Intercepted elements loop over each nodes with sz = 0
 *    (this will allow to check each nodes separately, you need to adjust coords though)
 *    to find In/Out nodes. All In marked nodes are boundary nodes.
 *
 */

  enum Partition { IN, OUT, INTERCEPTED };

  typedef std::function<Partition(const double *elemPhysCoords, double elemPhysSize)> DomainDecider;
}

#endif//DENDRO_KT_FILTER_FUNCTION_H
