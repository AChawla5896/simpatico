#ifndef MCMD_END_SWAP_MOVE_H
#define MCMD_END_SWAP_MOVE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>   // base class
#include <util/containers/DArray.h>    // member template
#include <util/space/Vector.h>          // member template parameter

namespace McMd
{

   using namespace Util;

   class McSystem;

   /**
   * A move that swaps the ends of a linear hetero-polymer.
   *
   * This move changes the energy if of a linear heteropolymer by
   * reversing the sequence of atom types. When applied to a diblock
   * copolymer, it swaps the A and B blocks.
   *
   * Limitation: This implementation only works with flexible chains,
   * with no angle or dihedral potentials.
   *
   * \sa \ref mcMd_mcMove_EndSwapMove_page "parameter file format"
   *
   * \ingroup McMd_McMove_Module
   */
   class EndSwapMove : public SystemMove
   {

   public:

      /**
      * Constructor.
      */
      EndSwapMove(McSystem& system);

      /**
      * Read species to which displacement is applied.
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Generate and accept or reject a move
      */
      virtual bool move();

   protected:

      /// Integer index for molecular species.
      int speciesId_;

      /// Array of atom type indices
      DArray<int> atomTypeIds_;

      /// Array of atomic positions (temporary).
      DArray<Vector> positions_;

   };

}
#endif
