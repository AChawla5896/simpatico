/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "PairPotential.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/storage/AtomStorage.h>
#include <ddMd/storage/AtomIterator.h>
#include <ddMd/storage/GhostIterator.h>
#include <ddMd/neighbor/PairIterator.h>
#include <ddMd/neighbor/PairList.h>
#include <ddMd/neighbor/CellList.h>
#include <ddMd/communicate/Domain.h>
#include <util/mpi/MpiLoader.h>
#include <util/space/Vector.h>
#include <util/global.h>

namespace DdMd
{
   using namespace Util;

   /*
   * Default constructor (for unit testing).
   */
   PairPotential::PairPotential()
    : skin_(0.0),
      cutoff_(0.0),
      pairCapacity_(0),
      domainPtr_(0),
      boundaryPtr_(0),
      storagePtr_(0),
      pairEnergies_(),
      methodId_(0),
      nPair_(0),
      hasMaxBoundary_(false)
   {  setClassName("PairPotential"); } 

   /*
   * Constructor.
   */
   PairPotential::PairPotential(Simulation& simulation)
    : skin_(0.0),
      cutoff_(0.0),
      pairCapacity_(0),
      domainPtr_(&simulation.domain()),
      boundaryPtr_(&simulation.boundary()),
      storagePtr_(&simulation.atomStorage()),
      pairEnergies_(),
      methodId_(0),
      nPair_(0),
      hasMaxBoundary_(false)
   {  setClassName("PairPotential"); } 

   /*
   * Associate with related objects. (for unit testing).
   */
   void PairPotential::associate(Domain& domain, Boundary& boundary, 
                                 AtomStorage& storage)
   {
      domainPtr_ = &domain;
      boundaryPtr_ = &boundary;
      storagePtr_ = &storage;
   } 

   /*
   * Destructor.
   */
   PairPotential::~PairPotential()
   {}

   /*
   * Allocate memory for the cell list.
   */
   void PairPotential::initialize(const Boundary& maxBoundary, double skin, 
                                  int pairCapacity)
   {
      maxBoundary_ = maxBoundary;
      skin_ = skin;
      pairCapacity_ = pairCapacity;
      cutoff_ = maxPairCutoff() + skin;
      hasMaxBoundary_ = true;
      pairList_.setCutoff(cutoff_);
      reserve();
   }

   /*
   * Read parameters for PairList and allocate memory.  
   */
   void PairPotential::readParameters(std::istream& in)
   {
      // Read skin
      read<double>(in, "skin", skin_);
      cutoff_ = maxPairCutoff() + skin_;
      pairList_.setCutoff(cutoff_);

      // Optionally read nCellCut
      nCellCut_ = 1; // Default value for optional parameter
      readOptional<int>(in, "nCellCut", nCellCut_); 

      // Optionally read pairCapacity
      pairCapacity_ = 0;
      readOptional<int>(in, "pairCapacity", pairCapacity_);

      // Optionally read maxBoundary
      hasMaxBoundary_ = 
        readOptional<Boundary>(in, "maxBoundary", maxBoundary_).isActive();

      // Reserve memory for cell and pair lists.
      reserve();
   }

   /*
   * Read parameters for PairList and allocate memory.  
   */
   void PairPotential::loadParameters(Serializable::IArchive& ar)
   {
      MpiLoader<Serializable::IArchive> loader(*this, ar);
  
      loadParameter<double>(ar, "skin", skin_);
      loadParameter<int>(ar, "nCellCut", nCellCut_, false);
      loadParameter<int>(ar, "pairCapacity", pairCapacity_);
      loader.load(hasMaxBoundary_);
      if (hasMaxBoundary_) {
         loadParameter<Boundary>(ar, "maxBoundary", maxBoundary_);
      }
      loader.load(cutoff_);
      loader.load(methodId_);

      pairList_.setCutoff(cutoff_);
      reserve();
   }

   /*
   * Save parameters to output/saving Archive.
   */
   void PairPotential::save(Serializable::OArchive& ar)
   {
      ar << skin_;
      Parameter::saveOptional(ar, nCellCut_, true);
      ar << pairCapacity_;
      ar << hasMaxBoundary_;
      if (hasMaxBoundary_) {
         ar << maxBoundary_;
      }
      ar << cutoff_;
      ar << methodId_;
   }

   /*
   * Reserve memory for the cell list and pair list.
   */
   void PairPotential::reserve()
   {

      // Set CellList atomCapacity
      int localCapacity = storage().atomCapacity();
      int totalCapacity = localCapacity + storage().ghostCapacity();
      cellList_.setAtomCapacity(totalCapacity);

      // Optionally make cell list grid
      if (hasMaxBoundary_) {
         UTIL_CHECK(cutoff_ > 0.0);
         Vector cutoffs;
         Vector lower;
         Vector upper;
         for (int i = 0; i < Dimension; ++i) {
            lower[i] = domain().domainBound(i, 0);
            upper[i] = domain().domainBound(i, 1);
            cutoffs[i] = cutoff_/maxBoundary_.length(i);
         }
         cellList_.makeGrid(lower, upper, cutoffs, nCellCut_);
      }

      // Optionally reserve memory for PairList
      if (pairCapacity_ > 0) {
         pairList_.reserve(localCapacity, pairCapacity_);
      }

   }

   /*
   * Build the cell list.
   */
   void PairPotential::buildCellList()
   {
      if (storage().isCartesian()) {
         UTIL_THROW("Coordinates are Cartesian entering buildCellList");
      }

      // Check CellList atomCapacity - resize if necessary.
      int totalCapacity = storage().atomCapacity() 
                        + storage().ghostCapacity();
      if (totalCapacity > cellList().atomCapacity()) {
         cellList().setAtomCapacity(totalCapacity);
      }

      // Set cutoff and domain bounds.
      Vector cutoffs;
      Vector lower;
      Vector upper;
      for (int i = 0; i < Dimension; ++i) {
         cutoffs[i] = cutoff_/boundaryPtr_->length(i);
         lower[i] = domain().domainBound(i, 0);
         upper[i] = domain().domainBound(i, 1);
      }

      // Make and clear cell list grid.
      cellList_.makeGrid(lower, upper, cutoffs, nCellCut_);
      cellList_.clear();

      // Compute cell indices for all local atoms.
      AtomIterator atomIter;
      storage().begin(atomIter);
      for ( ; atomIter.notEnd(); ++atomIter) {
         cellList_.placeAtom(*atomIter);
      }

      // Compute cell indices for all ghost atoms.
      GhostIterator ghostIter;
      storage().begin(ghostIter);
      for ( ; ghostIter.notEnd(); ++ghostIter) {
         cellList_.placeAtom(*ghostIter);
      }

      // Build cell list
      cellList_.build();
     
      // Postconditions
      assert(cellList_.isValid());
      assert(cellList_.nAtom() + cellList_.nReject() 
             == storage().nAtom() + storage().nGhost());
      if (storage().isCartesian()) {
         UTIL_THROW("Coordinates are Cartesian exiting buildCellList");
      }
   }

   /*
   * Build the pair list.
   */
   void PairPotential::buildPairList()
   {
      // Precondition
      UTIL_CHECK(storage().isCartesian());
      UTIL_CHECK(cellList_.isBuilt());
      pairList_.setCutoff(cutoff_);

      // Check arrays, resize if necessary.
      int localCapacity = storage().atomCapacity();
      if (pairList_.atomCapacity() == 0 || pairList_.pairCapacity() == 0) {
         int nPair = pairList_.countPairs(cellList_, reverseUpdateFlag());
         if (nPair > pairList_.atomCapacity()) {
            pairCapacity_ = 1.5*nPair;
         }
         pairList_.reserve(localCapacity, pairCapacity_);
      } else 
      if (localCapacity > pairList_.atomCapacity()) {
         UTIL_CHECK(pairCapacity_ == pairList_.pairCapacity()); 
         pairList_.reserve(localCapacity, pairCapacity_);
      }

      // Build the pair list.
      pairList_.build(cellList_, reverseUpdateFlag());
   }

   /*
   * Return value of pair energies.
   */
   DMatrix<double> PairPotential::pairEnergies() const
   {  return pairEnergies_.value(); }

   /*
   * Set a value for pair energies.
   */
   void PairPotential::setPairEnergies(DMatrix<double> pairEnergies)
   {  pairEnergies_.set(pairEnergies); }

   /*
   * Mark pair energies as unknown. 
   */
   void PairPotential::unsetPairEnergies()
   {  pairEnergies_.unset(); }

   /*
   * Compute total pair nPair on all processors.
   */
   #ifdef UTIL_MPI
   void PairPotential::computeNPair(MPI::Intracomm& communicator)
   #else
   void PairPotential::computeNPair()
   #endif
   { 
      double cut = maxPairCutoff();
      double cutSq = cut*cut;
      int localNPair = 0; 
      if (methodId() == 0) {
         localNPair = nPairList(cutSq); 
      } else 
      if (methodId() == 1) {
         localNPair = nPairCell(cutSq); 
      } else {
         localNPair = nPairNSq(cutSq); 
      }
      nPair_ = localNPair;
      #ifdef UTIL_MPI
      communicator.Reduce(&localNPair, &nPair_, 1, MPI::INT, MPI::SUM, 0);
      #else
      nPair_ = localNPair;
      #endif
   }

   /*
   * Return twice the number of pairs within the specified cutoff.
   * 
   * This method should only be called on the rank 0 processor. The
   * return value is computed by a previous call to computeNPair.
   */
   int PairPotential::nPair() const
   {  return nPair_; }

   /*
   * Count number pairs using pair list (private).
   */
   int PairPotential::nPairList(double cutoffSq)
   {
      Vector f;
      double rsq;
      PairIterator iter;
      Atom*  atom0Ptr;
      Atom*  atom1Ptr;
      int count = 0; // atom counter
      if (reverseUpdateFlag()) {
         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            assert(!atom0Ptr->isGhost());
            f.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = f.square();
            if (rsq < cutoffSq) {
               count += 2;
            }
         }
      } else {
         for (pairList_.begin(iter); iter.notEnd(); ++iter) {
            iter.getPair(atom0Ptr, atom1Ptr);
            assert(!atom0Ptr->isGhost());
            f.subtract(atom0Ptr->position(), atom1Ptr->position());
            rsq = f.square();
            if (rsq < cutoffSq) {
               if (!atom1Ptr->isGhost()) {
                  count += 2;
               } else {
                  count += 1;
               }
            }
         }
      }
      return count;
   }

   /*
   * Count number pairs using cell list (private).
   */
   int PairPotential::nPairCell(double cutoffSq)
   {
      // Find all neighbors (cell list)
      Cell::NeighborArray neighbors;
      Vector f;
      double rsq;
      Atom*  atomPtr0;
      Atom*  atomPtr1;
      const Cell*  cellPtr;
      int na, nn, i, j;
      int count = 0;

      // Iterate over linked list of local cells.
      cellPtr = cellList_.begin();
      while (cellPtr) {
         cellPtr->getNeighbors(neighbors, reverseUpdateFlag());
         na = cellPtr->nAtom();
         nn = neighbors.size();
         for (i = 0; i < na; ++i) {
            atomPtr0 = neighbors[i]->ptr();

            // Loop over atoms in this cell
            for (j = 0; j < na; ++j) {
               atomPtr1 = neighbors[j]->ptr();
               if (atomPtr1 > atomPtr0) {
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  if (rsq < cutoffSq) {
                     count += 2;
                  }
               }
            }

            // Loop over atoms in neighboring cells.
            if (reverseUpdateFlag()) {
               for (j = na; j < nn; ++j) {
                  atomPtr1 = neighbors[j]->ptr();
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  if (rsq < cutoffSq) {
                     count += 2;
                  }
               }
            } else {
               for (j = na; j < nn; ++j) {
                  atomPtr1 = neighbors[j]->ptr();
                  f.subtract(atomPtr0->position(), atomPtr1->position());
                  rsq = f.square();
                  if (rsq < cutoffSq) {
                     if (!atomPtr1->isGhost()) {
                        count += 2;
                     } else {
                        count += 1;
                     }
                  }
               }
            }

         }
         cellPtr = cellPtr->nextCellPtr();
      } // while (cellPtr) 
      return count;
   }

   /*
   * Count pairs using n-squared loop. (private)
   */
   int PairPotential::nPairNSq(double cutoffSq)
   {
      Vector f;
      double rsq;
      AtomIterator  atomIter0, atomIter1;
      GhostIterator ghostIter;
      int id0, id1;
      int count = 0;

      // Iterate over local atom 0
      storage().begin(atomIter0);
      for ( ; atomIter0.notEnd(); ++atomIter0) {
         id0 = atomIter0->id();

         // Iterate over local atom 1
         storage().begin(atomIter1);
         for ( ; atomIter1.notEnd(); ++atomIter1) {
            id1 = atomIter1->id();
            if (id0 < id1) {
               if (!atomIter0->mask().isMasked(id1)) {
                  f.subtract(atomIter0->position(), atomIter1->position());
                  rsq = f.square();
                  if (rsq < cutoffSq) {
                     count += 2;
                  }
               }
            }
         }

         // Iterate over ghost atoms
         storage().begin(ghostIter);
         if (reverseUpdateFlag()) {
            for ( ; ghostIter.notEnd(); ++ghostIter) {
               id1 = ghostIter->id();
               if (id0 < id1) {
                  if (!atomIter0->mask().isMasked(id1)) {
                     f.subtract(atomIter0->position(), ghostIter->position());
                     rsq = f.square();
                     if (rsq < cutoffSq) {
                        count += 2;
                     }
                  }
               }
            }
         } else {
            for ( ; ghostIter.notEnd(); ++ghostIter) {
               id1 = ghostIter->id();
               if (!atomIter0->mask().isMasked(id1)) {
                  f.subtract(atomIter0->position(), ghostIter->position());
                  rsq = f.square();
                  if (rsq < cutoffSq) {
                     count += 1;
                  }
               }
            }
         }

      }
      return count;
   }

}
