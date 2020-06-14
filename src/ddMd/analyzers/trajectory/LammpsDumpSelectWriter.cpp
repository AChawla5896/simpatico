/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LammpsDumpSelectWriter.h"
#include <ddMd/simulation/Simulation.h>
#include <ddMd/communicate/AtomCollector.h>
#include <ddMd/chemistry/Atom.h>
#include <util/archives/BinaryFileOArchive.h>
#include <util/format/Dbl.h>
#include <util/space/Vector.h>

namespace DdMd
{

   using namespace Util;

   /*  
   * Constructor.
   */  
   LammpsDumpSelectWriter::LammpsDumpSelectWriter(Simulation& simulation)
    : TrajectoryWriter(simulation)
   {  setClassName("LammpsDumpSelectWriter"); }

   /*  
   * Destructor.
   */  
   LammpsDumpSelectWriter::~LammpsDumpSelectWriter()
   {}  

   void LammpsDumpSelectWriter::readParameters(std::istream& in) 
   {   
      readInterval(in);
      readOutputFileName(in);
      read<int>(in, "numberOfTypes", nTypes_);
      desiredTypes_.allocate(nTypes_);
      readDArray<int>(in,"typesDesired", desiredTypes_, nTypes_);
      isInitialized_ = true;
   }   
   //Special Setup of lammpsSelect, counts the true number of atoms, and relabels them 
   void LammpsDumpSelectWriter::selectAtomSetup()
   {  
      atomStorage().computeNAtomTotal(domain().communicator());
      if (domain().isMaster()) {
         nAtom_ = atomStorage().nAtomTotal();
      }

      //Count up the total atoms in the desired subset
      if (domain().isMaster()) {
         atomSubsetIds_.allocate(nAtom_);
         int typeId;
         atomCollector().setup();
         // Count the atom subset
         AtomSubset_=0;
         Atom* batomPtr = atomCollector().nextPtr();
         while (batomPtr) {
            typeId = batomPtr->typeId();
            for (int i = 0; i < nTypes_; ++i)  {
               if (desiredTypes_[i] == typeId ) {
                  AtomSubset_ += 1;
               }
            }
         batomPtr=atomCollector().nextPtr();
         }
      } else {
       atomCollector().send();
      }

      // Now construct the array that converts real ids and subset ids
      if (domain().isMaster()) {
         int typeId;
         atomCollector().setup();
         Atom* catomPtr = atomCollector().nextPtr();
         DArray<int> sortingArray;
         sortingArray.allocate(AtomSubset_);
         int place=0;
         bool isDesired;
         //Get all the id's of all desired atoms
         while (catomPtr) {
            isDesired = false;
            typeId = catomPtr->typeId();
            for (int i = 0; i < nTypes_; ++i)  {
               if (typeId == desiredTypes_[i]) {
                  isDesired = true;
               }
            }
            if (isDesired) {
               sortingArray[place]=catomPtr->id();
               place += 1;
            } else {
               atomSubsetIds_[catomPtr->id()] = 0;
            }
         catomPtr=atomCollector().nextPtr();
         }
         //Sort all of the atoms and put them into the atomSubsetIds_ arrays
         int newId = 1;
         for (int i = 0; i < nAtom_; ++i) {
            for (int j = 0; j < AtomSubset_  ; ++j) {
               if (i == sortingArray[j]) {
                  atomSubsetIds_[i] = newId;
                  newId += 1;
               }
            }
         }
      } else {
         atomCollector().send();
      }
   }

   /*
   *  Write a configuration snapshot.
   */
   void LammpsDumpSelectWriter::writeFrame(std::ofstream &file, long iStep)
   {
      // Compute number of atoms in desired frames
      //atomStorage().computeNAtomTotal(domain().communicator());
      if (iStep == 0){
         selectAtomSetup();
      }
      int id;
      int typeId;
      if (domain().isMaster()) {
         Vector r;
         int id;
         int n = 0;
         bool isCartesian = atomStorage().isCartesian();
         file << "ITEM: TIMESTEP" << "\n";
         file << iStep << "\n";
         file << "ITEM: NUMBER OF ATOMS" << "\n";
         file << AtomSubset_ << "\n";
         file << "ITEM: BOX BOUNDS pp pp pp" << "\n";
         Vector lengths = boundary().lengths();
         file << Dbl(0.0) << Dbl(lengths[0]) << "\n";
      file << Dbl(0.0) << Dbl(lengths[1]) << "\n";
         file << Dbl(0.0) << Dbl(lengths[2]) << "\n";
         int shift = 0;
         int molId = 1;
         file << "ITEM: ATOMS id type mol x y z" << "\n";
         atomCollector().setup();
         Atom* atomPtr = atomCollector().nextPtr();
         while(atomPtr) {
            typeId = atomPtr->typeId();
            for (int i = 0; i < nTypes_; ++i)  {
               if (desiredTypes_[i] == typeId ) {
                  if (isCartesian) {
                     r = atomPtr->position(); }
                  else {
                     boundary().transformGenToCart(atomPtr->position(), r);
                  }
                  file << atomSubsetIds_[atomPtr->id()]  << " ";
                  file << typeId + 1 << " ";
                  file << molId << " ";
                  for (int i=0; i < Util::Dimension; ++i) {
                     file << Dbl(atomPtr->position()[i],13) << " ";
                  }
                  for (int i=0; i < Util::Dimension; ++i) {
                     file << shift << " ";
                  }
                  file << "\n";
               }
            }
         atomPtr = atomCollector().nextPtr();
         }
      } else {
        atomCollector().send();
      }
   }
}
