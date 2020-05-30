#ifndef DDMD_LAMMPS_DUMP_SELECT:_WRITER_H
#define DDMD_LAMMPS_DUMP_SELECT_WRITER_H

/*
 * Simpatico - Simulation Package for Polymeric and Molecular Liquids
 *
 * Copyright 2010 - 2014, The Regents of the University of Minnesota
 * Distributed under the terms of the GNU General Public License.
 */

#include <ddMd/analyzers/trajectory/TrajectoryWriter.h>   // base class

namespace DdMd
{

   using namespace Util;

   /**
   * Write a trajectory in the Lammps dump format.
   *
   * \ingroup McMd_TrajectoryWriter_Module
   */
   class LammpsDumpSelectWriter : public TrajectoryWriter
   {

   public:

      /**
      * Constructor.
      *
      * \param simulation parent Simulation object
      */
      LammpsDumpSelectWriter(Simulation& simulation);

      /**
      * Read parameters and initialize.
      *
      * \param in input parameter file
      */
      virtual void readParameters(std::istream& in);

      //Setup some of the molecule stuff
      /**
      * Destructor.
      */  
      virtual ~LammpsDumpSelectWriter();

      /**
      * Read a single frame. Frames are assumed to be read consecutively.
      *
      * \param file output file stream
      * \param iStep MD time step index
      */
      void writeFrame(std::ofstream &file, long iStep);

      void selectAtomSetup();

   private:

      struct IoAtom {
         Vector position;
         Vector velocity;
         int      typeId;
         int      id;

         IoAtom()
          : position(0.0),
            velocity(0.0),
            typeId(-1),
            id(-1)
            {}
       };
      std::vector<IoAtom> atoms_;
      /// Number of atoms in the file.
      int nAtom_;
      /// Number of atom types it is desired to record
      int nTypes_;
      int AtomSubset_;

      /// Atom types desired to save
      DArray<int> desiredTypes_;
      // Array of new atomIds
      DArray<int> atomSubsetIds_;
      long isInitialized_;
   };

}
#endif


