#ifdef  MCMD_PERTURB
#ifndef MCMD_MC_PERTURBATION_FACTORY_CPP
#define MCMD_MC_PERTURBATION_FACTORY_CPP

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include "McPerturbationFactory.h"  

// Subclasses of Perturbation
#include "McEnergyPerturbation.h"
#ifndef INTER_NOPAIR
#include "McPairPerturbation.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/pair/McPairPotential.h>
#include <inter/pair/LJPair.h>
#include <inter/pair/DpdPair.h>
#endif

#ifdef INTER_EXTERNAL
#include "McExternalPerturbation.h"
#include <mcMd/mcSimulation/McSystem.h>
#include <mcMd/potentials/external/ExternalPotential.h>
#include <inter/external/TanhCosineExternal.h>
#ifndef INTER_NOPAIR
#include "McPairExternalPerturbation.h"
#include <mcMd/potentials/pair/McPairPotential.h>
#include <inter/pair/LJPair.h>
#include <inter/pair/DpdPair.h>
#endif
#endif

namespace McMd
{

   using namespace Util;
   using namespace Inter;

   McPerturbationFactory::McPerturbationFactory(McSystem& system)
    : systemPtr_(&system)
   {}

   /* 
   * Return a pointer to a instance of Peturbation subclass className.
   */
   Perturbation* McPerturbationFactory::factory(const std::string &className) const
   {
      Perturbation *ptr = 0;

      // Try subfactories first
      ptr = trySubfactories(className);
      if (ptr) return ptr;

      if (className == "McEnergyPerturbation") {
         ptr = new McEnergyPerturbation(*systemPtr_);
      }
      #ifndef INTER_NOPAIR
      else if (className == "McPairPerturbation") {
         const std::string& interactionClassName = systemPtr_->pairPotential().
            interactionClassName();
         if (interactionClassName == "LJPair") {
            ptr = new McPairPerturbation<LJPair> (*systemPtr_);
         } else if (interactionClassName == "DpdPair") {
            ptr = new McPairPerturbation<DpdPair> (*systemPtr_);
         } else {
            UTIL_THROW("Unsupported pair potential.");
         }
      } 
      #endif
      #ifdef INTER_EXTERNAL
      #ifndef INTER_NOPAIR
      else if (className == "McPairExternalPerturbation") {
         const std::string& pairInteractionClassName = systemPtr_->pairPotential().interactionClassName();
         const std::string& externalInteractionClassName = systemPtr_->externalPotential().interactionClassName();
         if (pairInteractionClassName == "LJPair") {
            if (externalInteractionClassName == "TanhCosineExternal") {
               ptr = new McPairExternalPerturbation<LJPair,TanhCosineExternal> (*systemPtr_);
            } else {
               UTIL_THROW("Unsupported external potential.");
            }
         } else if (pairInteractionClassName == "DpdPair") {
            if (externalInteractionClassName == "TanhCosineExternal") {
               ptr = new McPairExternalPerturbation<DpdPair,TanhCosineExternal> (*systemPtr_);
            } else {
               UTIL_THROW("Unsupported external potential.");
            }
         } else {
            UTIL_THROW("Unsupported pair potential.");
         }
      }
      #endif 
      #endif
      #ifdef INTER_EXTERNAL
      else if (className == "McExternalPerturbation") {
         const std::string& interactionClassName = systemPtr_->externalPotential().
            interactionClassName();
         if (interactionClassName == "TanhCosineExternal") {
            ptr = new McExternalPerturbation<TanhCosineExternal> (*systemPtr_);
         } else {
            UTIL_THROW("Unsupported external potential.");
         }
      }
      #endif
      return ptr;
   }

}

#endif
#endif  // #ifdef  MCMD_PERTURB
