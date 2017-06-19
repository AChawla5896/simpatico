/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "Command.h"
#include <util/archives/Serializable_includes.h>

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   Command::Command() 
   {}

   /*
   * Default destructor.
   */
   Command::~Command()
   {}

   /*
   * Output at end - empty default implementation.
   */
   void Command::output()
   {}

}
