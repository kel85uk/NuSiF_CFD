/*    Copyright (C) 2013  kklloh

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/


#include "Debug.hh"

#include <iostream>
#include <cstdlib> // std::abort

//http://patrick.wagstrom.net/old/misc/dotfiles/dotbash_profile.html
const char * BLACK        = "\033[0;30m";
const char * RED          = "\033[0;31m";
const char * GREEN        = "\033[0;32m";
const char * BROWN        = "\033[0;33m";
const char * BLUE         = "\033[0;34m";
const char * MAGENTA      = "\033[0;35m";
const char * CYAN         = "\033[0;36m";
const char * WHITE        = "\033[0;37m";
const char * LIGHTBLACK   = "\033[1;30m";
const char * LIGHTRED     = "\033[1;31m";
const char * LIGHTGREEN   = "\033[1;32m";
const char * YELLOW       = "\033[1;33m";
const char * LIGHTBLUE    = "\033[1;34m";
const char * LIGHTMAGENTA = "\033[1;35m";
const char * LIGHTCYAN    = "\033[1;36m";
const char * LIGHTWHITE   = "\033[1;37m";


namespace internal
{

   void checkFct( bool b, const char * const expression, const std::string & message,
                  const char * const filename, int line )
   {
      if( !b )
      {
         std::cerr << "Check failed!\n\tFile:       " << filename << ":" << line << "\n"
                   << "\tExpression: " << expression << std::endl;
         if ( message.size() > 0 )
            std::cerr << "\tMessage: \n\t" << message << std::endl;
      }
      std::abort();
   }

   void assertFct( bool b, const char * const expression, const std::string & message,
                   const char * const filename, int line )
   {
      if( !b )
      {
         std::cerr << "Assertion failed!\n\tFile:       " << filename << ":" << line << "\n"
                   << "\tExpression: " << expression << std::endl;

         if ( message.size() > 0 )
            std::cerr << "\tMessage: \n\t" << message << std::endl;
      }
      std::abort();
   }

   void warnFct( const std::string & message, const char * const filename, int line )
   {
      std::cerr << "Warning!\n\tFile:       " << filename << ":" << line << "\n";
      std::cerr <<  "\t" << message << std::endl;
   }
   
   void progFct( const std::string & message, const char * const filename, int line ) /* Added progress routines */
   {
      std::cerr << "Currently in \n\tFile:       " << filename << ":" << line << "\n";
      std::cerr <<  "\t" << message << std::endl;
   }   


}



