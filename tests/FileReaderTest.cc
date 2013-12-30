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
#include "FileReader.hh"
#include "Debug.hh"

#include <cmath>
#include <iostream>


int main( )
{
   FileReader myReader;

   myReader.registerIntParameter   ("exampleInt" );
   myReader.registerRealParameter  ("exampleReal");
   myReader.registerStringParameter("exampleString" );
   PROGRESS("Doing stuff");
   bool res = myReader.readFile ("FileReaderTestInput.txt");
   CHECK_MSG(res, "Could not open file 'FileReaderTestInput.txt' which has to be in the current directory.");



   CHECK( myReader.getIntParameter   ("exampleInt")    == 42 );
   CHECK( myReader.getStringParameter("exampleString") == "someStringValue" );
   CHECK( std::abs( myReader.getRealParameter("exampleReal") - 42.4242 ) < 1e-5 );

   myReader.registerIntParameter("aNewInt");
   myReader.setParameter( "aNewInt",    43 ); // add new value ( no registration required )
   myReader.setParameter( "exampleInt", 44 ); // overwrite existing value

   CHECK( myReader.getIntParameter("aNewInt")     == 43 );
   CHECK( myReader.getIntParameter("exampleInt")  == 44 );
   
//   myReader.setParameter("TUD",10); //Deliberately cause an error, comment to pass
   
   myReader.printParameters();


   std::cout << "File Reader Test passed successfully" << std::endl;

   return 0;
}
