
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
