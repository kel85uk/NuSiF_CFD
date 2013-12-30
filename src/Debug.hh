#ifndef DEBUG_HH
#define DEBUG_HH

#include <sstream>

//===================================================================================================================
//
//  CHECK Macro used for tests, is activated in Debug and Release Mode
//
//===================================================================================================================

#define CHECK_MSG(X, MSG) \
   if( !(X) )  \
   { std::stringstream ss; \
     ss << MSG; \
     internal::checkFct ( (X), #X, ss.str(), __FILE__, __LINE__ );\
   }

#define CHECK(X) \
   if( !(X) ) {  internal::checkFct ( (X), #X, "", __FILE__, __LINE__ ); }

#define WARN(MSG) \
   { std::stringstream ss; \
     ss << MSG; \
     internal::warnFct ( ss.str(), __FILE__, __LINE__ );\
   }   




//===================================================================================================================
//
//  ASSERT Macro checks the given expression in Debug mode, disabled in Release mode
//
//===================================================================================================================


#ifndef NDEBUG

#define ASSERT_MSG(X, MSG) \
   if( !(X) )  \
   { std::stringstream ss; \
     ss << MSG; \
     internal::assertFct ( (X), #X, ss.str(), __FILE__, __LINE__ );\
   }

#define ASSERT(X) \
   if( !(X) ) {  internal::assertFct ( (X), #X, "", __FILE__, __LINE__ ); }
   
/* Added progress routines */   
#define PROGRESS(MSG) \
   { std::stringstream ss; \
     ss << MSG; \
     internal::progFct ( ss.str(), __FILE__, __LINE__ );\
   }   


#else

#define ASSERT_MSG(X, MSG)
#define ASSERT(X)
#define PROGRESS(MSG)

#endif //NDEBUG






namespace internal
{

   void checkFct ( bool b, const char * const expression, const std::string & message,
                   const char * const filename, int line );
   void assertFct( bool b, const char * const expression, const std::string & message,
                   const char * const filename, int line );
   void warnFct( const std::string & message,
                 const char * const filename, int line );
   void progFct( const std::string & message,
                 const char * const filename, int line );    /* Added progress routines */              

}



#endif // DEBUG_HH
