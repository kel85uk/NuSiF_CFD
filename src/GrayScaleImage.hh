#pragma once


#include "Types.hh"
#include "Debug.hh"

#include <string>
#include <vector>
#include <algorithm>


//*******************************************************************************************************************
/*! Class for loading and storing png images
*
* Very simple wrapper around lodepng library.
* You also have to compile and link lodepng.cc
*/
//*******************************************************************************************************************
class GrayScaleImage
{
public:
   GrayScaleImage() {} // required in getResizedImage
   
   GrayScaleImage(int i, int j);
   
   /// Loads a grayscale png image from the specified file
   GrayScaleImage( const std::string & pngFilename );

   void save( const std::string & pngFilename );

   GrayScaleImage getResizedImage( int newWidth, int newHeight ) const;

   int width()  const { return size_[0]; }
   int height() const { return size_[1]; }

   int size( int coord ) const;

   /// Returns a value between 0 and 1
   /// 0 means black - 1 means white
   real operator() ( int x, int y ) const;

   /// Returns the gray value of the specified pixel (between 0 and 255)
   unsigned char & getElement ( int x, int y );
   unsigned char   getElement ( int x, int y ) const;

   
protected:
   std::vector<unsigned char> image_; //< raw pixels   
   int size_[2];                      //< 0=width,  1=height
};







//===================================================================================================================
//
//  Implementation of inline functions
//
//===================================================================================================================


inline unsigned char &  GrayScaleImage::getElement ( int x, int y ) {
   ASSERT( x >= 0  && y >= 0 );
   ASSERT( x < size_[0] );
   ASSERT( y < size_[1] );
   return image_[ y * size_[0] + x ];
}

inline unsigned char GrayScaleImage::getElement ( int x, int y ) const {
   return const_cast<GrayScaleImage*> ( this )->getElement(x,y);
}


inline int GrayScaleImage::size( int coord ) const
{
   ASSERT(coord < 2);
   return size_[coord];
}




