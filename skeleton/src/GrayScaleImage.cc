//======================================================================================================================
/*!
 *  \file   GrayScaleImage.cpp
 *  \author Martin Bauer <martin.bauer@fau.de>
 */
//======================================================================================================================


// Local includes
#include "GrayScaleImage.hh"


// Extern
#include "lodepng.h"

#include <limits>


GrayScaleImage::GrayScaleImage( const std::string & pngFilename )
{
   unsigned int tmpWidth, tmpHeight;
   unsigned int error = lodepng::decode( image_, tmpWidth, tmpHeight, pngFilename, LCT_GREY, 8 );
   size_[0] = tmpWidth;
   size_[1] = tmpHeight;

   CHECK_MSG( error == 0 , "Error while loading PNG file: " << lodepng_error_text(error) );
}

GrayScaleImage::GrayScaleImage( int i, int j )
{
   size_[0] = i;
   size_[1] = j;
	 image_.resize(i*j);
	 static const real maxVal = real( std::numeric_limits<unsigned char>::max() );
	 std::fill (image_.begin(),image_.end(),maxVal);
}

void GrayScaleImage::save( const std::string & pngFilename )
{
   unsigned int error = lodepng::encode( pngFilename, image_,
                                         static_cast<unsigned int>( size_[0] ),
                                         static_cast<unsigned int>( size_[1] ),
                                         LCT_GREY, 8 );

   CHECK_MSG( error ==0 , "Error while loading PNG file: " << lodepng_error_text(error) );
}


real GrayScaleImage::operator() ( int x, int y ) const
{
   static const real maxVal = real( std::numeric_limits<unsigned char>::max() );

   const int yFlip = size_[1] - y - 1;

   return real ( getElement(x, yFlip) ) / maxVal;
}


GrayScaleImage GrayScaleImage::getResizedImage( int newWidth, int newHeight ) const
{
   if ( newWidth == size_[0]  && newHeight == size_[1] )
      return *this;


   GrayScaleImage resizedImage;

   resizedImage.size_[0] = newWidth;
   resizedImage.size_[1] = newHeight;

   resizedImage.image_.resize( newWidth * newHeight );

   real scaleX = real( size_[0] ) / real( newWidth );
   real scaleY = real( size_[1] ) / real( newHeight);

   for( int y = 0; y <  newHeight; ++y )
      for( int x = 0; x <  newWidth; ++x )
      {
         real oldX = x * scaleX;
         real oldY = y * scaleY;
         int oldXi = static_cast<int>( oldX );
         int oldYi = static_cast<int>( oldY );
         real xDiff = oldX - oldXi;
         real yDiff = oldY - oldYi;

         // bilinear interpolation

         resizedImage.getElement( x, y ) =
                  static_cast<unsigned char> (
                  (1 - xDiff) * (1 - yDiff ) * getElement( oldXi    , oldYi    ) +
                       xDiff  * (1 - yDiff ) * getElement( oldXi + 1, oldYi    ) +
                  (1 - xDiff) *      yDiff   * getElement( oldXi    , oldYi + 1) +
                       xDiff  *      yDiff   * getElement( oldXi + 1, oldYi + 1) );
      }

   return resizedImage;
}










