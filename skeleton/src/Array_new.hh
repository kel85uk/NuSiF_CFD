#ifndef ARRAY_HH
#define ARRAY_HH

#include <iostream>
#include <iomanip> //std::setw
#include <vector>
#include <algorithm>

#include "Types.hh"
#include "Debug.hh"

//*******************************************************************************************************************
/*!  Array class for 1,2 and 3 dimensions
*
*    - all elements are stored in a contiguous chunk of memory
*/
//*******************************************************************************************************************

template<typename T>
class Array
{

private:

   // For storing field data
   std::vector<T> Data;
   int Xsize,Ysize,Zsize;
   int dim;

public:
   Array():Data(1),Xsize(1),Ysize(1),Zsize(1),dim(0){}
   // Constructors for 1D,2D and 3D
   Array( int xSize );
   Array( int xSize, int ySize );
   Array( int xSize, int ySize, int zSize );

   // Destructor
   ~Array();

   // Access Operators for 1D, 2D and 3D
   inline T & operator () ( int i );
   inline T & operator () ( int i ,int j );
   inline T & operator () ( int i, int j, int k );

   // Access operators for const Arrays
   inline const T & operator () ( int i ) const;
   inline const T & operator () ( int i ,int j ) const;
   inline const T & operator () ( int i, int j, int k ) const;
   
   // Arithmetic operators
   Array<T>& operator+= (const Array &rhs);
   Array<T>& operator-= (const Array &rhs);   
   Array<T>& operator*= (const Array &rhs);   
   const Array<T> operator+ (const Array<T> &other) const;
   const Array<T> operator- (const Array<T> &other) const;
   const Array<T> operator* (const Array<T> &other) const;
   Array<T> operator- () const; 
   Array<T> operator* (const real value) const;
   
   // Minimum and maximum operators for Arrays
   inline T minimum();
   inline T maximum();
   inline T amax();
   inline T amin();

   // Minimum and maximum operators for const Arrays
   inline const T minimum() const;
   inline const T maximum() const;
   inline const T amin() const;
   inline const T amax() const;
      
   // Norm for Array
   inline real norm(std::string &s) const;
	
	// Sum for Array
   inline T sum() const;

   // initialize the whole array with a constant value
   void fill( T value );
   
   // Dot product of two Arrays (No dimension check)
   T dotNC( Array<T>& vec2);

   // return total size of the array
   int getSize() const;

   // return xSize for dimension==1, ySize for dimension==2 and zSize for dimension==3
   // other dimension values are not allowed
   int getSize(int dimension ) const;

   // Return dimension of Array
   int getDimen() const{return dim;}

   // Print the whole array ( for debugging purposes )
   const void print() const;   

};

//===================================================================================================================
//
//  Constructors
//
//===================================================================================================================

// For 1D Array
template <typename T>
Array<T>::Array( int xSize ):Data(xSize), Xsize(xSize), Ysize(1), Zsize(1), dim(0)
{}

// For 2D Array
template <typename T>
Array<T>::Array( int xSize, int ySize ):Data(xSize*ySize), Xsize(xSize), Ysize(ySize), Zsize(1), dim(1)
{}

// For 3D Array
template <typename T>
Array<T>::Array( int xSize, int ySize, int zSize ):Data(xSize*ySize*zSize), Xsize(xSize), Ysize(ySize),
                                                Zsize(zSize), dim(2)
{}

// Destructor
template <typename T>
Array<T>::~Array()
{}

//===================================================================================================================
//
//  Inline Access Operators and Sizes
//
//===================================================================================================================


// Operator() 1D
template <typename T>
inline T& Array<T>::operator ()(int i)
{
   ASSERT_MSG(dim >=0 && dim < 3, "Array is not 1D");
   ASSERT_MSG(i<Xsize*Ysize*Zsize, "Array out of bounds");
   return Data[i];
}

// Operator() 2D
template <typename T>
inline T& Array<T>::operator ()(int i,int j)
{
   ASSERT_MSG(1 == dim, "Array is not 2D");
   ASSERT_MSG(i<Xsize && j<Ysize, "Array out of bounds");
   return Data[i + j*(Xsize)];
}

// Operator() 3D
template <typename T>
inline T& Array<T>::operator ()(int i, int j, int k)
{
   ASSERT_MSG(2 == dim, "Array is not 3D");
   ASSERT_MSG(i<Xsize && j<Ysize && k<Zsize, "Array out of bounds");
   return Data[i + j*(Xsize) + k*(Xsize*Ysize)];
}

// Operator() const 1D
template <typename T>
inline const T& Array<T>::operator ()(int i) const
{
   ASSERT_MSG(0 == dim, "Array is not 1D");
   ASSERT_MSG(i<Xsize, "Array out of bounds");
   return Data[i];
}

// Operator() const 2D
template <typename T>
inline const T& Array<T>::operator ()(int i,int j) const
{
   ASSERT_MSG(1 == dim, "Array is not 2D");
   ASSERT_MSG(i<Xsize && j<Ysize, "Array out of bounds");
   return Data[i + j*(Xsize)];
}

// Operator() const 3D
template <typename T>
inline const T& Array<T>::operator ()(int i, int j, int k) const
{
   ASSERT_MSG(2 == dim, "Array is not 3D");
   ASSERT_MSG(i<Xsize && j<Ysize && k<Zsize, "Array out of bounds");
   return Data[i + j*(Xsize) + k*(Xsize*Ysize)];
}

template <typename T>
Array<T>& Array<T>::operator+= (const Array &rhs)
{
	ASSERT_MSG(Xsize*Ysize*Zsize == rhs.getSize(), "Array are of different sizes");
	for (auto it = Data.begin(); it!=Data.end();++it){
		*it += rhs.Data[it-Data.begin()];
	}
	return *this;
}

template <typename T>
Array<T>& Array<T>::operator-= (const Array &rhs)
{
	ASSERT_MSG(Xsize*Ysize*Zsize == rhs.getSize(), "Array are of different sizes");
	for (auto it = Data.begin(); it!=Data.end();++it){
		*it -= rhs.Data[it-Data.begin()];
	}
	return *this;
}

template <typename T>
Array<T>& Array<T>::operator*= (const Array &rhs)
{
	ASSERT_MSG(Xsize*Ysize*Zsize == rhs.getSize(), "Array are of different sizes");
	for (auto it = Data.begin(); it!=Data.end();++it){
		*it *= rhs.Data[it-Data.begin()];
	}
	return *this;
}

template <typename T>
const Array<T> Array<T>::operator+ (const Array<T> &other) const
{
	Array<T> result(*this);
	return result += other;
}
template <typename T>
const Array<T> Array<T>::operator- (const Array<T> &other) const
{
	Array<T> result(*this);
	return result -= other;
}

template <typename T>
const Array<T> Array<T>::operator* (const Array<T> &other) const
{
	Array<T> result(*this);
	return result *= other;
}

template <typename T>
Array<T> Array<T>::operator- () const
{
	Array<T> result((*this)*(-1));
	return result;
}
template <typename T>
Array<T> Array<T>::operator* (const real value) const
{
	Array<T> result(*this);
	for (auto i=result.Data.begin(); i!=result.Data.end();++i)
		result.Data[i-result.Data.begin()] *= value;
	return result;
}

//

/*
const Array<T> operator* (const Array<T> &other) const;
*/
// Minimum operators for Arrays
template <typename T>
inline T Array<T>::minimum()
{
   return (*min_element(Data.begin(), Data.end()));
}

// Maximum operators for Arrays
template <typename T>
inline T Array<T>::maximum()
{
   return (*max_element(Data.begin(), Data.end()));
}

// Absolute Min/Max operators for Arrays
template <typename T>
inline T Array<T>::amax()
{
   T result = 0;
	std::for_each(Data.begin(),Data.end(),[&](T val){result = (result < std::abs(val))? std::abs(val):result;});
	return result;
}

// Absolute Maximum operators for Arrays
template <typename T>
inline T Array<T>::amin()
{
   T result = 0;
	std::for_each(Data.begin(),Data.end(),[&](T val){result = (result > std::abs(val))? std::abs(val):result;});
	return result;
}

// Norm operator for Arrays
template <typename T>
inline real Array<T>::norm(std::string &s) const
{
	real result = 0.;
	if (s == "2"){
		std::for_each(Data.begin(),Data.end(),[&](T val){result += val*val;});
		result = std::sqrt(result);
	}
	else if (s == "1"){
		std::for_each(Data.begin(),Data.end(),[&](T val){result += std::abs(val);});
	}
	else if (s == "inf"){
		std::for_each(Data.begin(),Data.end(),[&](T val){result = (std::abs(val) > result)? std::abs(val):result;});
	}
   return result;
}

// Sum operator for Arrays
template <typename T>
inline T Array<T>::sum() const
{
	real result = 0.;
	std::for_each(Data.begin(),Data.end(),[&](T val){result += val;});
   return result;
}

// Minimum operators for const Arrays
template <typename T>
inline const T Array<T>::minimum() const
{
   return (*min_element(Data.begin(), Data.end()));
}

// Maximum operators for const Arrays
template <typename T>
inline const T Array<T>::maximum() const
{
   return (*max_element(Data.begin(), Data.end()));
}

template <typename T>
inline const T Array<T>::amax() const
{
   T result = 0;
	std::for_each(Data.begin(),Data.end(),[&](T val){result = (result < std::abs(val))? std::abs(val):result;});
	return result;
}

// Absolute Maximum operators for Arrays
template <typename T>
inline const T Array<T>::amin() const
{
   T result = 0;
	std::for_each(Data.begin(),Data.end(),[&](T val){result = (result > std::abs(val))? std::abs(val):result;});
	return result;
}

//===================================================================================================================
//
//  Convenience Functions
//
//===================================================================================================================


//initialize the whole array with a constant value
template <typename T>
void Array<T>::fill( T value )
{
   Data.assign(Data.size(),value);
}

// Dot product (no dimensions check)
template <typename T>
T Array<T>::dotNC( Array<T>& vec2)
{
	T result = 0;
	ASSERT_MSG(getSize() == vec2.getSize(), "Array of two different lengths!");
	for (auto it = Data.begin(); it!= Data.end(); ++it){
		result += (*it)*vec2(it - Data.begin());
	}
	return result;
}

// Print the whole array (for debugging purposes)
template <typename T>
const void Array<T>::print() const
{
   // For 2D Arrays the positive x-coordinate goes to the right
   //                   positive y-coordinate goes upwards
   //      -> the line with highest y-value is printed first
   switch(dim)
   {
      case 0: for (int i=0; i<Xsize; i++)
                 std::cout<< std::left<< std::setprecision(6)<< std::setw(12)<<Data[i];

              break;

      case 1: for (int j=Ysize-1; j>=0; j--) {
                 for (int i=0; i<Xsize; i++)
                     std::cout<< std::left<< std::setprecision(6)<< std::setw(12)<<Data[i + j*(Xsize)];
                 std::cout<<"\n";}
              break;
      case 2: for (int k=0; k<Zsize; k++){
                 std::cout<<"\nz = "<< k << "\n";
                 for (int j=Ysize-1; j>=0; j--) {
                    for (int i=0; i<Xsize; i++)
                        std::cout<< std::left<< std::setprecision(6)<<std::setw(12)<<Data[i + j*(Xsize) + k*(Xsize*Ysize)]<<" ";
                    std::cout<<"\n"; } }
              break;
      default:break;
   }
}

//returns required dimension of the array
template <typename T>
int Array<T>::getSize( int dimension ) const
{
   CHECK_MSG(dimension <= dim, "Invalid request for dimension");
   switch(dimension)
   {
      case 0: return(Xsize);
      case 1: return(Ysize);
      case 2: return(Zsize);
      default:return 0;
   }
}

//returns total size of the array
template <typename T>
int Array<T>::getSize() const
{
   switch(dim)
   {
      case 0: return(Xsize);
      case 1: return(Xsize*Ysize);
      case 2: return(Xsize*Ysize*Zsize);
      default:return 0;
   }
}

#endif //ARRAY_HH
