/**
  @file Point.h
  @brief A point class
  @author Hari Sundar
  */

/***************************************************************************
 *   Copyright (C) 2005 by Hari sundar   *
 *   hsundar@seas.upenn.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
***************************************************************************/
#ifndef __POINT_H
#define __POINT_H

#include  <cmath>
#include  <array>

/**
  @brief A point class
  @author Hari Sundar
  @author Masado Ishii
  */
template <unsigned int dim>
class Point{
  public:
    static constexpr unsigned int m_uiDim = (dim > 3 ? dim : 3);

    /** @name Constructors and Destructor */
    //@{
    Point();
    // virtual ~Point();

    Point(double scale);
    Point(const std::array<double, dim> &newCoords);
    Point(const double * newCoords);
    Point(double newx, double newy, double newz);
    Point(int newx, int newy, int newz);
    Point(unsigned int newx, unsigned int newy, unsigned int newz);
    Point(const Point &newpoint);
    //@}

    /** @name Getters */
    //@{
    const double& x() const {return _x; };
    const double& y() const {return _y; };
    const double& z() const {return _z; };
    const double& x(unsigned d) const { return _coords[d]; }

    int xint() const {return static_cast<int>(_x); };
    int yint() const {return static_cast<int>(_y); };
    int zint() const {return static_cast<int>(_z); };
    int xint(unsigned d) const { return static_cast<int>(_coords[d]); }
    //@}

    /** @name Overloaded Operators */
    //@{
    inline Point operator-() const;

    inline void operator += (const Point &other);
    inline void operator -= (const Point &other);
    inline void operator /= (const int divisor);
    inline void operator /= (const double divisor);
    inline void operator *= (const int factor);
    inline void operator *= (const double factor);

    inline Point& operator=(const Point &other);
    inline Point  operator+(const Point &other) const;
    inline Point  operator-(const Point &other) const;

    inline Point  operator/(const double divisor) const;
    inline Point  operator*(const double factor) const;
    
    double magnitude();

    inline bool operator != (const Point &other)///{ return ( (xint() != other.xint() ) || (yint() != other.yint()) || (zint() != other.zint())); };
    {
      unsigned int d = 0;
      while (d < dim && xint(d) == other.xint(d)) { d++; }
      return d != dim;
    }

    inline bool operator == (const Point &other)///{ return ( ( xint() == other.xint() ) && ( yint() == other.yint()) && ( zint() == other.zint())); };
    {
      unsigned int d = 0;
      while (d < dim && xint(d) != other.xint(d)) { d++; }
      return d == dim;
    }
    //@}

    inline double dot(Point other) { 
      double sum = 0.0;
      #pragma unroll(dim)
      for (unsigned int d = 0; d < dim; d++)
        sum += _coords[d] * other._coords[d];
      return sum;
    }

    inline double dot3(Point Other) {
      return  (_x*Other._x+_y*Other._y+_z*Other._z);
    };

    inline Point cross(Point  Other){
      return  Point(_y*Other._z-Other._y*_z, _z*Other._x-_x*Other._z, 
          _x*Other._y-_y*Other._x); 
    };

    inline double abs(){
      return sqrt(dot(*this));
    };

    void normalize();

    static Point TransMatMultiply( double* transMat, Point inPoint);
    static Point TransMatMultiply3( double* transMat, Point inPoint);
  protected:
    inline void initialize3(double newx, double newy, double newz);

    std::array<double,m_uiDim> _coords;
    double &_x = _coords[0];
    double &_y = _coords[1];
    double &_z = _coords[2];
};

template <unsigned int dim>
Point<dim>::Point()
{
  _coords.fill(0.0);
}

template <unsigned int dim>
Point<dim>::Point(double scale)
{
  _coords.fill(scale);
}

template <unsigned int dim>
Point<dim>::Point(const std::array<double, dim> &newCoords)
{
  std::copy(&newCoords[0], &newCoords[dim], &_coords[0]);
  std::fill(&_coords[dim], &_coords[m_uiDim], 0.0);   // 2d or 1d
}

template <unsigned int dim>
Point<dim>::Point(const double * newCoords)
{
  std::copy(newCoords, newCoords + dim, &_coords[0]);
  std::fill(&_coords[dim], &_coords[m_uiDim], 0.0);   // 2d or 1d
}

template <unsigned int dim>
Point<dim>::Point(double newx, double newy, double newz)
{
  initialize3(newx, newy, newz);
}

template <unsigned int dim>
Point<dim>::Point(int newx, int newy, int newz)
{ 
  initialize3(static_cast<double>(newx),
      static_cast<double>(newy),
      static_cast<double>(newz));
}

template <unsigned int dim>
Point<dim>::Point(unsigned int newx, unsigned int newy, unsigned int newz)
{ 
  initialize3(static_cast<double>(newx),
      static_cast<double>(newy),
      static_cast<double>(newz));
}

template <unsigned int dim>
Point<dim>::Point(const Point &newposition)
{
  _coords = newposition._coords;
}
/*
template <unsigned int dim>
Point<dim>::~Point()
{

}
*/

template <unsigned int dim>
inline void Point<dim>::initialize3(double newx, double newy, double newz)
{
  _x = newx;  _y = newy;  _z = newz;
}

template <unsigned int dim>
Point<dim> Point<dim>::operator - () const {
  Point ret(*this);
  #pragma unroll(dim)
  for (unsigned int d = 0; d < dim; d++)
    ret._coords[d] = -ret._coords[d];
  return ret;
}

template <unsigned int dim>
void Point<dim>::operator *= (const int factor){
  #pragma unroll(dim)
  for (unsigned int d = 0; d < dim; d++)
    _coords[d] *= factor;
}

template <unsigned int dim>
void Point<dim>::operator *= (const double factor){
  #pragma unroll(dim)
  for (unsigned int d = 0; d < dim; d++)
    _coords[d] *= factor;
}

template <unsigned int dim>
void Point<dim>::operator /= (const int divisor){
  if (divisor == 0) return;
  #pragma unroll(dim)
  for (unsigned int d = 0; d < dim; d++)
    _coords[d] /= static_cast<double>(divisor);
}

template <unsigned int dim>
void Point<dim>::operator /= (const double divisor){
  if (divisor == 0) return;
  #pragma unroll(dim)
  for (unsigned int d = 0; d < dim; d++)
    _coords[d] /= divisor;
}

template <unsigned int dim>
void Point<dim>::operator += (const Point& other){
  #pragma unroll(dim)
  for (unsigned int d = 0; d < dim; d++)
    _coords[d] += other._coords[d];
}

template <unsigned int dim>
void Point<dim>::operator -= (const Point& other){
  #pragma unroll(dim)
  for (unsigned int d = 0; d < dim; d++)
    _coords[d] -= other._coords[d];
}

template <unsigned int dim>
Point<dim> Point<dim>::operator - (const Point &other) const{
  Point ret(*this);
  ret -= other;
  return ret;
}

template <unsigned int dim>
Point<dim> Point<dim>::operator + (const Point &other) const{
  Point ret(*this);
  ret += other;
  return ret;
}

template <unsigned int dim>
Point<dim>& Point<dim>::operator=(const Point &other){
  #pragma unroll(dim)
  for (unsigned int d = 0; d < dim; d++)
    _coords[d] = other._coords[d];
  return *this;
}

template <unsigned int dim>
Point<dim> Point<dim>::operator /(const double divisor) const
{
  Point ret(*this);
  ret /= divisor;
  return ret;
}

template <unsigned int dim>
Point<dim> Point<dim>::operator *(const double factor) const
{
  Point ret(*this);
  ret *= factor;
  return ret;
}

template <unsigned int dim>
Point<dim> Point<dim>::TransMatMultiply3(double *transMat, Point inPoint)
{
  Point outPoint;

  outPoint._x = transMat[ 0]*inPoint._x +transMat[ 4]*inPoint._y +transMat[8]
    *inPoint._z +transMat[12];
  outPoint._y = transMat[ 1]*inPoint._x +transMat[ 5]*inPoint._y +transMat[9]
    *inPoint._z +transMat[13];
  outPoint._z = transMat[ 2]*inPoint._x +transMat[ 6]*inPoint._y
    +transMat[10]*inPoint._z +transMat[14];

  return outPoint;
}

template <unsigned int dim>
Point<dim> Point<dim>::TransMatMultiply(double *transMat, Point inPoint)
{
  if (dim == 3)
    return TransMatMultiply3(transMat, inPoint);

  Point outPoint;

  for (unsigned int i = 0; i < dim; i++)
  {
    outPoint._coords[i] = transMat[dim*(dim+1) + i];
    for (unsigned int j = 0; j < dim; j++)
      outPoint._coords[i] += transMat[j*(dim+1) + i] * inPoint._coords[j];
  }

  return outPoint;
}



template <unsigned int dim>
void Point<dim>::normalize() {
  operator/=(abs());
}

template <unsigned int dim>
double Point<dim>::magnitude()
{
  return abs();
}

// Template instantiations.
template class Point<2u>;
template class Point<3u>;
template class Point<4u>;

#endif // POINT_H
