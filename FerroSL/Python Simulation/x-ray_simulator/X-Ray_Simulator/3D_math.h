#ifndef MATH_3D_H
#define MATH_3D_H

#include <math.h>
#include <iostream>

struct Vec3D;

/////////////////
// Vec3D struc //
/////////////////

struct Vec3D {
private:
	double coords[3];
public:
	//Constructor
	Vec3D();
	Vec3D(double x, double y, double z);
	//Return coordinates
	double x(void);
	double y(void);
	double z(void);
	//Set coordinates
	void x(double x);
	void y(double y);
	void z(double z);
	//Vector operations
	double module();
	Vec3D cross3(Vec3D); //the cross product
	double dot3(Vec3D); //the dot product
	//Vector operators
	Vec3D&  operator =(Vec3D&);
	Vec3D operator +(Vec3D const &);
	Vec3D operator -(Vec3D const &);
	double operator *(Vec3D const &); //Dot product
	Vec3D operator %(Vec3D const &); //cross product
	Vec3D& operator +=(Vec3D const &);
	Vec3D& operator -=(Vec3D const &);
	//Vec3D& operator %=(Vec3D const &); //cross product
	bool operator ==(Vec3D const &);
	bool operator !=(Vec3D const &);
	double &  operator[](unsigned int i);
	double const &  operator[](unsigned int i) const;
	//double operators
	Vec3D operator *(double k); //multiplication of vector by a scalar k
	Vec3D& operator *=(double k); //multiplication of vector by a scalar k
};

std::ostream& operator<<( std::ostream& os, const Vec3D& V );

#endif
