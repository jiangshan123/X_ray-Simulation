#include "3D_math.h"

//////////////////
// Vec3D struct //
//////////////////

//Initialize Vec3D with 0
Vec3D::Vec3D() {
	coords[0] = 0;
	coords[1] = 0;
	coords[2] = 0;
}

//Initialize Vec3D with values
Vec3D::Vec3D(double x, double y, double z) {
	coords[0] = x;
	coords[1] = y;
	coords[2] = z;
}

//return coordinate 0 (x)
double Vec3D::x() {
	return coords[0];
}

//return coordinate 1 (y)
double Vec3D::y() {
	return coords[1];
}

//return coordinate 2 (z)
double Vec3D::z() {
	return coords[2];
}

//set coordinate 0 (x)
void Vec3D::x(double x) {
	coords[0] = x;
}

//set coordinate 1 (y)
void Vec3D::y(double y) {
	coords[1] = y;
}

//set coordinate 2 (z)
void Vec3D::z(double z) {
	coords[2] = z;
}

//returns the module of the vector
double Vec3D::module() {
	double sum = 0;
	for (int i = 0; i<3; i++) {
		sum += coords[i]*coords[i];
	}
	return sqrt(sum);
}

//returns the dot product of the vector with V2
double Vec3D::dot3(Vec3D V2) {
	double sum = 0;
	for (int i = 0; i<3; i++) {
		sum += coords[i]*V2[i];
	}
	return sum;
}

//returns the cross product of the vector with V2
Vec3D Vec3D::cross3(Vec3D V2) {
	Vec3D V3;
	V3[0] = coords[1]*V2[2]-coords[2]*V2[1];
	V3[1] = coords[2]*V2[0]-coords[0]*V2[2];
	V3[2] = coords[0]*V2[1]-coords[1]*V2[0];
	return V3;
}

//Sets the coordinates of the vector to the given vectore
Vec3D& Vec3D::operator =(Vec3D & V2) {
	for(int i= 0; i<3; i++){
        coords[i] = V2[i];
    }
	return *this;
}

//Adds two vectors
Vec3D Vec3D::operator +(Vec3D const & V2) {
	return Vec3D(coords[0]+V2[0],coords[1]+V2[1],coords[2]+V2[2]);
}

//Subtracts two vectors
Vec3D Vec3D::operator -(Vec3D const & V2) {
	return Vec3D(coords[0]-V2[0],coords[1]-V2[1],coords[2]-V2[2]);
}

//Dot product of two vectors
double Vec3D::operator *(Vec3D const & V2) {
	return this->dot3(V2);
}

//Cross product of two vectors
Vec3D Vec3D::operator %(Vec3D const & V2) {
	return this->cross3(V2);
}

//adds the coordinates of the vector to the given vector
Vec3D& Vec3D::operator +=(Vec3D const & V2) {
	for(int i= 0; i<3; i++){
        coords[i] += V2[i];
    }
	return *this;
}

//Subtracts the coordinates of the vector to the given vector
Vec3D& Vec3D::operator -=(Vec3D const & V2) {
	for(int i= 0; i<3; i++){
        coords[i] -= V2[i];
    }
	return *this;
}

//calculates the cross product with the given vector
/*Vec3D& Vec3D::operator %=(Vec3D const & V2) {
	Vec3D temp;
	temp = this->cross3(V2);
	coords[0] = temp[0];
	coords[1] = temp[1];
	coords[2] = temp[2];
	return *this;
}*/

//[] bracket operator to access coords
double & Vec3D::operator [](unsigned int i) {
	return coords[i];
}

//[] bracket operator to access coords and be able to set to constants
double const & Vec3D::operator [](unsigned int i) const{
	return coords[i];
}

//multiplies with a double
Vec3D Vec3D::operator *(double k) {
	return Vec3D(coords[0]*k,coords[1]*k,coords[2]*k);
}

//multiplies with a double
Vec3D& Vec3D::operator *=(double k) {
	for(int i= 0; i<3; i++){
		coords[i] *= k;
	}
	return *this;
}

// returns true if two vectors are the same
bool Vec3D::operator ==(Vec3D const & V2) {
	if (V2[0] != coords[0]) return false;
	if (V2[1] != coords[1]) return false;
	if (V2[2] != coords[2]) return false;
	return true;
}

//Returns True if two vectors aren't the same
bool Vec3D::operator !=(Vec3D const & V2) {
	if ( *this == V2) return false;
	return true;
}

//simple output for a vector
std::ostream& operator<<( std::ostream& os, const Vec3D& V ) {
	os << "(" << V[0] << "," << V[1] << "," << V[2] << ")";
	return os;
}
