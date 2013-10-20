/*  Copyright 2012 Daniel Wilson.
 *
 *  mydouble.h
 *  Part of the myutils library.
 *
 *  The myutils library is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *  
 *  The myutils library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *  
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with the myutils library. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef _MY_DOUBLE_H_
#define _MY_DOUBLE_H_

#include <limits>
#include <math.h>
#include <myerror.h>

using myutils::error;

/*	This class behaves to the user like a non-negative double, but
	is stored internally as the natural logarithm. Standard mathematical
	operations are performed on the logarithm of the number so that it
	should not underflow or overflow like a double. */
class mydouble {
protected:
	double _log;
	bool _zero;
public:
	/*Default constructor*/
	mydouble() {
		_zero = false;
	};
	/*Copy constructor*/
	mydouble(const double &_doub) {
		_zero = false;
		if(_doub<0.0) myutils::error("mydouble::mydouble(const double&): cannot initialize with negative number");
		if(_doub==0.0) setzero();
		else _log = log(_doub);
	};
	/*Copy constructor*/
	mydouble(const mydouble &_mydoub) {
		_zero = _mydoub._zero;
		_log = _mydoub._log;
	}
	// Construct a zero
	static mydouble zero() {
		mydouble z(0);
		return z;
	}
	/*Conversion operator
		THIS CONVERSION OPERATOR HAS BEEN DISABLED BECAUSE IT ALLOWED THE COMPILER TO
		IMPLICITLY MAKE MYDOUBLE->DOUBLE CONVERSIONS WHICH RESULTED IN LOSS OF PRECISION
		WHEN DOUBLE->MYDOUBLE CONVERSIONS WERE REQUIRED TO MAINTAIN PRECISION. IT HAS
		BEEN REPLACED BY THE SUBSEQUENT FUNCTION WHICH IS AN EXPLICIT CONVERSION TO TYPE
		DOUBLE WHICH THE COMPILER CANNOT CALL IMPLICITLY.
	operator double const() {
		return (_zero) ? 0.0 : exp(_log);
	};*/
	double todouble() const {
		return (_zero) ? 0.0 : exp(_log);
	}
	/*Assignment operator*/
	mydouble& operator=(const double &_doub) {
		_zero = false;
		if(_doub<0.0) myutils::error("mydouble::operator=(const double&): cannot assign a negative number");
		if(_doub==0.0) setzero();
		else _log = log(_doub);
		return *this;
	}
	/*Assignment operator*/
	mydouble& operator=(const mydouble &_mydoub) {
		_zero = _mydoub._zero;
		_log = _mydoub._log;
		return *this;
	}

	mydouble& setlog(const double &log) {
		_zero = false;
		_log = log;
		return *this;
	}
	mydouble& setzero() {
		_zero = true;
		_log = -std::numeric_limits<double>::max();
		return *this;
	}
	bool iszero() const {
		return _zero;
	}
	bool isinfinity() const {
		return !_zero && _log==std::numeric_limits<double>::infinity();
	}
	bool isbad() const {
		return !_zero && _log!=_log;
	}
	
	/*** MULTIPLICATION ***/
	mydouble operator*(const double &dbl) const {
		return operator*(mydouble(dbl));
	}
	mydouble operator*(const mydouble &mydbl) const {
		mydouble a;
		if(_zero || mydbl._zero) a.setzero();
		else a.setlog(_log + mydbl._log);
		return a;
	}
	mydouble& operator*=(const double &dbl) {
		if(_zero || dbl==0.0) setzero();
		else _log += mydouble(dbl)._log;
		return *this;
	}
	mydouble& operator*=(const mydouble &mydbl) {
		if(_zero || mydbl._zero) setzero();
		else _log += mydbl._log;
		return *this;
	}

	/*** DIVISION ***/
	mydouble operator/(const double &dbl) const {
		return operator/(mydouble(dbl));
	}
	mydouble operator/(const mydouble &mydbl) const {
		mydouble a;
		if(mydbl._zero) error("mydouble::operator/(const mydouble&): division by zero");
		else if(_zero) a.setzero();
		else a.setlog(_log - mydbl._log);
		return a;
	}
	mydouble& operator/=(const double &dbl) {
		if(dbl==0.0) error("mydouble::operator/=(const double&): division by zero");
		else if(!_zero) _log -= mydouble(dbl)._log;
		return *this;
	}
	mydouble& operator/=(const mydouble &mydbl) {
		if(mydbl._zero) error("mydouble::operator/=(const mydouble&): division by zero");
		else if(!_zero) _log -= mydbl._log;
		return *this;
	}

	/*** ADDITION ***/
	mydouble operator+(const double &dbl) const {
		if(dbl==0.0) return mydouble(*this);
		if(dbl<0.0) return operator-(mydouble(-dbl));
		return operator+(mydouble(dbl));
	}
	mydouble operator+(const mydouble &mydbl) const {
		mydouble a;
		if(_zero) a = mydouble(mydbl);
		else if(mydbl._zero) a = mydouble(*this);
		else {
			double diff = _log - mydbl._log;
			if(diff==0.0) a.setlog(log(2.0) + _log);
			else if(diff<0.0) a.setlog(mydbl._log + log(1.0 + exp(diff)));
			else a.setlog(_log + log(1.0 + exp(-diff)));
		}
		return a;
	}
	mydouble& operator+=(const double &dbl) {
		if(dbl==0.0) return *this;
		return operator+=(mydouble(dbl));
	}
	mydouble& operator+=(const mydouble &mydbl) {
		if(_zero) *this = mydbl;
		else if(!mydbl._zero) {
			double diff = _log - mydbl._log;
			if(diff==0.0) _log += log(2.0);
			else if(diff<0.0) _log = mydbl._log + log(1.0 + exp(diff));
			else _log += log(1.0 + exp(-diff));
		}
		return *this;
	}

	/*** SUBTRACTION - warning cannot have negative numbers ***/
	mydouble operator-(const double &dbl) const {
		if(dbl==0.0) return mydouble(*this);
		return operator-(mydouble(dbl));
	}
	mydouble operator-(const mydouble &mydbl) const {
		mydouble a;
		if(mydbl._zero) a = mydouble(*this);
		else if(_zero) error("mydouble::operator-(const mydouble&): subtracting a positive number from zero");
		else {
			/* diff must always be positive */
			double diff = _log - mydbl._log;
			if(diff==0.0) a.setzero();
			else if(diff<0.0) myutils::error("mydouble::operator-(const mydouble&) cannot handle negative numbers");
			else a.setlog(_log + log(1.0 - exp(-diff)));
		}
		return a;
	}
	mydouble& operator-=(const double &dbl) {
		if(dbl==0.0) return *this;
		return operator-=(mydouble(dbl));
	}
	mydouble& operator-=(const mydouble &mydbl) {
		if(!mydbl._zero) {
			if(_zero) error("mydouble::operator-=(const mydouble&): subtracting a positive number from zero");
			/* diff must always be positive */
			double diff = _log - mydbl._log;
			if(diff==0.0) setzero();
			else if(diff<0.0) myutils::error("mydouble::operator-=(const mydouble&) cannot handle negative numbers");
			else _log += log(1.0 - exp(-diff));
		}
		return *this;
	}

	/*** SPECIAL OPERATIONS ***/
	double LOG() const {
		return _log;
	}
	/* Caution: ^ has lower precedence than /*+- */
	mydouble operator^(const double &dbl) const {
		mydouble a;
		if(_zero) a.setzero();
		else a.setlog(_log * dbl);
		return a;
	}
	/* Caution: ^ has lower precedence than /*+- */
	mydouble operator^(const mydouble &mydbl) const {
		mydouble a;
		if(_zero) a.setzero();
		else a.setlog(_log * exp(mydbl._log));
		return a;
	}
	mydouble& operator^=(const double &dbl) {
		if(!_zero) _log *= dbl;
		return *this;
	}
	mydouble& operator^=(const mydouble &mydbl) {
		if(!_zero) _log *= exp(mydbl._log);
		return *this;
	}

	/*** COMPARISON OPERATORS ***/
	bool operator<(const double &dbl) const {
		return operator<(mydouble(dbl));
	}
	bool operator<(const mydouble &mydbl) const {
		return (_log < mydbl._log);
	}
	bool operator<=(const double &dbl) const {
		return operator<=(mydouble(dbl));
	}
	bool operator<=(const mydouble &mydbl) const {
		return (_log <= mydbl._log);
	}
	bool operator>(const double &dbl) const {
		return operator>(mydouble(dbl));
	}
	bool operator>(const mydouble &mydbl) const {
		return (_log > mydbl._log);
	}
	bool operator>=(const double &dbl) const {
		return operator>=(mydouble(dbl));
	}
	bool operator>=(const mydouble &mydbl) const {
		return (_log >= mydbl._log);
	}
	bool operator==(const double &dbl) const {
		return operator==(mydouble(dbl));
	}
	bool operator==(const mydouble &mydbl) const {
		return (_log == mydbl._log);
	}
	bool operator!=(const double &dbl) const {
		return operator!=(mydouble(dbl));
	}
	bool operator!=(const mydouble &mydbl) const {
		return (_log != mydbl._log);
	}
};

/*** MULTIPLICATION ***/
inline mydouble operator*(const double &dbl, const mydouble &mydbl) {
	mydouble a(dbl);
	return a *= mydbl;
}
/*** DIVISION ***/
inline mydouble operator/(const double &dbl, const mydouble &mydbl) {
	mydouble a(dbl);
	return a /= mydbl;
}
/*** ADDITION ***/
inline mydouble operator+(const double &dbl, const mydouble &mydbl) {
	mydouble a(dbl);
	return a += mydbl;
}
/*** SUBTRACTION - warning cannot have negative numbers ***/
inline mydouble operator-(const double &dbl, const mydouble &mydbl) {
	mydouble a(dbl);
	return a -= mydbl;
}
/*** SPECIAL OPERATIONS ***/
inline double log(const mydouble &mydbl) {
	return mydbl.LOG();
}
inline mydouble pow(const mydouble &_X, const mydouble &_Y) {
	return _X^_Y;
}
inline mydouble pow(const mydouble &_X, const double &_Y) {
	return _X^_Y;
}
/* Caution: ^ has lower precedence than /*+- */
inline mydouble operator^(const double dbl, const mydouble &mydbl) {
	mydouble a(dbl);
	return a ^= mydbl;
}
/*** COMPARISON OPERATORS ***/
inline bool operator<(const double &dbl, const mydouble &mydbl) {
	return (mydouble(dbl)<mydbl);
}
inline bool operator<=(const double &dbl, const mydouble &mydbl) {
	return (mydouble(dbl)<=mydbl);
}
inline bool operator>(const double &dbl, const mydouble &mydbl) {
	return (mydouble(dbl)>mydbl);
}
inline bool operator>=(const double &dbl, const mydouble &mydbl) {
	return (mydouble(dbl)>=mydbl);
}
inline bool operator==(const double &dbl, const mydouble &mydbl) {
	return (mydouble(dbl)==mydbl);
}
inline bool operator!=(const double &dbl, const mydouble &mydbl) {
	return (mydouble(dbl)!=mydbl);
}

#endif//_MY_DOUBLE_H_
