/*
 * Mat.hpp, part of LatAnalyze 3
 *
 * Copyright (C) 2013 - 2014 Antonin Portelli
 *
 * LatAnalyze 3 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LatAnalyze 3 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LatAnalyze 3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef Latan_Mat_hpp_
#define Latan_Mat_hpp_

#include <latan/Eigen/Dense>
#include <latan/Global.hpp>
#include <latan/IOObject.hpp>

BEGIN_NAMESPACE

// Eigen aliases
typedef Eigen::MatrixXd  DMatBase;
typedef Eigen::MatrixXcd CMatBase;
typedef Eigen::VectorXd  DVecBase;
typedef Eigen::VectorXcd CVecBase;

class DMat: public DMatBase, public IoObject
{
public:
    // constructors
    DMat(void);
    DMat(const DMatBase &m);
    DMat(const DMat &m);
    DMat(const unsigned int nRow, const unsigned int nCol);
    // IO
    virtual unsigned int getType(void) const;
};

END_NAMESPACE

#endif // Latan_Mat_hpp_
