/*
    Copyright 2011 Thibaut Paumard

    This file is part of Gyoto.

    Gyoto is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gyoto is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __YGYOTO_IDX_H
#define __YGYOTO_IDX_H

#include "yapi.h"

namespace YGyoto { class Idx; }

class YGyoto::Idx {
private:
  int _is_nuller;
  int _is_range;
  int _is_list;
  int _is_scalar;
  int _is_double;
  int _is_dlist;
  int _is_first;

  long _range[3];
  long _dims[Y_DIMSIZE];
  double _dval;
  double *_buf;
  long *_idx;
  long _nel;
  long _cur;
  int _valid;
  
public:
  Idx(int iarg, int res);
  int isNuller() const;
  int getNDims() const;
  long getNElements() const;
  long first();
  int  valid() const;
  long next();
  long current() const;
  double getDVal() const;
  int isDouble() const;
  int isDoubleArray() const;
  int isRangeOrScalar() const;
  int isFirst() const;
  int isLast() const;
  long range_min() const;
  long range_max() const;
  long range_dlt() const;
  long const * getBuffer() const;
  double const * getDoubleBuffer() const;
  long const * getDims() const;
};


#endif
