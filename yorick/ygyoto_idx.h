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

namespace YGyoto { class Idx; }

class YGyoto::Idx {
private:
  int _is_nuller;
  int _is_range;
  int _is_list;
  int _is_scalar;
  int _is_double;

  long _range[3];
  double _dval;
  long *_idx;
  long _nel;
  long _cur;
  int _valid;
  
public:
  Idx(int iarg, int res);
  int isNuller();
  int getNDims();
  long getNElements();
  long first();
  int  valid();
  long next();
  long current();
  double getDVal();
  int isDouble();
};


#endif
