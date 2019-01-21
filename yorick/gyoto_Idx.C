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

#include "ygyoto_idx.h"
#include "yapi.h"
#include "GyotoError.h"
#include "GyotoDefs.h"
#include "GyotoUtils.h"
#include <iostream>
using namespace YGyoto;
using namespace Gyoto;
using namespace std;
int YGyoto::Idx::isNuller() const {return _is_nuller;}
long YGyoto::Idx::getNElements() const {return _nel;}
long YGyoto::Idx::current() const {
  if (_is_list) return _idx[_cur];
  return _cur;
}
double YGyoto::Idx::getDVal() const {return _is_double?_dval:_range[0];}
int YGyoto::Idx::isDouble() const {return _is_double || _is_dlist;}
int YGyoto::Idx::isRangeOrScalar() const {return _is_range || _is_scalar;}
int YGyoto::Idx::isFirst() const {return _is_first;}

int YGyoto::Idx::isLast() const {
  if (_is_range) return _cur+_range[2] > _range[1];
  if (_is_scalar) return 1;
  if (_is_list) return _cur >= _nel;
  return 0;
}

long YGyoto::Idx::first() {
  _is_first=1;
  if (_is_range) return _cur=_range[0];
  if (_is_scalar) return _cur=_range[0];
  if (_is_list) return _idx[_cur=0];
  return 0;
}

int YGyoto::Idx::valid() const {
  if (_is_range) return _cur<=_range[1];
  if (_is_scalar) return _cur==_range[0];
  if (_is_list) return _cur<_nel;

  return 0;
}

long YGyoto::Idx::next() {
  _is_first=0;
  if (_is_range) return _cur+=_range[2];
  if (_is_scalar) return ++_cur;
  if (_is_list && (++_cur)<_nel) return _idx[_cur];

  return 0;
}

long YGyoto::Idx::range_min() const {
  if (!(_is_range || _is_scalar)) GYOTO_ERROR("BUG: not a range");
  return _range[0];
}

long YGyoto::Idx::range_max() const {
  if (!(_is_range || _is_scalar)) GYOTO_ERROR("BUG: not a range");
  return _range[1];
}

long YGyoto::Idx::range_dlt() const {
  if (!(_is_range || _is_scalar)) GYOTO_ERROR("BUG: not a range");
  return _range[2];
}

YGyoto::Idx::Idx(int iarg, int res) :
  _is_nuller(0), _is_range(0), _is_list(0), _is_scalar(0), _is_double(0),
  _is_dlist(0), _buf(NULL)
{
  int flags = yget_range(iarg, _range);
  if (flags) {
    _is_range=1;
    if (flags>=Y_MAX_DFLT) {
      flags-=Y_MAX_DFLT;
      _range[1]=res;
    }
    if (flags>=Y_MIN_DFLT) {
      flags-=Y_MIN_DFLT;
      _range[0]=1;
    }
    if (flags==Y_NULLER) {
      _is_nuller=1;
      flags=0;
      _nel=0;
    }
    if (flags>1)
      y_error("unsupported range syntax");

    if (_range[0]<=0) _range[0]+=res;
    if (_range[1]<=0) _range[1]+=res;
    if (_range[0]>res || _range[1]>res) y_error("max index too large");


    _nel=(_range[1]-_range[0]+_range[2])/_range[2];
    _dims[0]=1; _dims[1]=_nel;
    return;
  }
  if (yarg_number(iarg)==1){
    if (yarg_rank(iarg) > 0) {
      _is_list=1;
      _idx = ygeta_l(iarg, &_nel, _dims);
      return;
    }
    _is_scalar=1;
    long val=ygets_l(iarg);
    if (val>res) y_error("max index too large");
    if (val<=0) val+=res;
    _range[0]=_range[1]=val;
    _range[2]=1;
    _nel=1;
    _dims[0]=0;
    return;
  }
  if (yarg_number(iarg)==2) {
    _is_double=1;
    _buf=ygeta_d(iarg, &_nel, _dims);
    _dval=_buf[0];
    if (_dims[0]) _is_dlist=1;
    else _is_scalar=1;
    GYOTO_DEBUG_ARRAY(_dims, Y_DIMSIZE);
    GYOTO_DEBUG_EXPR(_is_scalar);
    GYOTO_DEBUG_EXPR(_is_dlist);
    return;
  }
  if (iarg<0 || yarg_nil(iarg)) {
      _is_range=1;
      _range[0]=1;
      _range[1]=res;
      _range[2]=1;
      _nel=res;
      _dims[0]=1;
      _dims[1]=_nel;
      return;
  }
  y_error("unsupported range syntax");
}

int YGyoto::Idx::getNDims() const {
  if (_is_range) return 1;
  if (_is_list) return 1;
  if (_is_scalar) return 0;
  GYOTO_ERROR("BUG: What does this YGyoto::Idx instance hold?");
  return 0;
}

long const * YGyoto::Idx::getDims() const {return _dims;}

long const * YGyoto::Idx::getBuffer() const {return _idx;}
double const * YGyoto::Idx::getDoubleBuffer() const {return _buf;}
