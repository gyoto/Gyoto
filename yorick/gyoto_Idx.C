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

using namespace YGyoto;

int YGyoto::Idx::isNuller() {return _is_nuller;}
long YGyoto::Idx::getNElements() {return _nel;}
long YGyoto::Idx::current() {
  if (_is_list) return _idx[_cur];
  return _cur;
}
double YGyoto::Idx::getDVal() {return _is_double?_dval:_range[0];}
int YGyoto::Idx::isDouble() {return _is_double;}
int YGyoto::Idx::isFirst() {return _is_first;}

int YGyoto::Idx::isLast() {
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

int YGyoto::Idx::valid() {
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


YGyoto::Idx::Idx(int iarg, int res) :
  _is_nuller(0), _is_range(0), _is_list(0), _is_scalar(0), _is_double(0)
{
  int flags = yget_range(iarg, _range), ndims=0;
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
    return;
  }
  if (yarg_rank(iarg) > 0) {
    _is_list=1;
    _nel=1;
    _idx = ygeta_l(iarg, &_nel, 0);
    return;
  }
  if (yarg_number(iarg)==1) {
    _is_scalar=1;
    long val=ygets_l(iarg);
    if (val>res) y_error("max index too large");
    if (val<=0) val+=res;
    _range[0]=_range[1]=val;
    _range[2]=1;
    _nel=1;
    return;
  }
  if (yarg_number(iarg)==2) {
    _is_scalar=1;
    _is_double=1;
    _dval=ygets_d(iarg);
    return;
  }
  if (iarg<0 || yarg_nil(iarg)) {
      _is_range=1;
      _range[0]=1;
      _range[1]=res;
      _range[2]=1;
      _nel=res;
      return;
  }
  y_error("unsupported range syntax");
}

int YGyoto::Idx::getNDims() {
  if (_is_range) return 1;
  if (_is_list) return 1;
  if (_is_scalar) return 0;
}
