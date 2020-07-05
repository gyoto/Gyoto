/*
    Copyright 2020 Thibaut Paumard & Frédéric Vincent

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

#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoShift.h"
#include "GyotoError.h"
#include "GyotoProperty.h"

using namespace std ; 
using namespace Gyoto ; 
using namespace Gyoto::Metric ; 

//// Property list:
//
// Note that none of those lines ends with punctation. "," and ";" are
// added by the macros where needed. Three steps:
//  1- GYOTO_PROPERTY_START(<classname>)
//  2- For each Property we want to support, a line such as:
//       GYOTO_PROPERTY_<type>(<classname>, <propertyname>, <accessorname>)
//     Note that the BOOL type is a bit special: the argument
//     <propertyname> is replaced by two arguments: <name_if_true> and
//     <name_if_false>.
//  3- GYOTO_PROPERTY_END(<classname>, <pointer to parent's property list>)
//
////
GYOTO_PROPERTY_START(Shift,
		     "Shift space-time.")
GYOTO_PROPERTY_METRIC(Shift, SubMetric, subMetric,
   "The underlying space-time to shift.")
GYOTO_PROPERTY_VECTOR_DOUBLE(Shift, Offset, offset,
			     "Amount by which to shift (4 components).")
GYOTO_PROPERTY_END(Shift, Generic::properties)

// This is the minimal constructor: it just sets the coordinate kind and
// the metric kind name.
Shift::Shift() :
Generic(GYOTO_COORDKIND_CARTESIAN, "Shift"),
  submet_(NULL)
{for (int i=0; i<4; ++i) offset_[i]=0.;}
Shift* Shift::clone() const { return new Shift(*this); }
Shift::~Shift(){
  if (submet_) submet_->unhook(this);
}

SmartPointer<Metric::Generic> Shift::subMetric() const { return submet_; }
void Shift::subMetric(SmartPointer<Metric::Generic> submet) {
  if (submet_) submet_->unhook(this);
  submet_=submet;
  if (submet_) {
    submet_->hook(this);
    mass(submet_->mass());
  }
}

void Shift::mass(double m) {
  submet_->mass(m);
  // hook system will also set mass_;
}

void Shift::tell(Hook::Teller* msg) {
  if (msg==submet_) Generic::mass(submet_->mass());
}

void Shift::offset(std::vector<double> const &v) {
  GYOTO_DEBUG_EXPR(v.size());
  if (v.size() !=4)
    GYOTO_ERROR("Shift offset needs exactly 4 tokens"); 
  for (int i=0; i<4; ++i) offset_[i]=v[i];
}

std::vector<double> Shift::offset() const {
  std::vector<double> res(4, 0.);
  for (int i=0; i<4; ++i) res[i]=offset_[i];
  return res;
}


void Shift::gmunu(double g[4][4], const double pos[4]) const
{
  double spos[4] = {pos[0]-offset_[0],
		    pos[1]-offset_[1],
		    pos[2]-offset_[2],
		    pos[3]-offset_[3]};
  submet_->gmunu(g,spos);
}

void Shift::gmunu_up(double gup[4][4], const double pos[4]) const
{
  double spos[4] = {pos[0]-offset_[0],
		    pos[1]-offset_[1],
		    pos[2]-offset_[2],
		    pos[3]-offset_[3]};
  submet_->gmunu_up(gup,spos);
}

void Shift::jacobian(double jac[4][4][4], const double pos[4]) const {
  double spos[4] = {pos[0]-offset_[0],
		    pos[1]-offset_[1],
		    pos[2]-offset_[2],
		    pos[3]-offset_[3]};
  submet_->jacobian(jac,spos);
}

int Shift::isStopCondition(double const coord[8]) const {
  double scoord[8] = {coord[0]-offset_[0],
		      coord[1]-offset_[1],
		      coord[2]-offset_[2],
		      coord[3]-offset_[3],
		      coord[4], coord[5], coord[6], coord[7]};
  return submet_->isStopCondition(scoord);
}


#ifdef GYOTO_USE_XERCES
void Shift::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
  if (p.type == Property::metric_t && submet_) {
    FactoryMessenger * childfmp = fmp -> makeChild ( "SubMetric" );
    submet_ -> fillElement(childfmp);
    delete childfmp;
  } else Object::fillProperty(fmp, p);
}

void Shift::setParameters(FactoryMessenger *fmp) {
  GYOTO_DEBUG << std::endl;

  std::string name="", content="", unit="";
  std::vector<std::string> plugin;
  FactoryMessenger * child = NULL;

  while (fmp->getNextParameter(&name, &content, &unit)) {
    GYOTO_DEBUG_EXPR(name);
    if (name=="SubMetric") {
      content = fmp -> getAttribute("kind");
      plugin  = split(fmp -> getAttribute("plugin"), ",");
      child   = fmp -> getChild();
      subMetric ((*Metric::getSubcontractor(content, plugin))(child, plugin));
      delete child;
    } else setParameter(name, content, unit);
  }

  GYOTO_DEBUG << "done" << std::endl;
}
#endif
