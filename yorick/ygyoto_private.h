/*
    Copyright 2011-2014, 2016 Thibaut Paumard

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

/*
  This header should be used only for the main gyoto.so plug-in, not
  in additional plug-ins such as gyoto_std.so.
 */

#ifndef __YGYOTO_PRIVATE_H
#define __YGYOTO_PRIVATE_H

#include "GyotoSmartPointer.h"

/*
  The following are to declare a new base (such as Metric or Astrobj)
  or independent (such as Photon) class. Seldom used in a Gyoto
  extension.
 */
#ifdef GYOTO_USE_XERCES
# define YGYOTO_PRINT_YUSEROBJ(NAME)				\
  void gyoto_##NAME##_print(void *obj) {				\
    string rest="", sub="";					\
    size_t pos=0, len=0;					\
    try {rest = Factory(((gyoto_##NAME*)obj)->smptr).format();}	\
    YGYOTO_STD_CATCH;						\
    while ( ( len=rest.length() ) )  {				\
      sub=rest.substr(0, pos=rest.find_first_of("\n",0));	\
      rest=rest.substr(pos+1, len-1);				\
      y_print( sub.c_str(),1 );					\
    }								\
  }
#else
# define YGYOTO_PRINT_YUSEROBJ(NAME)    \
  void gyoto_##NAME##_print(void *obj) { \
    y_print("GYOTO " #NAME,0);	       \
  }
#endif


// Two helper functions for the dot operator below;
// these are *not* public API, don't (ab)use them
/* return "__gyoto_var<id>" */
char const * const __ygyoto_var_name(size_t id);
/* return the index of "__gyoto_var<id>" */
long int __ygyoto_var_idx(size_t id);
 
#define YGYOTO_YUSEROBJ(NAME, CLASS)					\
  extern "C" {								\
    typedef struct gyoto_##NAME {					\
      Gyoto::SmartPointer<CLASS> smptr;					\
    } gyoto_##NAME;							\
    typedef struct gyoto_##NAME##_closure {				\
      Gyoto::SmartPointer<CLASS> smptr;					\
      char * member;						\
    } gyoto_##NAME##_closure;						\
    void gyoto_##NAME##_free(void *obj) {				\
      if (((gyoto_##NAME*)obj)->smptr) {				\
	((gyoto_##NAME*)obj)->smptr=NULL;				\
      } else printf("null pointer\n");					\
    }									\
    void gyoto_##NAME##_closure_free(void *obj) {			\
      if (((gyoto_##NAME##_closure*)obj)->smptr) {			\
	((gyoto_##NAME##_closure*)obj)->smptr=NULL;			\
	p_free(((gyoto_##NAME##_closure*)obj)->member);		\
      } else printf("null pointer\n");					\
    }									\
    YGYOTO_PRINT_YUSEROBJ(NAME)						\
    void gyoto_##NAME##_closure_print(void *obj) {			\
      std::string phrase = "Gyoto closure. Class: \"" #NAME "\", method: \"";	\
      phrase.append(((gyoto_##NAME##_closure*)obj)->member).append("\""); \
      y_print(phrase.c_str(),1);					\
      y_print("(Hint: I'm a functor, call me as a function)", 0);	\
    }									\
    void gyoto_##NAME##_eval(void *obj, int argc);			\
    void gyoto_##NAME##_closure_eval(void *obj, int argc) {		\
      Property const * prop =						\
	((Gyoto::Object*)((gyoto_##NAME##_closure*)obj)->smptr())	\
	->property(((gyoto_##NAME##_closure*)obj)->member);		\
      GYOTO_DEBUG_EXPR(prop);						\
	if (prop) {							\
	  std::string unit="";						\
	  std::string kwd="";						\
	  int parg=-1;							\
	  long kidx=0;							\
	  for (int iarg =argc-1; iarg>=0; --iarg) {			\
	    if ((kidx=yarg_key(iarg))>=0) {				\
	      /* this is a keyword */					\
	      if (strcmp(yfind_name(kidx),"unit"))			\
		y_error("Only the 'unit' keyword is supported");	\
	      unit=ygets_q(--iarg);					\
	    } else {							\
	      if (parg!=-1)						\
		y_error("Only one positional argument accepted");	\
	      parg=iarg;						\
	    }								\
	  }								\
	  if (yarg_nil(parg)) parg=-1;					\
	  if (parg==-1)							\
	    ypush_property(((gyoto_##NAME##_closure*)obj)->smptr,	\
			   *prop,					\
			   ((gyoto_##NAME##_closure*)obj)->member,	\
			   unit);					\
	  else	{							\
	    yget_property(((gyoto_##NAME##_closure*)obj)->smptr,	\
			  *prop, parg,					\
			   ((gyoto_##NAME##_closure*)obj)->member,	\
			  unit);					\
	    *ypush_##NAME()= ((gyoto_##NAME##_closure*)obj)->smptr;	\
	  }								\
	  return;							\
	}								\
      long used_var=0;							\
      /* we build a yorick statement into ss */				\
      stringstream ss;							\
      /* first variable (used_var=0) will be output */			\
      ss << "eq_nocopy, " <<  __ygyoto_var_name(used_var++) << ", ";	\
      /* push object second variable */					\
      *ypush_##NAME()= ((gyoto_##NAME##_closure*)obj)->smptr;		\
      yput_global(__ygyoto_var_idx(used_var), 0);			\
      yarg_drop(1);							\
      ss << __ygyoto_var_name(used_var++) << "(";			\
      /* use "member=" keyword */					\
      ss << ((gyoto_##NAME##_closure*)obj)->member << "=";		\
      bool coma=false;							\
      /* process arguments */						\
      long kidx=0;							\
      for (int iarg =argc-1; iarg>=0; --iarg) {				\
	if ((kidx=yarg_key(iarg))>=0) {					\
	  /* this is a keyword */					\
	  ss << ", " << yfind_name(kidx) << "=";			\
	  coma=false;							\
	} else {							\
	  /* add coma if preceding argument was not a keyword */	\
	  if (coma) ss << ", ";						\
	  ypush_use(yget_use(iarg));					\
	  yput_global(__ygyoto_var_idx(used_var), 0);			\
	  yarg_drop(1);							\
	  ss << __ygyoto_var_name(used_var++);				\
	  coma=true;							\
	}								\
      }									\
      /* terminate statement */						\
      ss << ");";							\
      /* interpret statement */						\
      long dims[Y_DIMSIZE]={1,1};					\
      *ypush_q(dims)=p_strcpy(ss.str().c_str());			\
      yexec_include(0, 1);						\
      yarg_drop(1);							\
      /* result is in var #0, push it on the stack */			\
      ypush_global(__ygyoto_var_idx(0));				\
      /* clean variables */						\
      ypush_nil();							\
      for (long k=0; k<used_var; ++k)					\
	yput_global(__ygyoto_var_idx(k), 0);				\
      yarg_drop(1);							\
    }									\
    void gyoto_##NAME##_closure_extract(void *obj, char*member) {	\
      long idxo = yget_global("__gyoto_obj", 0);			\
      long idxr = yget_global("__gyoto_res", 0);			\
      *ypush_##NAME()= ((gyoto_##NAME##_closure*)obj)->smptr;		\
      yput_global(idxo, 0);						\
      yarg_drop(1);							\
      long dims[Y_DIMSIZE]={1,1};					\
      string stmt = "eq_nocopy, __gyoto_res, __gyoto_obj(";		\
      stmt.append(((gyoto_##NAME##_closure*)obj)->member).append("=).")\
	.append(member);						\
      *ypush_q(dims)=p_strcpy(stmt.c_str());				\
      yexec_include(0, 1);						\
      yarg_drop(1);							\
      ypush_global(idxr);						\
    }									\
    static y_userobj_t gyoto_##NAME##_closure_obj =			\
    {const_cast<char*>("gyoto_" #NAME),					\
     &gyoto_##NAME##_closure_free,					\
     &gyoto_##NAME##_closure_print,					\
     &gyoto_##NAME##_closure_eval,					\
     &gyoto_##NAME##_closure_extract,					\
     0};								\
    void gyoto_##NAME##_extract(void *obj, char *member) {		\
      gyoto_##NAME##_closure * CLOSURE =				\
	(gyoto_##NAME##_closure *)ypush_obj(&gyoto_##NAME##_closure_obj, \
					  sizeof(gyoto_##NAME##_closure)); \
      CLOSURE -> smptr=((gyoto_##NAME*)obj)->smptr;			\
      CLOSURE -> member = p_strcpy(member);				\
    }								\
    static y_userobj_t gyoto_##NAME##_obj =				\
    {const_cast<char*>("gyoto_" #NAME),					\
     &gyoto_##NAME##_free,						\
     &gyoto_##NAME##_print,						\
     &gyoto_##NAME##_eval,						\
     &gyoto_##NAME##_extract,						\
     0};								\
  }									\
  Gyoto::SmartPointer<CLASS>* yget_##NAME(int iarg) {			\
    return &(((gyoto_##NAME*)yget_obj(iarg, &gyoto_##NAME##_obj))->smptr); \
  }									\
  Gyoto::SmartPointer<CLASS>* ypush_##NAME() {				\
  return &(((gyoto_##NAME*)ypush_obj(&gyoto_##NAME##_obj,		\
				     sizeof(gyoto_##NAME)))->smptr);	\
  }									\
  int yarg_##NAME(int iarg) {						\
    return yget_obj(iarg,0)==gyoto_##NAME##_obj.type_name;		\
  }									\
  extern "C" {								\
    void Y_is_gyoto_##NAME(int)					\
    {									\
      ypush_long(yarg_##NAME(0));					\
    }									\
  }

#ifdef GYOTO_USE_XERCES
# define YGYOTO_BASE_CONSTRUCTOR1_XML(GETBASE)	\
  GYOTO_DEBUG << "found no subcontractor for \"" << fname		\
  << "\", calling Factory now\n";					\
  *OBJ = Factory(fname).GETBASE();
#else
# define YGYOTO_BASE_CONSTRUCTOR1_XML(GETBASE)	\
  y_error("no XML support in this gyoto");
#endif

#define YGYOTO_BASE_CONSTRUCTOR1(BASE,GETBASE)					\
  extern "C" {								\
  void Y_gyoto_##BASE(int argc)	{					\
    Gyoto::SmartPointer<Gyoto::BASE::Generic> *OBJ = NULL;		\
    if (yarg_##BASE(argc-1)) {						\
      OBJ = yget_##BASE(argc);						\
    } else {								\
      if (!yarg_string(argc-1))						\
        y_error("Cannot allocate object of virtual class " #BASE);	\
      char * fname = ygets_q(argc-1);					\
      std::vector<std::string> plugin;					\
      if (argc >= 2 && yarg_string(argc-2)) {				\
	long ntot=0;							\
	ystring_t * plugs = ygeta_q(argc-2, &ntot, NULL);		\
	for (size_t i=0; i<ntot; ++i) plugin.push_back(plugs[i]);	\
      }									\
      OBJ = ypush_##BASE();						\
      Gyoto::BASE::Subcontractor_t * sub =				\
	Gyoto::BASE::getSubcontractor(fname, plugin, 1);		\
      if (sub) {							\
	GYOTO_DEBUG << "found a subcontractor for \"" << fname		\
		    << "\", calling it now\n";				\
	*OBJ = (*sub)(NULL, plugin);					\
      } else {								\
	YGYOTO_BASE_CONSTRUCTOR1_XML(GETBASE);				\
      }									\
      yarg_swap(0, argc);						\
      yarg_drop(1);							\
    }									\
    --argc;								\
    gyoto_##BASE##_eval(OBJ, argc);					\
  }									\
}
#define YGYOTO_BASE_CONSTRUCTOR(BASE)					\
  YGYOTO_BASE_CONSTRUCTOR1(BASE,get##BASE)


#endif
