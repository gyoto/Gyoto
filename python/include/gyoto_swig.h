#include "GyotoConfig.h"
#include "GyotoDefs.h"
#include "GyotoFactory.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoValue.h"
#include "GyotoProperty.h"
#include "GyotoObject.h"
#include "GyotoAstrobj.h"
#include "GyotoError.h"
#include "GyotoWorldline.h"
#include "GyotoPhoton.h"
#include "GyotoScreen.h"
#include "GyotoThinDisk.h"
#include "GyotoStandardAstrobj.h"
#include "GyotoSpectrometer.h"
#include "GyotoComplexSpectrometer.h"
#include "GyotoUniformSpectrometer.h"
#include "GyotoRegister.h"
#include "GyotoWIP.h"
#include "GyotoConverters.h"
#include "GyotoGridData2D.h"
#include "GyotoFitsRW.h"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

swig_type_info * __Gyoto_SWIGTYPE_p_Gyoto__Error() ;
