plug_in, "gyoto";
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


//////// SCENERY

extern gyoto_Scenery;
/* DOCUMENT scenery = gyoto_Scenery([filename,] [members=values ...])
             Create GYOTO Scenery object
         or scenery, [members=values]
             Set GYOTO Scenery member
         or res = scenery(member=)
             Get GYOTO Scenery member
         or scenery, xmlwrite=filename
             Save Scenery description to XML
         or data = scenery([ irange, jrange [, quantities ]])
             Ray-trace scenery

   PURPOSE:
    Create and use GYOTO Scenery objects.

    A GYOTO Scenery is used to render an image by relativistic
    ray-tracing. A GYOTO Scenery contains references to:
     - a GYOTO Metric ("the stage");
     - a GYOTO Astrobj ("the actors");
     - a GYOTO Screen ("the camera").

    See GYOTO for basic concepts and syntax.


   MEMBERS:

    metric=  see gyoto_Metric(): what "straight" means for light
             travel;
             
    screen=  see gyoto_Screen(), it specifies where the obseerver is
             located and the obseerving time;
             
    astrobj= see gyoto_Astrobj(): where the light comes from;
    
    delta=   a double scalar, the initial integration step for the
             Photons laucnched during ray-tracing;
             
    quantities= an array of strings giving the names of physical
             quantities that should be retrieved during
             ray-tracing. Beware that not all Astrobj kinds support
             all of those quantities. The order in which the
             quantities are listed when setting this member is not
             relevant, the following are equivalent:
                 data = scenery(quantities=["EmissionTime", "Intensity"])();
             and
                 data = scenery(quantities=["Intensity", "EmissionTime"])();
             Setting quantities here is not mandatory as the third
             positional argument used for ray-tracing permits to
             override it in an ordered fashion.

             Recognized quantities:
               "Intensity": apparent intensity of the Astrobj;
               "EmissionTime": time at which each photonreaching the
                    screen was emitted;
               "MinDistance": minimum distance ever reached between
                    each photon (whether comingfrom the object or not)
                    and the Astrobj);
               "FirstDistMin": First local minimum in the
                    Astrobj/Photon distance;
               "Redshift": ;
               "ImpactR", "ImpactX", "ImpactY" and "ImpactZ": R
                    (specrical), X, Y and Z (Cartsesian) coordinates
                    at which a photon was emitted by the Astrobj;
               "Spectrum": Inu spectrum of the Astrobj on this pixel,
                    the spectrometer is specified in the Screen
                    object;
               "BinSpectrum": spectrum of the Astrobj on this pixel
                    (the spectrometer is specified in the Screen
                    object), as would be detected by a real
                    spectrometer: the value in each spectral channel
                    is the integral of Inu over the spectral channel;
               "User1" to "User5": other specific scalar quantities an
                    Astrobj may be able to compute, refer to the
                    documentation for the Astrobj kind of your choice.

    nthreads=number of parallel threads to use in
             gyoto_Scenery_rayTrace. This has no effect when
             ray-tracing using the "data = scenery()" syntax below.
                    
    RAY-TRACING:
    
     Ray-traced data is retrieved calling the object like a function
     with no keyword:
        data = scenery ();
     or data = scenery (irange, jrange, quant )
     
    IRANGE and JRANGE are either scalars (i, j) or ranges in the usual
    form min:max:step. QUANT is an array of Yorick strings where each
    element selects one quantity to retrieve (see the QUANTITIES
    member above). QUANT may be void to use the quantities already set
    in the Scenery or the default for the Astrobj. If specifyng QUANT,
    DATA will be a MxNxP double array, where P=numberof(QUANT) and
    DATA(,,i) will contain the value of the quantity specifeid by
    QUANT(i).

    QUANTITIES may also be a scalar to retrieve a single
    quantity.

    The "Spectrum" quantity is a bit peculiar since it take more than
    one plane in data.
    
   SEE ALSO:
     gyoto_Metric, gyoto_Screen, gyoto_Astrobj, gyoto_Photon,
     gyoto_Spectrometer, gyoto_Scenery_rayTrace
*/

extern gyoto_Scenery_rayTrace
/* DOCUMENT res = gyoto_Scenery_rayTrace(scenery, imin, imax, jmin, jmax,
                                         impactcoords)

     if IMPACTCOORDS is an unadorned, nil variable it is output. If it
     is an expression or non-nil, it is input.
 */

func _gyoto_Scenery_adaptive_raytrace(sco, respmax, &computed) {
/* xDOCUMENT data = gyoto_Scenery_adaptive_raytrace(scenery, pmax, [computed])

BROKEN
   
     Ray-trace a GYOTO Scenery on an adaptive grid.

     For certain kinds of objects (in particular, Stars), this routine
     is much faster than the equivalent:
       data = scenery(resolution=3^pmax, raytrace=1);
     It is NOT guaranteed that the two methods yield the same result.

     The minimum distance between photon and object is first computed
     on a coarse grid which is then refined as required.

   SEE ALSO: gyoto_Scenery
 */
  write, format="%s\n",
    "WARNING: gyoto_Scenery_adaptive_raytrace() is under development";

  sc = sco(clone=); // don't modify calling object
  
  DBL_MAX=1e100;

  screen = sc(screen=);

  respmax=long(respmax);
  
  resp=1;
  resmax=3^respmax;
  step=long(3^(respmax-resp));
  first=step/2+1;

  data=array(double, resmax, resmax, 6);

  quantities = ["Intensity", "EmissionTime", "MinDistance", "ImpactX", "ImpactY", "ImpactZ"];
  
  screen, resolution=resmax;
  data(first::step, first::step, ) =
    sc(first::step, first::step, quantities);
  
  computed=array(long,resmax,resmax);

  nb=9;
  for (resp=2, res=3; resp<=respmax; ++resp) {

    // largest distance to neighbour
    dsub=data(first::step, first::step,);
    ind=where(dsub>=DBL_MAX);
    if (numberof(ind)) dsub(ind)=DBL_MAX;
    delta = array(double, res, res);
    d1=(dsub(dif,,4:6)^2)(,,sum);
    d2=(dsub(,dif,4:6)^2)(,,sum);
    delta(2:-1,2:-1)=[d1(:-1,2:-1), d1(2:,2:-1), d2(2:-1,:-1), d2(2:-1,2:)](,,max);
    delta(1,2:-1)=[d1(1,2:-1), d2(1,:-1), d2(1,2:)](,max);
    delta(0,2:-1)=[d1(0,2:-1), d2(0,:-1), d2(0,2:)](,max);
    delta(2:-1,1)=[d2(2:-1,1), d1(:-1,1), d1(2:,1)](,max);
    delta(2:-1,0)=[d2(2:-1,0), d1(:-1,0), d1(2:,0)](,max);
    delta(0,0)=[d1(0,0), d2(0,0)](max);
    delta(1,1)=[d1(1,1), d2(1,1)](max);
    delta(0,1)=[d1(0,1), d2(0,1)](max);
    delta(1,0)=[d1(1,0), d2(1,0)](max);


    // ! BEWARE : res is updated here
    res*=3;
    refine=array(int, res, res);
    refine(1::3, 1::3) =
    refine(2::3, 1::3) =
    refine(3::3, 1::3) =
    refine(1::3, 2::3) =
    refine(2::3, 2::3) =
    refine(3::3, 2::3) =
    refine(1::3, 3::3) =
    refine(2::3, 3::3) =
    refine(3::3, 3::3) =
      (dsub(,,3)<4*delta) | (dsub(,,3)<2);
    
    nstep=long(3^(respmax-resp));
    nfirst=first-nstep;
    data(nfirst     ::step,nfirst     ::step,)=
    data(first      ::step,nfirst     ::step,)=
    data(first+nstep::step,nfirst     ::step,)=
    data(nfirst     ::step,first      ::step,)=
    data(first+nstep::step,first      ::step,)=
    data(nfirst     ::step,first+nstep::step,)=
    data(first      ::step,first+nstep::step,)=
    data(first+nstep::step,first+nstep::step,)=
      dsub;

    step=nstep;
    first=nfirst;
    for (i=1; i<=res; ++i) {
      ibis=(i-1)*step+first;
      ind=where(refine(i,) & !computed(ibis,first::step));
      nb+=numberof(ind);
      if (numberof(ind)) {
        indbis=(ind-1)*step+first;
        data(ibis,indbis,)=sc(ibis, indbis, quantities);
        computed(ibis, indbis)=resp;
      }
    }

  }
  nb;
  return data;
}
