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


extern _gyoto_KerrBL_register_as_Metric;
_gyoto_KerrBL_register_as_Metric;
extern _gyoto_KerrKS_register_as_Metric;
_gyoto_KerrKS_register_as_Metric;
/* xDOCUMENT _gyoto_KerrBL_register_as_Metric
    should be called once so that gg(get_spin=1...) works.
 */

extern gyoto_KerrBL;
extern gyoto_KerrKS;
/* DOCUMENT gg = gyoto_KerrBL([filename,][members=values])
         or gg = gyoto_KerrKS([filename,][members=values])
         
         or gg, members=values
         or value = gg(member=)

         or coord = gg(makecoord=yinit, cst) // KerrBL only

         or coefs = gg(position, mu, nu)

   PURPOSE

     Instanciate and use a GYOTO KerrBL (for Kerr Boyer-Lindquist) or
     Kerr (for Kerr Kerr-Schild) Metric.

     For basics, first read "help, gyoto" and "help, gyoto_Metric".

   SPECIFIC MEMBER
   
     In addition to the generic Metric members, KerrBL and KerrKS
     metrics have a spin that can be set and retrieved as follows:
       gg, spin=value;
       value = gg(spin=)
       
   SPECIFIC METHOD

     In addition to the generic methods provided by gyoto_Metric,
     KerrBL metrics provide the following:
     
       coord=gg(get_coord=yinit, cst)
                get 8-coordinate (4-position & 4-velocity) COORD
                corresponding to the 6 coordinate (4-position & 2
                momenta) YINIT and the 3 motion constants in CST.

   SEE ALSO: gyoto, gyoto_Metric
 */
