/*
    Copyright 2014 Thibaut Paumard

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
  This file contains helpers for the various check*.i files in the
  Gyoto source.
*/

#include "gyoto.i"
#include "gyoto_std.i"
#include "matrix.i"
#include "linalg.i"


func doing(msg, dots) {
  if (is_void(dots)) dots=1;
  if (dots) format="%s...";
  else format="%s";
  write, format=format, msg;
}
func done {write, format="%s\n", " done.";}
func dot {write, format="%s", ".";}
func output(msg) {write, format="%s\n", msg;}

func hline(a,w)
{
  if (is_void(a)) a="_";
  if (is_void(w)) w=80;
  b="";
  for (i=1; i<=w; ++i) b+=a;
  output, b;
}

func begin_section(name, comment) {
  hline, "*";
  if (is_void(comment)) comment="";
  else comment=" ("+comment+")";
  output, "CHECKING "+strcase(1, name)+comment;
  hline, "*";
}

func end_section(name, comment) {
  hline, "*";
  if (is_void(comment)) comment="";
  else comment=" ("+comment+")";
  output, "ALL TESTS PASSED FOR "+strcase(1, name)+comment;
  hline, "*";
}

func jacobian(met, pos, eps=)
/* DOCUMENT jac=jacobian(met, pos[, eps=])

     Compute the Jacobian matrix of metric met, i.e.:

       jac(nu, mu, a) := d(met(pos)(nu, mu))/dx^a

     Derivatives are evaluated by finite difference using the
     infinitesimal step EPS (default: 1e-10).
   
 */
{
  grad=array(double, 4, 4, 4);
  if (is_void(eps)) eps=1e-6;
  
  for (i=1; i<=4; ++i) {
    delta=array(0., 4);
    delta(i)=eps;
    grad(,,i)=(met(pos+delta)-met(pos-delta))/(2.*eps);
  }

  return grad;
}

func gmunu_up(met, pos)
/* DOCUMENT gup=gmunu_up(met, pos)

   Return gup, the inverse matrix of the Metric:

     met(pos, nu, mu) := g_mu_nu
     gup(nu, mu):= g^mu^nu

   It is computed using LUsolve().
   
 */
{
  g=met(pos);
  return LUsolve(g);
}

func christoffel(met, pos, eps=)
/* DOCUMENT Gamma = christoffel(met, pos)

   Compute the 64 Christoffel symbols of metric MET:

     Gamma(nu, mu, a) := Gamma^a_mu_nu

   Uses gmunu_up() and jacobian(). The computation is done using only
   the method gmunu in the metric and can therefore be used for a
   consistency check with the christoffel method.
 */
{
  res=array(double, 4, 4, 4);
  
  gup=gmunu_up(met, pos);
  jac=jacobian(met, pos, eps=eps);
  
  for (a=1; a<=4; ++a) {
    for (mu=1; mu<=4; ++mu) {
      for (nu=1; nu<=4; ++nu) {
        for (i=1; i<=4; ++i) {
          res(nu, mu, a) += 0.5*gup(a, i)*( jac(nu, i, mu)
                                           +jac(i, mu, nu)
                                           -jac(nu, mu, i));
        }
      }
    }
  }
  return res;
}

func check_christoffels(gg, pos, tolerance=, eps=)
/* DOCUMENT check_christoffels, GG, positions

     Check that the two forms of the christoffel method are consistent
     with each other and that they are consistent with the gmunu
     method.
   
 */
{
  doing, "checking christoffel methods", 0;
  if (is_void(tolerance)) tolerance=1e-6;
  default_positions, gg, pos;
  d=dimsof(pos);

  for (n=1; n<=d(3); ++n) {
    Gamma1=gg(christoffel=pos(,n));
    Gamma2=array(double, 4, 4, 4);
    for (a=1; a<=4; ++a)
      for (i=1; i<=4; ++i)
        for (j=1; j<=4; ++j)
          Gamma2(j,i,a)=gg(christoffel=pos(,n), j, i, a);
    Gamma3=christoffel(gg, pos(,n), eps=eps);
    if (anyof(Gamma1!=Gamma2))
      error, "The two forms of the christoffel method don't yield the same result";
    if (max(abs(Gamma1-Gamma3))>tolerance)
      error, "The Christoffels don't agree with their numerical estimate";
    dot;
  }
  done;
}

func default_positions(gg, &pos)
{
  if (is_void(pos)) pos=[[0., 10., 1., 2.]];
  d=dimsof(pos);
  if (d(1)==1) {
    pos=[ pos ];
  }
}

func check_gmunu(gg, pos, tolerance=)
/* DOCUMENT check_gmunu, gg, positions

     Check that the two forms of the gmunu method are consistent
     with each other and that the matrix is symmetric.
   
 */
{
  if (is_void(tolerance)) tolerance=1e-15;
  default_positions, gg, pos;
  d=dimsof(pos);
  
  doing, "checking metric coefficients", 0;
  for (n=1; n<=d(3); ++n) {
    g=gg(pos);
    g2=array(double, 4, 4);
    for (i=1; i<=4; ++i) g2(i, )=gg(pos, i, );
    gm1=g*0.+1.;
    gm1(where(g))=g(where(g));
    if (max(abs((g2-g)*gm1)) > tolerance)
      error, "The two forms of the gmunu method don't yield the same result";
    if (anyof(g != transpose(g)))
      error, "The metric is not symmetric";
    dot;
  }
  done;
}

func check_gmunu_up(gg, pos, tolerance=)
/* DOCUMENT check_gmunu_up, gg, positions

     Check that the gmunu_up method yields the inverse of the gmunu
     method.

     Note that not all metric kinds implement the gmunu_up method.
   
 */
{
  doing, "checking gmunu_up method", 0;
  if (is_void(tolerance)) tolerance=1e-15;
  default_positions, gg, pos;
  d=dimsof(pos);

  for (n=1; n<=d(3); ++n) {
    g=gg(pos(,n));
    gup=gg(gmunu_up=pos(,n));
    gup2=array(double, 4, 4);
    for (i=1; i<=4; ++i) gup2(i, )=gg(gmunu_up=pos(,n), i, );
    prod=g(,+)*gup(+,);
    gm1=gup*0.+1.;
    gm1(where(gup))=gup(where(gup));
    if (max(abs((gup2-gup)*gm1)) > tolerance)
      error, "The two forms of the gmunu_up method don't yield the same result";
    if (max(abs((prod-diag([1., 1., 1., 1.])))) > tolerance)
      error, "gmunu_up is not the inverse of gmunu";
    dot;
  }
  done;
}
