/*
  Software from Michel & Stoitsov, Comput. Phys. Commun., 178, 535 (2008)
  Downloaded July 2018 from http://www.cpc.cs.qub.ac.uk/
 */

#define SIGN(a) (((a) < 0) ? (-1) : (1))

// Infinite norm of a std::complex number.
// ----------------------------------
// It is max(|Re[z]|,|Im[z]|)

inline double inf_norm (const std::complex<double> &z)
{
  return fmax (abs (real (z)),abs (imag (z)));
}


// Test of finiteness of a std::complex number
// --------------------------------------
// If real or imaginary parts are finite, true is returned.
// Otherwise, false is returned

inline bool isfinite (const std::complex<double> &z)
{
  const double x = real (z), y = imag (z);

  return (finite (x) && finite (y));
}




// Usual operator overloads of std::complex numbers with integers
// ---------------------------------------------------------
// Recent std::complex libraries do not accept for example z+n or z==n with n integer, signed or unsigned.
// The operator overload is done here, by simply putting a cast on double to the integer.

inline std::complex<double> operator + (const std::complex<double> &z,const int n)
{
  return (z+static_cast<double> (n));
}

inline std::complex<double> operator - (const std::complex<double> &z,const int n)
{
  return (z-static_cast<double> (n));
}

inline std::complex<double> operator * (const std::complex<double> &z,const int n)
{
  return (z*static_cast<double> (n));
}

inline std::complex<double> operator / (const std::complex<double> &z,const int n)
{
  return (z/static_cast<double> (n));
}

inline std::complex<double> operator + (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)+z);
}

inline std::complex<double> operator - (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)-z);
}

inline std::complex<double> operator * (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)*z);
}

inline std::complex<double> operator / (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)/z);
}








inline std::complex<double> operator + (const std::complex<double> &z,const unsigned int n)
{
  return (z+static_cast<double> (n));
}

inline std::complex<double> operator - (const std::complex<double> &z,const unsigned int n)
{
  return (z-static_cast<double> (n));
}

inline std::complex<double> operator * (const std::complex<double> &z,const unsigned int n)
{
  return (z*static_cast<double> (n));
}

inline std::complex<double> operator / (const std::complex<double> &z,const unsigned int n)
{
  return (z/static_cast<double> (n));
}

inline std::complex<double> operator + (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)+z);
}

inline std::complex<double> operator - (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)-z);
}

inline std::complex<double> operator * (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)*z);
}

inline std::complex<double> operator / (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n)/z);
}




inline bool operator == (const std::complex<double> &z,const int n)
{
  return (z == static_cast<double> (n));
}

inline bool operator != (const std::complex<double> &z,const int n)
{
  return (z != static_cast<double> (n));
}

inline bool operator == (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) == z);
}

inline bool operator != (const int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) != z);
}







inline bool operator == (const std::complex<double> &z,const unsigned int n)
{
  return (z == static_cast<double> (n));
}

inline bool operator != (const std::complex<double> &z,const unsigned int n)
{
  return (z != static_cast<double> (n));
}

inline bool operator == (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) == z);
}

inline bool operator != (const unsigned int n,const std::complex<double> &z)
{
  return (static_cast<double> (n) != z);
}






// Precise evaluation of exp[z]-1 for z std::complex
// --------------------------------------------
// When |Re[z]| >= 1 or |Im[z]| >= 1, one uses directly the standard exp function as it is precise.
// Otherwise, numerical cancellations can occur.
// So, one uses the always stable formula exp[z]-1 = expm1(x) - 2.exp(x).sin^2(y/2) + i.exp(x).sin(y) 
// with x = Re[z] and y = Im[z]. expm1(x) gives a precise evaluation of exp(x)-1 for x double.

std::complex<double> expm1 (const std::complex<double> &z)
{
  const double x = real (z),y = imag (z);

  if ((abs (x) >= 1.0) || (abs (y) >= 1.0)) return (exp (z) - 1.0);

  const double expm1_x = expm1 (x),exp_x = 1.0 + expm1_x,sin_y_over_two = sin (0.5*y),sin_y = sin (y);

  return std::complex<double> (expm1_x - 2.0*exp_x*sin_y_over_two*sin_y_over_two,exp_x*sin_y);
}



// Precise evaluation of log[1+z] for z std::complex
// --------------------------------------------
// When |Re[z]| >= 1 or |Im[z]| >= 1, one uses directly the standard log function as it is precise.
// Otherwise, numerical cancellations can occur.
// So, one uses the always stable formula log[1+z] = log1p(x) + log1p([y/(1+x)]^2)/2 + i.atan2(y,1+x)
// with x = Re[z] and y = Im[z]. log1p(x) gives a precise evaluation of log(1+x) for x double.
// atan2(x,y) gives the arc tangent of y/x so it is in ]-Pi:Pi]. 

std::complex<double> log1p (const std::complex<double> &z)
{
  const double x = real (z),y = imag (z); 

  const double xp1 = 1.0 + x,abs_x = abs (x), abs_y = abs (y);

  if ((abs_x >= 1.0) || (abs_y >= 1.0)) return log (1.0 + z);
 
  const double y_over_xp1 = y/xp1;

  return std::complex<double> (log1p (x) + 0.5*log1p (y_over_xp1*y_over_xp1),atan2 (y,xp1));
}




// Logarithm of Gamma[z] and Gamma inverse function
// ------------------------------------------------
// For log[Gamma[z]], if z is not finite or is a negative integer, the program returns an error message and stops.
// The Lanczos method is used. Precision : ~ 1E-15
// The method works for Re[z] >= 0.5 .
// If Re[z] <= 0.5, one uses the formula Gamma[z].Gamma[1-z] = Pi/sin (Pi.z).
// log[sin(Pi.z)] is calculated with the Kolbig method (K.S. Kolbig, Comp. Phys. Comm., Vol. 4, p.221 (1972)) : 
// If z = x+iy and y >= 0, log[sin(Pi.z)] = log[sin(Pi.eps)] - i.Pi.n, with z = n + eps so 0 <= Re[eps] < 1 and n integer.
// If y > 110, log[sin(Pi.z)] = -i.Pi.z + log[0.5] + i.Pi/2 numerically so that no overflow can occur.
// If z = x+iy and y < 0, log[Gamma(z)] = [log[Gamma(z*)]]*, so that one can use the previous formula with z*.
//
// For Gamma inverse, Lanczos method is also used with Euler reflection formula.
// sin (Pi.z) is calculated as sin (Pi.(z-n)) to avoid inaccuracy
// with z = n + eps with n integer and |eps| as small as possible.
//
// Variables:
// ----------
// x,y: Re[z], Im[z]
// log_sqrt_2Pi,log_Pi : log[sqrt(2.Pi)], log(Pi).
// sum : Rational function in the Lanczos method
// log_Gamma_z : log[Gamma(z)] value.
// c : table containing the fifteen coefficients in the expansion used in the Lanczos method.
// eps,n : z = n + eps so 0 <= Re[eps] < 1 and n integer for Log[Gamma].
//         z=n+eps and n integer so |eps| is as small as possible for Gamma_inv.
// log_const : log[0.5] + i.Pi/2
// g : coefficient used in the Lanczos formula. It is here 607/128.
// z,z_m_0p5,z_p_g_m0p5,zm1 : argument of the Gamma function, z-0.5, z-0.5+g, z-1 

std::complex<double> log_Gamma (const std::complex<double> &z)
{
  if (!isfinite (z)) std::cout<<"z is not finite in log_Gamma."<<std::endl, abort ();

  const double x = real (z),y = imag (z);

  if ((z == rint (x)) && (x <= 0)) std::cout<<"z is negative integer in log_Gamma."<<std::endl, abort ();

  if (x >= 0.5)
  {
    const double log_sqrt_2Pi = 0.91893853320467274177,g = 4.7421875;
    const std::complex<double> z_m_0p5 = z - 0.5, z_pg_m0p5 = z_m_0p5 + g, zm1 = z - 1.0;
    const double c[15] = {0.99999999999999709182,
			  57.156235665862923517,
			  -59.597960355475491248,
			  14.136097974741747174,
			  -0.49191381609762019978,
			  0.33994649984811888699E-4,
			  0.46523628927048575665E-4,
			  -0.98374475304879564677E-4,
			  0.15808870322491248884E-3,
			  -0.21026444172410488319E-3,
			  0.21743961811521264320E-3,
			 -0.16431810653676389022E-3,
			  0.84418223983852743293E-4,
			  -0.26190838401581408670E-4,
			  0.36899182659531622704E-5};
      
    std::complex<double> sum = c[0];
    for (int i = 1 ; i < 15 ; i++) sum += c[i]/(zm1 + i);

    const std::complex<double> log_Gamma_z = log_sqrt_2Pi + log (sum) + z_m_0p5*log (z_pg_m0p5) - z_pg_m0p5;

    return log_Gamma_z;
  }
  else if (y >= 0.0)
  {
    const int n = (x < rint (x)) ? (static_cast<int> (rint (x)) - 1) : (static_cast<int> (rint (x)));
    const double log_Pi = 1.1447298858494002;
    const std::complex<double> log_const(-M_LN2,M_PI_2),i_Pi(0.0,M_PI);
    const std::complex<double> eps = z - n,log_sin_Pi_z = (y > 110) ? (-i_Pi*z + log_const) : (log (sin (M_PI*eps)) - i_Pi*n);

    const std::complex<double> log_Gamma_z = log_Pi - log_sin_Pi_z - log_Gamma (1.0 - z);
 
    return log_Gamma_z;
  }
  else
    return conj (log_Gamma (conj (z)));
}




std::complex<double> Gamma_inv (const std::complex<double> &z)
{
  if (!isfinite (z)) std::cout<<"z is not finite in Gamma_inv."<<std::endl, abort ();

  const double x = real (z);

  if (x >= 0.5)
  {
    const double log_sqrt_2Pi = 0.91893853320467274177,g = 4.7421875;
    const std::complex<double> z_m_0p5 = z - 0.5, z_pg_m0p5 = z_m_0p5 + g, zm1 = z - 1.0;
    const double c[15] = {0.99999999999999709182,
			  57.156235665862923517,
			  -59.597960355475491248,
			  14.136097974741747174,
			  -0.49191381609762019978,
			  0.33994649984811888699E-4,
			  0.46523628927048575665E-4,
			  -0.98374475304879564677E-4,
			  0.15808870322491248884E-3,
			  -0.21026444172410488319E-3,
			  0.21743961811521264320E-3,
			 -0.16431810653676389022E-3,
			  0.84418223983852743293E-4,
			  -0.26190838401581408670E-4,
			  0.36899182659531622704E-5};
      
    std::complex<double> sum = c[0];
    for (int i = 1 ; i < 15 ; i++) sum += c[i]/(zm1 + i);

    const std::complex<double> Gamma_inv_z = exp (z_pg_m0p5 - z_m_0p5*log (z_pg_m0p5) - log_sqrt_2Pi)/sum;

    return Gamma_inv_z;
  }
  else
  {
    const int n = static_cast<int> (rint (x));

    const std::complex<double> eps = z - n;

    if (n%2 == 0)
      return (sin (M_PI*eps)*M_1_PI)/Gamma_inv (1.0 - z);
    else
      return (-sin (M_PI*eps)*M_1_PI)/Gamma_inv (1.0 - z);
  }
}



// Coulomb phase shift.
// --------------------
// It is given by the formula [Gamma[1+l+I.eta] - Gamma[1+l-I.eta]]/[2i].
// 0 is returned if 1+l+/-I.eta is a negative integer.
//
// Variables:
// ----------
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// Ieta,one_over_two_I : i.eta,1/[2.i] .
// arg_plus,arg_minus : 1+l+i.eta, 1+l-i.eta.
// log_Gamma_plus,log_Gamma_minus : logs of Gamma[1+l+I.eta], Gamma[1+l-I.eta].
// sigma_l : returned result.

std::complex<double> sigma_l_calc (const std::complex<double> &l,const std::complex<double> &eta)
{
  const std::complex<double> Ieta(-imag (eta),real (eta)),one_over_two_I(0,-0.5);
  const std::complex<double> arg_plus = 1.0 + l + Ieta, arg_minus = 1.0 + l - Ieta;

  if ((rint (real (arg_plus)) == arg_plus) && (rint (real (arg_plus)) <= 0.0)) return 0.0;
  if ((rint (real (arg_minus)) == arg_minus) && (rint (real (arg_minus)) <= 0.0)) return 0.0;

  const std::complex<double> log_Gamma_plus = log_Gamma (arg_plus),log_Gamma_minus = log_Gamma (arg_minus);  
  const std::complex<double> sigma_l = (log_Gamma_plus - log_Gamma_minus)*one_over_two_I;

  return sigma_l;
}






// log of C(l,eta)
// ---------------
// It is given by the formula l*log[2] - eta.Pi/2 + (log[Gamma[1+l+I.eta]] + log[Gamma[1+l-I.eta]])/2.0 - log[Gamma[2l+2]].
// 0 is returned if 1+l+/-I.eta is a negative integer.
// 2l+2 should not be a negative integer : one has to use -l-1 instead of l in this case.
//
// Variables:
// ----------
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// Ieta : i.eta .
// arg_plus,arg_minus : 1+l+i.eta, 1+l-i.eta.
// log_Gamma_plus,log_Gamma_minus,log_Gamma_2l_plus_2 : logs of Gamma[1+l+I.eta], Gamma[1+l-I.eta], Gamma[2l+2].
// log_Cl_eta : returned result.

std::complex<double> log_Cl_eta_calc (const std::complex<double> &l,const std::complex<double> &eta)
{
  const std::complex<double> Ieta(-imag (eta),real (eta));
  const std::complex<double> arg_plus = 1.0 + l + Ieta, arg_minus = 1.0 + l - Ieta; 

  if ((rint (real (arg_plus)) == arg_plus) && (rint (real (arg_plus)) <= 0.0)) return 0.0;
  if ((rint (real (arg_minus)) == arg_minus) && (rint (real (arg_minus)) <= 0.0)) return 0.0;

  const std::complex<double> log_Gamma_plus = log_Gamma (arg_plus),log_Gamma_minus = log_Gamma (arg_minus),log_Gamma_2l_plus_2 = log_Gamma (2.0*l + 2.0);
  const std::complex<double> log_Cl_eta = l*M_LN2 - M_PI_2*eta + 0.5*(log_Gamma_plus + log_Gamma_minus) - log_Gamma_2l_plus_2;

  return log_Cl_eta;
}



// Cut constant log for the asymptotic series.
// -------------------------------------------
// The asymptotic series and H[omega] behave differently near the negative real axis.
// Then, if one is in the bad quadrant of H[omega], one has to take into account the cut directly.
// One is in the bad quadrant of H[omega] if Re[z] < 0.0 and sign(Im[z]) = -omega.
// 
// H[omega] = [H[omega] from asymptotic function formula] + (1 - exp(2.i.Pi.(i.eta - l.omega))).[H[-omega] from asymptotic function formula]
//
//
// The cut constant log is then log [1 - exp(2.i.Pi.(i.eta - l.omega))].
// Its returned imaginary part is not necessarily in ]-Pi:Pi]
//
//
// Variables:
// ----------
// omega : 1 or -1.
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// Ieta : i.eta
// l_int,Ieta_int : closest integers to Re[l],Re[i.eta]
// eps :  (Ieta - Ieta_int) - (l.omega - l_int.omega)
// two_I_Pi, two_I_Pi_eps : 2.i.Pi, 2.i.Pi.eps .
// log_cut_constant : returned result.

std::complex<double> log_cut_constant_AS_calc (const int omega,const std::complex<double> &l,const std::complex<double> &eta)
{
  const std::complex<double> Ieta(-imag (eta),real (eta));
  const double l_int = rint (real (l)), Ieta_int = rint (real (Ieta));
  const std::complex<double> eps = (Ieta - Ieta_int) - omega*(l - l_int);

  const std::complex<double> two_I_Pi(0,2.0*M_PI),two_I_Pi_eps = two_I_Pi*eps;

  if (real (two_I_Pi_eps) > -0.1)
  {
    const std::complex<double> log_cut_constant = log (expm1 (-two_I_Pi_eps)) + two_I_Pi_eps;

    return log_cut_constant;
  }
  else
  {
    const std::complex<double> log_cut_constant = log1p (-exp (two_I_Pi_eps));

    return log_cut_constant;
  }
}





// Cut constant log for continued fractions : H[omega] from H[omega, not corrected] case.
// --------------------------------------------------------------------------------------
// The continued fraction has no cut on the negative real axis, whereas H[omega] has one.
// Then, if one is in the bad quadrant of H[omega], one has to take into account the cut directly.
// One is in the bad quadrant of H[omega] if Re[z] < 0.0 and sign(Im[z]) = -omega.
// 
// H[omega] = H[omega, not corrected] - cut_constant.F .
//
// The cut constant is 2i.omega.norm.(exp (2.i.Pi.[l.omega - i.eta]) - 1), and one takes its log.
// The imaginary part of the log is not necessarily in ]-Pi:Pi].
// Norm is 1.0 for normalized wave functions, C(l,eta)^2 for unnormalized wave functions.
//
//
// Variables:
// ----------  
// is_it_normalized : true if one wants normalized functions, i.e. the standard normalization,
//                    false if one wants F -> F/C(l,eta) and H+/H-/G -> H+/H-/G.C(l,eta), to avoid overflows for |eta| >> 1 and |z| small.
// omega : 1 or -1.
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// Ieta : i.eta .
// l_int,Ieta_int : closest integers to Re[l],Re[i.eta]
// eps : (l.omega - l_int.omega) - (Ieta - Ieta_int)
// log_norm : log[C(l,eta)^2] if is_it_normalized is false, 0.0 if it is true.
// two_I_Pi, two_I_Pi_eps : 2.i.Pi, 2.i.Pi.eps .
// log_two_I_omega : log[2.i.omega] = log[2] + i.omega.Pi/2 .
// log_cut_constant : returned result.

std::complex<double> log_cut_constant_CFa_calc (const bool is_it_normalized,const int omega,const std::complex<double> &l,const std::complex<double> &eta)
{
  const std::complex<double> Ieta(-imag (eta),real (eta));
  const double l_int = rint (real (l)), Ieta_int = rint (real (Ieta));
  const std::complex<double> eps = omega*(l - l_int) - (Ieta - Ieta_int);
  const std::complex<double> log_norm = (!is_it_normalized) ? (2.0*log_Cl_eta_calc (l,eta)) : (0.0);
  const std::complex<double> two_I_Pi(0,2.0*M_PI),two_I_Pi_eps = two_I_Pi*eps,log_two_I_omega(M_LN2,omega*M_PI_2);

  if (real (two_I_Pi_eps) < 0.1)
  {
    const std::complex<double> log_cut_constant = log_two_I_omega + log (expm1 (two_I_Pi_eps)) + log_norm;

    return log_cut_constant;
  }
  else
  {
    const std::complex<double> log_cut_constant = log_two_I_omega + log1p (-exp (-two_I_Pi_eps)) + two_I_Pi_eps + log_norm;

    return log_cut_constant; 
  }
}




// Cut constant log for continued fractions : H[omega] from H[-omega, not corrected] case.
// ---------------------------------------------------------------------------------------
// The continued fraction has no cut on the negative real axis, whereas H[-omega] has one.
// Then, if one is in the bad quadrant of H[-omega], one has to take into account the cut directly.
// One is in the bad quadrant of H[-omega] if Re[z] < 0.0 and sign(Im[z]) = omega.
// 
// H[omega] = H[-omega, not corrected] - cut_constant.F .
//
// The cut constant is 2i.omega.norm.exp (-2.i.Pi.[l.omega + i.eta]), and one takes its log.
// The returned imaginary part of the log is not necessarily in ]-Pi:Pi].
// Norm is 1.0 for normalized wave functions, C(l,eta)^2 for unnormalized wave functions.
//
//
// Variables:
// ----------
// is_it_normalized : true if one wants normalized functions, i.e. the standard normalization,
//                    false if one wants F -> F/C(l,eta) and H+/H-/G -> H+/H-/G.C(l,eta), to avoid overflows for |eta| >> 1 and |z| small.
// omega : 1 or -1.
// l : orbital angular momentum l.
// eta : Sommerfeld parameter = Coulomb_constant.Z.(2.mu/hbar^2)/(2k).
// Ieta : i.eta .
// l_int,Ieta_int : closest integers to Re[l],Re[i.eta]
// eps : (l.omega - l_int.omega) + (Ieta - Ieta_int)
// two_I_Pi : 2.i.Pi .
// log_norm : log[C(l,eta)^2] if is_it_normalized is false, 0.0 if it is true.
// log_two_I_omega : log[2.i.omega] = log[2] + i.omega.Pi/2 .
// log_cut_constant : returned result.

std::complex<double> log_cut_constant_CFb_calc (const bool is_it_normalized,const int omega,const std::complex<double> &l,const std::complex<double> &eta)
{
  const std::complex<double> Ieta(-imag (eta),real (eta));
  const double l_int = rint (real (l)), Ieta_int = rint (real (Ieta));
  const std::complex<double> eps = omega*(l - l_int) + (Ieta - Ieta_int);
  const std::complex<double> log_norm = (!is_it_normalized) ? (2.0*log_Cl_eta_calc (l,eta)) : (0.0);
  const std::complex<double> two_I_Pi(0,2.0*M_PI),log_two_I_omega(M_LN2,omega*M_PI_2),log_cut_constant = log_two_I_omega - two_I_Pi*eps + log_norm;

  return log_cut_constant;
}



















// Sin (chi) calculation
// ---------------------
// If 2l is integer, 0.0 is returned as chi is zero.
// If not, one calculates sin (chi) with chi = sigma(l,eta) - sigma(-l-1,eta) - (l+1/2).Pi .
// One uses the stable formula sin (chi) = -(2l+1).C(l,eta).C(-l-1,eta).
//
// Variables
// ---------
// l : orbital angular momentum l.
// eta : Sommerfeld parameter.
// sin_chi : sin (chi)

std::complex<double> sin_chi_calc (const std::complex<double> &l,const std::complex<double> &eta)
{
  if (rint (real (2.0*l)) == 2.0*l) return 0.0;

  const std::complex<double> sin_chi = -(2*l+1)*exp (log_Cl_eta_calc (l,eta) + log_Cl_eta_calc (-l-1,eta));
  
  return sin_chi;
}







// exp (i.omega.chi) calculation.
// ------------------------------
// One calculates exp (i.omega.chi), with chi = sigma(l,eta) - sigma(l,-eta) - (l+1/2).Pi .
// If 2l is integer, 1.0 is returned as chi is zero.
// If not, one first calculates sin (chi) with the previous routine.
// If |sin (chi)| > 0.5, chi obtained with the formula sigma(l,eta) - sigma(l,-eta) - (l+1/2).Pi is stable so exp[i.omega.chi] follows directly.
// If not, one uses exp[i.omega.chi] = cos (chi) + i.omega.sin (chi), with cos (chi) = sqrt [1 - sin (chi)*sin (chi)].sign[Re[cos (chi)]],
// with chi given by sigma(l,eta) - sigma(l,-eta) - (l+1/2).Pi .
//
// Variables
// ---------
// omega : 1 or -1
// l : orbital angular momentum l.
// eta : Sommerfeld parameter = Coulomb_constant.Z.(2.mu/hbar^2)/(2k).
// I_omega : i.omega
// sin_chi : sin (chi)
// chi : sigma(l,eta) - sigma(l,-eta) - (l+1/2).Pi . 
// cos_chi : sign[Re[cos (chi)]].sqrt[1 - [sin(chi)]^2]
// exp_I_omega_chi : exp[i.omega.chi], returned result.

std::complex<double> exp_I_omega_chi_calc (const int omega,const std::complex<double> &l,const std::complex<double> &eta)
{
  if (rint (real (2.0*l)) == 2.0*l) return 1.0;

  const std::complex<double> I_omega(0,omega),sin_chi = sin_chi_calc (l,eta);
  const std::complex<double> chi = sigma_l_calc (l,eta) - sigma_l_calc (-l-1,eta) - (l+0.5)*M_PI;

  if (abs (sin_chi) > 0.5)
  {
    const std::complex<double> exp_I_omega_chi = exp (I_omega*chi);

    return exp_I_omega_chi;
  }
  else
  {
    const std::complex<double> cos_chi = SIGN (real (cos (chi)))*sqrt (1.0 - sin_chi*sin_chi),exp_I_omega_chi = cos_chi + I_omega*sin_chi;

    return exp_I_omega_chi;
  }
}
