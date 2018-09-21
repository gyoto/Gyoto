// Gamma_inv denotes the entire inverse of the Gamma function, calculated in complex_functions.H .
// F(z) means 2F1(a,b,c,z) with the a, b, c and z given as inputs in the routine.

// Calculation of H(z,eps) = [Gamma(z+eps)/Gamma(z) - 1]/eps, with e and z std::complex so z,z+eps are not negative integers and 0 <= |eps|oo < 0.1
// -------------------------------------------------------------------------------------------------------------------------------------------
// The function H(z,eps) = [Gamma(z+eps)/Gamma(z) - 1]/e is calculated here with the Lanczos method.
// For the Lanczos method, the gamma parameter, denoted as g, is 4.7421875 and one uses a sum of 15 numbers with the table c[15], 
// so that it is precise up to machine accuracy.
// The H(z,eps) function is used in formulas occuring in1-z and 1/z transformations (see Comp. Phys. Comm. paper).
//
// One must have z and z+eps not negative integers as otherwise it is clearly not defined.
// As this function is meant to be precise for small |eps|oo, one has to have 0 <= |eps|oo < 0.1 .
// Indeed, a direct implementation of H(z,eps) with Gamma_inv or log_Gamma for |eps|oo >= 0.1 is numerically stable.
// The returned function has full numerical accuracy even if |eps|oo is very small.
//
// eps not equal to zero
// ---------------------
// If Re(z) >= 0.5 or Re(z+eps) >= 0.5, one clearly has Re(z) > 0.4 and Re(z+eps) > 0.4, 
// so that the Lanczos summation can be used for both Gamma(z) and Gamma(z+eps).
// One then has:
// log[Gamma(z+eps)/Gamma(z)] = (z-0.5) log1p[eps/(z+g-0.5)] + eps log(z+g-0.5+eps) - eps 
//                            + log1p[-eps \sum_{i=1}^{14} c[i]/((z-1+i)(z-1+i+eps)) / (c[0] + \sum_{i=1}^{14} c[i]/(z-1+i))]
// H(z,eps) = expm1[log[Gamma(z+eps)/Gamma(z)]]/eps .
//
// If Re(z) < 0.5 and Re(z+eps) < 0.5, Euler reflection formula is used for both Gamma(z) and Gamma(z+eps).
// One then has: 
// H(z+eps,-eps) = [cos(pi.eps) + sin(pi.eps)/tan(pi(z-n))].H(1-z,-eps) + (2/eps).sin^2(eps.pi/2) - sin(pi.eps)/(eps.tan(pi.(z-n)))
// H(1-z,-eps) is calculated with the Lanczos summation as Re(1-z) >= 0.5 and Re(1-z-eps) >= 0.5 .
// z-n is used in tan(pi.z) instead of z to avoid inaccuracies due the finite number of digits of pi.
// H(z,eps) = H(z+eps,-eps)/(1 - eps.H(z+eps,-eps)) provides the final result.
//
// eps equal to zero
// -----------------
// It is obtained with the previous case and eps -> 0 :
// If Re(z) >= 0.5, one has:
// H(z,eps) = (z-0.5)/(z+g-0.5) + log(z+g-0.5) - 1 - \sum_{i=1}^{14} c[i]/((z-1+i)^2) / (c[0] + \sum_{i=1}^{14} c[i]/(z-1+i))
//
// If Re(z) < 0.5, one has:
// H(z,0) = H(1-z,0) - pi/tan(pi.(z-n))
//
// Variables
// ---------
// z,eps: input variables of the function H(z,eps)
// g,c[15]: double and table of 15 doubles defining the Lanczos sum so that it provides the Gamma function precise up to machine accuracy.
// eps_pz,z_m_0p5,z_pg_m0p5,eps_pz_pg_m0p5,zm1,zm1_p_eps: z+eps,z-0.5,z+g-0.5,z+eps+g-0.5,z-1,z-1+eps
// x,eps_px: real parts of z and z+eps.
// n,m: closest integer ot the real part of z, same for z+eps.
// sum_num,sum_den: \sum_{i=1}^{14} c[i]/((z-1+i)(z-1+i+eps)) and (c[0] + \sum_{i=1}^{14} c[i]/(z-1+i)). 
//                  They appear respectively as numerator and denominator in formulas.
// Pi_eps,term,T1_eps_z: pi.eps, sin (pi.eps)/tan(pi.(z-n)), [cos(pi.eps) + sin(pi.eps)/tan(pi(z-n))].H(1-z,-eps)
// sin_Pi_2_eps,T2_eps_z,T_eps_z: sin^2(eps.pi/2), (2/eps).sin^2(eps.pi/2) - sin(pi.eps)/(eps.tan(pi.(z-n))), H(z+eps,-eps)

std::complex<double> Gamma_ratio_diff_small_eps (const std::complex<double> &z,const std::complex<double> &eps)
{
  const double g = 4.7421875;
  if (inf_norm (eps) > 0.1) std::cout<<"One must have |eps|oo < 0.1 in Gamma_ratio_diff_small_eps.", abort ();

  const std::complex<double> eps_pz = z + eps,z_m_0p5 = z - 0.5,z_pg_m0p5 = z_m_0p5 + g,eps_pz_pg_m0p5 = z_pg_m0p5 + eps, zm1 = z - 1.0,zm1_p_eps = zm1 + eps;
  const double x = real (z),eps_px = real (eps_pz);
  const int n = static_cast<int> (rint (x)),m = static_cast<int> (rint (eps_px));

  if ((z == n) && (n <= 0)) std::cout<<"z is negative integer in Gamma_ratio_diff_small_eps.", abort ();
  if ((eps_pz == m) && (m <= 0)) std::cout<<"z+eps is negative integer in Gamma_ratio_diff_small_eps.", abort ();

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

  if ((x >= 0.5) || (eps_px >= 0.5))
  {
    std::complex<double> sum_num = 0.0, sum_den = c[0];

    for (int i = 1 ; i < 15 ; i++) 
    {
      const std::complex<double> ci_zm1_pi_inv = c[i]/(zm1 + i);
      sum_num += ci_zm1_pi_inv/(zm1_p_eps + i), sum_den += ci_zm1_pi_inv;
    }

    if (eps != 0.0)
      return expm1 (z_m_0p5*log1p (eps/z_pg_m0p5) + eps*log (eps_pz_pg_m0p5) - eps + log1p (-eps*sum_num/sum_den))/eps;
    else
      return (z_m_0p5/z_pg_m0p5 + log (eps_pz_pg_m0p5) - 1.0 - sum_num/sum_den);
  }
  else
  {
    if (eps != 0.0)
    {
      const std::complex<double> Pi_eps = M_PI*eps,term = sin (Pi_eps)/tan (M_PI*(z-n)),T1_eps_z = (cos (Pi_eps) + term)*Gamma_ratio_diff_small_eps (1.0 - z,-eps);
      const std::complex<double> sin_Pi_2_eps = sin (M_PI_2*eps),T2_eps_z = (2.0*sin_Pi_2_eps*sin_Pi_2_eps - term)/eps, T_eps_z = T1_eps_z + T2_eps_z;
      
      return (T_eps_z/(1.0 - eps*T_eps_z));
    }
    else
      return (Gamma_ratio_diff_small_eps (1.0 - z,0.0) - M_PI/tan (M_PI*(z-n)));
  }
}




// Calculation of G(z,eps) = [Gamma_inv(z) - Gamma_inv(z+eps)]/eps, with e and z std::complex
// -------------------------------------------------------------------------------------
// The G(z,eps) function is used in formulas occuring in 1-z and 1/z transformations (see Comp. Phys. Comm. paper).
// Several case have to be considered for its evaluation. eps is considered equal to zero if z+eps and z are equal numerically.
//
// |eps|oo > 0.1
// -------------
// A direct evaluation with the values Gamma_inv(z) and Gamma_inv(z+eps) is stable and returned.
//
// |eps|oo <= 0.1 with z+eps and z numerically different
// -----------------------------------------------------
// If z is a negative integer, z+eps is not, so that G(z,eps) = -Gamma_inv(z+eps)/eps, for which a direct evaluation is precise and returned.
// If z+eps is a negative integer, z is not, so that G(z,eps) = Gamma_inv(z)/eps, for which a direct evaluation is precise and returned.
// If both of them are not negative integers, one looks for the one of z and z+eps which is the closest to a negative integer.
// If it is z, one returns H(z,eps).Gamma_inv(z+eps). If it is z+eps, one returns H(z+eps,-eps).Gamma_inv(z).
// Both values are equal, so that one chooses the one which makes the Gamma ratio Gamma(z+eps)/Gamma(z) in H(z,eps) the smallest in modulus.
//
// z+eps and z numerically equal
// -----------------------------
// If z is negative integer, G(z,0) = (-1)^(n+1) n!, where z = -n, n integer, which is returned.
// If z is not negative integer, one returns H(z,eps).Gamma_inv(z+eps) .
//
// Variables
// ---------
// z,eps: input variables of the function G(z,eps)
// eps_pz,x,eps_px: z+eps,real parts of z and z+eps.
// n,m: closest integer ot the real part of z, same for z+eps.
// fact: (-1)^(n+1) n!, returned when z = -n, n integer and z and z+eps identical numerically (eps ~ 0).
// is_z_negative_integer,is_eps_pz_negative_integer: true if z is a negative integer, false if not, same for z+eps.
// z_neg_int_distance, eps_pz_neg_int_distance: |z + |n||oo, |z + eps + |m||oo. 
//                                              If |z + |n||oo < |z + eps + |m||oo, z is closer to the set of negative integers than z+eps.
//                                              Gamma_inv(z+eps) is then of moderate modulus if Gamma_inv(z) is very small. 
//                                              If z ~ n, H(z,eps) ~ -1/eps, that so returning G(z,eps) = H(z,eps).Gamma_inv(z+eps) here is preferred.
//                                              Same for |z + |n||oo > |z + eps + |m||oo with z <-> z+eps.

std::complex<double> Gamma_inv_diff_eps (const std::complex<double> &z,const std::complex<double> &eps)
{
  const std::complex<double> eps_pz = z + eps;
  const double x = real (z),eps_px = real (eps_pz);
  const int n = static_cast<int> (rint (x)), m = static_cast<int> (rint (eps_px));

  const bool is_z_negative_integer = (z == n) && (n <= 0),is_eps_pz_negative_integer = (eps_pz == m) && (m <= 0);

  if (inf_norm (eps) > 0.1)
    return (Gamma_inv (z) - Gamma_inv (eps_pz))/eps;
  else if (eps_pz != z)
  { 
    if (is_z_negative_integer)
      return (-Gamma_inv (eps_pz)/eps);
    else if (is_eps_pz_negative_integer)
      return (Gamma_inv (z)/eps);
    else
    {
      const double z_neg_int_distance = inf_norm (z + abs (n)),eps_pz_neg_int_distance = inf_norm (eps_pz + abs (m));

      if (z_neg_int_distance < eps_pz_neg_int_distance)
	return Gamma_ratio_diff_small_eps (z,eps)*Gamma_inv (eps_pz);
      else
	return Gamma_ratio_diff_small_eps (eps_pz,-eps)*Gamma_inv (z);
    }
  }
  else if (is_z_negative_integer && is_eps_pz_negative_integer)
  {
    double fact = -1.0;
    for (int k = -1 ; k >= n ; k--) fact *= k;
    
    return fact;
  }
  else
    return Gamma_ratio_diff_small_eps (z,eps)*Gamma_inv (eps_pz);
}





// Calculation of Gamma_inv(1-m-eps)/eps of the A(z) polynomial in 1-z, and 1/z transformations
// --------------------------------------------------------------------------------------------
// This value occurs in A(z) in 1-z and 1/z transformations (see Comp. Phys. Comm. paper) for m > 0.
// Both cases of 1-m-eps numerically negative integer or not have to be considered
// 
// 1-eps-m and 1-m numerically different
// -------------------------------------
// One returns Gamma_inv(1-m-eps)/eps directly as its value is accurate.
// To calculate Gamma_inv(1-m-eps), one uses the value Gamma_inv(1-eps), needed in considered transformations,
// and one uses the equality Gamma_inv(1-m-eps) = Gamma_inv(1-eps) \prod_{i=1}^{m} (1-eps-i) for m > 0.
// It is trivially demonstrated from the equality Gamma(x+1) = x.Gamma(x). One Gamma function evaluation is removed this way from the calculation.
// 
// 1-eps-m and 1-m numerically equal
// ---------------------------------
// This implies that 1-m-eps is negative integer numerically.
// Here, eps ~ 0, so that one returns the limit of Gamma_inv(1-m-eps)/eps for eps -> 0, which is (-1)^m (m-1)!
//
// Variables
// ---------
// m,eps: variable inputs of the function (m,eps) -> Gamma_inv(1-m-eps)/eps
// Gamma_inv_one_meps: Gamma_inv(1-eps), previously calculated and here recycled to quickly calculate Gamma_inv(1-m-eps).
// one_meps: 1-eps


std::complex<double> A_sum_init (const int m,const std::complex<double> &eps,const std::complex<double> &Gamma_inv_one_meps)
{
  const std::complex<double> one_meps = 1.0 - eps;

  if (one_meps-m != 1-m)
  {
    std::complex<double> Gamma_inv_one_meps_mm = Gamma_inv_one_meps;
    for (int i = 1 ; i <= m ; i++) Gamma_inv_one_meps_mm *= one_meps - i;
    return Gamma_inv_one_meps_mm/eps;
  }
  else
  {
    double fact = 1.0;
    for (int n = 2 ; n < m ; n++) fact *= n;

    return (m%2 == 0) ? (fact) : (-fact);
  }
}




// Calculation of the log of Gamma_inv(1-m-eps)/eps
// ------------------------------------------------
// See previous function. It is used in case Gamma_inv(1-m-eps)/eps overflows.
//
// Variables
// ---------
// m,eps: variable inputs of the function (m,eps) -> log[Gamma_inv(1-m-eps)/eps]
// one_meps_mm: 1-eps-m
// i_Pi: i.Pi
// log_fact: logarithm of (-1)^m (m-1)!, here defined as log((m-1)!) + i.Pi if m is odd.

std::complex<double> log_A_sum_init (const int m,const std::complex<double> &eps)
{
  const std::complex<double> one_meps_mm = 1.0 - m - eps;

  if (one_meps_mm != 1-m)
    return (-log_Gamma (one_meps_mm) - log (eps));
  else
  {
    const std::complex<double> i_Pi(0,M_PI);

    double log_fact = 0.0;
    for (int n = 2 ; n < m ; n++) log_fact += log (static_cast<double>(n));

    return (m%2 == 0) ? (log_fact) : (log_fact + i_Pi);
  }
}







// Calculation of the first term of the B(z) power series in the 1-z transformation, divided by (1-z)^m
// ----------------------------------------------------------------------------------------------------
// In the 1-z transformation, the power series B(z) = \sum_{n=0}^{+oo} \beta_n (1-z)^n occurs (see Comp. Phys. Comm. paper).
// The first term \beta_0, divided by (1-z)^m, is calculated here. m is the closest integer to Re(c-a-b) >= 0 and eps = c-a-b-m.
//
// One has to consider |eps|oo > 0.1 and |eps|oo <= 0.1, where 1-m-eps and 1-m can be different or equal numerically, leading to some changes in this last case.
//
// |eps|oo > 0.1
// -------------
// One has \beta_0/(1-z)^m = [(a)_m (b)_m Gamma_inv(1-eps) Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) Gamma_inv(m+1)
//                         - (1-z)^eps Gamma_inv(a) Gamma_inv(b) Gamma_inv(1+m+eps)].[Gamma(c)/eps], stable in this regime for a direct evaluation.
//
// The values of Gamma(c), Gamma_inv(a+m+eps) and Gamma_inv(b+m+eps) were already calculated and recycled here.
// Gamma_inv(m+1) is calculated as 1/(m!).
//
// Gamma_inv(1+m+eps) is calculated from Gamma_inv(1-eps), using the equalities:
// Gamma_inv(1-m-eps) = Gamma_inv(1-eps) \prod_{i=1}^{m} (1-eps-i), where the product is 1 by definition if m = 0,
// Gamma_inv(1+m+eps) = (-1)^m sin (pi.eps)/[pi.(eps+m).Gamma_inv(1-m-eps)] from Euler reflection formula, Gamma(x+1) = x.Gamma(x) equality, and m+eps no zero.
// This scheme is much faster than to recalculate Gamma_inv(1+m+eps) directly.
// 
// |eps|oo <= 0.1
// --------------
// The \beta_0/(1-z)^m expression is rewritten so that it contains no instabilities:
// \beta_0/(1-z)^m = Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) [(G(1,-eps) Gamma_inv(m+1) + G(m+1,eps))
//                 - Gamma_inv(1+m+eps) (G(a+m,eps) Gamma_inv(b+m+eps) + G(b+m,eps) Gamma_inv(a+m)) 
//                 - E(log(1-z),eps) Gamma_inv(a+m) Gamma_inv(b+m) Gamma_inv(1+m+eps)] (a)_m (b)_m Gamma(c)
//
// E(log(1-z),eps) is [(1-z)^eps - 1]/eps if 1-m-eps and 1-m are different numerically, and log(1-z) otherwise (eps ~ 0).
// If 1-m-eps and 1-m are equal numerically, Gamma_inv(1+m+eps) is numerically equal to Gamma_inv(1+m), already calculated as 1/(m!).
// See |eps|oo > 0.1 case for data recycling of other values or for 1-m-eps and 1-m different numerically.
// In case it overflows, |eps|oo > 0.1 case is used as a last resort.
//
// Variables
// ---------
// a,b,c,one_minus_z: a,b,c and 1-z parameters and arguments of the 2F1(a,b,c,z) function.
// m,eps: closest integer to c-a-b, with Re(c-a-b) >= 0 and eps = c-a-b-m
// Gamma_c,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm,Gamma_inv_eps_pb_pm: recycled values of Gamma(c), Gamma_inv(1-eps), Gamma_inv(a+m+eps) and Gamma_inv(b+m+eps).
// inf_norm_eps,phase,a_pm,b_pm,one_meps,Pi_eps,Pi_eps_pm: |eps|oo,(-1)^m,a+m,b+m,1-eps,pi.eps,pi.(eps+m)
// Gamma_inv_one_meps_mm,Gamma_inv_eps_pm_p1: Gamma_inv(1-m-eps) and Gamma_inv(1+m+eps) calculated with the recycling scheme.
// prod1: (a)_m (b)_m Gamma_inv(1-eps) Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) Gamma_inv(m+1) in |eps|oo > 0.1 case.
// prod2: (1-z)^eps Gamma_inv(a) Gamma_inv(b) Gamma_inv(1+m+eps) in |eps|oo > 0.1 case.
// Gamma_inv_mp1,prod_ab: Gamma_inv(m+1) calculated as 1/(m!) and (a)_m (b)_m in |eps|oo <= 0.1 case.
// is_eps_non_zero: true if 1-m-eps and 1-m are different numerically, false if not.
// Gamma_inv_a_pm,Gamma_inv_b_pm,z_term: Gamma_inv(a+m),Gamma_inv(b+m),E(eps,log(1-z))
// prod1: Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) [(G(1,-eps) Gamma_inv(m+1) + G(m+1,eps)) in |eps|oo <= 0.1 case.
// prod2: Gamma_inv(1+m+eps) (G(a+m,eps) Gamma_inv(b+m+eps) + G(b+m,eps) Gamma_inv(a+m))
// prod3: E(eps,log(1-z)) Gamma_inv(a+m) Gamma_inv(b+m) Gamma_inv(1+m+eps) 
// res: returned \beta_0/(1-z)^m value in all cases.

std::complex<double> B_sum_init_PS_one (const std::complex<double> &a,const std::complex<double> &b,const std::complex<double> &c,
				   const std::complex<double> &Gamma_c,const std::complex<double> &Gamma_inv_one_meps,
				   const std::complex<double> &Gamma_inv_eps_pa_pm,const std::complex<double> &Gamma_inv_eps_pb_pm,
				   const std::complex<double> &one_minus_z,const int m,const std::complex<double> &eps)
{
  const double inf_norm_eps = inf_norm (eps),phase = (m%2 == 0) ? (1) : (-1);
  const std::complex<double> a_pm = a + m,b_pm = b + m,one_meps = 1.0 - eps,Pi_eps = M_PI*eps,Pi_eps_pm = M_PI*(eps + m);

  std::complex<double> Gamma_inv_one_meps_mm = Gamma_inv_one_meps;
  for (int i = 1 ; i <= m ; i++) Gamma_inv_one_meps_mm *= one_meps - i;
  if (inf_norm_eps > 0.1)
  {
    const std::complex<double> Gamma_inv_eps_pm_p1 = phase*sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm);
    std::complex<double> prod1 = Gamma_inv_one_meps*Gamma_inv_eps_pa_pm*Gamma_inv_eps_pb_pm;
    for (int n = 0 ; n < m ; n++) prod1 *= (a + n)*(b + n)/(n + 1.0);
    const std::complex<double> prod2 = Gamma_inv (a)*Gamma_inv (b)*Gamma_inv_eps_pm_p1*pow (one_minus_z,eps),res = Gamma_c*(prod1 - prod2)/eps;
    return res;
  }
  else
  {
    double Gamma_inv_mp1 = 1.0;
    std::complex<double> prod_ab = 1.0;
    for (int n = 0 ; n < m ; n++) Gamma_inv_mp1 /= n + 1.0, prod_ab *= (a + n)*(b + n);

    const bool is_eps_non_zero = (one_meps-m != 1-m);
    const std::complex<double> Gamma_inv_eps_pm_p1 = (is_eps_non_zero) ? (phase*sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm)) : (Gamma_inv_mp1);
    const std::complex<double> Gamma_inv_a_pm = Gamma_inv (a_pm),Gamma_inv_b_pm = Gamma_inv (b_pm);
    const std::complex<double> z_term = (is_eps_non_zero) ? (expm1 (eps*log (one_minus_z))/eps) : (log (one_minus_z));

    const std::complex<double> prod1 = Gamma_inv_eps_pa_pm*Gamma_inv_eps_pb_pm*(Gamma_inv_mp1*Gamma_inv_diff_eps (1.0,-eps) + Gamma_inv_diff_eps (m+1,eps));
    const std::complex<double> prod2 = Gamma_inv_eps_pm_p1*(Gamma_inv_eps_pb_pm*Gamma_inv_diff_eps (a_pm,eps) + Gamma_inv_a_pm*Gamma_inv_diff_eps (b_pm,eps));
    const std::complex<double> prod3 = Gamma_inv_a_pm*Gamma_inv_b_pm*Gamma_inv_eps_pm_p1*z_term;

    const std::complex<double> res = Gamma_c*prod_ab*(prod1 - prod2 - prod3); 
    if (isfinite (res))
      return res;
    else
    {
      const std::complex<double> Gamma_inv_eps_pm_p1 = phase*sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm);
      std::complex<double> prod1 = Gamma_inv_one_meps*Gamma_inv_eps_pa_pm*Gamma_inv_eps_pb_pm;
      for (int n = 0 ; n < m ; n++) prod1 *= (a + n)*(b + n)/(n + 1.0);
      const std::complex<double> prod2 = Gamma_inv (a)*Gamma_inv (b)*Gamma_inv_eps_pm_p1*pow (one_minus_z,eps),res_default = Gamma_c*(prod1 - prod2)/eps;
      return res_default;
    }
  }
}




// Calculation of the first term of the B(z) power series in the 1/z transformation, divided by z^{-m}
// ---------------------------------------------------------------------------------------------------
// In the 1/z transformation, the power series B(z) = \sum_{n=0}^{+oo} \beta_n z^{-n} occurs (see Comp. Phys. Comm. paper).
// The first term \beta_0, divided by z^{-m}, is calculated here. m is the closest integer to Re(b-a) >= 0 and eps = b-a-m.
//
// One has to consider |eps|oo > 0.1 and |eps|oo <= 0.1, where 1-m-eps and 1-m can be different or equal numerically, leading to some changes in this last case.
//
// |eps|oo > 0.1
// -------------
// One has \beta_0/z^{-m} = [(a)_m (1-c+a)_m Gamma_inv(1-eps) Gamma_inv(a+m+eps) Gamma_inv(c-a) Gamma_inv(m+1)
//          - (-z)^{-eps} (1-c+a+eps)_m Gamma_inv(a) Gamma_inv(c-a-eps) Gamma_inv(1+m+eps)].[Gamma(c)/eps], stable in this regime for a direct evaluation.
//
// The values of Gamma(c), Gamma_inv(c-a) and Gamma_inv(a+m+eps) were already calculated and recycled here.
// Gamma_inv(m+1) is calculated as 1/(m!). Gamma_inv(1+m+eps) is calculated from Gamma_inv(1-eps) as in the 1-z transformation routine.
// 
// |eps|oo <= 0.1
// --------------
// The \beta_0/z^{-m} expression is rewritten so that it contains no instabilities:
// \beta_0/z^{-m} = [((1-c+a+eps)_m G(1,-eps) - P(m,eps,1-c+a) Gamma_inv(1-eps)) Gamma_inv(c-a) Gamma_inv(a+m+eps) Gamma_inv(m+1)
//                + (1-c+a+eps)_m [G(m+1,eps) Gamma_inv(c-a) Gamma_inv(a+m+eps) - G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)]
//                - (G(c-a,-eps) - E(log(-z),-eps)) Gamma_inv(m+1+eps) Gamma_inv(a+m)]] (a)_m Gamma(c)
//
// Definitions and method are the same as in the 1-z transformation routine, except for P(m,eps,1-c+a).
// P(m,eps,s) = [(s+eps)_m - (s)_m]/eps for eps non zero and has a limit for eps -> 0.
// Let n0 be the closest integer to -Re(s) for s std::complex. A stable formula available for eps -> 0 for P(m,eps,s) is:
// P(m,eps,s) = (s)_m E(\sum_{n=0}^{m-1} L(1/(s+n),eps),eps) if n0 is not in [0:m-1],
// P(m,eps,s) = \prod_{n=0, n not equal to n0}^{m-1} (s+eps+n) + (s)_m E(\sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps),eps) if n0 is in [0:m-1].
// L(s,eps) is log1p(s eps)/eps if eps is not zero, and L(s,0) = s.
// This expression is used in the code.
//
// Variables
// ---------
// a,b,c,z: a,b,c and z parameters and arguments of the 2F1(a,b,c,z) function.
// m,eps: closest integer to b-a, with Re(b-a) >= 0 and eps = b-a-m.
// Gamma_c,Gamma_inv_cma,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm: recycled values of Gamma(c), Gamma_inv(c-a), Gamma_inv(1-eps) and Gamma_inv(a+m+eps).
// inf_norm_eps,phase,cma,a_mc_p1,a_mc_p1_pm,cma_eps,eps_pa_mc_p1,a_pm: |eps|oo,(-1)^m,c-a,1-c+a+m,c-a-eps,1-c+a+eps,a+m
// Gamma_inv_cma_meps,one_meps,Pi_eps,Pi_eps_pm: Gamma_inv(c-a-eps),1-eps,pi.eps,pi.(eps+m)
// Gamma_inv_one_meps_mm,Gamma_inv_eps_pm_p1: Gamma_inv(1-m-eps) and Gamma_inv(1+m+eps) calculated with the recycling scheme.
// prod1: (a)_m (1-c+a)_m Gamma_inv(1-eps) Gamma_inv(a+m+eps) Gamma_inv(c-a) Gamma_inv(m+1) in |eps|oo > 0.1 case.
// prod2: (-z)^{-eps} (1-c+a+eps)_m Gamma_inv(a) Gamma_inv(c-a-eps) Gamma_inv(1+m+eps) in |eps|oo > 0.1 case.
// n0: closest integer to -Re(1-c+a)
// is_n0_here: true is n0 belongs to [0:m-1], false if not.
// is_eps_non_zero: true if 1-m-eps and 1-m are different numerically, false if not.
// Gamma_inv_mp1,prod_a,prod_a_mc_p1: Gamma_inv(m+1) calculated as 1/(m!), (a)_m and (1-c+a)_m in |eps|oo <= 0.1 case.
// prod_eps_pa_mc_p1_n0: \prod_{n=0, n not equal to n0}^{m-1} (1-c+a+eps+n) if n0 belongs to [0:m-1], 0.0 if not, in |eps|oo <= 0.1 case.
// prod_eps_pa_mc_p1: (1-c+a+eps)_m in |eps|oo <= 0.1 case.
// sum: \sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps) if 1-m-eps and 1-m are different numerically, \sum_{n=0, n not equal to n0}^{m-1} 1/(s+n) if not.
// a_pn,a_mc_p1_pn,eps_pa_mc_p1_pn: a+n,1-c+a+n,1-c+a+eps+n values used in (a)_m, (1-c+a)_m and (1-c+a+eps)_m evaluations.
// sum_term,prod_diff_eps,z_term: E(\sum_{n=0, n not equal to n0}^{m-1} L(1/(s+n),eps),eps), P(m,eps,1-c+a), -E(-eps,log(-z))
// Gamma_inv_a_pm,Gamma_prod1: Gamma_inv(a+m), Gamma_inv(c-a).Gamma_inv(a+m+eps)
// prod1: ((1-c+a+eps)_m G(1,-eps) - P(m,eps,1-c+a) Gamma_inv(1-eps)) Gamma_inv(c-a) Gamma_inv(a+m+eps) Gamma_inv(m+1)
// prod_2a: Gamma_inv(c-a).Gamma_inv(a+m+eps).G(m+1,eps)
// prod_2b: G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)
// prod_2c: (G(c-a,-eps) - E(log(-z),-eps)) Gamma_inv(m+1+eps) Gamma_inv(a+m)
// prod2: (1-c+a+eps)_m [G(m+1,eps) Gamma_inv(c-a) Gamma_inv(a+m+eps) - G(a+m,eps) Gamma_inv(c-a) Gamma_inv(m+1+eps)] 
//        - (G(c-a,-eps) - E(log(-z),-eps)) Gamma_inv(m+1+eps) Gamma_inv(a+m)]]
// res: returned \beta_0/z^{-m} value in all cases.

std::complex<double> B_sum_init_PS_infinity (const std::complex<double> &a,const std::complex<double> &c,
					const std::complex<double> &Gamma_c,const std::complex<double> &Gamma_inv_cma,
					const std::complex<double> &Gamma_inv_one_meps,const std::complex<double> &Gamma_inv_eps_pa_pm,
					const std::complex<double> &z,const int m,const std::complex<double> &eps)
{
  const double inf_norm_eps = inf_norm (eps),phase = (m%2 == 0) ? (1) : (-1);
  const std::complex<double> cma = c - a,a_mc_p1 = 1.0 - c + a,a_mc_p1_pm = a_mc_p1 + m,cma_meps = cma - eps,eps_pa_mc_p1 = eps + a_mc_p1,a_pm = a + m;
  const std::complex<double> Gamma_inv_cma_meps = Gamma_inv (cma_meps),one_meps = 1.0 - eps,Pi_eps = M_PI*eps,Pi_eps_pm = M_PI*(eps + m);

  std::complex<double> Gamma_inv_one_meps_mm = Gamma_inv_one_meps;
  for (int i = 1 ; i <= m ; i++) Gamma_inv_one_meps_mm *= one_meps - i;
    
  if (inf_norm_eps > 0.1)
  { 
    const std::complex<double> Gamma_inv_eps_pm_p1 = phase*sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm);
    std::complex<double> prod1 = Gamma_inv_cma*Gamma_inv_eps_pa_pm*Gamma_inv_one_meps,prod2 = Gamma_inv (a)*Gamma_inv_cma_meps*Gamma_inv_eps_pm_p1*pow (-z,-eps);

    for (int n = 0 ; n < m ; n++)
    {
      const double n_p1 = n + 1.0;
      const std::complex<double> a_pn = a + n,a_mc_p1_pn = a_mc_p1 + n,eps_pa_mc_p1_pn = eps_pa_mc_p1 + n;
      prod1 *= a_pn*a_mc_p1_pn/n_p1, prod2 *= eps_pa_mc_p1_pn;
    }

    const std::complex<double> res = Gamma_c*(prod1 - prod2)/eps;
    return res;
  }
  else
  {
    const int n0 = -static_cast<int> (rint (real (a_mc_p1)));
    const bool is_eps_non_zero = (one_meps-m != 1-m),is_n0_here = (n0 >= 0) && (n0 < m);

    double Gamma_inv_mp1 = 1.0;
    std::complex<double> prod_a = 1.0,prod_a_mc_p1 = 1.0,prod_eps_pa_mc_p1_n0 = (is_n0_here) ? (1.0) : (0.0),prod_eps_pa_mc_p1 = 1.0,sum = 0.0;
    for (int n = 0 ; n < m ; n++)
    {
      const std::complex<double> a_pn = a + n,a_mc_p1_pn = a_mc_p1 + n,eps_pa_mc_p1_pn = eps_pa_mc_p1 + n;
      prod_a *= a_pn, prod_a_mc_p1 *= a_mc_p1_pn, prod_eps_pa_mc_p1 *= eps_pa_mc_p1_pn, Gamma_inv_mp1 /= n+1.0;

      if (n != n0) 
      {
	if (is_n0_here) prod_eps_pa_mc_p1_n0 *= eps_pa_mc_p1_pn;
	sum += (is_eps_non_zero) ? (log1p (eps/a_mc_p1_pn)) : (1.0/a_mc_p1_pn);
      }
    }

    const std::complex<double> Gamma_inv_eps_pm_p1 = (is_eps_non_zero) ? (phase*sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm)) : (Gamma_inv_mp1);
    const std::complex<double> sum_term = (is_eps_non_zero) ? (expm1 (sum)/eps) : (sum),prod_diff_eps = prod_eps_pa_mc_p1_n0 + prod_a_mc_p1*sum_term;

    const std::complex<double> z_term = (is_eps_non_zero) ? (expm1 (-eps*log (-z))/eps) : (-log (-z));
    const std::complex<double> Gamma_inv_a_pm = Gamma_inv (a_pm),Gamma_prod1 = Gamma_inv_cma*Gamma_inv_eps_pa_pm;
    const std::complex<double> prod1 = Gamma_prod1*Gamma_inv_mp1*(Gamma_inv_diff_eps (1.0,-eps)*prod_eps_pa_mc_p1 - Gamma_inv_one_meps*prod_diff_eps);
    const std::complex<double> prod_2a = Gamma_prod1*Gamma_inv_diff_eps (m+1,eps),prod_2b = Gamma_inv_cma*Gamma_inv_eps_pm_p1*Gamma_inv_diff_eps (a_pm,eps);
    const std::complex<double> prod_2c = Gamma_inv_eps_pm_p1*Gamma_inv_a_pm*(Gamma_inv_diff_eps (cma,-eps) + Gamma_inv_cma_meps*z_term);
    const std::complex<double> prod2 = prod_eps_pa_mc_p1*(prod_2a - prod_2b - prod_2c),res = Gamma_c*prod_a*(prod1 + prod2);

    if (isfinite (res))
      return res;
    else
    { 
      const std::complex<double> Gamma_inv_eps_pm_p1 = phase*sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm);
      std::complex<double> prod1 = Gamma_inv_cma*Gamma_inv_eps_pa_pm*Gamma_inv_one_meps,prod2 = Gamma_inv (a)*Gamma_inv_cma_meps*Gamma_inv_eps_pm_p1*pow (-z,-eps);
      
      for (int n = 0 ; n < m ; n++)
      {
	const double n_p1 = n + 1.0;
	const std::complex<double> a_pn = a + n,a_mc_p1_pn = a_mc_p1 + n,eps_pa_mc_p1_pn = eps_pa_mc_p1 + n;
	prod1 *= a_pn*a_mc_p1_pn/n_p1, prod2 *= eps_pa_mc_p1_pn;
      }

      const std::complex<double> res_default = Gamma_c*(prod1 - prod2)/eps;
      return res_default;
    }
  }
}



// Calculation of the derivative of the polynomial P(X) testing power series convergence
// -------------------------------------------------------------------------------------
// P(X) = |z(a+X)(b+X)|^2 - |(c+X)(X+1)|^2 = \sum_{i=0}^{4} c[i] X^{i}, for |z| < 1.
// It is positive when the power series term modulus increases and negative when it decreases, 
// so that its derivative provides information on its convergence (see Comp. Phys. Comm. paper).
// Its derivative components cv_poly_der_tab[i] = (i+1) c[i+1] for i in [0:3] so that P'(X) = \sum_{i=0}^{3} cv_poly_der_tab[i] X^{i} are calculated.
//
// Variables:
// ----------
// a,b,c,z: a,b,c and z parameters and arguments of the 2F1(a,b,c,z) function.
// cv_poly_der_tab[3]: table of four doubles containing the P'(X) components.
// mod_a2,mod_b2,mod_c2,mod_z2,R_a,Re_b,Re_c: |a|^2, |b|^2, |c|^2, |z|^2, Re(a), Re(b), Re(c), with which P(X) can be expressed.

void cv_poly_der_tab_calc (const std::complex<double> &a,const std::complex<double> &b,const std::complex<double> &c,const std::complex<double> &z,double cv_poly_der_tab[])
{
  const double mod_a2 = norm (a),mod_b2 = norm (b),mod_c2 = norm (c),mod_z2 = norm (z);
  const double Re_a = real (a),Re_b = real (b),Re_c = real (c);

  cv_poly_der_tab[0] = 2.0*((Re_a*mod_b2 + Re_b*mod_a2)*mod_z2 - Re_c - mod_c2);
  cv_poly_der_tab[1] = 2.0*((mod_a2 + mod_b2 + 4.0*Re_a*Re_b)*mod_z2 - 1.0 - 4.0*Re_c - mod_c2);
  cv_poly_der_tab[2] = 6.0*((Re_a + Re_b)*mod_z2 - Re_c - 1.0);
  cv_poly_der_tab[3] = 4.0*(mod_z2 - 1.0);
}




// Calculation of the derivative of the polynomial P(X) testing power series convergence at one x value
// ----------------------------------------------------------------------------------------------------
// P'(x) is calculated for a real x. See P'(X) components calculation routine for definitions.

double cv_poly_der_calc (const double cv_poly_der_tab[],const double x)
{
  const double Px = cv_poly_der_tab[0] + x*(cv_poly_der_tab[1] + x*(cv_poly_der_tab[2] + x*cv_poly_der_tab[3]));

  return Px;
}



// Calculation of an integer after which false convergence cannot occur
// --------------------------------------------------------------------
// See cv_poly_der_tab_calc routine for definitions.
// If P'(x) < 0 and P''(x) < 0 for x > xc, it will be so for all x > xc as P(x) -> -oo for x -> +oo and P(x) can have at most one maximum for x > xc. 
// It means that the 2F1 power series term modulus will increase or decrease to 0 for n > nc, with nc the smallest positive integer larger than xc.
//
// If P'(X) = C0 + C1.X + C2.X^2 + C3.X^3, the discriminant of P''(X) is Delta = C2^2 - 3 C1 C3.
//
// If Delta > 0, P''(X) has two different real roots and its largest root is -(C2 + sqrt(Delta))/(3 C3), because C3 = 4(|z|^2 - 1) < 0.
// One can take xc = -(C2 + sqrt(Delta))/(3 C3) and one returns its associated nc integer.
//
// If Delta <= 0, P''(X) has at most one real root, so that P'(X) has only one root and then P(X) only one maximum.
// In this case, one can choose xc = nc = 0, which is returned.
//
// Variables
// ---------
// cv_poly_der_tab: table of four doubles containing the P'(X) coefficients
// C1,C2,three_C3: cv_poly_der_tab[1], cv_poly_der_tab[2] and 3.0*cv_poly_der_tab[3], so that P''(X) = C1 + 2.C2.x + three_C3.x^2
// Delta: discriminant of P''(X), equal to C2^2 - 3 C1 C3.
// largest_root: if Delta > 0, P''(X) largest real root equal to -(C2 + sqrt(Delta))/(3 C3).

int min_n_calc (const double cv_poly_der_tab[])
{
  const double C1 = cv_poly_der_tab[1],C2 = cv_poly_der_tab[2],three_C3 = 3.0*cv_poly_der_tab[3];
  const double Delta = C2*C2 - three_C3*C1;

  if (Delta <= 0.0) 
    return 0;
  else
  {
    const double largest_root = -(C2 + sqrt (Delta))/three_C3;
    return fmax (static_cast<int> (ceil (largest_root)),0);
  }
}





// Calculation of the 2F1 power series converging for |z| < 1
// ----------------------------------------------------------
// One has 2F1(a,b,c,z) = \sum_{n = 0}^{+oo} (a)_n (b)_n / ((c)_n n!) z^n,
// so that 2F1(a,b,c,z) = \sum_{n = 0}^{+oo} t[n] z^n, with t[0] = 1 and t[n+1] = (a+n)(b+n)/((c+n)(n+1)) t[n] for n >= 0.
// If a or b are negative integers, F(z) is a polynomial of degree -a or -b, evaluated directly.
// If not, one uses the test of convergence |t[n] z^n|oo < 1E-15 to truncate the series after it was checked that false convergence cannot occur.
//
// Variables:
// ----------
// a,b,c,z: a,b,c and z parameters and arguments of the 2F1(a,b,c,z) function. One must have here |z| < 1.
// term,sum: term of the 2F1 power series equal to t[n] z^n, truncated sum at given n of the 2F1 power series.
// na,nb: absolute values of the closest integers to Re(a) and Re(b). a = -na or b = -nb means one is in the polynomial case.
// cv_poly_der_tab: coefficients of the derivative of the polynomial P(X) = |z(a+X)(b+X)|^2 - |(c+X)(X+1)|^2
// min_n: smallest integer after which false convergence cannot occur. It is calculated in min_n_calc.
// possible_false_cv: always true if n < min_n. If n >= min_n, it is true if P'(n) > 0. 
//                    If n >= min_n and P'(n) < 0, it becomes false and remains as such for the rest of the calculation. 
//                    One can then check if |t[n] z^n|oo < 1E-15 to truncate the series.

std::complex<double> hyp_PS_zero (const std::complex<double> &a,const std::complex<double> &b,const std::complex<double> &c,const std::complex<double> &z)
{  
  std::complex<double> term = 1.0,sum = 1.0;
  
  const int na = abs (static_cast<int> (rint (real (a)))),nb = abs (static_cast<int> (rint (real (b))));

  if (a == -na)
    for (int n = 0 ; n < na ; n++)
    {
      term *= z*(a+n)*(b+n)/((n+1.0)*(c+n));
      sum += term;
    }
  else if (b == -nb)
    for (int n = 0 ; n < nb ; n++)
    {
      term *= z*(a+n)*(b+n)/((n+1.0)*(c+n));
      sum += term;
    }
  else
  {
    double cv_poly_der_tab[4];
    cv_poly_der_tab_calc (a,b,c,z,cv_poly_der_tab);

    const int min_n = min_n_calc (cv_poly_der_tab);
    bool possible_false_cv = true;
    int n = 0;

    while (possible_false_cv || (inf_norm (term) > 1E-15))
    {
      term *= z*(a+n)*(b+n)/((n+1.0)*(c+n));
      sum += term;
      if (possible_false_cv && (n >= min_n)) possible_false_cv = (cv_poly_der_calc (cv_poly_der_tab,n) > 0);
      n++;
    }
  }

  return sum;
}




// Calculation of the 2F1 power series converging with the 1-z transformation
// --------------------------------------------------------------------------
// The formula for F(z) in the 1-z transformation holds:
// F(z) = (-1)^m (pi.eps)/sin (pi.eps) [A(z) + B(z)] for eps not equal to zero, F(z) = (-1)^m [A(z) + B(z)] for eps = 0
// where m = |Re(c-a-b)], eps = c-a-b-m, A(z) = \sum_{n=0}^{m-1} alpha[n] (1-z)^n, B(z) = \sum_{n=0}^{+oo} beta[n] (1-z)^n, and:
//
// alpha[0] = [Gamma_inv(1-m-eps)/eps] Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) Gamma(c)
// [Gamma_inv(1-m-eps)/eps] is calculated in A_sum_init. 
// alpha[0] is recalculated with log[Gamma] if the previous expression overflows, and its imaginary part removed if a, b and c are real.
// alpha[n+1] = (a+n)(b+n)/[(n+1)(1-m-eps+n)] alpha[n], n in [0:m-2].
//
// beta[0] is defined in B_sum_init_PS_one function comments.
// gamma[0] = Gamma(c) (a)_m (b)_m (1-z)^m Gamma_inv(a+m+eps) Gamma_inv(b+m+eps) Gamma_inv(m+1) Gamma_inv(1-eps)
//
// beta[n+1] = (a+m+n+eps)(b+m+n+eps)/[(m+n+1+eps)(n+1)] beta[n] + [(a+m+n)(b+m+n)/(m+n+1) - (a+m+n) - (b+m+n) - eps + (a+m+n+eps)(b+m+n+eps)/(n+1)]
//             x gamma[n]/[(n+m+1+eps)(n+1+eps)], n >= 0.
// gamma[n+1] = (a+m+n)(b+m+n)/[(m+n+1)(n+1-eps)] gamma[n], n >= 0.
//
// B(z) converges <=> |1-z| < 1
// The test of convergence is |beta[n] (1-z)^n|oo < 1E-15 |beta[0]|oo, for n large enough so that false convergence cannot occur.
//
// Variables
// ---------
// a,b,c,one_minus_z: a,b,c parameters and 1-z from z argument of 2F1(a,b,c,z)
// m,phase,m_m1,m_p1,eps,eps_pm,eps_pm_p1,a_pm,b_pm,one_meps,one_meps_pm: |Re(c-a-b)], (-1)^m, m-1, m+1, c-a-b-m, eps+m, eps+m+1, a+m, b+m, 1-eps, 1-eps-m
// eps_pa,eps_pb,eps_pa_pm,eps_pb_pm,Pi_eps,Gamma_c: eps+a, eps+b, eps+a+m, eps+b+m, pi.eps, Gamma(c)
// Gamma_inv_eps_pa_pm,Gamma_inv_eps_pb_pm,Gamma_prod: Gamma_inv(eps+a+m), Gamma_inv(eps+b+m), Gamma(c).Gamma_inv(eps+a+m).Gamma_inv(eps+b+m)
// Gamma_inv_one_meps,A_first_term,A_sum,A_term: Gamma_inv(1-eps), alpha[0], A(z), alpha[n] (1-z)^n
// pow_mzp1_m,B_first_term,prod_B,ratio: (1-z)^m, beta[0], (a)_m (b)_m (1-z)^m, (a+n)(b+n)/(n+1) for n in [0:m-2].
// B_extra_term,B_term,B_sum,B_prec: gamma[n], beta[n] (1-z)^n, B(z), 1E-15 |beta[0|oo
// cv_poly1_der_tab,cv_poly2_der_tab: P1'(X) and P2'(X) coefficients of the potentials derivatives of P1(X) and P2(X) defined in cv_poly_der_tab_calc 
//                                    with parameters a1 = a, b1 = b, c1 = 1-m-eps, z1 = 1-z and a2 = eps+b+m, b2 = eps+a+m,c2 = eps+m+1, z2 = 1-z.
// min_n: smallest integer after which false convergence cannot occur. 
//        It is calculated in min_n_calc with both P1'(X) and P2'(X), so one takes the largest integer coming from both calculations.
// possible_false_cv: always true if n < min_n. If n >= min_n, it is true if P1'(n) > 0 or P2'(n) > 0. 
//                    If n >= min_n and P1'(n) < 0 and P2'(n) < 0, it becomes false and remains as such for the rest of the calculation. 
//                    One can then check if |beta[n] z^n|oo < 1E-15 to truncate the series.
// n,n_pm_p1,n_p1,a_pm_pn,b_pm_pn: index of power series, n+m+1, n+1, a+m+n, b+m+n
// eps_pm_p1_pn,n_p1_meps,eps_pa_pm_pn,eps_pb_pm_pn,eps_pm_pn: eps+m+n+1, n+1-eps, eps+a+m+n, eps+b+m+n, eps+m+n
// prod1,prod2,prod3: (eps+a+m+n)(eps+b+m+n), (eps+m+1+n)(n+1), (a+m+n)(b+m+n)

std::complex<double> hyp_PS_one (const std::complex<double> &a,const std::complex<double> &b,const std::complex<double> &c,const std::complex<double> &one_minus_z)
{
  const int m = static_cast<int> (rint (real (c - a - b))),phase = (m%2 == 0) ? (1) : (-1),m_m1 = m - 1, m_p1 = m + 1;
  const std::complex<double> eps = c - a - b - m,eps_pm = eps + m,eps_pm_p1 = eps_pm + 1.0,a_pm = a + m,b_pm = b + m,one_meps = 1.0 - eps,one_meps_mm = one_meps - m;
  const std::complex<double> eps_pa = a + eps,eps_pb = b + eps,eps_pa_pm = eps_pa + m,eps_pb_pm = eps_pb + m,Pi_eps = M_PI*eps,Gamma_c = 1.0/Gamma_inv (c);
  const std::complex<double> Gamma_inv_eps_pa_pm = Gamma_inv (eps_pa_pm),Gamma_inv_eps_pb_pm = Gamma_inv (eps_pb_pm);
  const std::complex<double> Gamma_prod = Gamma_c*Gamma_inv_eps_pa_pm*Gamma_inv_eps_pb_pm,Gamma_inv_one_meps = Gamma_inv (one_meps);
  
  const std::complex<double> A_first_term = (m > 0) ? (Gamma_prod*A_sum_init (m,eps,Gamma_inv_one_meps)) : (0.0);
  std::complex<double> A_sum = A_first_term,A_term = A_first_term;

  if (!isfinite (A_first_term)) 
  {
    A_sum = A_term = exp (log_Gamma (c) - log_Gamma (eps_pa_pm) - log_Gamma (eps_pb_pm) + log_A_sum_init (m,eps));
    if ((imag (a) == 0.0) && (imag (b) == 0.0) && (imag (c) == 0.0)) A_sum = A_term = real (A_term);
  }

  const std::complex<double> pow_mzp1_m = pow (one_minus_z,m);
  const std::complex<double> B_first_term = B_sum_init_PS_one (a,b,c,Gamma_c,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm,Gamma_inv_eps_pb_pm,one_minus_z,m,eps)*pow_mzp1_m;

  std::complex<double> prod_B = pow_mzp1_m;
  for (int n = 0 ; n < m_m1 ; n++) 
  {
    const std::complex<double> ratio = (a + n)*(b + n)/(n + 1.0);
    A_term *= one_minus_z*ratio/(n + one_meps_mm), A_sum += A_term, prod_B *= ratio;
  }

  if (m > 0) prod_B *= (a + m - 1.0)*(b + m - 1.0)/m;

  std::complex<double> B_extra_term = prod_B*Gamma_prod*Gamma_inv_one_meps,B_term = B_first_term,B_sum = B_first_term;
  const double B_prec = 1E-15*inf_norm (B_first_term);

  double cv_poly1_der_tab[4],cv_poly2_der_tab[4];
  cv_poly_der_tab_calc (a,b,one_meps_mm,one_minus_z,cv_poly1_der_tab), cv_poly_der_tab_calc (eps_pb_pm,eps_pa_pm,eps_pm_p1,one_minus_z,cv_poly2_der_tab);

  const int min_n = fmax (min_n_calc (cv_poly1_der_tab),min_n_calc (cv_poly2_der_tab));
  bool possible_false_cv = true;

  int n = 0;
  while (possible_false_cv || (inf_norm (B_term) > B_prec))
  {
    const int n_pm_p1 = n + m_p1,n_p1 = n + 1;
    const std::complex<double> a_pm_pn = a_pm + n,b_pm_pn = b_pm + n,eps_pm_p1_pn = eps_pm_p1 + n,n_p1_meps = one_meps + n;
    const std::complex<double> eps_pa_pm_pn = eps_pa_pm + n,eps_pb_pm_pn = eps_pb_pm + n,eps_pm_pn = eps_pm + n;
    const std::complex<double> prod1 = eps_pa_pm_pn*eps_pb_pm_pn,prod2 = eps_pm_p1_pn*n_p1,prod3 = a_pm_pn*b_pm_pn;

    B_term = one_minus_z*(B_term*prod1/prod2 + B_extra_term*(prod3/n_pm_p1 - a_pm_pn - b_pm_pn - eps + prod1/n_p1)/(eps_pm_p1_pn*n_p1_meps));
    B_sum += B_term;
    B_extra_term *= one_minus_z*prod3/(n_pm_p1*n_p1_meps);
    if (possible_false_cv && (n >= min_n)) possible_false_cv = (cv_poly_der_calc (cv_poly1_der_tab,n) > 0) || (cv_poly_der_calc (cv_poly2_der_tab,n) > 0);
    n++;
  }

  return (eps == 0.0) ? (phase*(A_sum + B_sum)) : ((A_sum + B_sum)*phase*Pi_eps/sin (Pi_eps));
} 



// Calculation of the 2F1 power series converging with the 1/z transformation
// --------------------------------------------------------------------------
// The formula for F(z) in the 1/z transformation holds:
// F(z) = (-1)^m (pi.eps)/sin (pi.eps) [A(z) + B(z)] for eps not equal to zero, F(z) = (-1)^m [A(z) + B(z)] for eps = 0
// where m = |Re(b-a)], eps = b-a-m, A(z) = \sum_{n=0}^{m-1} alpha[n] z^{-n}, B(z) = \sum_{n=0}^{+oo} beta[n] z^{-n}, and:
//
// alpha[0] = [Gamma_inv(1-m-eps)/eps] Gamma_inv(c-a) Gamma_inv(a+m+eps) Gamma(c)
// [Gamma_inv(1-m-eps)/eps] is calculated in A_sum_init. 
// alpha[0] is recalculated with log[Gamma] if the previous expression overflows, and its imaginary part removed if a, b and c are real.
// alpha[n+1] = (a+n)(1-c+a+n)/[(n+1)(1-m-eps+n)] alpha[n], n in [0:m-2].
//
// beta[0] is defined in B_sum_init_PS_infinity function comments.
// gamma[0] = Gamma(c) (a)_m (1-c+a)_m z^{-m} Gamma_inv(a+m+eps) Gamma_inv(c-a) Gamma_inv(m+1) Gamma_inv(1-eps)
//
// beta[n+1] = (a+m+n+eps)(1-c+a+m+n+eps)/[(m+n+1+eps)(n+1)] beta[n] 
//           + [(a+m+n)(1-c+a+m+n)/(m+n+1) - (a+m+n) - (1-c+a+m+n) - eps + (a+m+n+eps)(1-c+a+m+n+eps)/(n+1)]
//           x gamma[n]/[(n+m+1+eps)(n+1+eps)], n >= 0.
// gamma[n+1] = (a+m+n)(b+m+n)/[(m+n+1)(n+1-eps)] gamma[n], n >= 0.
//
// B(z) converges <=> |z| > 1
// The test of convergence is |beta[n] z^{-n}|oo < 1E-15 |beta[0]|oo, for n large enough so that false convergence cannot occur.
//
// Variables
// ---------
// a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
// m,phase,m_m1,m_p1,eps,a_mc_p1,one_meps,one_meps_pm,a_pm,a_mc_p1_pm,cma: |Re(b-a)], (-1)^m, m-1, m+1, b-a-m, 1-c+a, 1-eps, 1-eps-m, a+m, 1-c+a+m, c-a
// eps_pa,eps_pm_p1,eps_pa_mc_p1_pm,Pi_eps,eps_pa_pm,eps_pm,Gamma_c: eps+a, eps+m+1, eps+1-c+a+m, pi.eps, eps+a+m, eps+m, Gamma(c)
// Gamma_inv_eps_pa_pm,Gamma_inv_cma,z_inv,pow_mz_ma: Gamma_inv(eps+a+m), Gamma_inv(c-a), 1/z, (-z)^(-a)
// Gamma_inv_one_meps,Gamma_prod: Gamma_inv(1-eps), Gamma(c) Gamma_inv(c-a) Gamma_inv(eps+a+m)
// A_first_term,A_sum,A_term: alpha[0], A(z), alpha[n] z^{-n}
// pow_z_inv_m,B_first_term,prod_B,ratio: z^{-m}, beta[0], (a)_m (1-c+a)_m z^{-m}, (a+n)(1-c+a+n)/(n+1) for n in [0:m-2].
// B_extra_term,B_term,B_sum,B_prec: gamma[n], beta[n] z^{-n}, B(z), 1E-15 |beta[0|oo
// cv_poly1_der_tab,cv_poly2_der_tab: P1'(X) and P2'(X) coefficients of the potentials derivatives of P1(X) and P2(X) defined in cv_poly_der_tab_calc 
//                                    with parameters a1 = a, b1 = 1-c+a, c1 = 1-m-eps, z1 = 1/z and a2 = b, b2 = eps+1-c+a+m,c2 = eps+m+1, z2 = 1/z.
// min_n: smallest integer after which false convergence cannot occur. 
//        It is calculated in min_n_calc with both P1'(X) and P2'(X), so one takes the largest integer coming from both calculations.
// possible_false_cv: always true if n < min_n. If n >= min_n, it is true if P1'(n) > 0 or P2'(n) > 0. 
//                    If n >= min_n and P1'(n) < 0 and P2'(n) < 0, it becomes false and remains as such for the rest of the calculation. 
//                    One can then check if |beta[n] z^n|oo < 1E-15 to truncate the series.
// n,n_pm_p1,n_p1,a_pm_pn,a_mc_p1_pm_pn,eps_pm_p1_pn,n_p1_meps,eps_pa_pm_pn,eps_pa_mc_p1_pm_pn,eps_pm_pn: 
// index of power series, n+m+1, n+1, a+m+n, 1-c+a+m+n, eps+m+n+1,n+1-eps, eps+a+m+n, eps+1-c+a+m+n, eps+m+n
// prod1,prod2,prod3: (eps+a+m+n)(eps+1-c+a+m+n), (eps+m+1+n)(n+1), (a+m+n)(1-c+a+m+n)

std::complex<double> hyp_PS_infinity (const std::complex<double> &a,const std::complex<double> &b,const std::complex<double> &c,const std::complex<double> &z)
{
  const int m = static_cast<int> (rint (real (b - a))),phase = (m%2 == 0) ? (1) : (-1), m_m1 = m - 1,m_p1 = m + 1;
  const std::complex<double> eps = b - a - m,a_mc_p1 = 1.0 - c + a,one_meps = 1.0 - eps,one_meps_mm = one_meps - m,a_pm = a + m,a_mc_p1_pm = a_mc_p1 + m,cma = c - a;
  const std::complex<double> eps_pa = eps + a, eps_pm_p1 =  eps + m + 1,eps_pa_mc_p1_pm = eps + a_mc_p1_pm,Pi_eps = M_PI*eps;
  const std::complex<double> eps_pa_pm = eps_pa + m,eps_pm = eps + m,Gamma_c = 1.0/Gamma_inv (c),Gamma_inv_eps_pa_pm = Gamma_inv (eps_pa_pm);
  const std::complex<double> Gamma_inv_cma = Gamma_inv (cma),z_inv = 1.0/z,pow_mz_ma = pow (-z,-a),Gamma_inv_one_meps = Gamma_inv (one_meps);
  const std::complex<double> Gamma_prod = Gamma_c*Gamma_inv_cma*Gamma_inv_eps_pa_pm;

  const std::complex<double> A_first_term = (m > 0) ? (Gamma_prod*A_sum_init (m,eps,Gamma_inv_one_meps)) : (0.0);
  std::complex<double> A_sum = A_first_term,A_term = A_first_term;

  if (!isfinite (A_first_term)) 
  {
    A_sum = A_term = exp (log_Gamma (c) - log_Gamma (cma) - log_Gamma (b) + log_A_sum_init (m,eps));
    if ((imag (a) == 0.0) && (imag (b) == 0.0) && (imag (c) == 0.0)) A_sum = A_term = real (A_term);
  }

  const std::complex<double> pow_z_inv_m = pow (z_inv,m);
  const std::complex<double> B_first_term = B_sum_init_PS_infinity (a,c,Gamma_c,Gamma_inv_cma,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm,z,m,eps)*pow_z_inv_m;

  std::complex<double> prod_B = pow_z_inv_m;
  for (int n = 0 ; n < m_m1 ; n++) 
  {
    const std::complex<double> ratio = (a + n)*(a_mc_p1 + n)/(n + 1.0);
    A_term *= z_inv*ratio/(n + one_meps_mm), A_sum += A_term, prod_B *= ratio;
  }

  if (m > 0) prod_B *= (a + m - 1.0)*(a_mc_p1 + m - 1.0)/m;

  std::complex<double> B_extra_term = prod_B*Gamma_prod*Gamma_inv_one_meps,B_term = B_first_term,B_sum = B_first_term;
  const double B_prec = 1E-15*inf_norm (B_first_term);
  
  double cv_poly1_der_tab[4],cv_poly2_der_tab[4];
  cv_poly_der_tab_calc (a,a_mc_p1,one_meps_mm,z_inv,cv_poly1_der_tab), cv_poly_der_tab_calc (b,eps_pa_mc_p1_pm,eps_pm_p1,z_inv,cv_poly2_der_tab);

  const int min_n = fmax (min_n_calc (cv_poly1_der_tab),min_n_calc (cv_poly2_der_tab));
  bool possible_false_cv = true;

  int n = 0;
  while (possible_false_cv || (inf_norm (B_term) > B_prec))
  {
    const int n_pm_p1 = n + m_p1,n_p1 = n + 1;
    const std::complex<double> a_pm_pn = a_pm + n,a_mc_p1_pm_pn = a_mc_p1_pm + n,eps_pm_p1_pn = eps_pm_p1 + n,n_p1_meps = one_meps + n;
    const std::complex<double> eps_pa_pm_pn = eps_pa_pm + n,eps_pa_mc_p1_pm_pn = eps_pa_mc_p1_pm + n,eps_pm_pn = eps_pm + n;
    const std::complex<double> prod1 = eps_pa_pm_pn*eps_pa_mc_p1_pm_pn,prod2 = eps_pm_p1_pn*n_p1,prod3 = a_pm_pn*a_mc_p1_pm_pn;

    B_term = z_inv*(B_term*prod1/prod2 + B_extra_term*(prod3/n_pm_p1 - a_pm_pn - a_mc_p1_pm_pn - eps + prod1/n_p1)/(eps_pm_p1_pn*n_p1_meps));
    B_sum += B_term;
    B_extra_term *= z_inv*prod3/(n_pm_p1*n_p1_meps);
    if (possible_false_cv && (n >= min_n)) possible_false_cv = (cv_poly_der_calc (cv_poly1_der_tab,n) > 0) || (cv_poly_der_calc (cv_poly2_der_tab,n) > 0);
    n++;
  }

  return (eps == 0.0) ? (phase*pow_mz_ma*(A_sum + B_sum)) : ((A_sum + B_sum)*phase*pow_mz_ma*Pi_eps/sin (Pi_eps));
}





// Calculation of F(z) in transformation theory missing zones of the std::complex plane with a Taylor series
// ----------------------------------------------------------------------------------------------------
// If z is close to exp(+/- i.pi/3), no transformation in 1-z, z, z/(z-1) 
// or combination of them can transform z in a std::complex number of modulus smaller than a given Rmax < 1 .
// Rmax is a radius for which one considers power series summation for |z| > Rmax is too slow to be processed. One takes Rmax = 0.9 .
// Nevertheless, for Rmax = 0.9, these zones are small enough to be handled 
// with a Taylor series expansion around a point z0 close to z where transformation theory can be used to calculate F(z).
// One then chooses z0 to be 0.9 z/|z| if |z| < 1, and 1.1 z/|z| if |z| > 1, 
// so that hyp_PS_zero or hyp_PS_infinity can be used (see comments of these functions above).
// For this z0, F(z) = \sum_{n=0}^{+oo} q[n] (z-z0)^n, with:
// q[0] = F(z0), q[1] = F'(z0) = (a b/c) 2F1(a+1,b+1,c+1,z0)
// q[n+2] = [q[n+1] (n (2 z0 - 1) - c + (a+b+c+1) z0) + q[n] (a+n)(b+n)/(n+1)]/(z0(1-z0)(n+2))
// As |z-z0| < 0.1, it converges with around 15 terms, so that no instability can occur for moderate a, b and c.
// Convergence is tested with |q[n] (z-z0)^n|oo + |q[n+1] (z-z0)^{n+1}|oo. Series is truncated when this test is smaller than 1E-15 (|q[0]|oo + |q[1] (z-z0)|oo).
// No false convergence can happen here as q[n] behaves smoothly for n -> +oo.
//
// Variables
// ---------
// a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
// abs_z,is_abs_z_small: |z|, true if |z| < 1 and false if not.
// z0,zc_z0_ratio,z0_term1,z0_term2: 0.9 z/|z| if |z| < 1, and 1.1 z/|z| if |z| > 1, (z-z0)/(z0 (1-z0)), 2 z0 - 1, c - (a+b+c+1) z0
// hyp_PS_z0,dhyp_PS_z0,prec: F(z0), F'(z0) calculated with 2F1 as F'(z0) = (a b/c) 2F1(a+1,b+1,c+1,z0),
// precision demanded for series truncation equal to 1E-15 (|q[0]|oo + |q[1] (z-z0)|oo).
// n,an,anp1,anp2,sum: index of the series, q[n] (z-z0)^n, q[n+1] (z-z0)^{n+1}, q[n+2] (z-z0)^{n+2}, truncated sum of the power series.

std::complex<double> hyp_PS_complex_plane_rest (const std::complex<double> &a,const std::complex<double> &b,const std::complex<double> &c,const std::complex<double> &z)
{
  const double abs_z = abs (z);
  const bool is_abs_z_small = abs_z < 1.0;
  
  const std::complex<double> z0 = (is_abs_z_small) ? (0.9*z/abs_z) : (1.1*z/abs_z),zc = z - z0;
  const std::complex<double> zc_z0_ratio = zc/(z0*(1.0 - z0)),z0_term1 = 2.0*z0 - 1.0,z0_term2 = c - (a + b + 1.0)*z0;
  
  const std::complex<double> hyp_PS_z0 = (is_abs_z_small) ? (hyp_PS_zero (a,b,c,z0)) : (hyp_PS_infinity (a,b,c,z0));
  const std::complex<double> dhyp_PS_z0 = (is_abs_z_small) ? (hyp_PS_zero (a + 1.0,b + 1.0,c + 1.0,z0)*a*b/c) : (hyp_PS_infinity (a + 1.0,b + 1.0,c + 1.0,z0)*a*b/c);

  int n = 0;
  std::complex<double> an = hyp_PS_z0, anp1 = zc*dhyp_PS_z0, sum = an + anp1;

  const double prec = 1E-15*(inf_norm (an) + inf_norm (anp1));

  while (inf_norm (an) + inf_norm (anp1) > prec)
  {
    const std::complex<double> anp2 = zc_z0_ratio*(anp1*(n*z0_term1 - z0_term2) + an*zc*(a + n)*(b + n)/(n + 1))/(n + 2);
    sum += anp2;
    n++;
    an = anp1;
    anp1 = anp2;
  }

  return sum;
}




// Calculation of F(z) for arbitrary z using previous routines
// -----------------------------------------------------------
// Firstly, it is checked if a,b and c are negative integers.
// If neither a nor b is negative integer but c is, F(z) is undefined so that the program stops with an error message.
// If a and c are negative integers with c < a, or b and c are negative integers with b < a, 
// or c is not negative integer integer but a or b is, one is in the polynomial case.
// In this case, if |z| < |z/(z-1)| or z = 1, hyp_PS_zero is used directly, as then |z| <= 2 
// and no instability arises with hyp_PS_zero as long the degree of the polynomial is small (<= 10 typically).
// If not, one uses the transformation F(z) = (1-z)^{-a} 2F1(a,c-b,c,z/(z-1)) if a is negative integer 
// or F(z) = (1-z)^{-b} 2F1(b,c-a,c,z/(z-1)) if b is negative integer along with hyp_PS_zero.
// Indeed, 2F1(a,c-b,c,X) is a polynomial if a is negative integer, and so is 2F1(b,c-a,c,X) if b is negative integer, 
// so that one has here |z/(z-1)| <= 2 and the stability of the method is the same as for the |z| < |z/(z-1)| case.
// If one is in the non-polynomial case, one checks if z >= 1. If it is, one is the cut of F(z) so that z is replaced by z - 10^{-307}i.
// Then, using F(z) = 2F1(b,a,c,z) and F(z) = (1-z)^{c-a-b} 2F1(c-a,c-b,c,z), 
// one replaces a,b,c parameters by combinations of them so that Re(b-a) >= 0 and Re(c-a-b) >= 0.
// Exchanging a and b does not change convergence properties, while having Re(c-a-b) >= 0 accelerates it (In hyp_PS_zero, t[n] z^n ~ z^n/(n^{c-a-b}) for n -> +oo).
// If |1-z| < 1E-5, one uses hyp_PS_one as the vicinity of the singular point z = 1 is treated properly.
// After that, one compares |z| and |z/(z-1)| to R in {0.5,0.6,0.7,0.8,0.9}. 
// If one of them is smaller than R, one uses hyp_PS_zero without transformation or with the transformation F(z) = (1-z)^{-a} 2F1(a,c-b,c,z/(z-1)).
// Then, if both of them are larger than 0.9, one compares |1/z|, |(z-1)/z|, |1-z| and |1/(1-z)| to R in {0.5,0.6,0.7,0.8,0.9}. 
// If one of them is found smaller than R, 
// with the condition that |c-b|oo < 5 for (z-1)/z transformation, |a,b,c|oo < 5 for |1-z| transformation and |a,c-b,c|oo < 5 for |1/(1-z)| transformation,
// the corresponding transformation is used. If none of them was smaller than 0.9, 
// one is in the missing zones of transformation theory so that the Taylor series of hyp_PS_complex_plane_rest is used.
//
// Variables
// ---------
// a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
// Re_a,Re_b,Re_c,na,nb,nc: real parts of a,b,c, closest integers to a,b,c.
// is_a_neg_int,is_b_neg_int,is_c_neg_int: true if a,b,c is negative integers and false if not.
// zm1,z_over_zm1,z_shift: z-1, z/(z-1), z - 10^{-307}i in case z >= 1.
// ab_condition, cab_condition: true if Re(b-a) >= 0 and false if not, true if Re(c-a-b) >= 0 and false if not.
// abs_zm1,abz_z,abs_z_inv,abs_z_over_zm1,abs_zm1_inv,abs_zm1_over_z: |z-1|, |z|, |1/z|, |z/(z-1)|, |1/(z-1)|, |(z-1)/z|
// are_ac_small: true if |a|oo < 5 and |c|oo < 5, false if not.
// is_cmb_small: true if |c-b|oo < 5, false if not.
// are_abc_small: true if |a|oo < 5, |b|oo < 5 and |c|oo < 5, false if not.
// are_a_cmb_c_small: true if |a|oo < 5, |c-b|oo < 5 and |c|oo < 5, false if not.
// R_tab,R: table of radii {0.5,0.6,0.7,0.8,0.9}, one of these radii.

std::complex<double> hyp_2F1 (const std::complex<double> &a,const std::complex<double> &b,const std::complex<double> &c,const std::complex<double> &z)
{
  const double Re_a = real (a), Re_b = real (b), Re_c = real (c);

  const int na = static_cast<int> (rint (Re_a)),nb = static_cast<int> (rint (Re_b)),nc = static_cast<int> (rint (Re_c));
  const bool is_a_neg_int = (a == na) && (na <= 0),is_b_neg_int = (b == nb) && (nb <= 0),is_c_neg_int = (c == nc) && (nc <= 0);
  const std::complex<double> zm1 = z-1.0;

  if (is_c_neg_int)
  {
    if (is_a_neg_int && (nc < na))
    {
      const std::complex<double> z_over_zm1 = z/zm1;
      return ((z == 1.0) || (abs (z) < abs (z_over_zm1))) ? (hyp_PS_zero (a,b,c,z)) : (pow (-zm1,-a)*hyp_PS_zero (a,c-b,c,z_over_zm1));
    }
    else if (is_b_neg_int && (nc < nb))
    {
      const std::complex<double> z_over_zm1 = z/zm1;
      return ((z == 1.0) || (abs (z) < abs (z_over_zm1))) ? (hyp_PS_zero (a,b,c,z)) : (pow (-zm1,-b)*hyp_PS_zero (b,c-a,c,z_over_zm1));
    }
    else 
      std::cout<<"2F1 undefined"<<std::endl, abort ();
  }

  if (is_a_neg_int)
  {
    const std::complex<double> z_over_zm1 = z/zm1;
    return ((z == 1.0) || (abs (z) < abs (z_over_zm1))) ? (hyp_PS_zero (a,b,c,z)) : (pow (-zm1,-a)*hyp_PS_zero (a,c-b,c,z_over_zm1));
  }
  else if (is_b_neg_int)
  {
    const std::complex<double> z_over_zm1 = z/zm1;
    return ((z == 1.0) || (abs (z) < abs (z_over_zm1))) ? (hyp_PS_zero (a,b,c,z)) : (pow (-zm1,-b)*hyp_PS_zero (b,c-a,c,z_over_zm1));
  }

  const std::complex<double> z_shift(real (z),-1E-307);
  if ((real (z) >= 1.0) && (imag (z) == 0.0)) return hyp_2F1 (a,b,c,z_shift);

  const bool ab_condition = (Re_b >= Re_a), cab_condition = (Re_c >= Re_a + Re_b);
  if (!ab_condition || !cab_condition)
  {
    if (!ab_condition && cab_condition) return hyp_2F1 (b,a,c,z);
    else if (!cab_condition && ab_condition) return pow (-zm1,c-a-b)*hyp_2F1 (c-b,c-a,c,z);
    else return pow (-zm1,c-a-b)*hyp_2F1 (c-a,c-b,c,z);
  }

  const double abs_zm1 = abs (zm1);
  if (abs_zm1 < 1E-5) return hyp_PS_one (a,b,c,-zm1);

  const double abs_z = abs (z),abs_z_inv = 1.0/abs_z,abs_z_over_zm1 = abs_z/abs_zm1,abs_zm1_inv = 1.0/abs_zm1,abs_zm1_over_z = 1.0/abs_z_over_zm1;
  const bool are_ac_small = (inf_norm (a) < 5.0) && (inf_norm (c) < 5.0),is_cmb_small = (inf_norm (c-b) < 5.0);
  const bool are_abc_small = are_ac_small && (inf_norm (b) < 5.0),are_a_cmb_c_small = are_ac_small && is_cmb_small;
  
  const double R_tab[5] = {0.5,0.6,0.7,0.8,0.9};
  
  for (unsigned int i = 0 ; i < 5 ; i++)
  {
    const double R = R_tab[i];

    if (abs_z <= R) return hyp_PS_zero (a,b,c,z);
    if (is_cmb_small && (abs_z_over_zm1 <= R)) return pow (-zm1,-a)*hyp_PS_zero (a,c-b,c,z/zm1);
  }
  
  for (unsigned int i = 0 ; i < 5 ; i++)
  {
    const double R = R_tab[i];

    if (abs_z_inv <= R) return hyp_PS_infinity (a,b,c,z);
    if (is_cmb_small && (abs_zm1_over_z <= R)) return pow (-zm1,-a)*hyp_PS_infinity (a,c-b,c,z/zm1);

    if (are_abc_small && (abs_zm1 <= R)) return hyp_PS_one (a,b,c,-zm1);
    if (are_a_cmb_c_small && (abs_zm1_inv <= R)) return pow (-zm1,-a)*hyp_PS_one (a,c-b,c,-1.0/zm1);
  }

  return hyp_PS_complex_plane_rest (a,b,c,z);  
}



// Test of 2F1 numerical accuracy using hypergeometric differential equation
// -------------------------------------------------------------------------
// If z = 0, F(z) = 1 so that this value is trivially tested.
// To test otherwise if the value of F(z) is accurate, one uses the fact that z(z-1) F''(z) + (c - (a+b+1) z) F'(z) - a b F(z) = 0.
// If z is not equal to one, a relative precision test is provided by |F''(z) + [(c - (a+b+1) z) F'(z) - a b F(z)]/[z(z-1)]|oo/(|F(z)|oo + F'(z)|oo + |F''(z)|oo).
// If z is equal to one, one uses |(c - (a+b+1)) F'(z) - a b F(z)|oo/(|F(z)|oo + F'(z)|oo + 1E-307).
// F'(z) and F''(z) are calculated using equalities F'(z) = (a b/c) 2F1(a+1,b+1,c+1,z) and F'(z) = ((a+1)(b+1)/(c+1)) (a b/c) 2F1(a+2,b+2,c+2,z).
//
// Variables
// ---------
// a,b,c,z: a,b,c parameters and z argument of 2F1(a,b,c,z)
// F,dF,d2F: F(z), F'(z) and F''(z) calculated with hyp_2F1 using F'(z) = (a b/c) 2F1(a+1,b+1,c+1,z) and F'(z) = ((a+1)(b+1)/(c+1)) (a b/c) 2F1(a+2,b+2,c+2,z).

double test_2F1 (const std::complex<double> &a,const std::complex<double> &b,const std::complex<double> &c,const std::complex<double> &z,const std::complex<double> &F)
{
  if (z == 0.0)
    return inf_norm (F - 1.0);
  else if (z == 1.0)
  {
    const std::complex<double> dF = hyp_2F1(a+1.0,b+1.0,c+1.0,z)*a*b/c;

    return inf_norm ((c - (a+b+1.0))*dF - a*b*F)/(inf_norm (F) + inf_norm (dF) + 1E-307);
  }
  else 
  {
    const std::complex<double> dF = hyp_2F1(a+1.0,b+1.0,c+1.0,z)*a*b/c;
    const std::complex<double> d2F = hyp_2F1(a+2.0,b+2.0,c+2.0,z)*a*(a+1.0)*b*(b+1.0)/(c*(c+1.0));

    return inf_norm (d2F + ((c - (a+b+1.0)*z)*dF - a*b*F)/(z*(1.0-z)))/(inf_norm (F) + inf_norm (dF) + inf_norm (d2F));
  }
}
