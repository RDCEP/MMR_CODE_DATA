$ontext
Solve G-M model with certain kappa 

Authors: Yongyang Cai, Hoover Institution and University of Chicago

If using material from this code, the user should cite the following paper:
Cai, Yongyang, and Alan H. Sanstad (2016). Model uncertainty and energy technology policy: 
the example of induced technical change. Computers & Operations Research 66, 362-373. 

$offtext

SETS	t		stages				/ 1*491 /;

PARAMETERS
  E0(t)			emission rate
  Emissions(t)
  kappa			
;

$GDXIN EmissionRCP85.gdx
$load Emissions

E0(t) = Emissions(t) / 2.12;

SCALARS
  deltat		time unit			/ 1 /
  rho 			discount rate			/ 0.05 /
  M_C			abatement parameter		/ 83 /
  alpha_C1		abatement parameter		/ 3 /
  alpha_C2		abatement parameter		/ 2 /
  alpha			rate of autonomous change	/ 0.005 /
  gamma			technological change parameter	/ 0.5 /
  phi			technological change parameter 	/ 0.5 /
  H0			initial knowledge		/ 1 /
  kappa0		purely autonomous technical change	/ 0 /
  kappa1		both autonomous and induced technical change	/ 1 /
  beta			parameter for equation fo motion	/ 0.64 /
  delta			parameter for equation fo motion	/ 0.008 /
  S0			initial CO2 concentration in ppm	/ 392 /
  S_PRE			pre-industrial CO2 concentration in ppm	/ 278 /
  M_D			damage parameter			/ 0.0012 /
  alpha_D		damage parameter			/ 2 /
  MPsi			technological change parameter 	/ 0.0022 /
  deg			degree for cost function 		/ 1 /
;


*******************************
*******************************

* Define variables
POSITIVE VARIABLES
  A(t) 		abatement
  H(t)	 	stock of knowledge
  S(t)		CO2 concentration
;

VARIABLES
  TotalCostDamage		total cost and damage
;


EQUATIONS
  TotalCostDamageFun
  TransitionH(t)
  TransitionS(t)
;

TotalCostDamageFun..	TotalCostDamage =e= sum( t$(ord(t)<card(t)), exp(-deltat*rho*(ord(t)-1)) * deltat * 
	( M_C*A(t)**alpha_C1/((E0(t)-A(t))**alpha_C2*H(t)**deg ) + M_D*S(t)**alpha_D ) );

TransitionH(t)$(ord(t)<card(t)).. 	(H(t+1)-H(t))/deltat =e= alpha*H(t) + kappa*MPsi*H(t)**phi*A(t)**gamma;
TransitionS(t)$(ord(t)<card(t)).. 	(S(t+1)-S(t))/deltat =e= beta*(E0(t)-A(t)) - delta*(S(t)-S_PRE);

*******************************
**  Upper and Lower Bounds

A.lo(t) = 0.001;
A.up(t) = E0(t)-0.001;
H.lo(t) = 0.001;

A.l(t) = E0(t)/5;
H.l('1') = H0;
S.l('1') = S0;
kappa = (kappa0+kappa1)/2;
loop(t$(ord(t)<card(t)),
  H.l(t+1) = H.l(t) + deltat*( alpha*H.l(t) + kappa*MPsi*H.l(t)**phi*A.l(t)**gamma );
  S.l(t+1) = S.l(t) + deltat*( beta*(E0(t)-A.l(t)) - delta*(S.l(t)-S_PRE) );
);

H.fx('1') = H0;
S.fx('1') = S0;

*******************************

** Solution options
option iterlim = 99900;
option reslim = 99999;
option solprint = on;
option limrow = 0;
option limcol = 0;
option nlp = conopt;

model TotalCostDamageProblem / all /;

kappa = (kappa0+kappa1)/2;
solve TotalCostDamageProblem minimizing TotalCostDamage using nlp ;

display A.l;

file AnnualConstantKappa_mean_sol;
put AnnualConstantKappa_mean_sol;
AnnualConstantKappa_mean_sol.nw = 12;
AnnualConstantKappa_mean_sol.nr = 2;
AnnualConstantKappa_mean_sol.nz = 1e-15;
loop(t,
    put t.tl:4:0;
    put S.l(t):14:6;
    put A.l(t):14:6;
    put E0(t):14:6;
    put /;
);


kappa = kappa0;
solve TotalCostDamageProblem minimizing TotalCostDamage using nlp ;

display A.l;

file AnnualConstantKappa_0_sol;
put AnnualConstantKappa_0_sol;
AnnualConstantKappa_0_sol.nw = 12;
AnnualConstantKappa_0_sol.nr = 2;
AnnualConstantKappa_0_sol.nz = 1e-15;
loop(t,
    put t.tl:4:0;
    put S.l(t):14:6;
    put A.l(t):14:6;
    put E0(t):14:6;
    put /;
);


kappa = kappa1;
solve TotalCostDamageProblem minimizing TotalCostDamage using nlp ;

display A.l;

file AnnualConstantKappa_1_sol;
put AnnualConstantKappa_1_sol;
AnnualConstantKappa_1_sol.nw = 12;
AnnualConstantKappa_1_sol.nr = 2;
AnnualConstantKappa_1_sol.nz = 1e-15;
loop(t,
    put t.tl:4:0;
    put S.l(t):14:6;
    put A.l(t):14:6;
    put E0(t):14:6;
    put /;
);

