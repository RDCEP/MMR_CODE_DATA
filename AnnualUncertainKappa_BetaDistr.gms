$ontext
Solve G-M model with uncertain kappa with Beta distribution

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

set        n		nodes for integration		/ 1*101 /;
alias(n, evenn);

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

PARAMETERS
  kappas(n)        	nodes for integration
  weights(n)		weights for integration
;

kappas(n) = kappa0 + kappa1*(ord(n)-1)/(card(n)-1);

*******************************
*******************************

* Define variables
POSITIVE VARIABLES
  A(t) 		abatement
  H(t,n) 	stock of knowledge
  S(t)		CO2 concentration
;

VARIABLES
  TotalCostDamage		total cost and damage
;


EQUATIONS
  TotalCostDamageFun
  TransitionH(t,n)
  TransitionS(t)
;

TotalCostDamageFun..	TotalCostDamage =e= sum(n, weights(n) * sum( t$(ord(t)<card(t)), exp(-deltat*rho*(ord(t)-1)) * deltat * 
	( M_C*A(t)**alpha_C1/((E0(t)-A(t))**alpha_C2*H(t,n)**deg) + M_D*S(t)**alpha_D ) ));

TransitionH(t,n)$(ord(t)<card(t)).. 	(H(t+1,n)-H(t,n))/deltat =e= alpha*H(t,n) + kappas(n)*MPsi*H(t,n)**phi*A(t)**gamma;
TransitionS(t)$(ord(t)<card(t)).. 	(S(t+1)-S(t))/deltat =e= beta*(E0(t)-A(t)) - delta*(S(t)-S_PRE);

*******************************
**  Upper and Lower Bounds

A.lo(t) = 0.001;
A.up(t) = E0(t)-0.001;
H.lo(t,n) = 0.001;

A.l(t) = E0(t)/5;
H.l('1',n) = H0;
S.l('1') = S0;
loop(t$(ord(t)<card(t)),
  loop(n,
    H.l(t+1,n) = H.l(t,n) + deltat*( alpha*H.l(t,n) + kappas(n)*MPsi*H.l(t,n)**phi*A.l(t)**gamma );
  );
  S.l(t+1) = S.l(t) + deltat*( beta*(E0(t)-A.l(t)) - delta*(S.l(t)-S_PRE) );
);

H.fx('1',n) = H0;
S.fx('1') = S0;

*******************************

** Solution options
option iterlim = 99900;
option reslim = 99999;
option solprint = off;
option limrow = 0;
option limcol = 0;
option nlp = conopt;

model TotalCostDamageProblem / all /;

parameters beta_alpha
	beta_beta;

beta_alpha = 1;
beta_beta = 3;
weights(n) = (kappa1-kappa0)/(card(n)-1)*2/3 * sqr(1-kappas(n)) * 3;
weights('1') = (kappa1-kappa0)/(card(n)-1)/3 * sqr(1-kappas('1')) * 3;
weights(n)$(ord(n)=card(n)) = (kappa1-kappa0)/(card(n)-1)/3 * sqr(1-kappas(n)) * 3;
loop(evenn$(ord(evenn)<=card(n)/2),
  weights(n)$(ord(n)=2*ord(evenn)) = (kappa1-kappa0)/(card(n)-1)*4/3 * sqr(1-kappas(n)) * 3;
);

solve TotalCostDamageProblem minimizing TotalCostDamage using nlp ;
solve TotalCostDamageProblem minimizing TotalCostDamage using nlp ;


file AnnualUncertainKappa_BetaDistr1_sol;
put AnnualUncertainKappa_BetaDistr1_sol;
AnnualUncertainKappa_BetaDistr1_sol.nw = 12;
AnnualUncertainKappa_BetaDistr1_sol.nr = 2;
AnnualUncertainKappa_BetaDistr1_sol.nz = 1e-15;
loop(t,
    put t.tl:4:0;
    put S.l(t):14:6;
    put A.l(t):14:6;
    put E0(t):14:6;
    put /;
);


beta_alpha = 3;
beta_beta = 1;
weights(n) = (kappa1-kappa0)/(card(n)-1)*2/3 * sqr(kappas(n)) * 3;
weights('1') = (kappa1-kappa0)/(card(n)-1)/3 * sqr(kappas('1')) * 3;
weights(n)$(ord(n)=card(n)) = (kappa1-kappa0)/(card(n)-1)/3 * sqr(kappas(n)) * 3;
loop(evenn$(ord(evenn)<=card(n)/2),
  weights(n)$(ord(n)=2*ord(evenn)) = (kappa1-kappa0)/(card(n)-1)*4/3 * sqr(kappas(n)) * 3;
);

solve TotalCostDamageProblem minimizing TotalCostDamage using nlp ;
solve TotalCostDamageProblem minimizing TotalCostDamage using nlp ;


file AnnualUncertainKappa_BetaDistr2_sol;
put AnnualUncertainKappa_BetaDistr2_sol;
AnnualUncertainKappa_BetaDistr2_sol.nw = 12;
AnnualUncertainKappa_BetaDistr2_sol.nr = 2;
AnnualUncertainKappa_BetaDistr2_sol.nz = 1e-15;
loop(t,
    put t.tl:4:0;
    put S.l(t):14:6;
    put A.l(t):14:6;
    put E0(t):14:6;
    put /;
);
