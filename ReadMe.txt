There are four GAMS code files in this package, for the examples in the following paper 
Cai, Yongyang, and Alan H. Sanstad (2016). Model uncertainty and energy technology policy: 
the example of induced technical change. Computers & Operations Research 66, 362-373. 

The code file AnnualConstantKappa.gms is solving G-M model with a fixed kappa (which 
could be mean of uncertain kappa); AnnualUncertainKappa_uniform.gms solve solves the 
G-M problem with uniformly distributed kappa using the expected cost minimization method 
discussed in the paper; The code AnnualUncertainKappa_BetaDistr.gms is similar to 
AnnualUncertainKappa_uniform.gms except that kappa is assumed to have Beta distribution. 
The code file MinmaxRegret.gms solves the robust decision making problem with uncertain 
kappa but unknown distributions using the computational min-max regret method introduced 
in the paper. 

Running these code files require a full version of GAMS and the data file EmissionRCP85.gdx 
(which is provided in the package). If you do not have the full version of GAMS, you can 
submit the code files to the online NEOS server (using the GAMS version of CONOPT): 
https://neos-server.org/neos/ for a free run. 
 