

*******get simulated data********;
proc import datafile = '1-data\simdata.csv'
 out = WORK._dat dbms = CSV ;
run;


***sort by id***;
proc sort data=_dat; by id visit; 
run;

***curate data***;
data _dat; set _dat; by id;
_last= last.id; /*for non time-varying survival models, just use last or first obs;*/
_first=first.id;
run;

data _dat; set _dat; 
  time= years/20;
  brainloss_time = brainloss*time; /*create brainloss*time interactions;*/
  run;


**************************************************************************;
**************************************************************************;
**step1: run separate LDA model to get initial estimates for joint model***;
**************************************************************************;
***************************************************************************;

proc mixed data=_dat noclprint method=ml cl covtest empirical ;  
  class id ;
  model globz =  brainloss time brainloss_time age0 male   / s cl;
  random intercept time  / s subject=id type=un; 
  ods output solutionF=modelbased;
  ods output CovParms=covp;
  ods output SolutionR=sr(keep= effect id estimate); /*output predicted random effet*/
run;


***format predicted random effect for merging***;
proc transpose data=sr out=re prefix=estimate;
  by id;
  id effect;
  var estimate;
run;

***merge predicted random effect to _dat**;
data _dati;  merge _dat  Re(drop= _name_ rename=(estimateIntercept=u0i estimatetime=u1i)); by id; run;


**************************************************************************;
**************************************************************************;
*step2: run separate EVENT model to get initial estimates for joint model*;
**************************************************************************;
***************************************************************************;

proc lifereg data=_dati(where=(_last=1));
model demyears*dementia(0)= brainloss age0 male u0i u1i/dist=weibull;
ods output ParameterEstimates=Pesurv;
run;

*reparameterize output from proc lifereg to change estimates to same scale as proc nlmixed*******;
%let scale = 0.3126; *0.3126 is the scale from proc lifereg;
data Pesurv;
set Pesurv(drop =DF);
 if _N_ <7 then Estimate = -estimate/&scale; 
run;


**************************************************************************;
**************************************************************************;
**********step3: run joint model, initializing from step 1****************;
**************************************************************************;
**************************************************************************;

***curate data***;
data _dat; set _dat(where=(not missing(globz))); by id;
  _event1= dementia;
  _time1= demyears;
  _last = last.id; 
run;



***run first joint model***;
  proc nlmixed data=_dat  empirical gconv=0; 
  parms  a0=1.3162 a1=-0.1103 a2=-0.9009 a3=-0.2134 a4=-0.02718 a5=-0.1673 
         tau11=1.0273 tau21=-0.08817 tau22=1.0706 sigma=0.4042  /*LMM : initialize from step 1*/
         b0=-18.0332 b1=0.5600 b2=0.09530 b3=0.05985  scale1=0.3126 
         r10=-0.9467 r11=-0.9212; /*event model:initialize from reparametized estimates of  step 2*/

  eta0  = a0 + a1*brainloss + a2*time + a3*brainloss_time 
        + a4*age0 + a5*male ;   /*LMM:linear predictors;*/
  eta1 = b0 + b1*brainloss + b2*age0 + b3*male ;/*event:linear predictors;*/

  linpredLDA = eta0 + u0i + time*u1i; /*LMM submodel with RE;*/

  resid = (globz-linpredLDA); /*Gaussian residual term;*/

  s2 = sigma**2; /*Joint LMM submodels;*/
  if (abs(resid) > 1.3E100) or (s2 < 1e-12) then do;
       loglikeLDA = -1e20;
        end; 
  else do;
       loglikeLDA = -0.5*(1.837876 + resid**2 / s2  + log(s2));  /* Joint LDA subModel Likelihood specification;*/
       end;
  
   if (_last) then do;  /*Joint Event submodels;*/

        linpredEvent1 = eta1 + r10*u0i + r11*u1i ; /*event model with random ints & slopes;*/
        *Weibull survival;
         gam1   = 1/scale1;
         alph1  = exp(linpredEvent1);
         loglikeEvent1 = (_event1=1)*(-alph1*(_time1**gam1) + log(gam1) + linpredEvent1 + (gam1-1)*log(_time1))
                       + (_event1=0)*(-alph1*(_time1**gam1));  /* Dementia Submodel likelihood;*/
        
         end; 
    else do;
         loglikeEvent1=0;
         end;

    model _last ~ general(loglikeLDA + loglikeEvent1); /*general log-likelihood specification;*/

    random u0i u1i ~ normal([0, 0],[tau11,tau21,tau22]) subject=id; /*RE distb;*/

    ods output ParameterEstimates=newests0 (drop = Gradient rename =(StandardError = StdErr Probt = pvalue));

    run;
	title"first joint model";
 


**************************************************************************;
**************************************************************************;
**********step4: final  joint model, initializing from step 2**************;
**********to reach final estimates****************************************;
**************************************************************************;

proc nlmixed data=_dat    empirical gconv=0; 
  parms  / data = newests0; /*parameter initializing from step 3 ;*/
  eta0  = a0 + a1*brainloss + a2*time + a3*brainloss_time 
        + a4*age0 + a5*male ;   /*linear predictors;*/
  eta1 = b0 + b1*brainloss + b2*age0 + b3*male ;

  linpredLDA = eta0 + u0i + time*u1i; /* *LMM submodel with random ints & slopes;*/

      resid = (globz-linpredLDA); /*Gaussian residual term;*/
    s2 = sigma**2;

    if (abs(resid) > 1.3E100) or (s2 < 1e-12) then do;
       loglikeLDA = -1e20;
       end;
    else do;
       loglikeLDA = -0.5*(1.837876 + resid**2 / s2  + log(s2));  /* Joint LDA subModel Likelihood specification;*/
       end;

    if (_last) then do;  /*Joint Event submodels;*/

        linpredEvent1 = eta1 + r10*u0i + r11*u1i ; /*event model with random ints & slopes;*/
        *Weibull survival;
         gam1   = 1/scale1;
         alph1  = exp(linpredEvent1);
         loglikeEvent1 = (_event1=1)*(-alph1*(_time1**gam1) + log(gam1) + linpredEvent1 + (gam1-1)*log(_time1))
                       + (_event1=0)*(-alph1*(_time1**gam1));  /* Dementia Submodel likelihood;*/
         end; 
    else do;
      loglikeEvent1=0;
      end;

    model _last ~ general(loglikeLDA + loglikeEvent1); /* general log-likelihood specification;*/

     random u0i u1i ~ normal([0, 0],[tau11,tau21,tau22]) subject=id; /*random effects distb;*/

     ods output ParameterEstimates=newests (drop = Gradient rename =(StandardError = StdErr Probt = pvalue));

    run;
	title"last joint model";


