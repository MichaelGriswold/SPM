**************************************************************************
***** Step 0: Data curation
**************************************************************************

** get data;
import delimited "1-data\simdata.csv", clear
  destring globz ,replace force  // globaz was read as string
  gen time = years/20 // scale time to 20 year effect
  global adj male age0 // macro variable with adjustment vars for all models

	
**************************************************************************
***** Step 1: Initial (separate) LDA estimates
**************************************************************************
gsem (globz  <- c.brainloss##c.time $adj U0[id]@1 c.time#U1[id]@1) 
  matrix lda = e(b) // store LDA beta & var estimates
  capture drop _b0i _b1i
  predict _b0i, latent(U0[id]) // EB random intercept ests
  predict _b1i, latent(U1[id]) // EB random slope ests
	
	
**************************************************************************
***** Step 2: Initial (separate) Event estimates
**************************************************************************
** Run gsem, get initial Event parameter ests (random intercept only)
gsem (demyears <- brainloss $adj _b0i _b1i, family(weibull, failure(dementia))) 
  matrix event = e(b) // store Event lambda, alpha & rho estimates

	
**************************************************************************
***** Step 3: Improve initial estimates
**************************************************************************	
** Run submodels together, but use EB estimated latent effects in EVENT submodel
matrix init0 = lda, event // stack initial ests   
gsem ///
	(globz  <- c.brainloss##c.time $adj U0[id]@1 c.time#U1[id]@1) /// LMM submodel
	(demyears <- brainloss $adj _b0i _b1i, family(weibull, failure(dementia))) /// Dementia Event submodel
	 , from(init0, skip)				
  matrix init1 = e(b) // store estimates
	** IMPORTANT: names in the initialization vector must match what gsem is looking for
	global oldnames : colfullnames init1
	global newnames = regexr("$oldnames","demyears:_b0i","demyears:U0[id]")
	global newnames = regexr("$newnames","demyears:_b1i","demyears:U1[id]")
	matrix colnames init1 = $newnames

**************************************************************************
***** Step 4: Final SPM estimates
**************************************************************************
** 4.1 Model Based SEs
gsem ///
  (globz  <- c.brainloss##c.time $adj U0[id]@1 c.time#U1[id]@1) ///
  (demyears <- brainloss $adj U0[id]@r0 U1[id]@r1, fam(weib, fail(dementia))) ///
   , from(init1, skip)		
   matrix spm_modelse = e(b) // store estimates

** 4.1 Huber-White semi-Robust Based SEs
gsem ///
  (globz  <- c.brainloss##c.time $adj U0[id]@1 c.time#U1[id]@1) ///
  (demyears <- brainloss $adj U0[id]@r0 U1[id]@r1, fam(weib, fail(dementia))) ///
   , from(spm_modelse, skip) vce(robust)
	 
**************************************************************************	 
** For comparison: full/complete data	 estimates
gsem (globz_complete  <- c.brainloss##c.time $adj U0[id]@1 c.time#U1[id]@1) 
 	 
	 
	 
	 
