// 	$Id: bestParamUSSRV.cpp 920 2006-09-21 09:26:12Z ninio $	
#include "bestParamUSSRV.h"

/* structure of this method:
(1) checks of the number of parameters to optimize, and decide how many parameters optimizations iteration,
and how many parameters+bbl iterations will be done.
(2) A loop over the parameters+bbl iterations
	(2.1) A loop over the parameters optimization iterations
		(2.1.1) Optimize alpha
		(2.1.2) Optimize nu
		(2.1.3) Optimize f
		if the likelihood wasn't changed during this loop --> parameters converged --> break
	(2.2) BBL
	if the likelihood wasn't changed during this loop --> parameters+bbl converged --> break
(3) return likelihood
*/

// @@@@ should use different parameter for maxIterations and maxBblIterations
MDOUBLE bestParamUSSRV::operator() (tree& et,
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights /* =NULL */,
					   const MDOUBLE AlphaUpperBound /* = 15 */, 
					   const MDOUBLE NuUpperBound /* = 15 */, 
					   const MDOUBLE FUpperBound /* = 1 */, 
					   const MDOUBLE epsilonParamOptimization /* = 0.01 */,
					   const MDOUBLE epsilonLikelihoodImprovment /* = 0.05 */,
					   const int maxIterations /* = 10 */)
{
	_bestL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;	

	bestAlphaFixedTreeSSRV alphaOptimization;
	bestNuFixedTreeSSRV nuOptimization;
	bestFFixedTreeSSRV fOptimization;
	
	int it, bblIt;
	int numberOfIterations(maxIterations);
	int numberOfParametersAndBblIterations(maxIterations);
	
	// if only one parameter is optimize (only Alpha or only Nu or only F) then we need only one iteration.
	// if we only do bbl, without any optimization of the parameters, then we don't need iterations at all.
	int countParameters2Optimize(0);
	if (_AlphaOptimizationFlag) countParameters2Optimize++;
	if (_NuOptimizationFlag) countParameters2Optimize++;
	if (_FOptimizationFlag) countParameters2Optimize++;

	if (countParameters2Optimize==0)
	{
		numberOfIterations=0;
		numberOfParametersAndBblIterations=1;
	}
	else if (countParameters2Optimize==1)
		numberOfIterations=1;
	
	if (_bblOptimizationFlag == false)
		numberOfParametersAndBblIterations = 1;
	
	_bestAlpha = model.getAlpha();
	_bestNu = model.getNu();
	_bestF = model.getF();

	bool changes(false);
	bool bblChanges(false);
	for (bblIt=0; bblIt < numberOfParametersAndBblIterations; ++bblIt)
	{
		bblChanges = false;

		// parameters optimizations (without bbl)
		// in each iteration : optimization of Alpha and then optimization of Nu, and then of F.
		for (it=0; it < numberOfIterations; ++it)
		{
			changes = false;	
			// Alpha optimization
			if (_AlphaOptimizationFlag)
			{
				newL = alphaOptimization(et,sc,baseSc,model,weights,AlphaUpperBound,epsilonParamOptimization);

				//the improvemnt in Likelihood is smaller than epsilon
				if (newL < _bestL)
				{				
					LOG(5,<<"likelihood went down in LS! (Alpha optimization)"<<endl<<"oldL = "<<_bestL<<" newL= "<<newL<<endl);
					//go back to previous alpha
					alphaOptimization.setAlpha(_bestAlpha,model);
					alphaOptimization.setBestL(_bestL); // @@@@ maybe this is unnecessary 
					//break;
				}
				else 
				{// update of likelihood and model.
					if (newL > _bestL+epsilonLikelihoodImprovment) 
					{
						changes = true;
						bblChanges = true;
					}
					_bestL = newL;
					_bestAlpha = alphaOptimization.getBestAlpha();
					LOG(5,<<"new L = " << _bestL<<"  new Alpha = " << _bestAlpha<<endl);		
				}
			}
		
			// Nu optimization
			if (_NuOptimizationFlag)
			{
				newL = nuOptimization(et,sc,baseSc,model,weights,NuUpperBound,epsilonParamOptimization);
			
				//the improvemnt in Likelihood is smaller than epsilon
				if (newL < _bestL)
				{
					LOG(5,<<"likelihood went down in LS! (Nu optimization)"<<endl<<"oldL = "<<_bestL<<" newL= "<<newL<<endl);
					//go back to previous Nu
					nuOptimization.setNu(_bestNu,model);
					nuOptimization.setBestL(_bestL); // @@@@ maybe this is unnecessary 
					//break;
				}
				else
				{// update of likelihood and model.
					if (newL > _bestL+epsilonLikelihoodImprovment) 
					{
						changes = true;
						bblChanges = true;
					}
					_bestL = newL;
					_bestNu = nuOptimization.getBestNu();
					LOG(5,<<"new L = " << _bestL<<"  new Nu = " << _bestNu<<endl);		
				}
			}

			// F optimization
			if (_FOptimizationFlag)
			{
				newL = fOptimization(et,sc,baseSc,model,weights,FUpperBound,epsilonParamOptimization);

				//the improvemnt in Likelihood is smaller than epsilon
				if (newL < _bestL)
				{
					LOG(5,<<"likelihood went down in LS! (F optimization)"<<endl<<"oldL = "<<_bestL<<" newL= "<<newL<<endl);
					//go back to previous F
					fOptimization.setF(_bestF,model);
					fOptimization.setBestL(_bestL); // @@@@ maybe this is unnecessary 
					//break;
				}
				else 
				{// update of likelihood and model.
					if (newL > _bestL+epsilonLikelihoodImprovment ) 
					{
						changes = true;
						bblChanges = true;
					}

					_bestL = newL;
					_bestF = fOptimization.getBestF();
					LOG(5,<<"new L = " << _bestL<<"  new F = " << _bestF<<endl);						}
			}
			if (changes == false)
			{
				LOG(5,<<"bestParamUSSRV parameters alpha,nu,f converged!"<<endl);
				break;
			}
		}

		if (it == numberOfIterations)
			LOG(5,<<"bestParamUSSRV parameters alpha, nu, f, did not converge after " << numberOfIterations << "iterations"<<endl);


		// BBL
		if (_bblOptimizationFlag == true)
		{
			// debug OZ
			LOG(5,<<"the like before = "<<likelihoodComputation2USSRV::getTreeLikelihoodAllPosAlphTheSame(et,sc,baseSc,model));
			// end of debug
			bblEM2USSRV bbl(et,sc,baseSc,model,NULL,maxIterations);
			newL = bbl.getTreeLikelihood();
			LOG(5,<<"current best L= "<<_bestL<<endl);
			LOG(5,<<"new L After BBL = " << newL<< " = "<<likelihoodComputation2USSRV::getTreeLikelihoodAllPosAlphTheSame(et,sc,baseSc,model)<<endl);
			if (newL > _bestL+epsilonLikelihoodImprovment)
				bblChanges = true;
			if (newL < _bestL){
				LOG(5,<<"likelihood went down in LS! (BBL)"<<endl<<"oldL = "<<_bestL);
				LOG(5,<<" newL= "<<newL<<endl) ;
			}
			else
				_bestL = newL;
		}

		if (bblChanges == false)
		{
			LOG(5,<<"bestParamUSSRV bbl and parameters converged!"<<endl);
			break;
		}
	}

	if (bblIt == numberOfParametersAndBblIterations)
		LOG(5,<<"bestParamUSSRV bbl and parameters alpha did not converge after " << numberOfParametersAndBblIterations << "iterations"<<endl);	

	return _bestL;
}

