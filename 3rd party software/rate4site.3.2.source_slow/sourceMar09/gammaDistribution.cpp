// $Id: gammaDistribution.cpp 391 2005-06-12 12:44:57Z ninio $

 #include "definitions.h"
#include "gammaDistribution.h"
#include "gammaUtilities.h"
#include "errorMsg.h"
#include "logFile.h"
#include <cmath>


gammaDistribution::gammaDistribution(MDOUBLE alpha,int in_number_of_categories) :distribution(){
	_globalRate=1.0;
	setGammaParameters(in_number_of_categories,alpha);
}

//copy constructor
gammaDistribution::gammaDistribution(const gammaDistribution& other) : 
	_alpha(other._alpha), 
	_bonderi(other._bonderi),
	_rates(other._rates),
	_ratesProb(other._ratesProb),
	_globalRate(other._globalRate) {
}

void gammaDistribution::setAlpha(MDOUBLE in_alpha) {
	if (in_alpha == _alpha) return;
	setGammaParameters( categories(), in_alpha);
}

void gammaDistribution::change_number_of_categories(int in_number_of_categories) {
	setGammaParameters( in_number_of_categories, _alpha);
}

//this function builds the gamma distribution
void gammaDistribution::setGammaParameters(int in_number_of_categories, MDOUBLE in_alpha) {
	if (in_alpha < MINIMUM_ALPHA_PARAM)	in_alpha = MINIMUM_ALPHA_PARAM;// when alpha is very small there are underflow problems
	_alpha = in_alpha;
	_rates.clear();
	_rates.resize(in_number_of_categories);
	_ratesProb.erase(_ratesProb.begin(),_ratesProb.end());
	_ratesProb.resize(in_number_of_categories,1.0/in_number_of_categories);
	_bonderi.clear();
	_bonderi.resize(in_number_of_categories+1);
	if (in_number_of_categories==1) {
		_rates[0] = 1.0;
		return;
	}
	if (categories()>1) {	
		// fills _rates with the mean values for each category,and creates the boundaries 
		// for the categories
		fill_mean(); 
		return ;
	}
	
}

int gammaDistribution::fill_mean() {
	fill_bonderi();
	int i;
	//for (i=0; i<=categories(); ++i) cout<<endl<<bonderi[i];
	//LOG(5,<<"\n====== the r categories are =====\n");
	for (i=0; i<categories(); ++i) {
		_rates[i]=the_avarage_r_in_category_between_a_and_b(_bonderi[i],_bonderi[i+1],_alpha,_alpha,categories());
		//LOG(5,<<meanG[i]<<endl);
	}
	//LOG(5,<<endl<<alpha<<endl);
	return 0;
}

int gammaDistribution::fill_bonderi() {
	int i;
	for (i=1; i<categories(); ++i) {
		_bonderi[i]=search_for_z_in_dis_with_any_beta(_alpha,_alpha,(MDOUBLE)i/categories());
	}
	_bonderi[0]=0;
	_bonderi[i]=VERYBIG/10000.0;// this is becuase we multiply bondei[i] by alpha or beta, and 
	// by this manipulation we avoid overflows...
	
	return 0;
}



//gets cumulative probability till a certain point
const MDOUBLE gammaDistribution::getCumulativeProb(const MDOUBLE x) const
{
	//since r~gamma(alpha, beta) then beta*r~ gamma(alpha,1)=gammp
	//here we assume alpha=beta
	return gammp(_alpha, x*_alpha);
}







