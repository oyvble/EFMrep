//This script contains functions for calculating the likelihoods function for replicates.
//AUTHOR: Oyvind Bleka, July 2021
/*ABOUT:

- No structuring of data is needed (already done in prepareC)
- All markers are looped outside the inner large-sum loop.
- Possible to traverse certain genotype combinations (tremendous speedup for 3 and more contributors)

*/


#include <vector> //vector storage
#include <cmath> //includes lgamma
#include <thread> //used to obtain number of logical processes
#include <Rmath.h> //includes pgamma
#include "genoProbFuns.h" //helpfunction for genotype probs
#ifdef _OPENMP
#include <omp.h> //parallelization
#endif

using namespace std;

/*MAIN IMPORTS DATA AN STRUCTURES THE DATA FOR MARKER BASED FUNCTIONS*/
extern "C" {
		

//The P(E|theta)= P(E|g,theta)*P(g) expression is returned
void loglikGamma_allcomb2( double *loglikVEC, int *nJointCombs, int *NOC, int *NOK, 
	double *mixProp, double *PHexp, double *PHvar, double *DEG, double *stutt,  //stutt is stutter proportion vector (1 element for each kind)
	double *AT, double *fst, double *dropinProb, double *dropinWeight, 
   	int *nReps, int *nMarkers, int *nRepMarkers, int *nAllelesVEC, int *nAlleles2VEC, int *startIndMarker_nAlleles, int *startIndMarker_nAllelesReps,	
	double *peaksLong, double *freqsLong, double *nTypedLong, double *maTypedLong, double *basepairLong,  int *repIDvec, int *startIndMarker_nRepMarkers,
	int *nGenos, int *outG1allele, int *outG1contr, int *startIndMarker_outG1allele, int *startIndMarker_outG1contr,  int *startIndMarker_nJointCombs,
	int *nStutters, int *stuttFromIndVEC, int *stuttToIndVEC, int *stuttParamIndVEC, int *startIndMarker_nStutters, 	
	int *knownGind, int *maxThreads, int *relGind, double *ibdLong) {	
	
	#ifdef _OPENMP
	int numThreads = thread::hardware_concurrency();
	int useThreads = min(numThreads,*maxThreads);
	if(*maxThreads==0) useThreads=numThreads; //use all available threads
	omp_set_num_threads(useThreads);  //set number of threads to use 
	#endif
	
	int aa, kk, rr; //indices for alleles, contributors and replicates
	//Prepare transformation of parameters (per replicate) to save computation:
	vector<double> PHvarSq(*nReps,0.0);
	vector<double> shape0(*nReps,0.0);
	vector<double> scale0(*nReps,0.0);
	vector<double> const1(*nReps,0.0);
	vector<double> const2(*nReps,0.0);
	for(rr=0; rr< *nReps; rr++) { //Traverse each replicate vars (param potential different for each)
		PHvarSq[rr] = PHvar[rr]*PHvar[rr]; //Square param
		shape0[rr] = 1/PHvarSq[rr]; // //theta_Am[locind]/theta_omegasq; //shape for 'full het allele'    
		scale0[rr] = PHexp[rr]*PHvarSq[rr]; //obtain scale param 
		const1[rr] = 1/scale0[rr]; //constant 1 (used in pgamma,dgamma)
		const2[rr] = log(scale0[rr]); //constant 2 (used in dgamma)		
	}
	
	//Calculating over all markers
	const double smalltol = 1.0e-30; //a tiny number > 0 (avoiding zero roundoff errors)
	for(int locind=0; locind< *nMarkers; locind++) {	//for each marker:
	
		//OBTAIN CONSTANTS FOR SPECIFIC MARKER:
	
		//default settings
		double fst0 = fst[locind]; //theta-correction param	(marker specific)
				
		//Prepare dimensions and data vectors
		int NOK0 = NOK[locind]; //obtain number of contributors (may be different for different markers)
		int NOU = *NOC - NOK0; //number of unknowns (may be different for markers)
		int nReps0 = nRepMarkers[locind]; //number of replicates for specific marker (may be different for different markers)
		
		//Prepare dimen
		int numGenos1p = nGenos[locind]; //Number of genotypes (1 contributor)
		int nAlleles = nAllelesVEC[locind]; //number of alleles
		int nAlleles2 = nAlleles2VEC[locind]; //number of alleles (including potential stutters)
		int nPS = nAlleles2 - nAlleles; //number of potential stutters
													
		int SI_nAlleles0 = startIndMarker_nAlleles[locind]; //index start number of alleles
		int SI_nAllelesReps0 = startIndMarker_nAllelesReps[locind]; //index start number of alleles (taking into account Reps)
		int SI_nRepMarkers0 = startIndMarker_nRepMarkers[locind]; //index start number of alleles (taking into account Reps)
		
		
		int SI_nJointCombs0 = startIndMarker_nJointCombs[locind];
		int SI_nStutters0 = startIndMarker_nStutters[locind];
		int SI_outG1allele0 = startIndMarker_outG1allele[locind];
		int SI_outG1contr0 = startIndMarker_outG1contr[locind];
				
		//Prepare vector for known contributors (and also unknown): Need to know positions!
		vector<int> GindKnown(NOK0,0); //genotype index in vector 
		vector<int> kindKnown(NOK0,0); //contributor index in vector
		vector<int> kindUnknown(NOU,0); //contributor index in vector		
		vector<int> kindRel(*NOC,-1); //genotype index in vector, one index for each contributors  (-1 means no related)
	
		int cc = 0; //counters for known
		int jj = 0; //counter for unknowns
		int knownGind0; //obtain genotype index of known 
		for(kk=0; kk< *NOC; kk++) { //for each contributors:
			kindRel[kk] = relGind[ (*NOC)*locind + kk ]; //copy genotype index 		
			knownGind0 = knownGind[ (*NOC)*locind + kk ]; //copy genotype index  (KNOWN)			
			if(knownGind0>=0) { //If contributor is known (genotype given)
				GindKnown[cc] = knownGind0; //copy genotype index
				kindKnown[cc] = kk; //insert contributor index
				cc++; //update counter for knowns
			} else { //if contributor is unknown (genotype not given)
				kindUnknown[jj] = kk; //insert contributor index
				jj++; //update counter for unknowns
			}
		}
				
		//Prepare variables (before sum-iterations): Save computations
		double nTyped = nTypedLong[locind]; //copy total number of typed alleles
		vector<double> maTypedvec(nAlleles,0.0); //copy number of typed alleles (each type)
		vector<double> shapevK(nAlleles2*nReps0, 0.0); //Precalculation for known contributors (vectorized over all replicates of particular marker)
		vector<double> shapev0(nAlleles2*nReps0, 0.0); //degrad scaling of shape (init vector) (vectorized over all replicates of particular marker)	
		int aaind; //vectorization index (see info below)
		int repIDmarker; //Used to indicate correct replicate (some markers may have fewer replicates!!!)
		int mixPropContrInd; //used to indicate correct mixture proportion vector (of particular replicate)
		for (aa = 0; aa < nAlleles; aa++) { //traverse each observed alleles (indices), also the Q-allele. Potential not necessary!
			maTypedvec[aa] = maTypedLong[ SI_nAlleles0 + aa ]; //copy previously typed alleles
			for (rr = 0; rr < nReps0; rr++) { //traverse each replicate	
				aaind = aa*nReps0 + rr; //get index vectorized vector which includes (allele,rep) indices as [(1,1), (1,2), (1,3), (2,1), (2,2) etc]
				repIDmarker = repIDvec[ SI_nRepMarkers0 + rr ]; //obtain correct repID for marker
				shapev0[aaind] = exp( log( shape0[repIDmarker] ) + basepairLong[ SI_nAllelesReps0 + aaind ]*log( DEG[repIDmarker] ) ); //scaling with degradation model (assumed already scaled)
			
				//Sum up contribution for each alleles (Taking into account mix proportions): ONLY CALCULATED FOR KNOWN CONTRIBUTORS INITIALLY
				for (kk = 0; kk < NOK0; kk++) { //for each known contributors 	
					mixPropContrInd = kindKnown[kk] + repIDmarker*(*NOC); //obtain correct index of mixture proportion (vectorized across replicates)
					shapevK[aaind] += outG1contr[ SI_outG1contr0 + nAlleles*GindKnown[kk] + aa] * mixProp[ mixPropContrInd ]; //contr from contr k to allele aa
				} //end for each known contributor
			} //end for each reps
		}
		
		#pragma omp parallel for //shared(shapev0, dropin0,constY0) //perform parallel on summing up bigsum
		for (int iter = 0; iter < nJointCombs[locind]; iter++) { //for each combined/joint genotype outcome
		
			//int iter = combUseVEC[ SI_nJointCombs0 + iter2]; //obtain correct iteration index to consider

			//Following variables PHexpst be declared for each iteration (non-shared):
			int a,k,r; //used to traverse alleles(a), contributors(k) and replicates(r)
			int aind; //used as index for alleles in genotypes
			int repID; //Used to indicate correct replicate (some markers may have fewer replicates!!!)
			int mixPropInd; //used to indicate correct mixture proportion vector (of particular replicate)
			vector<int> jointGind(NOU, 0); //This will be the contribution index for the unknown inds (directly corresponds to indices of Gmarg)
			vector<double> shapev = shapevK; //make a copy of existing shapevector
			vector<double> maTypedvec2 = maTypedvec; //Creating copy of counter (per allele)
			double nTyped2 = nTyped; //creating copy of total counter

			//jointGind = digits(iter, base = numGenos1p, pad = NOU); #Equivalent operaion
			int startInd_alleles; //init for allele index
			double genoProd = 1.0; //calculating the genotype probability of the unknowns
			int modrest = iter; //used to keep remained after modulo (init as iter number)
			bool inserted = false; //boolean of whether all digits are inserted (Rest is zero padded)
			int allele_contr; //allele contribution (decided by genotype combinations)
			for (k = 0; k < NOU; k++) { //for each unknown contributors (summing up wrt both contr (outG1contr) and mx (mixProp)): Need each contr to derive shapev
				if (!inserted) { //if not all digits inserted
					if ( k>0 ) {
						modrest = int((modrest - jointGind[k - 1]) / numGenos1p); //extract remaining, divide to get to next digit (necessary with int converion?)
					}
					jointGind[k] = modrest % numGenos1p; //INSERT NUMBER: convert number to "numGenos1p" basis 	
					if (modrest < numGenos1p) { //check if rest is smaller than base (only run if not inserted)
						inserted = true;  //then all digits are inserted
					}	
				} //else { jointGind[k] = 0; //zero pad	}
		
				//Sum up contribution for each allele:rep (Taking into account mix proportions for each replicate)
				for (a = 0; a < nAlleles; a++) { //traverse each alleles (indices)				
					allele_contr = outG1contr[  SI_outG1contr0 + nAlleles*jointGind[k] + a]; //obtain allele contribution
					if( allele_contr==0 ) continue; //skip allele if no contribution from genotypes (SPEEDUP???)
					
					for (r = 0; r < nReps0; r++) { //traverse each replicate
						aind = a*nReps0 + r; //get index vectorized vector which includes (allele,rep) indices as [(1,1), (1,2), (1,3), (2,1), (2,2) etc]
						repID = repIDvec[ SI_nRepMarkers0 + r ]; //obtain correct repID for marker
		
						mixPropInd = kindUnknown[k] + repID*(*NOC); //obtain correct index of mixture proportion (vectorized across replicates)
						shapev[aind] +=  allele_contr* mixProp[ mixPropInd ]; //contr from contr k to allele a. NOTICE THE k+*NOK0 shift!
					}
				}
		
				//CALCULATE GENOTYPE PROBS OF UNKNOWNS (MAY BE RELATED:  kindRel[kindUnknown[k]]!= -1)
				//Sending pointer where marker indices are starting for outG1 and freqLong vector
				genoProd *= prob_relUnknown( outG1allele + SI_outG1allele0, freqsLong + SI_nAlleles0,  jointGind[k], kindRel[kindUnknown[k]], ibdLong + 3*kindUnknown[k], fst0, maTypedvec2, nTyped2 );	 //Note: scale with 3 because ibdLong is a '3-long vector' per contributor						
					
				//LAST: UPDATE COUNTERS FOR ALLELES (a and b)
				startInd_alleles = SI_outG1allele0 + 2*jointGind[k];
				maTypedvec2[ outG1allele[startInd_alleles  ] ] +=1; //update allele count for particular allele (1)
				maTypedvec2[ outG1allele[startInd_alleles+1] ] +=1; //update allele count for particular allele (2)					
				nTyped2 += 2; //update total count
			} //end for each unknown contributor

			//////////////////////////////////////////////
			//Calculating the inner sum-part begins here//
			//////////////////////////////////////////////

			//Scaling shape parameter with degrad model (could be done elsewhere?  done last in prev operatoin)
			for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (indices), having PH>0
				for (r = 0; r < nReps0; r++) { //traverse each replicate
					aind = a*nReps0 + r; //get index vectorized vector which includes (allele,rep) indices as [(1,1), (1,2), (1,3), (2,1), (2,2) etc]
					shapev[aind] *= shapev0[aind]; //scaling with degradation model (assumed already scaled)					
				}
			}
			
			//Scaling shape parameter with stutter model:			
			vector<double> shapev2 = shapev;  //make a copy of existing shapevector
			int stuttind,stuttind2,stuttFromInd,stuttToInd,stuttParamInd; //init vars
			for (r = 0; r < nReps0; r++) { //traverse each replicate (necessary if the stutter params are different for each replicate)
				repID = repIDvec[ SI_nRepMarkers0 + r ]; //obtain correct repID for marker
				for(stuttind=0; stuttind < nStutters[locind]; stuttind++) { 
					stuttind2 = SI_nStutters0 + stuttind; //obtain index of stutter
					stuttFromInd = stuttFromIndVEC[stuttind2]*nReps0 + r;  //NB: CAREFUL WITH INDEX
					
					if( shapev[stuttFromInd]>smalltol) { //ONLY NECESSARY TO PROVIDE MODIFICATION IF shapeval>0
						stuttToInd = stuttToIndVEC[stuttind2]*nReps0 + r;  //NB: CAREFUL WITH INDEX
						stuttParamInd = stuttParamIndVEC[stuttind2] + 2*repID;  //Obtain parameter index (NOTE: ALWAYS using 2 elements per replicate). NB: CAREFUL WITH INDEX					
						shapev2[stuttToInd] += stutt[stuttParamInd] * shapev[stuttFromInd]; // #OBTAINED stutters
						shapev2[stuttFromInd] -= stutt[stuttParamInd] * shapev[stuttFromInd]; // #SUBTRACTED stutters
					}
				}
			}
			
			//Summing up likelihood contribution of each alleles:reps:
			vector<int> nDropin(nReps0,0); //count number of drop-in for each replicate
			double logevidProb = 0.0; //evidence probability (weight)
			
			double peak; //obtaining observeed peak height 
			double AT0; //analytical threshold (depends on both marker and replicate)
			for(r = 0; r < nReps0; r++) { //traverse each replicates (observed alleles indicated by PH)
				repID = repIDvec[ SI_nRepMarkers0 + r ]; //obtain correct repID for marker
				AT0 = AT[SI_nRepMarkers0 + r]; //obtain AT to use (follows same as nReps per marker)
				
				for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)				
					aind = a*nReps0 + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
					peak = peaksLong[SI_nAllelesReps0 + aind]; //Notice the cumulative locus shift (Long vector)
					
					//IF NOT DROPOUT
					if( peak >= AT0 ){ //if PH>0 (this is same as PH>=AT since threshold has already been applied)
						if(shapev2[aind] > smalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)
							logevidProb += -lgamma(shapev2[aind]) - shapev2[aind]*const2[repID] + (shapev2[aind]-1)*log(peak) - peak*const1[repID]; //gamma/lgamma is found in cmath							
						} else { //If contribution and PH>0 	//DROPIN SET (C)							
							logevidProb += dropinWeight[ SI_nAllelesReps0 + aind ]; //likelihood for observed dropin (both PH and lambda are included)							
							nDropin[r] += 1; //count dropin for particular replicate															   
						}							
						
					//OTHERWISE IT IS DROPOUT
					} else if(shapev2[aind] > smalltol) { //IF PH missing  and contribution(it't a dropout //CONTRIBUTION SET (B)
						logevidProb += pgamma(AT0, shapev2[aind],scale0[repID], 1, 1); //Add log(dropout-probability)
						//log( gamma_p(shapev2[aind],AT0*const1[repID]) );  
					}
					
				} //end for each replicates (r)
			} //end for each observed alleles
			
			//Weight potential drop-outs (stutters)
			 if(nPS>0) {
				for(r = 0; r < nReps0; r++) { //traverse each replicates (observed alleles indicated by PH)
					repID = repIDvec[ SI_nRepMarkers0 + r ]; //obtain correct repID for marker
					AT0 = AT[SI_nRepMarkers0 + r]; //obtain AT to use (follows same as nReps per marker)
					for(a=nAlleles; a<nAlleles2; a++) { //traverse remaining alleles
					  aind = a*nReps0 + r; //get index vectorized vector which includes (allele,rep) indices as [(1,1), (1,2), (1,3), (2,1), (2,2) etc]
					  if(shapev2[aind]>smalltol) { //IF DROPOUT						
						logevidProb += pgamma(AT0, shapev2[aind],scale0[repID], 1, 1); //Add log(dropout-probability)
					  }
					}
				}
			 }
						
			//Liklihood OF NUMBER OF NOISE/DROP-IN: Traverse each replicates (observed alleles indicated by PH)
			double dropinProb0; //dropin prob
			for(r = 0; r<nReps0; r++) {
				dropinProb0 = dropinProb[SI_nRepMarkers0 + r]; //obtain dropin prob (follows same as nReps per marker)
				if(nDropin[r]==0) { //in case of no droping
					logevidProb += log(1-dropinProb0);
				} else {
					logevidProb += nDropin[r]*log(dropinProb0);
				}
			}
			
			//FINAL INSERTION OF VALUES:
  		   loglikVEC[SI_nJointCombs0 + iter] = logevidProb + log(genoProd); //calculate P(E|gj)P(gj)			   		   
		} //end for each combination iterations (bigsum)
		//logLik[0] += log(bigsum); //calculate logLik by adding log(innerSUM)	
	} //end for each marker
} //end main function

} //end external