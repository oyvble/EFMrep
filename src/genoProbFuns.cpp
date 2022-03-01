//helpfunction to get allele probability:
#include <vector> //vector storage

double prob_a(double Pa, double mm, double nn, double fst2) {
	return( (mm*(fst2) + (1-(fst2))*Pa)/(1+(nn-1)*(fst2)) ); 
}

double prob_relUnknown(int *outG1vec, double *Fvec, int Ugind, int Rgind, double *ibd2, double fst2, std::vector<double> maTypedvec2, double nTyped2) {
	//Ugind = genotype index of unknown indiviudal 
	//Rgind = genotype index of related indiviudal 		
	//ibd2 = ibd vector (NOC long)
	
	double genoSum; //used to sum the genotype probability tfor unknowns		
	int aindU = outG1vec[2*Ugind  ]; //get allele index of genotype g_1 (unknown)
	int bindU = outG1vec[2*Ugind+1]; //get allele index of genotype g_2 (unknown)
	int aindR = outG1vec[2*Rgind  ]; //get allele index of genotype g_1 (related)
	int bindR = outG1vec[2*Rgind+1]; //get allele index of genotype g_2 (related)
	
	bool Uhom = aindU==bindU; //boolean of whether unknown genotype is homozygote
	
	//First step: Calculate random match probability of unrelated situation:
	genoSum = prob_a(Fvec[aindU],maTypedvec2[aindU],nTyped2,fst2); //init with prob 1st allele  
	if(Uhom) { //if unknown is homozygote					
		genoSum *= prob_a(Fvec[aindU],maTypedvec2[aindU]+1,nTyped2+1,fst2); //calculate random match prob (always used) 					
	} else { //if unknown is heterozygote variant
		genoSum *= 2*prob_a(Fvec[bindU],maTypedvec2[bindU],nTyped2+1,fst2); //calculate prob 2st allele  (and scale with 2)
	}
		
	//Extension with kappa-coefficient: SEE FORMULAS IN TABLE A.3 in Book "A forensic practicioners guide...."
	if( Rgind != -1 ) { //if related is specified (not -1 index)
		genoSum *= ibd2[0]; //multiply with kappa0
		if( Ugind==Rgind ) { //if Unknown and Related are same genotype
			genoSum+=ibd2[2]; //sum with kappa2
			
			if(Uhom) { //if unknown is homozygote
				genoSum += prob_a(Fvec[aindU],maTypedvec2[aindU],nTyped2,fst2)*ibd2[1] ; //multiply with kappa1 
			} else { //if unknown is heterozygote variant
				genoSum += (prob_a(Fvec[aindU],maTypedvec2[aindU],nTyped2,fst2)+prob_a(Fvec[bindU],maTypedvec2[bindU],nTyped2,fst2))*ibd2[1]/2 ; //multiply with kappa1 						
			}	
			
		} else { //if not the same genotype we need to check overlap (a,b)~(c,d)
			bool A1eq = aindU==aindR; //check if a=c 
			bool A2eq = aindU==bindR; //check if a=d 
			bool B1eq = bindU==aindR; //check if b=c 
			bool B2eq = bindU==bindR; //check if b=d 

			if( A1eq || A2eq || B1eq || B2eq) { //if any overlap (one shared allele): THERE ARE 3 OUTCOME!!
				if(Uhom) { //if Unknown is homozygote we know what allele to use
					genoSum += prob_a(Fvec[aindU],maTypedvec2[aindU],nTyped2,fst2)*ibd2[1]/2; //multiply with kappa1 	
				} else { //if Unknown is heterozygote we need to found the non-overlapping allele 
					bool Rhom = aindR==bindR; //check if related genotype is homozygote
					int cind; //temporary variable to select non-overlapping allele
					if(A1eq || A2eq) { //Allele bindU is non-overlapping
						cind = bindU;
					} else {
						cind = aindU; //Allele aindU is non-overlapping
					}
					if(Rhom) { //should not divide by 2 if related geno is homozygote
						genoSum += prob_a(Fvec[cind],maTypedvec2[cind],nTyped2,fst2)*ibd2[1]; //multiply with kappa1 																			
					} else { //should divide by 2 if related geno is heterozygote
						genoSum += prob_a(Fvec[cind],maTypedvec2[cind],nTyped2,fst2)*ibd2[1]/2; //multiply with kappa1 																												
					}
				}
			}	//otherwise none shared allele (and we are done (kappa0 already included)												
		}
	}
	return( genoSum ); //return genotype prob
}
