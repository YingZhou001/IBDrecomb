#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#include "head.h"
#include "check.h"
#include "read.h"
#include "estimateRec.h"

int main(int arg, char* argv[])
{
    srand((unsigned int)time(NULL));
    clock_t checkpt, begin0, begin;
    begin0 = clock();
    int i, j, k, tag=0;
    // input checking
    begin = clock(); tag=checkinput(arg, argv);if(tag==-1)return tag;

    IBDNUMtot=check_ibd(ibdinp);
    fprintf(stderr, "  Estimation range = (%d, %d)\n", POSmin + TRIM, POSmax - TRIM);
    POSdrift = POSmin + TRIM;
    megL = (int)floor((POSmax - POSmin - 2*TRIM - 1)*1.0 / BIN) + 1;

    checkpt=clock();
    fprintf(stderr,"  Done! Time = %lfs\n", (checkpt-begin)*1.0/CLOCKS_PER_SEC);

    //load data
    begin = clock();
    fprintf(stderr, "\nData loading\n");
    int **sourceIBD;
    sourceIBD = (int **)malloc(IBDNUMtot*sizeof(int*));
    for(i=0;i<IBDNUMtot;i++){sourceIBD[i] = (int *)malloc(2*sizeof(int));}

    read_source_ibd(ibdinp, sourceIBD);
    checkpt = clock();
    fprintf(stderr,"  Done! Time = %lfs\n", (checkpt-begin)*1.0/CLOCKS_PER_SEC);

    // Estimation on large scale
    begin = clock();
    struct map outMegMap[megL]; // store map in large scale
    struct mapPar outMegMapPar; //record the parameters for constructing outMegMap
    meg_Rec(sourceIBD, POSmin + TRIM, POSmax - TRIM, BIN, outMegMap, &outMegMapPar);
    if(fBIN < 0 || fBIN >= BIN ){
	//skip to output if fBIN is smaller than BIN
	if(ohapmap == 1)print_map(outMegMap, outMegMapPar);
        else print_plinkmap(outMegMap, outMegMapPar);
	checkpt = clock();
	fprintf(stderr,"Succuss! Time = %lfs\n", (checkpt-begin0)*1.0/CLOCKS_PER_SEC);
    }
    else {
	// Estimation on fine scale
	//begin = clock();
	double newbin;
	fineL = (int)floor(( outMegMapPar.bin -0.1 ) / fBIN) + 1;
	newbin = outMegMapPar.bin / fineL;
	fprintf(stderr, "\n\nFine scale estimation:\n\tScale = %.1lf\n", newbin);
	struct map outfineMapOld[ megL ][ fineL ]; // store map in fine scale
	struct map outfineMapNew[ megL ][ fineL ]; // output map in fine scale
	struct map outfineMap[ megL ][ fineL ]; 
	struct mapPar outfineMapParOld[ megL ]; //record the parameters for constructing outFineMap
	struct mapPar outfineMapParNew[ megL ]; 
	struct mapPar outfineMapPar[megL]; 


	init_fine_map(outfineMap, outfineMapPar, outMegMap, outMegMapPar);//initiate fine maps from mega maps
	cp_map(outfineMap, outfineMapPar, outfineMapOld, outfineMapParOld);

	double G, sumX0, sumX1;
	for( i = 1; i <= ITER; i++){
	    //for( j = 1; j <= 1; j++){
	    fprintf(stderr, "\n  Iter %d:\t", i);
	    sumX=0;
	    for( j = 0; j < megL; j++){
		G = outMegMap[ j ].r * outMegMapPar.bin / 1e6; //centiMorgan
		if(j>=LENDindex && j <=RENDindex)fineRec(sourceIBD, j, G, outfineMapOld, outfineMapParOld, outfineMap, outfineMapPar);
		if(j<LENDindex || j >RENDindex)fineSideRec(sourceIBD, j, G, outfineMapOld, outfineMapParOld, outfineMap, outfineMapPar);
	    }
	    if( i == 1){
		cp_map(outfineMap, outfineMapPar, outfineMapNew, outfineMapParNew);
		sumX1 = sumX;
	    }
	    else {
		cp_map(outfineMap, outfineMapPar, outfineMapNew, outfineMapParNew);
		add_map(outfineMapNew, outfineMapParNew, outfineMapOld, outfineMapParOld);
		sumX1 = (sumX + sumX0)/2;
	    }
	    cp_map(outfineMap, outfineMapPar, outfineMapOld, outfineMapParOld);
	    sumX0 = sumX;
	    fprintf(stderr, "IBD Ends num = %.0lf", sumX1);
	}
	if(ohapmap == 1)print_fine_map(outfineMapNew, outfineMapParNew);
	else print_fine_plinkmap(outfineMapNew, outfineMapParNew);
	free_matrix_d(sourceIBD, IBDNUMtot);
	checkpt = clock();
	fprintf(stderr,"\n\nDone! Total Time = %lfs\n", (checkpt-begin0)*1.0/CLOCKS_PER_SEC);
	}
	return 0;
    }
