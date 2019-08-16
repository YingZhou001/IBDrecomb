#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#include "head.h"
#include "tools.h"
#include "end.h"


void meg_Rec(int **sourceIBD, int from, int to, double bin, struct map megmap[], struct mapPar *par)
{
    clock_t checkpt, begin;
    int binnum;
    double newbin;
    int i, j, k;
    binnum  = megL;
    newbin = (double)(to - from) / binnum;//renew bins to make each bin has the same width

    fprintf(stderr, "\n\nLarge scale estimation:\n");
    fprintf(stderr, "\tScale = %.1lf, ", newbin);
    fprintf(stderr, "Estimation ranges = (%.1lf, %.1lf)\n", newbin*0 + POSdrift, newbin*binnum + POSdrift);

    double *coverage; // IBD coverage for each bin
    double *X, *X0, *X1; //
    coverage=(double*)malloc(binnum*sizeof(double));
    X=(double*)malloc(binnum*sizeof(double));
    X0=(double*)malloc(binnum*sizeof(double));
    X1=(double*)malloc(binnum*sizeof(double));
    for(i=0;i<binnum;i++){X[i]=0;X0[i]=0;}

    int **end; //IBD ends with update genetic position and length
    end=(int **)malloc(IBDNUMtot*sizeof(int*));
    for(i=0;i<IBDNUMtot;i++){end[i]=(int *)malloc(5*sizeof(int));}


    par->L = binnum;
    par->bin = newbin;
    par->minibd = minIBD;
    constoutmap(megmap, *par); // initiate recombination map with constant rate
    //print_map(megmap, par);
    if(mute == 0)fprintf(stderr, "\n\nIteration = 0\n");
    begin = clock();
    par->ibdcut = update_IBDCUT(megmap, *par);

    IBDCUT = par->ibdcut;//fprintf(stderr, "IBDCUT: %lf\n", IBDCUT);
    ////read and sort ibd
    par->ibdnum = read_ibd(sourceIBD, end, megmap, *par);
    sort_ibd(end, STEP, par->ibdnum);
    cal_end_region_legnth(end, par->ibdnum);
    update_TRIM(megmap, *par);
    maxIBDnum = check_max_IBDnum_megbin(sourceIBD, *par);
    par->mindepth = cal_bin_depth(end, coverage, megmap, *par);

    end_counts_in_bin(X, end, sourceIBD, *par);
    for(j=0;j<par->L;j++)X0[j]=X[j];
    //for(i=0;i<par->L;i++)fprintf(stderr,"%lf\n",X[i]);
    update_map(X, megmap, *par);
    checkpt=clock();
    if(mute ==0 )fprintf(stderr, "\tTotal Time for iter 0 = %lfs\n", (checkpt-begin)*1.0/CLOCKS_PER_SEC);
    //print_map(megmap, *par);


    for(i=1;i<ITER;i++){
	sumX=0;
	begin = clock();
	//fprintf(stderr, "\n\nIteration = %d\n", i);

	fprintf(stderr, "  Iter %d:", i+1);

	par->ibdcut = update_IBDCUT(megmap, *par);
	par->ibdnum = read_ibd(sourceIBD, end, megmap, *par);
	sort_ibd(end, STEP, par->ibdnum);
	cal_end_region_legnth(end, par->ibdnum);//out: ENDlength 
	update_TRIM(megmap, *par);//in: ENDlength; out: end-region index, LENDindex & RENDindex
	par->mindepth = cal_bin_depth(end, coverage, megmap, *par);
	end_counts_in_bin(X, end, sourceIBD, *par);
	adjust_end_counts(X, sourceIBD, megmap, *par);
	if(i == 1){
	    for(j=0;j<par->L;j++){X1[j]=X[j];sumX+=X1[j];}
	}
	else {
	    for(j=0;j<par->L;j++){X1[j]=X[j]+X0[j];sumX+=X1[j]/2;}
	}
	update_map(X1, megmap, *par);
	for(j=0;j<par->L;j++)X0[j]=X[j];

	fprintf(stderr, "\tIBD Ends num = %.0lf\n", sumX);

	checkpt=clock();
	if(mute ==0 )fprintf(stderr, "\tTotal Time for iter %d = %lfs\n", i, (checkpt-begin)*1.0/CLOCKS_PER_SEC);
    }
    sumX=0;
    j=0;while(j<par->L){if(j<LENDindex||j>RENDindex)sumX+=X[j]/2;j++;}
    //print_map(megmap, *par);
    free_array(coverage);
    free_array(X);
    free_array(X0);
    free_array(X1);
    free_matrix_d(end, IBDNUMtot);

    if(mute == 0)fprintf(stderr, "Done!\n");
}


void fineRec(int **sourceIBD, int I, double G, struct map mOld[megL][fineL], struct mapPar pOld[megL], struct map m[megL][fineL], struct mapPar p[megL])
{
    // range of meg bins [from, to) and the genetic length G
    clock_t checkpt, begin;
    int i, j, k;
    int binnum = fineL;
    double bin = mOld[I][1].p - mOld[I][0].p;
    //fprintf(stderr, "%d %lf\n", binnum, bin);
    if(mute == 0)fprintf(stderr, "\tBin %d: [%.1lf, %.1lf)\n", I, mOld[I][0].p + POSdrift, mOld[I][fineL-1].p + POSdrift + bin);

    double *coverage; // IBD coverage for each bin
    double *X; //
    coverage=(double*)malloc(binnum*sizeof(double));
    X=(double*)malloc(binnum*sizeof(double));
    for(i=0;i<binnum;i++){X[i]=0;}

    int **end; //IBD ends with update genetic position and length

    end=(int **)malloc(maxIBDnum*sizeof(int*)); //###can be improved later

    for(i=0;i<maxIBDnum;i++){end[i]=(int *)malloc(5*sizeof(int));}

    begin = clock();
    p[I].ibdnum = read_local_ibd(sourceIBD, end, I, mOld, pOld); // sourceIBD to end
    //fprintf(stderr, "%d\n", p[I].ibdnum);
    sort_ibd(end, STEP, p[I].ibdnum);
    p[I].mindepth = cal_fine_bin_depth(end, coverage, mOld[I], p[I]);
    end_counts_in_fine_bin(X, end, sourceIBD, m[I], p[I]); //JOE!!!!
    //for(i=0;i<binnum;i++){sumX+=X[i];}
    for(i=0;i<binnum;i++){sumX+=X[i];
	//fprintf(stderr,"%d-%d %.0lf\n",I,i, X[i]);
    }
    //fprintf(stderr, "\n%d %.0lf\n",I, sumX);
    update_fine_map(X, m[I], p[I], G);
    //for(i=0; i < fineL; i++){fprintf(stderr, "%lf %lf %lf %lf\n", m[I][i].p, m[I][i].r,m[I][i].g, G);}
    free_array(coverage);
    free_array(X);
    free_matrix_d(end, maxIBDnum);

}


void fineSideRec(int **sourceIBD, int I, double G, struct map mOld[megL][fineL], struct mapPar pOld[megL], struct map m[megL][fineL], struct mapPar p[megL])
{
    clock_t checkpt, begin;
    int i, j, k;
    int binnum = fineL;
    double bin = mOld[I][1].p - mOld[I][0].p;
    if(mute == 0)fprintf(stderr, "\tBin %d: [%.1lf, %.1lf)\n", I, mOld[I][0].p + POSdrift, mOld[I][fineL-1].p + POSdrift + bin);

    double *coverage; // IBD coverage for each bin
    double *X; //
    coverage=(double*)malloc(binnum*sizeof(double));
    X=(double*)malloc(binnum*sizeof(double));
    for(i=0;i<binnum;i++){X[i]=0;}

    int **end; //IBD ends with update genetic position and length

    end=(int **)malloc(maxIBDnum*sizeof(int*)); //###can be improved later

    for(i=0;i<maxIBDnum;i++){end[i]=(int *)malloc(5*sizeof(int));}

    begin = clock();
    p[I].ibdnum = read_local_ibd(sourceIBD, end, I, mOld, pOld); // sourceIBD to end
    p[I].mindepth = cal_side_fine_bin_depth(I, end, coverage, mOld[I], p[I]);
    mute=0;
    end_counts_in_side_fine_bin(I, X, end, m, p);
    mute=1;
    for(i=0;i<binnum;i++){sumX+=X[i];
	//fprintf(stderr,"%.0lf,", X[i]);
    }
    //fprintf(stderr, "\n%d %.0lf\n",I, sumX);
    update_fine_map(X, m[I], p[I], G);
    free_array(coverage);
    free_array(X);
    free_matrix_d(end, maxIBDnum);

}

