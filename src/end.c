#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <limits.h>

#include "head.h"
#include "end.h"
#include "tools.h"

int adjust_end_counts(double *X, int **sourceIBD, struct map m[], struct mapPar p)
{
	clock_t checkpt, begin = clock();
	int i, j, k, c;
	int endl_i, endr_i;
	int mgl_i, mgr_i;
	int start_i, end_i;
	int L;
	double fold;
	fold=FOLD;
	endl_i = LENDindex;
	endr_i = RENDindex;
	mgl_i = (int)floor(LENDindex * (fold + 1)) + 1;
	mgr_i = p.L - (int)floor((p.L - RENDindex) * (fold + 1)) - 1;
	while(mgl_i > RENDindex || mgr_i < LENDindex){
		fold=fold/2;
		mgl_i = (int)floor(LENDindex * (fold + 1)) + 1;
		mgr_i = p.L - (int)floor((p.L - RENDindex) * (fold + 1)) - 1;
		fprintf(stderr, "Warning: temperally adjust fold to %lf in this round\n", fold);
		if(fold < 0.1){
			fprintf(stderr, "Warning: Adjunct region is too small, bias may be brought in if you see this in the last iteration\n");
			break;}
	}
	//fprintf(stderr, "%d %d %d %d\n", endl_i, mgl_i, mgr_i, endr_i);
	int E_r, E_l;
	double l, ibdmax, sumx, sumr;
	double *onesideends, *onesidecoverage, lmincoverage, rmincoverage;
	onesideends=(double *)malloc(p.L*sizeof(double));
	onesidecoverage=(double *)malloc(p.L*sizeof(double));

	//adjust the left ends

	int **sourceibd, sourceibd_size_l;
	sourceibd=(int **)malloc(IBDNUMtot*sizeof(int*));
	for(i=0;i<IBDNUMtot;i++){sourceibd[i]=(int *)malloc(4*sizeof(int));}
	sourceibd_size_l=0;
	for(i=0;i<IBDNUMtot;i++){
		k=sourceibd_size_l;
		E_l=sourceIBD[i][0];E_r=sourceIBD[i][1];
		if( E_l < (mgl_i+1)*p.bin ){
			sourceibd[k][0]=E_l;sourceibd[k][1]=E_r;
			sourceibd[k][2]=(int)floor(fastphy2gen(E_l, m, p)*1e6);
			sourceibd[k][3]=(int)floor(fastphy2gen(E_r, m, p)*1e6);
//printf("joe=%d %d %d %d\n", sourceibd[k][0], sourceibd[k][1],sourceibd[k][2],sourceibd[k][3]);
			sourceibd_size_l++;
		}
	}

	ibdmax=1e6*(m[p.L-1].g-m[mgl_i+1].g);
	//ibdmax=1e6*30;
	lmincoverage=INT_MAX;
	for(i=0;i<mgl_i;i++){
		onesidecoverage[i]=cal_left_coverage(i, sourceibd, sourceibd_size_l, m, ibdmax, p);
		if(lmincoverage>=onesidecoverage[i])lmincoverage=onesidecoverage[i];
		//fprintf(stderr, "%lf %lf %lf\n", ratio[i], X[i]/X[endl_i], m[i].r/m[endl_i].r);
		//fprintf(stderr, "left %d %lf\n", i, onesidecoverage[i]);
	}
	for(i=0;i<mgl_i;i++){        
		onesideends[i]=count_left_end(i, lmincoverage, sourceibd, sourceibd_size_l, m, ibdmax, p);
		//fprintf(stderr, "joe- left %d %lf %lf\n", i, onesideends[i], onesidecoverage[i]);
	}

	sumX=0;//recording the total number of IBD ends used for estiamtion
	sumx=0;sumr=0;
	for(i=endl_i;i<mgl_i;i++){sumx+=X[i];sumr+=onesideends[i];}
	for(i=0;i<endl_i;i++){
		if(sumr>0)X[i]=sumx*onesideends[i]/sumr;else X[i]=0;
		sumX=sumX+onesideends[i];}

	for(i=endl_i;i<endr_i;i++)sumX=sumX+X[i];
	//adjust the right ends
	////find the central region for 7 bins
	int sourceibd_size_r;
	sourceibd_size_r=0;
	for(i=0;i<IBDNUMtot;i++){
		k=sourceibd_size_r;
		//IBD overlapping with the target region endr_i-p.L-1, c-c+7
		E_l=sourceIBD[i][0];E_r=sourceIBD[i][1];
		if( E_r > (mgr_i-1)*p.bin ){
			sourceibd[k][0]=E_l;sourceibd[k][1]=E_r;
			sourceibd[k][2]=(int)floor(fastphy2gen(E_l, m, p)*1e6);
			sourceibd[k][3]=(int)floor(fastphy2gen(E_r, m, p)*1e6);
			sourceibd_size_r++;
		}
	}

	ibdmax=1e6*(m[mgr_i-1].g-m[0].g);
	//ibdmax=1e6*30;

	//fprintf(stderr, "%d %lf\n", c, ibdmax/1e6);
	//for(j=0;j<7;j++){fprintf(stderr, "%d %lf\n", c+j, X[c+j]);}fprintf(stderr, "\n");

	rmincoverage=INT_MAX;
	for(i=mgr_i;i<p.L;i++){
		onesidecoverage[i]=cal_right_coverage(i, sourceibd, sourceibd_size_r, m, ibdmax, p);
		if(rmincoverage>=onesidecoverage[i])rmincoverage=onesidecoverage[i];
		//fprintf(stderr, "right %d %lf\n",i, onesidecoverage[i]);
	}
	if(mute == 0)fprintf(stderr, "\n\tCoverage: left = %lf, right = %lf\n", lmincoverage, rmincoverage);
	for(i=mgr_i;i<p.L;i++){
		onesideends[i]=count_right_end(i, rmincoverage, sourceibd, sourceibd_size_r, m, ibdmax, p);
		//fprintf(stderr, "joe- right %d %lf %lf\n", i, onesideends[i], onesidecoverage[i]);
	}
	//abort();
	sumx=0;sumr=0;

	for(i=mgr_i;i<endr_i;i++){sumx+=X[i];sumr+=onesideends[i];}
	for(i=endr_i;i<p.L;i++){
		if(sumr>0)X[i]=sumx*onesideends[i]/sumr;else X[i]=0;
		sumX=sumX+onesideends[i];}

	for(i=0;i<IBDNUMtot;i++){free(sourceibd[i]);}free(sourceibd);
	free(onesideends);free(onesidecoverage);
	checkpt=clock();
	//fprintf(stderr,"adjust_end_counts time:%lfs\n", (checkpt-begin)*1.0/CLOCKS_PER_SEC);
	return 0;

}

double count_left_end(int b1, double bcoverage, int **sourceIBD, int sourceIBD_L, struct map m[], double ibdmax, struct mapPar p)
{
	//bin b1 is at the left end while bin b2 is at the central
	int i, j, k, L;
	int E_l, E_r, G_l, G_r;
	int **b1end;
	double a, b;
	b1end=(int **)malloc(sourceIBD_L*sizeof(int*));
	for(i=0;i<sourceIBD_L;i++){
		b1end[i]=(int *)malloc(5*sizeof(int));
	}
	int b1endsize;
	b1endsize=0;
	//calculate the coverage for these two bins
	for(i=0;i<sourceIBD_L;i++){
		E_l=sourceIBD[i][0];E_r=sourceIBD[i][1];
		G_l=sourceIBD[i][2];G_r=sourceIBD[i][3];
		//for b1
		a=p.bin*b1;b=p.bin*(b1+1);
		//fprintf(stderr, "%d %lf %lf %d %d %d\n",i, a, b, E_l, E_r, minIBD);
		if(E_l<b&&E_r>=b){
			k=b1endsize;
			if(E_l<=a){b1end[k][1]=(int)floor(a);b1end[k][3]=m[b1].g*1e6;}
			else {b1end[k][1]=E_l;b1end[k][3]=G_l;}
			b1end[k][2]=E_r;b1end[k][4]=G_r;
			b1end[k][0]=b1end[k][4]-b1end[k][3];
			if((b1end[k][2]-b1end[k][1])>=minIBD&&b1end[k][0]>(p.ibdcut*1e6*LENthreshold)&&b1end[k][0]<ibdmax){
				b1endsize++;
			}
		}
	}
	sort_ibd(b1end, STEP, b1endsize);
	//count the left ends in each bin
	////for bin b1
	double X1, l;
	a=p.bin*b1;b=p.bin*(b1+1);l=0;X1=0;
	for(i=0;i<b1endsize;i++){
		E_l=b1end[i][1];E_r=b1end[i][2];
		if(E_l>a&&E_l<b){X1++;l+=cal_overlap_rate(a, b, E_l, E_r, p.bin);}
		if(E_l<=a&&E_r>b)l+=1;
		if(l>=bcoverage)break;
//fprintf(stderr, "%lf %lf %d %lf\n", l, bcoverage, b1endsize, cal_overlap_rate(a, b, E_l, E_r, p.bin));
	}
	for(i=0;i<sourceIBD_L;i++){free(b1end[i]);}
	free(b1end);
	return X1;
}


double count_right_end(int b1, double bcoverage, int **sourceIBD, int sourceIBD_L, struct map m[], double ibdmax, struct mapPar p)
{
	int i, j, k;
	int a, b, E_l, E_r, G_l, G_r;
	int **b1end;
	b1end=(int **)malloc(sourceIBD_L*sizeof(int*));
	for(i=0;i<sourceIBD_L;i++){
		b1end[i]=(int *)malloc(5*sizeof(int));
	}
	int b1endsize;
	b1endsize=0;
	//calculate the coverage for these two bins
	for(i=0;i<sourceIBD_L;i++){
		E_l=sourceIBD[i][0];E_r=sourceIBD[i][1];
		G_l=sourceIBD[i][2];G_r=sourceIBD[i][3];
		//for b1
		a=p.bin*b1;b=p.bin*(b1+1);
		if(E_l<a&&E_r>=a){
			k=b1endsize;
			b1end[k][1]=E_l;b1end[k][3]=G_l;
			//if(b1+1 >= p.L)fprintf(stderr, "%d %lf\n", b1+1, m[b1+1].g);
			if(E_r>=b){
				b1end[k][2]=b;
				if(b1<p.L-1)b1end[k][4]=m[b1+1].g*1e6;
				else b1end[k][4]=(int)floor(m[b1].g*1e6+m[b1].r*p.bin);
			}
			else {b1end[k][2]=E_r;b1end[k][4]=G_r;}
			b1end[k][0]=b1end[k][4]-b1end[k][3];
			if((b1end[k][2]-b1end[k][1])>=minIBD&&b1end[k][0]>(p.ibdcut*LENthreshold*1e6)&&b1end[k][0]<ibdmax){b1endsize++;}
		}
	}
	sort_ibd(b1end, STEP, b1endsize);

	//count the left ends in each bin
	double type1, type2;
	//////for bin b1
	double X1, l;
	a=p.bin*b1;b=p.bin*(b1+1);l=0;X1=0;type1=0;type2=0;
	for(i=0;i<b1endsize;i++){
		E_l=b1end[i][1];E_r=b1end[i][2];
		if(E_r>a&&E_r<b){X1++;l+=cal_overlap_rate(a, b, E_l, E_r, p.bin);type1++;}
		//if(E_r>a&&E_r<b){X1++;l+=1;type1++;}
		if(E_r>=b&&E_l<a){l+=1;type2++;}
		if(l>=bcoverage)break;
	}
	for(i=0;i<sourceIBD_L;i++){free(b1end[i]);}
	free(b1end);
	return X1;
}

double cal_left_coverage(int bin, int **sourceIBD, int sourceIBD_L, struct map m[], double ibdmax, struct mapPar p)
{
	int i;
	double  a, b;
	int E_r, E_l, G_r, G_l;
	double x, bcoverage=0;
	double bend[5];
	for(i=0;i<sourceIBD_L;i++){
		E_l=sourceIBD[i][0];E_r=sourceIBD[i][1];
		G_l=sourceIBD[i][2];G_r=sourceIBD[i][3];

		a=p.bin*bin;b=p.bin*(bin+1);
		if(E_l<b&&E_r>=b){
			if(E_l<=a){bend[1]=a;bend[3]=m[bin].g*1e6;}
			else {bend[1]=E_l;bend[3]=G_l;}
			bend[2]=E_r;bend[4]=G_r;
			bend[0]=bend[4]-bend[3];
			if((bend[2]-bend[1])>=minIBD&&bend[0]>(p.ibdcut*LENthreshold*1e6)&&bend[0]<ibdmax){
				x=cal_overlap_rate(a, b, E_l, E_r, p.bin);
				bcoverage+=x;
			}
		}
	}
	return bcoverage;
}

double cal_right_coverage(int bin, int **sourceIBD, int sourceIBD_L, struct map m[], double ibdmax, struct mapPar p)
{
	int i;
	double  a, b;
	int E_r, E_l, G_r, G_l;
	double x, bcoverage=0;
	double bend[5];
	for(i=0;i<sourceIBD_L;i++){
		E_l=sourceIBD[i][0];E_r=sourceIBD[i][1];
		G_l=sourceIBD[i][2];G_r=sourceIBD[i][3];

		a=p.bin*bin;b=p.bin*(bin+1);
		if(E_l<a&&E_r>=a){
			bend[1]=E_l;bend[3]=G_l;
			if(E_r>=b){bend[2]=b;bend[4]=fastphy2gen(b, m, p)*1e6;}
			else {bend[2]=E_r;bend[4]=G_r;}
			bend[0]=bend[4]-bend[3];
			if((bend[2]-bend[1])>=minIBD&&bend[0]>(p.ibdcut*LENthreshold*1e6)&&bend[0]<ibdmax){
				x=cal_overlap_rate(a, b, E_l, E_r, p.bin);
				bcoverage+=x;
			}
		}
	}
	return bcoverage;
}

double cal_side_fine_bin_depth(int I, int **end, double *coverage, struct map m[], struct mapPar p)
{
	clock_t checkpt, begin = clock();
	int i, j;
	double from,to;
	int E_l, E_r;
	int L;
	double bin;
	L = fineL;
	bin = p.bin;
	int k= 0;
	for(i=0;i<L;i++){coverage[i]=0;}
	for(i=0;i<L;i++){
		from = floor(m[i].p); to = floor(m[i].p + bin);
		for(j=0;j<p.ibdnum;j++){
			E_l = end[j][1]; E_r = end[j][2];
			if(I < LENDindex && E_l < from) E_l = (int)from - 1;
			if(I > RENDindex && E_r > to) E_r = (int)to + 1;
			if((E_r - E_l) > fBIN || E_r > E_l){
				if(I < LENDindex && E_r > to)coverage[i] += cal_overlap_rate(from, to, E_l, E_r, bin);
				if(I > RENDindex && E_l < from)coverage[i] += cal_overlap_rate(from, to, E_l, E_r, bin);
			}
		}
	}
	double mindepth=INT_MAX;
	for(i=0;i<L;i++){
		if(mindepth>=coverage[i])mindepth=coverage[i];
	}
	checkpt=clock();
	if(mute == 0)fprintf(stderr, "\t%d-Minimum depth = %lf\n",I, mindepth);
	return mindepth;
}

int end_counts_in_side_fine_bin(int I, double *X, int **end, struct map m[megL][fineL], struct mapPar p[megL])
{
	int i, j, k;
	int E_l, E_r, G_l, G_r;
	int **bend;
	double from, to;
	int bendsize;
	int L = p[I].ibdnum;
	double l;


double midibdcm;
double minibdcm1;
double minibdcm2;
	
	bend=(int **)malloc(L*sizeof(int*));
	for(i=0;i<L;i++){
		bend[i]=(int *)malloc(5*sizeof(int));
	}
	i = 0;
	for(i=0;i<fineL;i++){
		from = floor(m[I][i].p); 
		to = floor(m[I][i].p + p[I].bin);
		l=0;
		bendsize=0;
		minibdcm1 = INT_MAX;
		minibdcm2 = INT_MAX;
		for(j=0;j<L;j++){
			E_l = end[j][1]; E_r = end[j][2];
			if(I < LENDindex && E_l < from) E_l = (int)from - 1;
			if(I > RENDindex && E_r > to) E_r = (int)to + 1;
			if((E_r - E_l) > fBIN || E_r > E_l){
				k=bendsize;
				bend[k][1]=E_l;bend[k][2]=E_r;
				bend[k][3]=(int)floor(finephy2gen(E_l, m, p)*1e6);;
				bend[k][4]=(int)floor(finephy2gen(E_r, m, p)*1e6);;
				bend[k][0]=bend[k][4]-bend[k][3];
				bendsize++;
			}
		}
		sort_ibd(bend, STEP, bendsize);
		for(j=0;j<bendsize;j++){
			E_l=bend[j][1];E_r=bend[j][2];
			if(cal_overlap_rate(from, to, E_l, E_r, p[I].bin) > 0.00001 \
					&& bend[j][0] < minibdcm2)minibdcm2=bend[j][0];
			if((E_l>from&&E_l<to) && I < LENDindex){
				X[i]++;
				if(bend[j][0]<minibdcm1)minibdcm1=bend[j][0];
				l+=cal_overlap_rate(from, to, E_l, E_r, p[I].bin);}
			if((E_r>from&&E_r<to) && I > RENDindex){
				X[i]++;
				if(bend[j][0]<minibdcm1)minibdcm1=bend[j][0];
				l+=cal_overlap_rate(from, to, E_l, E_r, p[I].bin);}
			if(E_l<=from&&E_r>=to)l+=1;

			if(l>=p[I].mindepth)break;
		}
		midibdcm = bend[j/2][0];
//fprintf(stderr, "%d %d %.0lf %.0lf %.0lf\n", I, i, midibdcm, minibdcm1, minibdcm2);
if(minibdcm1 >= INT_MAX-1 )minibdcm1 = -9e6;
		m[I][i].midibdcm = midibdcm/1e6;
		m[I][i].minibdcm1 = minibdcm1/1e6;
		m[I][i].minibdcm2 = minibdcm2/1e6;

	}
	for(i=0;i<L;i++){free(bend[i]);}
	free(bend);
	return 0;
}

