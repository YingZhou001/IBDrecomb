#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include<unistd.h>
#include<math.h>
#include "head.h"
#include <limits.h>


int checkinput(int arg, char**argv)
{
    fprintf(stderr, "####################\n");
    fprintf(stderr, "#Welcome to IBDrec!#\n");
    fprintf(stderr, "####################\n");
    fprintf(stderr, "Copyright (C) 2019 Ying Zhou\n\n");
    int i, j , k;
    BIN=500000;
    fBIN=-1;
    TRIM=1;
    mute=1;STEP=500;STAT=0;
    ITER=20;FOLD=1;RECrate=1;ENDrate=1;
    ibd=0;simibd=0;refinedibd=0;
    LENthreshold=0;
    oplink=1;
    ohapmap=0;
    sprintf(CHR, ".");
    sprintf(out, "./");
    for(i=1;i<arg;i++){
	if(strcmp(argv[i], "-ibd")==0){ibd=1;sprintf(ibdinp, "%s", argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-simibd")==0){simibd=1;sprintf(ibdinp, "%s", argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-refinedibd")==0){refinedibd=1;sprintf(ibdinp, "%s", argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-rate")==0){RECrate=atof(argv[i+1]);i=i+1;}
//	else if(strcmp(argv[i], "-endrate")==0){ENDrate=atof(argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-bin")==0){BIN=atoi(argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-fbin")==0){fBIN=atoi(argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-fold")==0){FOLD=atof(argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-trim")==0){TRIM=atoi(argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-chr")==0){sprintf(CHR, "%s", argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-iter")==0){ITER=atoi(argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-lcut")==0){LENthreshold=atoi(argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-mute")==0){mute=atoi(argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-stat")==0){STAT=atoi(argv[i+1]);i=i+1;}
	else if(strcmp(argv[i], "-ohapmap")==0){ohapmap=1;i=i+1;}
	else if(strcmp(argv[i], "-help")==0){print_help();}
	else {fprintf(stderr, "Error!, can't recognize %s!\n", argv[i]);print_help();}
    }
    fprintf(stderr, "\nCheck parameters:\n");
    fprintf(stderr, "  IBD inputs: %s ", ibdinp);
    if(access(ibdinp, R_OK) != -1)fprintf(stderr, "exist!\n"); 
    else {fprintf(stderr, "Does not exist\n");
	print_help();
    }
    fprintf(stderr, "  Trim length: %d\n", TRIM);
    fprintf(stderr, "  Bin size: %.0lf\n", BIN);
    if(fBIN>0){
	if(fBIN<=BIN)fprintf(stderr, "  Fine bin size: %.0lf\n", fBIN);
	else {
	    BIN=fBIN;
	    fprintf(stderr, "  Resize fine bin to: %.0lf\n", fBIN);
	}
    }
    else {fprintf(stderr, "output scale is not specified, use -fbin\n");
	print_help();}
    fprintf(stderr, "  Scale rate: %lf\n", RECrate);
    fprintf(stderr, "  Iterations: %d\n", ITER);
    if(ohapmap == 1){
	fprintf(stderr, "  Output format: hapmap\n");
	oplink = 0;
    }
    else {
	fprintf(stderr, "  Output format: plink\n");
    }

    //fprintf(stderr, "  \n");
    return 0;
}

/*\t-fold 1.0\t# Relative size of adjunct-region to the end-region (see Methods of the published paper)\n \
*/
int print_help(void)
{
    fprintf(stderr, "Usage: IBDrec [options] parameters\n\n \
(Required inputs:)\n \
\t-refinedibd/-ibd <filename>\t# IBD segments input, '-refinedibd' for the refined-IBD's output format and '-ibd' for two-column generic format\n \
\t-fbin <integer>\t# Bin size(bp) for outputing recombination rates\n\n \
(Optional parameters:)\n \
\t-bin 500000\t# Bin size(bp) for first step of estimation\n \
\t-rate 1.0\t# Average cM/Mb to normalize the estimated recombination rates\n \
\t-trim 1\t\t# Length(bp) to trim at chromosome ending\n \
\t-iter 20\t# Iteration number\n \
\t-chr .\t\t# Chromosome to estimate\n \
\t-help \t\t# Print this help file\n \
\n");
    exit(-1);
}

int check_gz(char *inpfile)
{
    unsigned char buffer[2]; // note: 1 byte
    FILE *ifp;
    ifp=fopen(inpfile, "rb");
    fread(buffer, 2, 1, ifp);
    fclose(ifp);
    char tmp[32];
    sprintf(tmp, "%x%x",buffer[0], buffer[1]);
    if(strcmp(tmp, "1f8b")==0)return 1; //the input is gz file
    else return 0;
}

int check_ibd(char *ibdinp)
{
    FILE *ifp;
    double a, b;
    int gztag, length, chr_tag;
    char chr[32], chr_old[32];
    fprintf(stderr, "Check IBD input:\n");
    gztag=check_gz(ibdinp);
    if(gztag==1){char cmd[10240];sprintf(cmd, "zcat %s", ibdinp);ifp=popen(cmd, "r");}
    else {ifp=fopen(ibdinp, "r");}
    length=0;
    minIBD=INT_MAX;
    POSmin=INT_MAX;POSmax=0;
    chr_tag=0;
    while(!feof(ifp)){
	if(ibd)fscanf(ifp, "%lf %lf", &a, &b);
	if(simibd)fscanf(ifp, "%*s %*s %*s %*s %*s %lf %lf %*s", &a, &b);
	if(refinedibd){
	    fscanf(ifp, "%*s %*s %*s %*s %s %lf %lf %*s %*s", chr, &a, &b);
	    if(strcmp(CHR, ".")==0){
		if(strcmp(chr,chr_old)!=0)chr_tag++;
		if(chr_tag==2){
		    fprintf(stderr, "Error, IBD segments from multiple chromosomes are detected, please use '-chr' to specify the one to be estimated.\n");
		    exit(-1);
		}
		sprintf(chr_old, "%s", chr);
		sprintf(CHR, "%s", chr);
	    }
	    else {
		if(strcmp(chr,chr)!=0)continue;
	    }
	}
	length++;
	if(a<POSmin)POSmin=a;
	if(b>POSmax)POSmax=b;
	if((b-a)<minIBD){
	    minIBD=(int)floor(b-a);
	}
    }
    if(gztag==1)pclose(ifp);
    else fclose(ifp);
    fprintf(stderr, "  %d IBD segments, cover region = (%d %d)\n",length, POSmin, POSmax);
    return length;
}

int check_max_IBDnum_megbin(int **sourceIBD, struct mapPar p)
{
    int i, j;
    int a,b;
    double *cov;
    cov=(double*)malloc(p.L*sizeof(double));
    for(i=0;i<p.L;i++){cov[i]=0;}
    for(j=0;j<IBDNUMtot;j++){
	a=(int)floor(sourceIBD[j][0]/p.bin);
	b=(int)floor(sourceIBD[j][1]/p.bin);
	for(i=a;i<=b;i++){if(i>=0&&i<p.L)cov[i]++;}
    }
    double maxdepth=0;
    for(i=0;i<p.L;i++){
	if(maxdepth<=cov[i])maxdepth=cov[i];
    }
    free_array(cov);
    if(mute == 0)fprintf(stderr, "\tMaximum depth = %.0lf\n", maxdepth);
    if(maxdepth > INT_MAX)return INT_MAX;
    else return (int)floor(maxdepth)+1;
}
