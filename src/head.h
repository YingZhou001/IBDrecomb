#define FLEN 10240


char ibdinp[FLEN], posinp[FLEN], vcfinp[FLEN], out[FLEN];
char CHR[32];

double BIN, fBIN;
double sumX; //total number of IBD ends being used
int ITER;
int STEP;
int STAT; //whether print out statistics of IBD length for each bin
int LENthreshold;//length threshold for step 1, automo


int IBDNUM, IBDNUMtot;
int maxIBDnum; //for each mega bins
int minIBD; //minimum ibd physical length in bp
//int SNPNUM;
int POSmin, POSmax;
int POSdrift;// move the starting position to 1
//int BINNUM;
int TRIM;
//int TRIMl, TRIMr;
int LENDindex, RENDindex;
int END;
double RECrate;
double gMID_l,gMID_r;///min/max-imum IBD genetic length 
double gMID;///midian IBD genetic length 
double ENDrate, LENDlength, RENDlength;//to define the range of end-region
double IBDCUT;
//int REFSNPNUM;
//double MINdepth;
double FOLD;

int megL, fineL;

int simibd;
int refinedibd;
int ibd;
int mute;
int oplink, ohapmap;


struct map {
	double p;
	double r;
	double g;
	double midibdcm;
	double minibdcm1;
	double minibdcm2;
};

struct mapPar {
	int L;
	double bin;
	double minibd;
	double ibdcut;
	double mindepth;
	int ibdnum;
	int trimL;
	int trimR;
};


