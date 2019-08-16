#include "head.h"
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<limits.h>

double fastphy2gen(double x0, struct map m[], struct mapPar p)
{
    int i, L;
    double bin;
    double x1, y1, x2,y2;
    double y0;

    L = p.L;
    bin = p.bin;

    if(x0<=m[0].p){
	x1=m[0].p;y1=m[0].g;
	x2=m[3].p;y2=m[3].g;
	y0=y1+(x0-x1)*(y2-y1)/(x2-x1);
    }
    else if(x0>=m[L-1].p){
	x1=m[L-4].p;y1=m[L-4].g;
	x2=m[L-1].p;y2=m[L-1].g;
	y0=y2+(x0-x2)*(y2-y1)/(x2-x1);
    }
    else {
	i=(int)floor(x0/bin);
	x1=m[i].p;y1=m[i].g;
	x2=m[i+1].p;y2=m[i+1].g;
	y0=y1+(x0-x1)*(y2-y1)/(x2-x1);
    }
    //fprintf(stderr, "L:%d bin:%lf in:%lf out:%lf\n",L, bin, x0, y0);
    return y0;
}

double finephy2gen(double x0, struct map m[megL][fineL], struct mapPar p[megL])
{
    int i, j, L;
    double bin, megbin;
    double x1, y1, x2,y2;
    double y0;
    L = megL;
    megbin = m[1][0].p - m[0][0].p;
    bin = p[0].bin;
    // search the mega map
    if(x0<=m[0][0].p){
	x1=m[0][0].p;y1=m[0][0].g;
	x2=m[2][0].p;y2=m[2][0].g;
	y0=y1+(x0-x1)*(y2-y1)/(x2-x1);
    }
    else if(x0>=(m[L-1][fineL-1].p)){
	x1=m[L-3][0].p;y1=m[L-3][0].g;
	x2=m[L-1][fineL-1].p;y2=m[L-1][fineL-1].g;
	y0=y2+(x0-x2)*(y2-y1)/(x2-x1);
    }
    else {
	i=(int)floor(x0/megbin);
	j=(int)floor((x0 - m[i][0].p)/bin);
	x1=m[i][j].p;y1=m[i][j].g;
	if((j+1) < fineL){
	    x2=m[i][j+1].p;y2=m[i][j+1].g;
	    y0=y1+(x0-x1)*(y2-y1)/(x2-x1);
	}
	else {
	    y0=y1+(x0-x1)*m[i][j].r/1e6;
	}
	//fprintf(stderr, "**L:%d bin:%lf megbin: %lf in:%lf out:%lf %lf %lf\n",L, bin, megbin, x0, y0, x1, y1);
	//fprintf(stderr, "i: %d, j: %d, %d %d\n", i, j, megL, fineL);
	//fprintf(stderr, "p:%lf g:%lf\n", m[i][j].p, m[i][j].g);
    }
    return y0;
}



void constoutmap(struct map m[], struct mapPar p)
{
    int i;
    int L, bin;
    for ( i = 0; i < p.L; i++ ){
	m[i].p = i * p.bin;
	m[i].r = 1;
    }
    m[0].g = 0;
    for( i = 1; i < p.L; i++ ){
	m[i].g = m[i-1].g + m[i].r * p.bin/1e6;
    }
}

void update_TRIM(struct map m[], struct mapPar p)
{
    // LENDindex and RENDindex record the ranges of the mid area
    int i, p0;
    int begin_L, end_L;
    if((m[p.L-1].g/4) < LENDlength || (m[p.L-1].g/4) < RENDlength){
	fprintf(stderr, "Error! -endrate %lf is too large, you need to set a smaller value.\n", ENDrate);
	exit(-1);
    }

    for(i=0;i<p.L;i++){
	if((m[i].g-m[0].g)>LENDlength){
	    LENDindex = i; //[0, LEND)
	    begin_L = m[i].p;
	    break;}}
    for(i=p.L-1;i>=0;i--){
	if((m[p.L-1].g + m[p.L-1].r*p.bin/1e6 - m[i].g)>RENDlength){
	    RENDindex = i; //[REND, BINNUM)
	    end_L = (p.L - i) * p.bin;
	    break;}}
    if(mute==0)fprintf(stderr, "\tEnd region length: left = %d,  right = %d\n", begin_L, end_L);
    if(mute==0)fprintf(stderr, "\tEnd region index: left = %d,  right = %d\n", LENDindex, RENDindex);

}

double update_IBDCUT(struct map m[], struct mapPar p)
{
    int i;
    double ibdcut, newibdcut;
    ibdcut = 0;
    for(i = 0; i < p.L; i ++){
	newibdcut = fastphy2gen(m[i].p + minIBD, m, p) - m[i].g;
	if(newibdcut > ibdcut)ibdcut = newibdcut;
    }
    if (p.ibdcut > m[p.L-1].g/3){
	fprintf(stderr, "Warning: reset IBDCUT from %lf to %lf, bias may be brought in if you see this in the last iteration\n", ibdcut, m[p.L-1].g/4);
	ibdcut = m[p.L-1].g/4;
    }
    return(ibdcut);
}

//sort_ibd
//this is a roughly way to sort integers
////results may not so accurate but enough for the IBD sorting
////numbers in [a-step, a+step) are in the bin a
void sort_ibd(int **end, int step, int ibdnum)
{
    clock_t checkpt, begin = clock();
    int i, j, k, l, max_L;
    int **end_sort;
    end_sort=(int **)malloc(ibdnum*sizeof(int*));
    for(i=0;i<ibdnum;i++){end_sort[i]=(int *)malloc(5*sizeof(int));}
    max_L=0;
    for(i=0;i<ibdnum;i++){if(max_L<=end[i][0])max_L=end[i][0];}
    int binnum;
    binnum=(int)floor((max_L+1)*1.0/step)+1;

    int **bin;;
    int *bin_ibd_count;
    bin_ibd_count=(int *)malloc(binnum*sizeof(int));
    for(i=0;i<binnum;i++){bin_ibd_count[i]=0;}
    bin=(int **)malloc(binnum*sizeof(int*));
    for(i=0;i<ibdnum;i++){k=binnum-(int)floor(end[i][0]*1.0/step)-1;bin_ibd_count[k]++;}
    for(i=0;i<binnum;i++){
	bin[i]=(int *)malloc((bin_ibd_count[i]+1)*sizeof(int));
	bin_ibd_count[i]=0;
    }
    for(i=0;i<ibdnum;i++){
	k=binnum-(int)floor(end[i][0]*1.0/step)-1;
	if(k>binnum){fprintf(stderr, "2test%d %d %d\n",k, binnum, end[i][0]);}
	bin[k][bin_ibd_count[k]]=i;
	bin_ibd_count[k]++;
    }
    k=0;
    for(i=0;i<binnum;i++){
	if(bin_ibd_count[i]>0){
	    for(j=0;j<bin_ibd_count[i];j++){
		memcpy(end_sort[k], end[bin[i][j]], sizeof(int)*5);
		k++;
	    }
	}
    }
    for(i=0;i<ibdnum;i++)memcpy(end[i], end_sort[i], sizeof(int)*5);
    checkpt=clock();
    for(i=0;i<ibdnum;i++){free(end_sort[i]);}free(end_sort);
    for(i=0;i<binnum;i++){free(bin[i]);}free(bin);
    free(bin_ibd_count);
    //if(mute == 0)fprintf(stderr,"\tsort_ibd time:%lfs\n", (checkpt-begin)*1.0/CLOCKS_PER_SEC);
}


//cal_bin_depth
double cal_bin_depth(int **end, double *coverage, struct map m[], struct mapPar p)
{
    clock_t checkpt, begin = clock();
    int i, j;
    int a,b;
    for(i=0;i<p.L;i++){coverage[i]=0;}
    for(j=0;j<p.ibdnum;j++){
	a=(int)floor(end[j][1]/p.bin);
	b=(int)floor(end[j][2]/p.bin);
	if(a<b){//dont allow entire IBD seg staying in the bin
	    for(i=a+1;i<b;i++){if(i>=0&&i<p.L)coverage[i]++;}
	    if(a>=0&&b<p.L){
		coverage[a]+=(a+1)-end[j][1]/p.bin;
		coverage[b]+=end[j][2]/p.bin-b;
	    }
	}
	if(a == b ){
	    if(a>=0&&a<p.L)coverage[a]+=(end[j][2]-end[j][1])/p.bin;
	}
    }
    double mindepth=INT_MAX;
    for(i=LENDindex;i<RENDindex;i++){
	if(mindepth>=coverage[i])mindepth=coverage[i];
	//fprintf(stderr, "%lf %lf\n", mindepth, coverage[i]);
    }
    //for(i=0;i<p.L;i++)fprintf(stderr, "joe-%lf %lf\n", mindepth, coverage[i]);
    checkpt=clock();
    if(mute == 0)fprintf(stderr, "\tMinimum depth = %lf\n", mindepth);
    return mindepth;
}

double cal_fine_bin_depth(int **end, double *coverage, struct map m[], struct mapPar p)
{
    //fprintf(stderr, "\ttest\n");
    clock_t checkpt, begin = clock();
    int i, j;
    int a,b;
    int L;
    double bin;
    L = fineL;
    bin = p.bin;
    int k= 0;
    for(i=0;i<L;i++){coverage[i]=0;}
    for(j=0;j<p.ibdnum;j++){
	a=(int)floor((end[j][1]-m[0].p)/bin);
	b=(int)floor((end[j][2]-m[0].p)/bin);
	if(a<b){
	    for(i=0;i<L;i++){if(i>a&&i<b)coverage[i]++;}
	    if(a>=0&&a<L)coverage[a]+=(a+1)-(end[j][1]-m[0].p)/bin;
	    if(b>=0&&b<L)coverage[b]+=(end[j][2]-m[0].p)/bin-b;
	}
	if(a==b){
	    if(a>=0&&a<L)coverage[a]+=(end[j][2]-end[j][1])/bin;
	}
    }
    //for(i=0;i<L;i++){fprintf(stderr, "%d %d %lf\n",i, p.ibdnum,coverage[i]);}
    double mindepth=INT_MAX;
    for(i=0;i<L;i++){
	if(mindepth>=coverage[i])mindepth=coverage[i];
    }
    checkpt=clock();
    if(mute == 0)fprintf(stderr, "\tMinimum depth = %lf\n", mindepth);
    return mindepth;
}


double cal_overlap_rate(double a, double b, int E_l, int E_r, double bin)
{
    double overlap, out;
    overlap=0;
    if(E_l>=a&&E_r<b)overlap=E_r-E_l;
    if(E_l<a&&E_r>=b)overlap=bin;
    if(E_l>=a&&E_l<b&&E_r>=b)overlap=b-E_l;
    if(E_l<a&&E_r<b&&E_r>=a)overlap=E_r-a;
    out = overlap/bin;
    if(out>1.01)fprintf(stderr, "joe**%lf %lf %d %d %lf %lf\n", a,b,E_l,E_r, overlap, out);
    if(out>1)out=1;
    return out;
}

int end_counts_in_bin(double *X, int **end, int **sourceIBD, struct mapPar p)
{
    clock_t checkpt, begin = clock();
    int i, j;
    double a, b;
    double l;
    int E_r, E_l;
    for(i=0;i<p.L;i++){
	a=i*p.bin;b=(i+1)*p.bin;l=0;
	X[i]=0;
	for(j=0;j<p.ibdnum;j++){
	    E_l=end[j][1];
	    E_r=end[j][2];
	    if(E_l>a&&E_r<b){
		X[i]+=2;
		l+=cal_overlap_rate(a, b, E_l, E_r, p.bin);
	    }
	    else {
		if(E_l>=a&&E_l<b){X[i]++;}
		if(E_r>=a&&E_r<b){X[i]++;}
		if((a<=E_l&&b>E_l)||(a>E_l&&b<=E_r)||(a<=E_r&&b>E_r)){
		    l+=cal_overlap_rate(a, b, E_l, E_r, p.bin);
		}
	    }
	    if(l>=p.mindepth){break;}
	}		
	//fprintf(stderr, "%lf, %lf %lf\n", X[i], l, p.mindepth);
    }
    checkpt=clock();
    return 0;
}

int end_counts_in_fine_bin(double *X, int **end, int **sourceIBD, struct map m[], struct mapPar p)
{
    clock_t checkpt, begin = clock();
    int i, j;
    double a, b;
    double l;
    int E_r, E_l;

    double midibdcm;
    double minibdcm1;
    double minibdcm2;


    for(i=0;i<p.L;i++){
	//a=i*p.bin;b=(i+1)*p.bin;l=0;
	a=m[i].p;b=m[i].p+p.bin;l=0;

	X[i]=0;
	minibdcm1 = INT_MAX;
	minibdcm2 = INT_MAX;
	for(j=0;j<p.ibdnum;j++){
	    //E_l=end[j][1]-(int)floor(m[0].p);
	    //E_r=end[j][2]-(int)floor(m[0].p);
	    E_l=end[j][1];
	    E_r=end[j][2];

	    //fprintf(stderr, "%lf %lf %d %d\n", a, b, E_l, E_r);
	    if(E_l>=a&&E_r<=b){
		X[i]+=2;
		if(end[j][0]<minibdcm1)minibdcm1=end[j][0];
		l+=cal_overlap_rate(a, b, E_l, E_r, p.bin);
		if(end[j][0]<minibdcm2)minibdcm2=end[j][0];
	    }
	    else {
		if(E_l>=a&&E_l<b){
		    X[i]++;
		    if(end[j][0]<minibdcm1)minibdcm1=end[j][0];
		}
		if(E_r>=a&&E_r<b){
		    X[i]++;
		    if(end[j][0]<minibdcm1)minibdcm1=end[j][0];
		}
		if((a<=E_l&&b>E_l)||(a>E_l&&b<=E_r)||(a<=E_r&&b>E_r)){
		    l+=cal_overlap_rate(a, b, E_l, E_r, p.bin);
		    if(end[j][0]<minibdcm2)minibdcm2=end[j][0];
		}
	    }
	    if(l>=p.mindepth){break;}
	}
	//if(tmp == 1)fprintf(stderr, "%.0lf %.0lf %.0lf %d\n", a, b, X[i], p.ibdnum);
	midibdcm = end[j/2][0];
	if(minibdcm1 >= INT_MAX-1 )minibdcm1 = -9e6;
	m[i].midibdcm = midibdcm/1e6;
	m[i].minibdcm1 = minibdcm1/1e6;
	m[i].minibdcm2 = minibdcm2/1e6;

	//fprintf(stderr, "%lf, %lf %lf\n", X[i], l, p.mindepth);
    }
    checkpt=clock();
    return 0;
}

void cal_end_region_legnth(int **end, int ibd_num)
{
    double *lend, *rend;
    lend = (double *) calloc (IBDNUMtot, sizeof(double));
    rend = (double *) calloc (IBDNUMtot, sizeof(double));
    long int Ll, Lr;
    Ll = 0; Lr = 0;
    long int i;

    for( i = 0; i <IBDNUMtot; i++) {
	if(end[i][1] < (POSmin - POSdrift + 10)){lend[Ll]=end[i][0]; Ll++;}
	if(end[i][2] > (POSmax - POSdrift - 10)){rend[Lr]=end[i][0]; Lr++;}
    }

    int a, b;
    a = (int) (Ll / 2.0);
    b = (int) (Lr / 2.0);

    LENDlength = lend[a] * ENDrate /1e6;
    RENDlength = rend[b] * ENDrate /1e6;
    if(mute==0)fprintf(stderr, "Left=%lf Right=%lf", LENDlength, RENDlength);
}



void init_fine_map(struct map outfineMap[megL][fineL], struct mapPar outfineMapPar[megL], struct map outMegMap[megL], struct mapPar outMegMapPar)
{
    int i, j, fineL;
    double newbin;
    fineL = (int)floor(( outMegMapPar.bin -0.1 ) / fBIN) + 1;
    newbin = outMegMapPar.bin / fineL;

    for( i = 0; i < megL; i++){
	outfineMap[ i ][ 0 ].p = outMegMap[ i ].p;
//fprintf(stderr, "%lf\n", outfineMap[ i ][ 0 ].p+POSdrift);
	outfineMap[ i ][ 0 ].r = outMegMap[ i ].r;
	outfineMap[ i ][ 0 ].g = outMegMap[ i ].g;
	outfineMapPar[ i ].L = fineL;
	outfineMapPar[ i ].bin = newbin;
	outfineMapPar[ i ].minibd = outMegMapPar.minibd;
	outfineMapPar[ i ].ibdcut = outMegMapPar.ibdcut;
	outfineMapPar[ i ].mindepth = outMegMapPar.mindepth;
	outfineMapPar[ i ].ibdnum = outMegMapPar.ibdnum;
	for( j = 1; j < fineL; j++){
	    outfineMap[ i ][ j ].p = outfineMap[ i ][ j-1 ].p + newbin;
//fprintf(stderr, "%lf\n", outfineMap[ i ][ j ].p+POSdrift);
	    outfineMap[ i ][ j ].r = outMegMap[ i ].r;
	    outfineMap[ i ][ j ].g = outfineMap[ i ][ j-1 ].g + newbin * outfineMap[ i ][ j-1 ].r/1e6;
	}
    }

}


void update_map(double *X, struct map m[], struct mapPar p)
{
    clock_t checkpt, begin = clock();
    int i;
    double sum_X;
    double sum_r;
    sum_X=0;sum_r=0;
    for(i=0;i<p.L;i++){
	sum_X+=X[i];
	sum_r+=1;
    }
    i=0;
    while(i<p.L){
	if(sum_X>0)m[i].r=X[i]*sum_r/sum_X;
	else m[i].r=0;
	i++;}

    for(i=1;i<p.L;i++){
	m[i].g=m[i-1].g+m[i-1].r*p.bin/1e6;
    }
    checkpt=clock();
    if(mute == 0)fprintf(stderr,"\tupdate_map time:%lfs\n", (checkpt-begin)*1.0/CLOCKS_PER_SEC);
}

void update_fine_map(double *X, struct map m[], struct mapPar p, double G)
{
    int i;
    double sum_X, sum_R, r;
    sum_X=0;sum_R=0;
    r = G/p.L*1e6/p.bin;
    for(i=0;i<p.L;i++){
	sum_X+=X[i];
	sum_R+=r;
    }
    i=0;
    while(i<p.L){
	if(sum_X>0)m[i].r=X[i]*sum_R/sum_X;
	else m[i].r=0;
	i++;}


    for(i=1;i<p.L;i++){
	m[i].g=m[i-1].g+m[i-1].r*p.bin/1e6;
    }
}


void cp_map(struct map m[megL][fineL], struct mapPar p[megL], struct map mOld[megL][fineL], struct mapPar pOld[megL])
{
    int i, j;
    for( i = 0; i < megL; i++){
	pOld[i].L = p[i].L;
	pOld[i].bin = p[i].bin;
	pOld[i].minibd = p[i].minibd;
	pOld[i].ibdcut = p[i].ibdcut;
	pOld[i].mindepth = p[i].mindepth;
	pOld[i].ibdnum = p[i].ibdnum;
	for( j = 0; j < fineL; j++){
	    mOld[i][j].p = m[i][j].p;
	    mOld[i][j].r = m[i][j].r;
	    mOld[i][j].g = m[i][j].g;
	    mOld[i][j].midibdcm = m[i][j].midibdcm;
	    mOld[i][j].minibdcm1 = m[i][j].minibdcm1;
	    mOld[i][j].minibdcm2 = m[i][j].minibdcm2;

	}
    }
}

void add_map(struct map m[megL][fineL], struct mapPar p[megL], struct map mOld[megL][fineL], struct mapPar pOld[megL])
{
    int i, j;
    for( i = 0; i < megL; i++){
	pOld[i].L = p[i].L;
	pOld[i].bin = p[i].bin;
	pOld[i].minibd = p[i].minibd;
	pOld[i].ibdcut = p[i].ibdcut;
	pOld[i].mindepth = p[i].mindepth;
	pOld[i].ibdnum = p[i].ibdnum;
	for( j = 0; j < fineL; j++){
	    mOld[i][j].p = m[i][j].p;
	    mOld[i][j].midibdcm = m[i][j].midibdcm;
	    mOld[i][j].minibdcm1 = m[i][j].minibdcm1;
	    mOld[i][j].minibdcm2 = m[i][j].minibdcm2;
	    mOld[i][j].r = (mOld[i][j].r + m[i][j].r)/2;
	    mOld[i][j].g = (mOld[i][j].g + m[i][j].g)/2;
	}
    }
}



void print_map(struct map m[], struct mapPar p)
{
    int i;
    printf("Pos rate GenPos\n");
    for ( i = 0; i < p.L; i++){
	printf("%.0lf %lf %lf\n", m[i].p + POSdrift, m[i].r*RECrate, m[i].g*RECrate);
    }
}

void print_plinkmap(struct map m[], struct mapPar p)
{
    int i;
    for ( i = 0; i < p.L; i++){
        printf("%s . %lf %.0lf\n", CHR, m[i].g*RECrate, m[i].p + POSdrift);
    }    
}



void print_fine_map(struct map m[megL][fineL], struct mapPar p[])
{
    int i, j;
    int from, to;
    from = 0; to = megL;
    if(strcmp(CHR, ".")!=0)printf("Chr Pos_bp rate_cM.Mp Pos_cM");
    else printf("Pos_bp rate_cM.Mp Pos_cM");
    if(STAT == 1)printf(" midIBD minUSEIBD1 minALLIBD2");
    printf("\n");
    for ( i = from; i < to; i++){
	for (j = 0; j < fineL; j++){
	    if(strcmp(CHR, ".")!=0)printf("%s %.0lf %lf %lf",CHR, m[i][j].p + POSdrift, m[i][j].r*RECrate, m[i][j].g*RECrate);
	    else printf("%.0lf %lf %lf",m[i][j].p + POSdrift, m[i][j].r*RECrate, m[i][j].g*RECrate);
	    if(STAT == 1)printf(" %lf %lf %lf", m[i][j].midibdcm*RECrate, m[i][j].minibdcm1*RECrate, m[i][j].minibdcm2*RECrate);
	    printf("\n");
	}
    }
}

void print_fine_plinkmap(struct map m[megL][fineL], struct mapPar p[])
{
    int i, j;
    int from, to;
    from = 0; to = megL;
    //printf("Pos rate GenPos");
    if(STAT == 1)printf(" midIBD minUSEIBD1 minALLIBD2");
    printf("\n");
    for ( i = from; i < to; i++){
        for (j = 0; j < fineL; j++){
            printf("%s . %lf %.0lf",CHR, m[i][j].g*RECrate, m[i][j].p + POSdrift);
            printf("\n");
        }
    }
    i = to -1; j = fineL -1; 
    printf("%s . %lf %.0lf",CHR, (m[i][j].g + m[i][j].r * (m[i][j].p - m[i][j-1].p) /1e6)*RECrate, m[i][j].p * 2 - m[i][j-1].p + POSdrift);
    printf("\n");

}


void free_array(void *p){free(p);p=NULL;}

void free_matrix_d(int **p, int length){
    int i;
    for(i=0;i<length;i++){free(p[i]);p[i]=NULL;}
    free(p);p=NULL;
}

void free_matrix_lf(double **p, int length){
    int i;
    for(i=0;i<length;i++){free(p[i]);p[i]=NULL;}
    free(p);p=NULL;
}
