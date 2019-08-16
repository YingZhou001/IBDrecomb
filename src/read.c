#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

#include "head.h"
#include "tools.h"

void read_source_ibd(char *ibdinp, int **sourceIBD)
{
    clock_t checkpt, begin = clock();
    double a, b;
    FILE *ifp;
    int gztag, i;
    char chr[32];
    gztag=check_gz(ibdinp);
    if(gztag==1){char cmd[10240];sprintf(cmd, "zcat %s", ibdinp);ifp=popen(cmd, "r");}
    else {ifp=fopen(ibdinp, "r");}
    i=0;
    while(!feof(ifp)){
	if(ibd)fscanf(ifp, "%lf %lf", &a, &b);
	if(simibd)fscanf(ifp, "%*s %*s %*s %*s %*s %lf %lf %*s", &a, &b);
	if(refinedibd){
	    fscanf(ifp, "%*s %*s %*s %*s %s %lf %lf %*s %*s", chr, &a, &b);
	    if(strcmp(chr, CHR)!=0)continue;
	}
	sourceIBD[i][0]=(int)floor(a)-POSdrift;
	sourceIBD[i][1]=(int)floor(b)-POSdrift;
	i++;
    }
    if(gztag==1)pclose(ifp);
    else fclose(ifp);
    checkpt = clock();
    if(mute==0)fprintf(stderr,"\tread_source_ibd time:%lfs\n", (checkpt-begin)*1.0/CLOCKS_PER_SEC);
}

int read_ibd(int **sourceIBD, int **end, struct map m[], struct mapPar p)
{
    clock_t checkpt, begin = clock();
    int ibdnum;
    int j, i=0;


    for(j=0;j<IBDNUMtot;j++){
	end[i][1]=sourceIBD[j][0];end[i][2]=sourceIBD[j][1];
	end[i][3]=(int)floor(fastphy2gen(end[i][1], m, p)*1e6);
	end[i][4]=(int)floor(fastphy2gen(end[i][2], m, p)*1e6);
	end[i][0]=end[i][4]-end[i][3];
	if(end[i][0]>p.ibdcut*LENthreshold*1e6){i++;}
    }
    ibdnum = i;
    checkpt=clock();
    if(mute == 0)fprintf(stderr, "\tIBDCUT = %lf, %d IBDs pass the filtration\n",p.ibdcut*LENthreshold, ibdnum);
    if(mute == 0)fprintf(stderr,"\tread_ibd time:%lfs\n", (checkpt-begin)*1.0/CLOCKS_PER_SEC);
    return(ibdnum);
}

int read_local_ibd(int **sourceIBD, int **end, int I, struct map m[megL][fineL], struct mapPar p[megL])
{
    clock_t checkpt, begin = clock();
    double from, to;
    double ovlap;
    double tmp1, tmp2;
    double tmpout=999999999;
    int ibdnum;
    int j, i=0;
    int k=0;
    from = m[I][0].p; to = m[I][fineL-1].p + p[I].bin;
    for(j=0;j<IBDNUMtot;j++){
	end[i][1]=sourceIBD[j][0];end[i][2]=sourceIBD[j][1];
	if((end[i][2]-end[i][1])>fBIN || end[i][2]>end[i][1]){
	    ovlap = cal_overlap_rate(from, to, end[i][1], end[i][2], m[1][0].p-m[0][0].p);
	    tmp1=0;tmp2=0;
	    if(ovlap > 1e-7){
		tmp1=finephy2gen(end[i][1], m, p);
		tmp2=finephy2gen(end[i][2], m, p);
		end[i][3]=(int)floor(tmp1*1e6);
		end[i][4]=(int)floor(tmp2*1e6);
		end[i][0]=end[i][4]-end[i][3];
		if(end[i][0] < tmpout)tmpout = end[i][0];
		//if(end[i][0]<0)
		//fprintf(stderr, "joe--%d %d %d %lf %lf %lf\n",end[i][1], end[i][2], end[i][0], tmp1, tmp2, ovlap);
		//if(end[i][0]>p[I].ibdcut*1e6){i++;}
		i++;
	    }
	}
    }
    ibdnum = i;
    checkpt=clock();
    if(mute == 0)fprintf(stderr, "\t%d-%d IBDs pass the filtration %lf\n",I, ibdnum, tmpout/1e6);
    return(ibdnum);
}
