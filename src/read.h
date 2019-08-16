
void read_source_ibd(char *ibdinp, int **sourceIBD);
int read_ibd(int **sourceIBD, int **end, struct map m[], struct mapPar p);
int read_local_ibd(int **sourceIBD, int **end, int I, struct map m[megL][fineL], struct mapPar p[megL]);
