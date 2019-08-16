int adjust_end_counts(double *X, int **sourceIBD, struct map m[], struct mapPar p);
double count_left_end(int b1, double bcoverage, int **sourceIBD, int sourceIBD_L, struct map m[], double ibdmax, struct mapPar p);
double count_right_end(int b1, double bcoverage, int **sourceIBD, int sourceIBD_L, struct map m[], double ibdmax, struct mapPar p);
double cal_left_coverage(int bin, int **sourceIBD, int sourceIBD_L, struct map m[], double ibdmax, struct mapPar p);
double cal_right_coverage(int bin, int **sourceIBD, int sourceIBD_L, struct map m[], double ibdmax, struct mapPar p);

double cal_side_fine_bin_depth(int I, int **end, double *coverage, struct map m[], struct mapPar p);
int end_counts_in_side_fine_bin(int I, double *X, int **end, struct map m[megL][fineL], struct mapPar p[megL]);
