
double fastphy2gen(double x0, struct map m[], struct mapPar p);
double finephy2gen(double x0, struct map m[megL][fineL], struct mapPar p[megL]);

void constoutmap(struct map m[], struct mapPar p);

void update_TRIM(struct map m[], struct mapPar p);
double update_IBDCUT(struct map m[], struct mapPar p);
void sort_ibd(int **end, int step, int ibdnum);
double cal_bin_depth(int **end, double *coverage, struct map m[], struct mapPar p);
double cal_overlap_rate(double a, double b, int E_l, int E_r, double bin);

void cal_end_region_legnth(int **end, int ibd_num);

int end_counts_in_bin(double *X, int **end, int **sourceIBD, struct mapPar p);
void update_map(double *X, struct map m[], struct mapPar p);
void init_fine_map(struct map outfineMap[megL][fineL], struct mapPar outfineMapPar[megL], struct map outMegMap[megL], struct mapPar outMegMapPar);
void cp_map(struct map mm[megL][fineL], struct mapPar p[megL], struct map mOldm[megL][fineL], struct mapPar pOld[megL]);
void add_map(struct map mm[megL][fineL], struct mapPar p[megL], struct map mOldm[megL][fineL], struct mapPar pOld[megL]);

double cal_fine_bin_depth(int **end, double *coverage, struct map m[], struct mapPar p);
int end_counts_in_fine_bin(double *X, int **end, int **sourceIBD, struct map m[], struct mapPar p);
void update_fine_map(double *X, struct map m[], struct mapPar p, double G);



void print_map(struct map m[], struct mapPar p);
void print_plinkmap(struct map m[], struct mapPar p);
void print_fine_map(struct map m[megL][fineL], struct mapPar p[]);
void print_fine_plinkmap(struct map m[megL][fineL], struct mapPar p[]);


void free_array(void *p);
void free_matrix_d(int **p, int length);
void free_matrix_lf(double **p, int length);
