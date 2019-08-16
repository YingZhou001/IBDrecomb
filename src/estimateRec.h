void meg_Rec(int **sourceIBD, int from, int to, double bin, struct map megmap[], struct mapPar *par);
void fineRec(int **sourceIBD, int I, double G, struct map mOld[megL][fineL], struct mapPar pOld[megL], struct map m[megL][fineL], struct mapPar p[megL]);
void fineSideRec(int **sourceIBD, int I, double G, struct map mOld[megL][fineL], struct mapPar pOld[megL], struct map m[megL][fineL], struct mapPar p[megL]);
