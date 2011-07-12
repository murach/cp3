/*
   frng.h                                       B Bunk 12/2010
*/
#ifdef LONG64
#define INT64 long
#else
#define INT64 long long
#endif

float frng();
double dfrng();
int ifrng();
INT64 lfrng();

void frngv(int n, float *rvec);
void dfrngv(int n, double *dvec);
void ifrngv(int n, int *ivec);
void lfrngv(int n, INT64 *lvec);

void frnget(INT64 *iseed);
void frnset(INT64 *iseed);
void frnini(int num);
