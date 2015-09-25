#include <stdio.h>

/* Period parameters */  
#define MERSENNE_N 624
#define MERSENNE_M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)


static unsigned long mt[MERSENNE_N]; /* the array for the state vector  */
static int mti=MERSENNE_N+1; /* mti==MERSENNE_N+1 means mt[MERSENNE_N] is not initialized */

/* initializing the array with NONZERO seed */
void setmt19937( unsigned long seed)
{
    /* setting initial seeds to mt[MERSENNE_N] using         */
    /* the generator Line 25 of Table 1 in          */
    /* [KNUTH 1981, The Art of Computer Programming */
    /*    Vol. 2 (2nd Ed.), pp102]                  */
    mt[0]= seed & 0xffffffff;
    for (mti=1; mti<MERSENNE_N; mti++)
        mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

double mt19937()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= MERSENNE_N) { /* generate MERSENNE_N words at one time */
        int kk;

        if (mti == MERSENNE_N+1)   /* if sgenrand() has not been called, */
            setmt19937(4357); /* a default initial seed is used   */

        for (kk=0;kk<MERSENNE_N-MERSENNE_M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MERSENNE_M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<MERSENNE_N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MERSENNE_M-MERSENNE_N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[MERSENNE_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[MERSENNE_N-1] = mt[MERSENNE_M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return ( (double)y * 2.3283064365386963e-10 );
}

/* -- added routines for checkpointing -- h. katzgraber -------------- */
#include <stdlib.h>
void mt_restore(char *fn)
{
    int i,icnt;
    FILE *fp;
    fp = fopen(fn,"r");
    if(fp == NULL) {
	printf("\n *** Cannot open %s for write\n",fn);
	exit(1);
    }
    icnt = fscanf(fp,"%d\n",&mti);
    if(icnt != 1){
	printf("\n *** Read %d indices, not 1 \n",icnt);
	exit(1);
    }
    for(i = 0; i < MERSENNE_N; i++){
	fscanf(fp,"%lu\n",&mt[i]);
    }
    fclose(fp);
}

void mt_save(char *fn)
{
    int i;
    FILE *fp;
    fp = fopen(fn,"w");
    if(fp == NULL) {
	printf("\n *** Cannot open %s for write\n",fn);
	exit(1);
    }
    fprintf(fp,"%d\n",mti);
    for(i=0;i<MERSENNE_N;i++){
	fprintf(fp,"%lu\n",mt[i]);
    }
    fclose(fp);
}
/* -- end added routines --------------------------------------------- */


unsigned long imt19937()
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    
    if (mti >= MERSENNE_N) { /* generate MERSENNE_N words at one time */
	int kk;
	
	if (mti == MERSENNE_N+1)   /* if sgenrand() has not been called, */
	    setmt19937(4357); /* a default initial seed is used   */
	
	for (kk=0;kk<MERSENNE_N-MERSENNE_M;kk++) {
	    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	    mt[kk] = mt[kk+MERSENNE_M] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	for (;kk<MERSENNE_N-1;kk++) {
	    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	    mt[kk] = mt[kk+(MERSENNE_M-MERSENNE_N)] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	y = (mt[MERSENNE_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
	mt[MERSENNE_N-1] = mt[MERSENNE_M-1] ^ (y >> 1) ^ mag01[y & 0x1];
	
	mti = 0;
    }
    
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    
    return y; 
}

unsigned long imt19937range(int imin, int imax)
{
    unsigned long y, tmpval;
    double range;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */
    
    if (mti >= MERSENNE_N) { /* generate MERSENNE_N words at one time */
	int kk;
	
	if (mti == MERSENNE_N+1)   /* if sgenrand() has not been called, */
	    setmt19937(4357); /* a default initial seed is used   */
	
	for (kk=0;kk<MERSENNE_N-MERSENNE_M;kk++) {
	    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	    mt[kk] = mt[kk+MERSENNE_M] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	for (;kk<MERSENNE_N-1;kk++) {
	    y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
	    mt[kk] = mt[kk+(MERSENNE_M-MERSENNE_N)] ^ (y >> 1) ^ mag01[y & 0x1];
	}
	y = (mt[MERSENNE_N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
	mt[MERSENNE_N-1] = mt[MERSENNE_M-1] ^ (y >> 1) ^ mag01[y & 0x1];
	
	mti = 0;
    }
    
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);
    
    range = (imax - imin) * 2.3283064365386963e-10;
    tmpval = imin + (int)( y * range);
    if (tmpval > imax) tmpval = imax;
    
    return tmpval;
}
