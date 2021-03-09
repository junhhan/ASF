#define MTn 624
#define MTm 397
#define MATRIX_A 0x9908b0dfUL
#define UMASK 0x80000000UL
#define LMASK 0x7fffffffUL
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

void _MTinit(unsigned long MTstate[], int *MTleft, int *MTinitf, unsigned long seed) {
	int i;
	MTstate[0]= seed & 0xffffffffUL;
	for (i= 1; i < MTn; ++i) {
		MTstate[i] = (1812433253UL * (MTstate[i-1] ^ (MTstate[i-1] >> 30)) + i); 
		MTstate[i] &= 0xffffffffUL;
	}
	(*MTleft)= 1; (*MTinitf)= 1;
}

void _MTnext(unsigned long MTstate[], unsigned long MTnext[], int *MTleft , int *MTcnt) {
	unsigned long *p= MTstate;
	
	(*MTcnt)= 0;
	(*MTleft)= MTn;
	
	int i;
	for (i= MTn - MTm + 1; --i; ++p) {
		*p= p[MTm] ^ TWIST(p[0], p[1]);
	} 
	
	for (i= MTm; --i; ++p) {
		*p = p[MTm - MTn] ^ TWIST(p[0], p[1]);
	} 
	
	*p= p[MTm - MTn] ^ TWIST(p[0], MTstate[0]);
	
	for (i= 0; i < MTn; ++i) {
		MTnext[i]= MTstate[i];
	}
}

double _Rand(unsigned long MTstate[], unsigned long MTnext[], int *MTleft , int *MTcnt) {
    unsigned long y;

    (*MTleft) -=1;
	if ((*MTleft) == 0) {_MTnext(MTstate, MTnext, MTleft, MTcnt);}
    y= MTnext[(*MTcnt)];
    (*MTcnt) += 1;

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return ((double)y * (1.0/4294967295.0));
}
