/*
Copyright (c) 2020 deletecode

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Also, let me know if you use this algorithm or the code.

*/


#include <stdlib.h>
#include "time.h"
#include "util.h"
#include "conio.h"
#include "math.h"
#include "stdio.h"
#include "windows.h"

#include <cassert>

/*
This is a working demo of an idea I call skip chaining, which is described before the class SkipChaining.

The code might not be used by other projects (most of which are in Java), but this is more about the algorithm.
This can also serve as a benchmark to compare Java.  For the specific case of finding 9 mushroom island chunks, this scans a seed space at about 2.5 exaseeds/second on 1 cpu.

Note that this may not actually generate because there is another check that needs to pass for mushroom.

As this is a demo, it's not fully set up to be modular, since it's not known how (or if) it'll be used.
*/

// Tested in VS 2005 and g++ 4.8.3-2 (as part of win-builds)

// For g++ use this:
// g++ skipchains.cpp -O1
// Don't use -O2 or higher.  Something in the math is wrong.





// Some options for ease of use for others for demoing.

// Print each step as it traverses through the seed space.
const int config_print_skipchain_details=false;

// Check mod this ==0 (100 for mushroom islands, 10 to simulate MapIsland, but it should handle any mod).
// Mod should be 3 or more, and *not* a power of 2
// The speed will vary tremendously depending on the mod.

// Modulus to use for all chunks
const int config_mod=100;
// Number of chunks to test.  If using a smaller mod for testing, you want this higher.
const int config_numspots=9;
#define config_SEEDFUNC GetShroomChunkSeed

#define config_staticmod 1


// example for finding a lot of land in MapIsland 
// more for testing speed for a different modulus
//const int config_mod=10;
//const int config_numspots=15;
//#define config_SEEDFUNC GetMapIslandSeed

const int config_numskiparrays=config_numspots-1; 
const int config_do_inv_space=1;

// With more skip arrays, this may need to be higher
// But it is allocated on the stack.
const int config_skiparraysize=1000;  

const int config_print_seeds=1;

// Uncomment for linking/including cuBiomes code
//#define CHECKCUBIOMESGEN




#ifdef __GNUC__

#include "stdint.h"
#define __int64 int64_t

#endif



// I wanted to keep this stand-alone since the primary algorithm doesn't depend on it
// 
#ifdef CHECKCUBIOMESGEN

#include "finders.h"
#include "generator.h"

int CheckFullGenForMushrooms(__int64 seed,int baseX,int baseY,int sizeX,int sizeY) {
   	initBiomes();
    // Allocate and initialize a stack of biome layers.
    LayerStack g = setupGenerator(MC_1_16);
    // Extract the desired layer.
    Layer *layer = &g.layers[L_ADD_MUSHROOM_256];
	
    // Allocate a sufficient buffer for the biomes and for the image pixels.
    int *biomes = allocCache(layer, sizeX, sizeY);
    // Apply the seed only for the required layers and generate the area.
    setWorldSeed(layer, seed);
	genArea(layer, biomes, baseX, baseY, sizeX, sizeY);

	int count=0;
	for(int i=0;i<sizeX*sizeY;i++) {
		if(biomes[i]==mushroom_fields)
			count++;
	}

	// Clean up.
    freeGenerator(g);
    free(biomes);
	return count;
    
}    

#else

// Taken from cubiomes
static inline __int64 processWorldSeed(register __int64 ws, const __int64 bs)
{
    ws *= ws * 6364136223846793005LL + 1442695040888963407LL;
    ws += bs;
    ws *= ws * 6364136223846793005LL + 1442695040888963407LL;
    ws += bs;
    ws *= ws * 6364136223846793005LL + 1442695040888963407LL;
    ws += bs;
    ws *= ws * 6364136223846793005LL + 1442695040888963407LL;
    return ws;
}

static inline __int64 getChunkSeed(register __int64 ss, const __int64 x, const __int64 z)
{
    ss += x;
    ss *= ss * 6364136223846793005LL + 1442695040888963407LL;
    ss += z;
    ss *= ss * 6364136223846793005LL + 1442695040888963407LL;
    ss += x;
    ss *= ss * 6364136223846793005LL + 1442695040888963407LL;
    ss += z;
    return ss;
}
#endif


// Convenience function for returning time in seconds using high precision counter
inline double Time()
{
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	
	LARGE_INTEGER time;
	QueryPerformanceCounter(&time);

	static LARGE_INTEGER lasttime;
	if(lasttime.QuadPart==0)
		QueryPerformanceCounter(&lasttime);

	return (time.QuadPart-lasttime.QuadPart)/(double)freq.QuadPart;
}




inline __int64 GetShroomChunkSeed(__int64 seed,int x, int z) {
	__int64 ss=processWorldSeed(seed,-7479281634960481323LL);
	return getChunkSeed(ss, x, z);
}

// For MapIsland
inline __int64 GetMapIslandSeed(__int64 ws, int x, int z) {
	__int64 shroom_layer_seed=3107951898966440229LL;
	ws *= ws * 6364136223846793005LL + 1442695040888963407LL;
    ws += shroom_layer_seed;
    ws *= ws * 6364136223846793005LL + 1442695040888963407LL;
    ws += shroom_layer_seed;
    ws *= ws * 6364136223846793005LL + 1442695040888963407LL;
    ws += shroom_layer_seed;

	ws *= ws * 6364136223846793005LL + 1442695040888963407LL;
	
	__int64 ss=ws;
	ss += x;
    ss *= ss * 6364136223846793005LL + 1442695040888963407LL;
    ss += z;
    ss *= ss * 6364136223846793005LL + 1442695040888963407LL;
    ss += x;
    ss *= ss * 6364136223846793005LL + 1442695040888963407LL;
    ss += z;
	return ss;
}


// modular inverse of a, mod b
// taken off rosettacode
// credit really goes to Euler
__int64 ModInverse(__int64 a, __int64 b)
{
	if(a<0)
		a+=b;
	__int64 b0 = b, t, q;
	__int64 x0 = 0, x1 = 1;
	if (b == 1) return 1;
	int steps=0;
	while (a > 1) {
		if(b==0)
			return -1;
		q = a / b;
		t = b, b = a % b, a = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
		steps++;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}

__int64 GCF(__int64 x,__int64 y) {
	__int64 a=x;
	__int64 b=y;

	while (b != 0)
	{
		__int64 t = b;
		b = a % b;
		a = t;
	}

	__int64 gcd = a;
	__int64 lcm = (x*y)/gcd;
	return gcd;
}

// testing/experimentation
void MathTestingModInv() {
	for(int i=0;i<9256;i++) {
		int fwd=i;
		int inv=(int)ModInverse(i,1LL<<32);
		printf("%d %d\n",fwd,inv);
		if(inv!=-1) {
			int t=3;
			printf("%d ",(int)t);
			for(int k=0;k<3;k++) {
				t*=fwd;
				printf("%d ",(int)t);
			}
			for(int k=0;k<3;k++) {
				t*=inv;
				printf("%d ",(int)t);
			}
			printf("\n");
			int g=1;
		}
	}
}
// testing/experimentation
void MathTestingModInv2() {
	int bs=5; 
	__int64 seed0=GetMapIslandSeed(bs+(1LL<<32)*0,0,0);
	__int64 seed1=GetMapIslandSeed(bs+(1LL<<32)*1,0,0);
	int diff=(seed1-seed0)>>32;
	int diff_inv=ModInverse(diff,1LL<<32);
	printf("diff: %d inv: %d\n",diff,diff_inv);

	for(int k=0;k<100;k++) {
		int j=k*diff_inv*5;
		__int64 g=GetMapIslandSeed(bs+(1LL<<32)*j,0,0);
		//int highg=g>>32;
		printf("%lld %lld %lld\n",g,g>>32,(g>>24)%10);
	}
}


// The upper bits of a output seed are linearly-dependent on the upper bits of the input seed, when the lower bits are locked.
// This exploits that fact.
class SeedAcclerator2 {
private:
	__int64 base_cs;

public:
	int interval;
	void Init(__int64 seed0,__int64 seed1) {
		base_cs =  seed0;
		interval = (seed1-seed0)>>32;
	}
	__int64 GetCS(int j) {
		// Uses a 32 bit multiply instead of 64 bit
		// Specific to little-endian.  On big endian switch 1 to 0
		__int64 ret=base_cs;
		((int*)(&ret))[1]+=interval*j;
		return ret;
	}
};


int MaxLift(int mod) {
	int bit=1;
	while(bit!=256) {
		if(mod&bit)
			return bit;
		bit<<=1;
	};
	return bit;
}


// This is like SeedAcclerator2, but the inverse.
// It allows stepping through a seed space for *only* passing seeds.
// It only works for a non power of 2 mod, which shouldn't be passed in in the first place.

// It uses the modular inverse of the difference of the upper bits of two seeds 2^32 apart, to explicitly control the upper bits of the output seed
// Using Get(i) on this function returns upper bits of a seed, which when passed into e.g. GetShroomChunkSeed, will have upper bits
// that correspond linearly to i.
// It does the added steps of chopping up the seed space and finding the first seed.

// Note the mod passed into this is values like 10,100,5,3 etc.  The modular inverse space is respect to 2^32.
// This takes both into account.

// (Credit to Matthew Bolan to showing me this could be done, though I had to figure out how)

class ModInverseSpace {
private:
	int diff_inv;  // modular inverse of difference between 2 elements when adding 2^32
	int shift;     // a shift such that diff*Get(0) returns a value near INT_MIN that passes the mod test
	int chop;      // mod/(highest power of 2 that divides mod)  (if this equals 1, this class will not help speed)
	void CalcChop(const int mod) {
		// hopefully this gets precomputed.  otherwise, template parameters can be used.
		//_BitScanReverse()
		int bit=1;
		while(bit!=256) {
			if(mod&bit) {
				chop=mod/bit;
				return;
			}
			bit<<=1;
		};
		chop=mod/bit;
	}

public:
	const int GetChop() {
		return chop;
	}
	// slower method
	//2687
	/*bool Init(__int64 seed0, __int64 seed1, int mod) {
		int diff = (seed1-seed0)>>32;  // calculate the differential
		diff_inv = ModInverse(diff, 1LL<<32); // get the inverse
		shift = -(seed0>>32) + INT_MIN;  // shift back to the start
		CalcChop(mod);

		for(int i=0;i<chop;i++) {  // the 25 would be (mod/(highest power of 2 that divides mod)) (has to do with the hensel lifting stuff)
			                       // for mushroom that's 100/4 (25), for MapIsland it's 10/2 (5)
			
			int m = ((seed0+(diff*Get(0)*(1LL<<32)))>>24) % mod;  
			if(m == 0) {
				break;
			}
			shift++;
		}
		return true;
	}
	// could probably be made into one multiply and one add, and might be worth it.
	int Get(int k) {
		return (k*chop+shift)*diff_inv;
	}*/
	
	//2748
	bool Init(__int64 seed0, __int64 seed1, int mod) {
		int diff = (seed1-seed0)>>32;  // calculate the differential
		diff_inv = ModInverse(diff, 1LL<<32); // get the inverse
		CalcChop(mod);
		shift=(-(seed0>>32) + INT_MIN)*diff_inv;

		for(int i=0;i<chop;i++) {
			int m = ((seed0+(diff*Get(0)*(1LL<<32)))>>24) % mod;  
			if(m == 0) {
				break;
			}
			shift+=diff_inv;
		}
		diff_inv*=chop;
		return true;
	}
	// faster calculation
	int Get(int k) {
		return shift+diff_inv*k;
	}
};
 
inline int FloorMod(__int64 a, int m,int test) {
	int ret=a%m;
	if(test==0)
		return ret==0;
	if(ret<0)
		ret+=m;
	return ret==test;
}

// A basic utility array storing distances and sorting.
struct Skip2 {
	int distance;
	//int count;
};
class SkipArrayBasic {
public:
	static const int maxcount=config_skiparraysize;
	Skip2 a[maxcount];
	int count;
	SkipArrayBasic() {
		count=0;
		memset(a,0,sizeof(a));
	}
	void Add(int distance) {
		// Make sure it isn't already added
		bool found=0;
		for(int i=0;i<count;i++) {
			if(a[i].distance==distance) {
				found=1;
				break;
			}
		}
		// Add, if it isn't.
		if(!found) {
			if(count<maxcount) {
				//printf("New link: %d\n",distance);
				a[count].distance=distance;
				count++;
				Sort();
			} else {
				printf("Error: too many in skip array ");
			}
		}
	}

	void Init() {
		count=0;
	}
	void Sort() {
		// a dumb sort, cause it won't happen often
		for(int i=0;i<count;i++) {
			for(int j=i+1;j<count;j++) {
				if(a[i].distance>a[j].distance)  // sort lowest to highest so all (or most) seeds are found
				//	if(a[i].count<a[j].count)   // sort by frequency, faster but misses like 30% of them.
				{
					Skip2 tmp=a[i];
					a[i]=a[j];
					a[j]=tmp;
				}
			}
		}
	}
};

/*
I call this multi-level skip chaining.  It's remotely similar to a skip-list, in the sense that it jumps between "towers".
It exploits some mathematical properties of the RNG used in minecraft's map/biome generation.
It's approximately 70,000 times faster than brute forcing through the upper bits, in this specific case.
(for mushroom seed finding, requiring 9 values return true - this example is used in this explanation because it allows for such extreme speedups),
and it's fairly accurate, finding like 80-90%.


General idea:

It works because iterating only the upper 32 bits of an input seed results in regular patterns in a typical mod test,
and in fact it is highly exploitable - the upper bits in fact become linear.
For instance, the expression T(u)=GetShroomChunkSeed(seed+(1LL<<32)*u,X,0)>>24%100==0 would only pass at distances in the upper bits of ONLY 5, 19, 60 (for instance; these values depend on seed)
For example, if T(2) passes, then one or more of T(2+5), T(2+19), T(2+60) *must also* pass.

But these regular distances also appear when requiring multiple X to pass the test.  A bit like looking at interference pattern.




Math background:

Mathematically, this is similar to how "Hensel lifting" allows to lock output bits below a certain bit, by locking the respective input bits.
This is a special case, locking the lower 32 bits (half or more of the total bits).  The upper bits of the output
become linearly-dependent on the the upper bits of the input. (SeedAcclerator2 exploits this for quick calculation, for example)
This is because the non-linear component drops completely out of the equation.
(Side note: if you wish to confirm it's linearly-dependent, solve the function e.g. GetShroomChunkSeed in terms
 of upper and lower bits. you'll see a squared component drop out because it's multiplied by 2^64 then implicitly modded by 2^64)

The linear dependency also allows explicitly iterating over *passing* values for one X value.  This is used in ModInverseSpace, which
uses the modular inverse of this linear dependency.  This gives a speedup by shrinking the space 25x, but since this
algorithm is fairly efficient to begin with ( O(log2(seed space)) ?), the actual speed gains are only about 10x.

Because of these math facts, this works mostly with regular 32 bit int values (upper bits).
It still has to do a 40 bit modulus (done as a 64 bit), and that may be unavoidable.


The algorithm:

Now the algorithm.  It stores "skip arrays" at multiple levels.
skiparray[0] is the list of possible distances between where X=0 passes 
skiparray[1] is the list of possible distances between where X=0 AND X=1 pass
skiparray[2] is the list of possible distances between where X=0 AND X=1 AND X=2 pass
etc.



What the algorithm effectively does:
The first thing it will do is populate skiparray[0], by basically brute forcing until it finds the pattern (as mentioned before, this is always 3 distances).
Then it will utilize skiparray[0] to populate skiparray[1]
Then it will utilize skiparray[1] to populate skiparray[2]
and so on.
If a higher skip array fails to find the next match, it will fall back on
lower skip arrays, until it can find a new distance for the next higher one.


It might seem excessively complex to build the skip arrays like this, rather than just store a single one of them (this was my previous attempt), but it actually 
gives a speedup of about 200x over that attempt.  This is because the 32 bit seed space being searched turns out to be kinda small, so
the effort to build up the skip array was actually taking most of the time, so it's far better to use the accelerator to build a faster accelerator.


Since most of the seed finding utilities are in Java, this might remain a c++ demo, unless.
I have tentative plans to put this on a GPU as well.


The speedup for 9 mushroom chunks goes as follows:
250,000x speedup by pre-checking base seeds and lifting (.04^9 vs .01^9)
70,000x speedup using this multi-level skip chaining method
4x speedup (guess) by using SeedAccelerator2
10x speedup using ModInverseSpace
That's a 700,000,000,000x speedup over a pure brute force approach.


Big thanks to Matthew Bolan for teaching me a lot about this modular arithmetic stuff and Hensel lifting.
I don't know if this has been invented before, but if it has been, credit to them too.
*/


class SkipChaining {
public:
	// Number to test. With a smaller mod, you probably want this higher.
	static const int numspots=config_numspots; 
	// Make sure this is less than numspots if using do_inv_space is true
	static const int numskiparrays=config_numskiparrays; 
	// Set to 1 to use ModInverseSpace to shrink seed space.
	static const int do_inv_space=config_do_inv_space;
	// Modulus.  100 would be for mushroom chunks.  change it to estimate the speedup this method would give.
#ifdef config_staticmod
	static const int mod=config_mod;  
#else
	int mod;
#endif
	
		
	SkipChaining() {
#ifndef config_staticmod
		mod=0;
#endif
	}
	SeedAcclerator2 accelerators[numspots];

	ModInverseSpace inv_space;

	int NumMatches(int testpos,int numtest) {
		if(do_inv_space)
			testpos=inv_space.Get(testpos);
		int matched;
		for(matched=0;matched<numtest;matched++) {
			if(!FloorMod(accelerators[matched].GetCS(testpos)>>24,mod,0))
			//if((accelerators[matched].GetCS(testpos)>>24)%mod)
				break;
		}
		return matched;
	}
	int Do(int *outputseeds,int *outputseeds_count,int outputseeds_count_max) {
		
		
		HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE); // For pretty color printing only.

		inv_space.Init(
			accelerators[numspots-1].GetCS(0),
			accelerators[numspots-1].GetCS(1),mod);

		
		double starttime=Time();

		// It's a bit faster to put these arrays on the stack.
		SkipArrayBasic skiparray[numskiparrays];
		//SkipArrayBasic *skiparray=new SkipArrayBasic[numskiparrays];


		int lastseen[numskiparrays];
		for(int i=0;i<numskiparrays;i++)
			lastseen[i]=-1;
		
		int iteration=0;

		int num_fullmatch=0;
		int j=0;
		int last_num_matched=0;

		while(true) {
			int num_matched=0;
			int arraynum;
			int testpos=-1;

			// Test the skip in the skiparrays from highest to lowest, starting from the match count from the last iteration.
			for(arraynum=last_num_matched-1;arraynum>=0;arraynum--) {
				SkipArrayBasic &s=skiparray[arraynum];
				
				// Hueristics to guess when a skiparray is has enough values to be considered accurate enough
				// I believe technically it might be 3^arraynum but I wouldn't be able to prove it.
				// If it's too high, it'll be slower than necessary.  If it's too low, it'll miss more values.
				if(s.count<
					3+arraynum*2
					//(2<<arraynum)
					//2*(int)(pow(2.1,arraynum))
					)
					continue;

				// Test the skips
				for(int i=0;i<s.count;i++) {
					testpos = j+s.a[i].distance;
					num_matched = NumMatches(testpos,numskiparrays);
					if(num_matched > arraynum) {
						goto foundmatch;
					}
				}
			}
			// All known skips failed.  Brute force it.  This should only happen at the beginning.
			for(int i=1;i<1000000;i++) {
				testpos = j+i;
				num_matched = NumMatches(testpos,numskiparrays);
				if(num_matched>=1)
					goto foundmatch;
			}
			printf("Failed to do basic iteration %d ",skiparray[0].count);
			break;
foundmatch:;

			if(do_inv_space) {
				if(testpos>(1LL<<32)/inv_space.GetChop()) break;  // test in mod inverse space
			} else {
				if(testpos>=0&&j<0) break;  // test full range
			}

			


			for(int i=arraynum+1;i<num_matched;i++) {
 				if(lastseen[i]!=-1)
					skiparray[i].Add(testpos-lastseen[i]);
			}

			for(int i=0;i<num_matched;i++)
				lastseen[i]=testpos;

			if(num_matched==numskiparrays) {
				if(NumMatches(testpos,numspots)==numspots) {
					if((*outputseeds_count)<outputseeds_count_max) {
						if(do_inv_space)
							outputseeds[(*outputseeds_count)++]=inv_space.Get(testpos);
						else
							outputseeds[(*outputseeds_count)++]=testpos;
					} else {
						printf("too many seeds");
					}
				}
				num_fullmatch++;
			}


			// This is just used to visualize/debug the process
			if(config_print_skipchain_details&&(iteration%1)==0) {
				for(int i=0;i<numskiparrays+1;i++) {
					if(i<arraynum+1)
						printf("|");
					else if(i<num_matched)
						printf("+");
					else
						printf(" ");
				}
				printf("%d "/*Match %d %d*/"%*d",
					testpos,
					/*num_matched,
					arraynum,*/
					8,j?(testpos-j):0);
				for(int i=0;i<numskiparrays;i++) {
					if(skiparray[i].count==0) {
						//printf(".");
						break;
					} else {
						printf(" (");
						for(int k=0;k<skiparray[i].count;k++) {
							if(k>0)
								printf(",");
								
							if(i<4&&k<5&&(skiparray[i].a[k].distance<500000||k<1)/*||true*/) {
								SetConsoleTextAttribute(hConsole, arraynum==i&&(skiparray[i].a[k].distance==testpos-j)?10:7);
								printf("%d",skiparray[i].a[k].distance);
								SetConsoleTextAttribute(hConsole, 7);
							} else {
								for(int s=k;s<skiparray[i].count;s++) 
									if(arraynum==i&&skiparray[i].a[s].distance==testpos-j)
										SetConsoleTextAttribute(hConsole, 10);
								printf("%d more",skiparray[i].count-k);
								SetConsoleTextAttribute(hConsole, 7);
								break;
							}
						}
						
						printf(")");
					}
				}
				printf("\n");

				//Sleep(10);
			}

			iteration++;

			last_num_matched=num_matched;
			j=testpos;
		};


		//printf("%d",skiparray[0].a[0].distance,skiparray[0].a[1].distance);
		/*for(int i=0;i<numskiparrays;i++)
			printf("%d, ",skiparray[i].count);
		printf("\n");*/

		//delete []skiparray;

		return num_fullmatch;
	}
};






void MathTestingModInv4() {
	int bs=58;
	int mod=100;
	for(;bs<200;bs++) {
		if((GetMapIslandSeed(bs,0,0)>>24)%4==0)
			break;
	}
	ModInverseSpace r;

	r.Init(
		GetMapIslandSeed(bs,0,0),
		GetMapIslandSeed(bs+(1LL<<32),0,0),mod);

	for(int k=0;k<(1LL<<32)/25;k+=1000) {
		int j=r.Get(k);
		int m=(GetMapIslandSeed(bs+j*(1LL<<32),0,0)>>24)%100;
		printf("%d",m);
	}


}

void MathTestingModInv3() {
	int bs=8;
	int mod=100;
	for(;bs<200;bs++) {
		if((GetMapIslandSeed(bs,0,0)>>24)%4==0)
			break;
	}

	int diff[2],diff_inv[2];
	for(int i=0;i<2;i++) {
		diff[i]=(
			GetMapIslandSeed(bs+(1LL<<32)*1,i,0)-
			GetMapIslandSeed(bs+(1LL<<32)*0,i,0))>>32;

		diff_inv[i]=ModInverse(diff[i],1LL<<32);
	}

	printf("diff_inv: %d\n",diff_inv[0]);
	//int gcf=GCF(diff_inv[0],diff_inv[1]);

	__int64 seed0=GetMapIslandSeed(bs+(1LL<<32)*0,0,0);
	int highpart0=seed0>>32;
	int shift=-highpart0+INT_MIN;//highpart0*diff[0];

	
	for(int i=0;i<25;i++) {
		int m0=(GetMapIslandSeed(bs+(shift*diff_inv[0])*(1LL<<32),0,0)>>24)%mod;
		if(m0<0)
			m0+=mod;
		if(m0<4)
			break;
		shift++;
	}
	printf("highpart: %d shift: %d\n",highpart0,shift);


	int matches=0;
	int notmatchcount=0;
	for(int k=0;k<(1LL<<32)/25;k+=1) {
		int j=(k*25+shift)*diff_inv[0]/**diff_inv[1]*1*/;
		__int64 g[2];
		g[0]=GetMapIslandSeed(bs+(1LL<<32)*j,0,0);
		g[1]=GetMapIslandSeed(bs+(1LL<<32)*j,4,0);
		//int highg=g>>32;
		//printf("%lld %lld %lld\n",g,g>>32,(g>>24)%10);

		int m=((g[1]>>24)%mod);
		if(m<0)
			m+=mod;
		bool match=m==0;
		notmatchcount++;
		if(match) {
			/*printf("%lld %lld %d ",g[0]>>32,(g[0]>>24)%mod,m);
		//if(match)
			printf("dist: %d",notmatchcount);
			printf("\n");*/
			matches++;
		}
		if(match)
			notmatchcount=0;
	}
	printf("Matches: %d\n",matches);
}

inline __int64 GetSeed(__int64 ws,int count) {
	__int64 shroom_layer_seed=-7479281634960481323LL;
	for(int i=0;i<count;i++) {
		ws *= ws * 6364136223846793005LL + 1442695040888963407LL;
		ws += shroom_layer_seed;
	}
	return ws;
}

void MathTesting() {
	printf("MathTesting\n");
	for(int k=0;k<33;k++) {
		printf("k: %d\n",k);
		for(int j=0;j<3;j++) {
			int bs=5;
			printf("%lld\n",
				GetSeed(bs+(1LL<<29)*(j+1),k)-
				GetSeed(bs+(1LL<<29)*(j+0),k));
		}
	}
}


void MathTestingPattern() {

	int bs=21;
	int diff=59294991;
		//a
	while(true) {
		__int64 seed0=GetSeed(bs,1); // just some arbitrary 64bit value

		for(int q=0;q<5;q++) {
			for(int k=00;k<190;k++) {
				__int64 calcseed=seed0;
				calcseed+=__int64(k*(diff+q*2))<<32;
				int res=((calcseed>>24)%10)==0;
				printf("%c",res?'1':' ');
				/*if(k==100)
					printf("%lld",calcseed>>24);*/
			}
			printf("\n");
		}
		
		int key=getch();
		diff+=((key=='d')-(key=='a'))*2;
		diff+=((key=='q')-(key=='e'))*200000;
		bs+=(key=='w')-(key=='s');
		
	};
}



// to aid in a test
class SkipArrayHelper : public SkipArrayBasic {
public:
	int lastmatch;
	SkipArrayHelper() {
		lastmatch=-1;
	}
	void AddEasy(int j) {
		if(lastmatch!=-1) {
			Add(j-lastmatch);
		}
		lastmatch=j;
	}
};



void TestPredictPatterns() {

	//printf("%d",GCF(170,92));
	const int numtest=2;
	int baseseed=0;
	while (1) {
		SkipArrayHelper sa[3];

		SeedAcclerator2 s[2];
		for(int x=0;x<numtest;x++)
			s[x].Init(GetShroomChunkSeed(baseseed,x,0),GetShroomChunkSeed(baseseed+(1LL<<32),x,0));

		int mods[2];
		mods[0]=9;
		mods[1]=11;

		for(int j=0;j<100000000;j++) {
			int passcount=0;
			for(int x=0;x<numtest;x++) {
				__int64 cs=s[x].GetCS(j);//GetShroomChunkSeed(baseseed+j*(1LL<<32),x,0);
				int pass=((cs>>24)%mods[x])==0;
				passcount+=pass;
				if(pass)
					sa[x].AddEasy(j);
			}
			if(passcount==2)
				sa[2].AddEasy(j);
		}
		printf("%d",GCF(s[0].interval,s[1].interval));
			
		for(int x=0;x<3;x++) {
			printf("\n%d: ",sa[x].count);
			for(int i=0;i<sa[x].count;i++)
				printf("%d ",sa[x].a[i].distance);
			printf("\n");
		}
		baseseed=getch();
	};

	
}




void UnitTests() {
	int baseseed=5;
	__int64 s0=GetShroomChunkSeed(baseseed+0*(1LL<<32),0,0);
	__int64 s1=GetShroomChunkSeed(baseseed+1*(1LL<<32),0,0);
	assert(s0==-7472044254813591856LL);
	assert(s1==5491871351961970384LL);

	ModInverseSpace mi;
	mi.Init(s0,s1,100);
	int j=mi.Get(0);
	assert(j==-1920341363);
	int u=GetShroomChunkSeed(baseseed+j*(1LL<<32),0,0)>>32;
	assert(u-INT_MIN<=25);

	assert((GetShroomChunkSeed(baseseed+mi.Get(1)*(1LL<<32),0,0)>>32)-INT_MIN<=50);
	assert((GetShroomChunkSeed(baseseed+mi.Get(2)*(1LL<<32),0,0)>>32)-INT_MIN<=75);
	assert((GetShroomChunkSeed(baseseed+mi.Get(3)*(1LL<<32),0,0)>>32)-INT_MIN>=75);

	int g=1;
}




void TestDistances() {
	int baseseed=15;
	while(1) {
		int lastmatch=-1;
		SkipArrayBasic sa;
		int totalfullmatches=0;

		int diffs[10];
			
		
		int numtest=2;
		for(int x=0;x<numtest;x++)
			diffs[x]=(GetShroomChunkSeed(baseseed,x,0)-GetShroomChunkSeed(baseseed+(1LL<<32),x,0))>>32;
		
			

		for(int j=0;j<1000000;j++) {
			int count=0;
			for(int x=0;x<numtest;x++) {
				__int64 cs=GetShroomChunkSeed(baseseed+j*(1LL<<32),x,0);
				count+=FloorMod(cs>>24,5,0);//&&cs>0;
			}
			if(j>100000&&totalfullmatches==0)
				break;
			int pass=count==numtest;
			totalfullmatches+=pass;
			if(j<10000) {
				if(pass) {
					if(lastmatch!=-1) {
						printf("%c",(j-lastmatch)+'0');
						sa.Add(j-lastmatch);
					}
					lastmatch=j;
				} else {
					printf(" ");
				}
			}
		}

		printf("\n%d: ",sa.count);
		for(int i=0;i<sa.count;i++)
			printf("%d ",sa.a[i].distance);
		printf("\n");
		int gcf=GCF(diffs[0],diffs[1]);
		//printf("gcf: %d",gcf);
		baseseed=getch();
	}
}



// This does the steps necessary to successfully use skip chaining, and benchmarking.
// First is preparing the base seeds (or more precisely, the lower bits), that will be lifted later.  If you aren't familiar with this,
// it's making the list of seeds which can possibly pass the full test - e.g. 9 of 9 output seeds can possibly all pass the %100==0 test.
// In the case of %100==0, check that %4==0 passes for each one, and do this for 0 to 2^26.  The lower 26 bits are locked, so when >>24 is applied before %100, the 2 lowest bits
// are predetermined to both be 0.
// Since this is a common technique, it's not going to be explained more

// Next it does a benchmark vs a brute force test on the upper bits.

// Next it calls the skipchain class

// I would 

double DoSkipChaining(int basex, int basey, int mod,int loud,int benchmark) {
	SkipChaining chain;

	int numspots=chain.numspots;
	//__int64 outputseeds_full[50000];
	int baseseeds[50000];
	int baseseeds_count=0;

	int baseseeds_total_count=0;

#ifndef config_staticmod
	chain.mod=mod;
#endif



	if(loud)
		printf("testing for %d chunks at %d,%d\n",numspots,basex,basey);

	
	int lift=MaxLift(chain.mod);
	if(loud)
		printf("mod: %d lift: %d\n",chain.mod,lift);
		
	// Technically you can add 2 here (doing only evens, because each seed has a sister seed which is the opposite oddness)
	// but prefer leaving that out.
	for(int bs=0;bs<(1<<26);bs+=1) {
		int x=0;
		for(;x<numspots;x++)
			if((config_SEEDFUNC(bs,basex+x,basey)>>24)%lift!=0)
				break;
		if(x==numspots) {
			baseseeds[baseseeds_count++]=bs;
		}
		// break early, more for testing
		if(baseseeds_count>=5000)
			break;
		if(benchmark&&baseseeds_count>5)
			break;
		/*if(bs>10000000)
			break;*/
	}
	if(loud)
		printf("%d base seeds found\n",baseseeds_count);


	// calculate brute force speed for comparison
	double bruteforcespeed=0;
	{
		double starttime=Time();
		int bs=baseseeds[0];
		for(int x=0;x<numspots;x++)
			chain.accelerators[x].Init(
				config_SEEDFUNC(bs          ,basex+x,basey),
				config_SEEDFUNC(bs+(1LL<<32),basex+x,basey));
			
		int num=1<<23;
		int count=0;
		for(int j=0;j<num;j++) {
			int x;
			for(x=0;x<chain.numspots;x++)
				if(!FloorMod(chain.accelerators[x].GetCS(j),chain.mod,0))
					break;
			if(x==chain.numspots)
				count++;
		}
		bruteforcespeed=num/(Time()-starttime);
		if(loud)
			printf("brute force: %f Megaseed/sec\n",bruteforcespeed/pow(10.0,6));
	}


	// start timer here 

	double slicesize=pow(2.0,64-26);

	double starttime=Time();
	
	int bsnum;
	for(bsnum=0;bsnum<baseseeds_count;bsnum++) {
		if(benchmark&&Time()>starttime+5)
			break;

		int bs=baseseeds[bsnum];

		int outputseeds[50000];
		int outputseedscount_slice=0;
		double slice_starttime=Time();
		
		// this is where the lifting occurs.
		// It's hard coded to 26 bits because that covers common cases
		// It's not hard to make it more general, but trying to keep it simple.
		for(int k=0;k<(1<<(32-26));k++) {
			int bs2=bs+(1<<26)*k;
			for(int x=0;x<numspots;x++)
				chain.accelerators[x].Init(
					config_SEEDFUNC(bs2          ,basex+x,basey),
					config_SEEDFUNC(bs2+(1LL<<32),basex+x,basey));
			
			// Do the skip chaining
			int outputseeds_count=0;
			chain.Do(outputseeds,&outputseeds_count,50000);

			baseseeds_total_count+=outputseeds_count;
			outputseedscount_slice+=outputseeds_count;
			
			
			// There's a bunch of options for printing stuff out here.  Just uncomment to use them.
			if(1) {
				/*if(outputseeds_count)
					printf(" %d ",outputseeds_count);*/

				
				for(int i=0;i<outputseeds_count;i++) {
					__int64 bs3=bs2+(1LL<<32)*outputseeds[i];
					
					
					//printf(" %lld ",bs3);
					
					// re-confirm this actually passes
					int match=0;
					for(int x=0;x<numspots;x++)
						match+=((config_SEEDFUNC(bs3,basex+x,basey)>>24)%chain.mod)==0;
					//printf("%d",match);
					if(match==numspots) {
						//printf(".",match);
						if(config_print_seeds)
							printf(" %lld ",bs3);
					}


					
#ifdef CHECKCUBIOMESGEN
					// check if this will actually generate mushroom islands using cuBiomes code
					// to actually generate the map

					int actualshroom=CheckFullGenForMushrooms(bs3,basex,basey,numspots,1);
					//printf(" actual generated: %d\n ",actualshroom);
					if(actualshroom==numspots) {
						printf("SEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEED at %d %d",basex*256,basey*256);
						printf(" %d : %lld ",actualshroom,bs3);
						//getch();
					}
#endif

					//printf("%d",match);

				}
				//printf(".");
				//printf("%d",chain.
			}
		
		}
		if(loud) {
			double slicetime=(Time()-slice_starttime);
			double progress=(double)bs/(1<<26);
			printf("\n%0.2f %d %.0f %f petaseed/sec slice: %05.2f teraseed/sec, %.2f s speedup %.1f | ",
				progress,
				baseseeds_total_count,
				(baseseeds_total_count/progress),
				progress*pow(2.0,64)/pow(10.0,15)/(Time()-starttime),
				slicesize/pow(10.0,12)/slicetime,
				slicetime,
				slicesize/slicetime/bruteforcespeed);
		}
	}
	if(loud) {
		printf("\ndone, %d found, %f petaseed/sec %f total time %.1f teraseed/sec (slices)\n",
			baseseeds_total_count,
			pow(2.0,64)/pow(10.0,15)/(Time()-starttime),
			(Time()-starttime),
			baseseeds_count*slicesize/pow(10.0,12)/(Time()-starttime));
	}
	return bsnum*(1LL<<(32+6))/(Time()-starttime);
}


//int thread_data[10];

int numthreads=0;
int totalseeds_x32=0;

static DWORD WINAPI AdvancedAdvancedChainsThread(LPVOID threadnumv) {
	int threadnum=(int)threadnumv;
	//int x=(int)thread_data[threadnum];

	static int pos=39400;

	while(true) {
		while(threadnum>=numthreads) {
			Sleep(1000);
		};
		int height=10000;
		int y=(pos%height)-5000;
		int x=pos/height;
		pos++;
		if(pos%200==0)
			printf("%d(%d,%d) ",pos,x,y);
		DoSkipChaining(6+x,y,100,false,false);
		totalseeds_x32++;;
	};
	/*int x=(int)data;
	for(int y=2590;y<5800;y++) {
		if(x==0&&(y%10)==0)
			printf("pos:%d\n",y);
		DoSkipChaining(-14+x,y,100,false);
	}*/
	return 0;
}

void DoModSpeedTest() {
	for(int m=80;m<100;m++) {
		if(MaxLift(m)!=m&&MaxLift(m)<8)
			DoSkipChaining(0,2,m,true,true);
	}
}

void DoMultithreadedFind() {

	for(int x=0;x<10;x++)
		CreateThread(0,0,AdvancedAdvancedChainsThread,(LPVOID)x,0,0);
	numthreads=4;
	double starttime=Time();
	int c=0;
	while(true) {
		if(kbhit()) {
			numthreads=getch()-'0';
			printf("NumThreads:%d ",numthreads);
		};
		Sleep(100);
		//printf("%d\n",totalseeds_x32);
		if(c%1000==0)
			printf("%f petaseed/sec\n",(totalseeds_x32*pow(2.0,64))/(Time()-starttime)/pow(10.0,15));
		c++;
	};
}

void DoModSpeedTest2() {
	for(int mod=12;mod<40;mod++) {
		int lift=MaxLift(mod);
		//if((mod%2)==1)
		if(lift<=4&&mod!=lift)
		{
			for(int x=0;x<50;x++) {
				double speed=DoSkipChaining(x,1,mod,true,true);
				//int num32spaces=0;

				printf("\nMod: %d Speed (upper bits): %f Gseed/sec\n\n",mod,speed/pow(10.0,9));
			}
		}
	}
}

int main(int argc, char *argv[]) {

	//DoMultithreadedFind();
	DoSkipChaining(0,2,100,true,0);
	//DoModSpeedTest2();

	getch();
}