# SkipChains



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


It might seem excessively complex to build the skip arrays like this, rather than just store one of them (this was my previous attempt), but it actually 
gives a speedup of about 200x over that attempt.  This is because the 32 bit seed space being searched turns out to be kinda small, so
the effort to build up the skip array was actually taking most of the time.


I intend to use template arguments to configure this and stay the highest performance.
However, since most of the seed finding utilities are in Java, this might remain a c++ demo.
I have tentative plans to put this on a GPU.


The speedup for 9 mushroom chunks goes as follows:
250,000x speedup by pre-checking base seeds (.04^9 vs .01^9)
70,000x speedup using this multi-level skip chaining method
4x speedup (guess) by using SeedAccelerator2
10x speedup using ModInverseSpace
That's a 700,000,000,000x speedup over a pure brute force approach.


Big thanks to Matthew Bolan for teaching me a lot about this modular arithmetic stuff and Hensel lifting.
I don't know if this has been invented before, but if it has been, credit to them too.
