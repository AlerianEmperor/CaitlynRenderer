#ifndef _RND_H_
#define _RND_H_
#include <random>

using namespace std;

thread_local uint32_t s_RndState = 1;
static const float imax = 1.0f / UINT32_MAX;
static const float irand_max = 1.0f / RAND_MAX;
float randf()
{
	uint32_t x = s_RndState;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 15;
	s_RndState = x;
	return x * imax;
	//return rand() * irand_max;
}

static uint32_t PCG_Hash(uint32_t input)
{
	uint32_t state = input * 747796405u + 2891336453u;
	uint32_t word = ((state >> ((state >> 28u) + 4u)) ^ state) * 277803737u;
	return (word >> 22u) ^ word;
}

/*
static float randf2(uint32_t& seed)
{
	seed = PCG_Hash(seed);
	return seed * imax;//(float)numeric_limits<uint32_t>::max();
}
*/

static float randf2()
{
	s_RndState = PCG_Hash(s_RndState);
	return s_RndState * imax;//(float)numeric_limits<uint32_t>::max();
}

#endif // !_RND_H_
