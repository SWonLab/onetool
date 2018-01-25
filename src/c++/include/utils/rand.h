#pragma once
#ifndef __WISARD_RAND_H__
#define __WISARD_RAND_H__

#include <stdlib.h>
#include <stdio.h>
#include "global/common.h"

namespace ONETOOL {

double wsUnifrand();

/******************************************************************************
 * @file	CmRandomNumberGenerator.h
 * @brief	Random Number Generator H file
 * @author  Takayuki HARUKI (University of Toyama, Japan)
 * @since	Nov. 2005
 *
 ******************************************************************************/

//!< Period Parameter for Mersenne Twister
#define	NNNN (624)
#define	MMMM (397)

/******************************************************************************
 * @class	CmRandomNumberGenerator
 * @brief	Class (of mathematical module) for generating random number
 *			by Mersenne Twister
 * @author  Takayuki HARUKI (University of Toyama, Japan)
 * @data	Nov. 2005
 * @par	how to use
 *			-# get singleton object by using getInstance
 *			-# change seed by changeSeed if you want
 *			-# get random number by calling getFloat or getDouble functions
 *			-# call release finally only at once
 *
 ******************************************************************************/
class cWisardRandom
{
private:
	
/**
 * @name Constructor and Destructor (private for SINGLETON)
 */
//@{
	cWisardRandom();
	virtual ~cWisardRandom();
//@}

public:
	
/**
 * @name Singleton
 */
//@{
	static cWisardRandom* getInstance();
	void release();
//@}

/**
 * @name Random number
 */
//@{
	float			getFloat();
	double			getDouble();
	unsigned long	getInt();
//@}

/**
 * @name Seed
 */
//@{
	void changeSeed(unsigned long a_ulSeed);
//@}

private:
	
	void init_genrand(unsigned long s);	// initialize mt[N] with a seed
	unsigned long genrand_int32(void);	// [0, 0xffffffff]
	long genrand_int31(void);			// [0, 0x7fffffff]
	double genrand_real1(void);			// [0, 1]
	double genrand_real2(void);			// [0, 1)
	double genrand_real3(void);			// (0, 1)
	
private:
	
	static cWisardRandom* s_pInstance;
	
	static unsigned long	mt[NNNN];	//!< array for the state vector
	static int				mti;	//!< mti==N+1 means mt[N] is not initialized
};

int wsAdvrand(wsUint N_range);

#if RAND_MAX == INT_MAX
#	define wsRand rand
#else
inline int wsRand()
{
#if RAND_MAX == INT_MAX
	return rand();
#else
	//	LOG("X\n");
	int const l = rand() & 0x7FFF ;
	int const h = rand() & 0x7FFF ;
	return ( h << 15 ) | l ;
#endif
}
#endif

} // End namespace ONETOOL

#endif
