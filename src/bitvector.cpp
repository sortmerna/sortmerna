/*
 * @file bivector.cpp
 * @brief Functions for manipulating the Levenshein bitvector table.
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2013-2015 Bonsai Bioinformatics Research Group
 * 2015 Knight Lab, Department of Pediatrics, UCSD, La Jolla
 *
 * SortMeRNA is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SortMeRNA is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 *
 * file: bitvector.cpp
 * contact: jenya.kopylov@gmail.com, laurent.noe@lifl.fr, helene.touzet@lifl.fr
 *
 */

/** @file */ 

#include <iostream>
#include <iterator>

#include "../include/bitvector.hpp"

using namespace std;


/// mask_4 = 15 = 00001111; this is to keep the bitvectors of length 4
MYBITSET mask_4 = 15; 
/// mask the first four bitvectors of substrings $xxx
MYBITSET mask_3 = 7;



/*
 *
 * FUNCTION 	: void init_win_r()
 *		  (see bitvector.hpp for a description)                  
 *
 ************************************************************************************************/
void 
init_win_f ( char* ptrf, 
			MYBITSET* bittable_000, 
			MYBITSET* bittable_010,
            int numbvs)
{
	/// [w_1] forward

	MYBITSET *reset  = bittable_000;

	/// set manually the bitvectors at position i = 0 of [w_1] reverse
	for ( int bitn = 2; bitn >= 0; bitn-- )
	{
		*(reset+*ptrf++) |= (1<<bitn);
	}

	/// set the bitvectors for positions i > 0 of [w_2] forward
	MYBITSET *setbit = bittable_010;
	MYBITSET *win_ptr1 = bittable_000;
	MYBITSET *win_ptr2 = bittable_010;

	for ( int i = 1; i <= numbvs; i++ )
	{
	 	*win_ptr2 = (*win_ptr1++)<<1;
	  	*win_ptr2++ &= mask_4;
	  	/// if i%4 == 0
		if ( !(i&3) )
	  	{
			/// set the LSB of candidate nt bitvector to 1
			*(setbit+*ptrf++) |= 1;
			/// reset the setbit pointer to subsequent bitvector
			setbit = win_ptr2;
	  	}	
	}

}//~init_win_f()



/*
 *
 * FUNCTION 	: void init_win_r()
 *		  (see bitvector.hpp for a description)                  
 *
 ************************************************************************************************/
void 
init_win_r ( char* ptrr, 
			MYBITSET* bittable_000, 
			MYBITSET* bittable_010,
            int numbvs)
{
 	/// [w_1] reverse

	MYBITSET *reset  = bittable_000;

	/// set manually the bitvectors at position i = 0 of [w_1] reverse
	for ( int bitn = 2; bitn >= 0; bitn-- )
	{
		*(reset+*ptrr--) |= (1<<bitn);
	}

	/// set the bitvectors for positions i > 0 of [w_2] forward (note: mask = 15 = 00001111; this is to keep the bitvectors of length 4)

	MYBITSET *setbit = bittable_010;
	MYBITSET *win_ptr1 = bittable_000;
	MYBITSET *win_ptr2 = bittable_010;

	for ( int i = 1; i <= numbvs; i++ )
	{
	 	*win_ptr2 = (*win_ptr1++)<<1;
	  	*win_ptr2++ &= mask_4;
		if ( !(i&3) )
	  	{
			*(setbit+*ptrr--) |= 1;
			setbit = win_ptr2;
	  	}
	}

}//~init_win_r()



/*
 *
 * FUNCTION 	: void offset_win_k1()
 *		  (see bitvector.hpp for a description)
 *
 ***************************************************************************************/
void 
offset_win_k1 ( char *FW2, 
		char *RW1,
		MYBITSET* bittable_0P0,
		MYBITSET* bittable_100,
		MYBITSET* bittable_110,
        int numbvs)
{
	/// [w_1] reverse

	/// compute bitvectors for shifted window from previous window
	MYBITSET *win_ptr1 = bittable_0P0;
	MYBITSET *win_ptr2 = bittable_100;

	for ( int i = 0; i < numbvs; i++ )
	{
	 	*--win_ptr2 = *--win_ptr1;
	}/// for every bitvector in [w_1] reverse

	/// compute bitvectors for depth = 0
	MYBITSET *setbit = win_ptr1;

	for ( int i = 4; i > 0; --i )
	{
	 	*win_ptr1++ >>= 1;
	}
	*(setbit+*RW1) |= MSB4; /// set second highest MSB to 1

	/// set bit for depth = 1
	*(win_ptr2+*RW1) |= MSB8; /// set MSB to 1


	/// [w_2] forward

	/// 1. offset all bitvectors for [w_2] forward
	win_ptr1 = bittable_100;
	win_ptr2 = bittable_110;

	for ( int j = 0; j < numbvs; j++ )
	{
	 	*win_ptr1++ = *win_ptr2++;
	}//~for every bitvector in [w_2] forward

	/// mask the first four bitvectors of substrings $xxx
	win_ptr2 = bittable_100;

	for ( int j = 4; j > 0; --j )
	{
	 	*win_ptr2++&=mask_3;
	}

	/// compute the bitvectors for newly added letter from window shift
	setbit = win_ptr1;

	for ( int j = 4; j > 0; --j )
	{
		(*win_ptr1++<<=1)&=15;
	}
	*(setbit+*FW2) |= 1;

}//~offset_win_k1()


/*
 *
 * FUNCTION 	: void output_win_k1()
 *		  (see bitvector.hpp for a description)
 *
 ***************************************************************************************/
void 
output_win_k1 ( MYBITSET* bittable_000, bool w, int partialwin )
{
		/// [w_2] forward
		if ( w )
			cout << "forward" << endl;
		else 
			cout << "reverse" << endl;

		cout << "\n\t\t\tA\tC\tG\tT\n\t\t";
		int ceil = partialwin-2;

		/// for each nt letter
		for ( int i = 0; i < ceil; i++ )
		{
			cout << i << "\t";
	       
	       	/// for each depth
	       	for ( int j = 0; j < 4; j++ )
	       	{
		 		cout << (int)*(bittable_000++) << "\t";
	       	}
	       	cout << "\n\t\t";
	   }

	   cout << endl;

}//~output_win_k1()

