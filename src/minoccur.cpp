/*
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * Copyright (C) 2014 Bonsai Bioinformatics Research Group
 *
 * This file is part of SortMeRNA.
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
 * file: minoccur.cpp
 * contact: jenya.kopylov@gmail.com, laurent.noe@lifl.fr, helene.touzet@lifl.fr
 *
 */

#include "../include/minoccur.hpp"

using namespace std;


/* to sort the kmers_occurrences array during search for MINOCCUR value using RATIO */
int _qsort(const void *x, const void *y) 
{
	return (*(int*)x - *(int*)y);
}

#ifdef WINDOWS
int isnewline( char *ch )
{
	/* '\n' */
	if ( *ch == 0x0A )
	{
		/* '\r ' */
		if ( *(ch-1) == 0x0D )
		{
			return 1;	
		}
		/* '\n' */
		return 0;
	}
	/* '\r' */
	else if ( *ch == 0x0D )
	{
		return 0;
	}
	/* other */
	else return -1;
}
#endif

size_t strs = 0;
char **rrnas = NULL;
char* rrnatext = NULL;
off_t filesize = 0;
int maxlen = 0;
int minlen = 100000;

unsigned int MINOCCUR = 0;




/*
 *
 * FUNCTION 	: void loadrrnaraw   ( void )
 *		  load the rRNA original fasta file to compute minoccur
 **************************************************************************************************************/
void loadrrnaraw ( void )
{
	int   fd       = -1;
	/* rRNA database input file */
	string fname = rrnadbfile;

	if ((fd = open(fname.c_str(), O_RDONLY)) == -1) 
	{
		printf("  Could not open rRNA file!\n");
		printf("  Check that your database filename has suffix .fasta\n\n");
		exit(EXIT_FAILURE);
	}
	if ((filesize = lseek(fd, 0, SEEK_END)) == -1) 
	{
		printf("  Could not seek rRNA file!\n\n");
		exit(EXIT_FAILURE);
	}
	if (lseek(fd, 0, SEEK_SET) == -1) 
	{
		printf("  Could not seek rRNA file!\n\n");
		exit(EXIT_FAILURE);
	}
	if ((rrnatext = (char*)mmap(0, filesize, PROT_READ, MAP_SHARED, fd, 0))
			== MAP_FAILED) 
	{
		printf("  Could not mmap rRNA file !\n\n");
		exit(EXIT_FAILURE);
	}
        if (close(fd) == -1) 
	{
		printf("  Could not close rRNA file!\n\n");
		exit(EXIT_FAILURE);
	}

	strs = 0;

	/* count the number of strings in the file  */
	if ( rrnatext[0] == '>' ) 
		strs+=2;

	for ( size_t i = 1; i < (size_t)filesize; i++ )
	{
		if ( rrnatext[i] == '>' && rrnatext[i-1] == '\n' ) 
			strs+=2;
	}

	/* create a char* array for pointers to each string in mmap */
	rrnas = (char**)malloc ( strs*sizeof(char*) );
	if ( rrnas == NULL )
	{
	      fprintf(stderr, "out of memory\n");
	      exit(EXIT_FAILURE);
	}

	char* line = &rrnatext[0];
	int rrnalen = 0;

	int j = 0;

	while ( *line != '\0' )
	{
		if ( *line == '>' )
		{
			/* the rrna tag */
			rrnas[j++] = line;	
			while ( *line++ != '\n' );

			/* the read */
			rrnas[j++] = line;			
		}

		rrnalen = 0;

		while ( *line != '\0' && *line != '>' )
	  {
			if ( *line != '\n' || *line != ' ' ) 
				rrnalen++;
			line++;
		}

		/* compute the minimum length rrna */
	  rrnalen > maxlen ? maxlen = rrnalen : maxlen;
		rrnalen < minlen ? minlen = rrnalen : minlen;
	}

	return ;
}//~loadrrnaraw()



/*
 *
 * FUNCTION 	: void find_minoccur   ( int minlenread )
 *		  (see `../include/minoccur.hpp' for a description)
 **************************************************************************************************************/
void find_minoccur ( int minlenread )
{
	/* create the rrnas array */
	loadrrnaraw();

	if ( minlenread > minlen ) minlenread = minlen;

	/* the minimum occurrence of an s-mer possible in an minlenread-nt window,
           initialized to a high number for the first rRNA */
	unsigned int minoccur = 1<<30;

	/* the total number of s-mers in a minlenread-nt window */
	int number_of_kmers = (minlenread - LNWIN + 1)*2;

	/* the number of possible s-mer `hits' on a minlenread-nt window, initialized to
           the minimum number possible */
    unsigned int min_num_win = ceil(RATIO*number_of_kmers);
	unsigned int minoccur_win = number_of_kmers-min_num_win;

	char* rrnaf;

	/* go through every rrna in the database */ 
	for ( int i = 1; i < (int)strs; i+=2 )
	{
		char* mynt = rrnas[i];
		int _j = 0;
		int rrnalen = 0;
		rrnaf = (char*)malloc( sizeof(char)*maxlen );

		/* an ordered set of number of occurrences of s-mers which are less
           	than minoccur */
		multiset<unsigned int> minoccurrence_tree;

		/* the actual count of active s-mers in a minlenread-nt window */
		unsigned int num_win = min_num_win;

	  	/* encode forward rRNA sequence into tertiary integer alphabet {0,1,2,3} */
    while ( *mynt != '\0' && *mynt != '>' ) 
	  {
			/* skip line feed, carriage return or empty space in rRNA sequence */
			if ( *mynt == '\n' or *mynt == '\r' or *mynt == ' ' )
			{
				mynt++; 
				continue;
			}
			/* exact character */
			else if ( *mynt == 'A' or *mynt == 'a' 
                   or *mynt == 'C' or *mynt == 'c' 
				   or *mynt == 'G' or *mynt == 'g' 
				   or *mynt == 'T' or *mynt == 't' 
				   or *mynt == 'U' or *mynt == 'u')
			{
	    		rrnaf[_j++] = map_nt[(int)*mynt++];
			}
			/* iupac ambigiuty code */
			else
			{
				int ran = rand();
				switch( *mynt++ )
				{
					case 'B':
					{
							      /*C  G  T*/
						char amb[3] = { 1, 2, 3 };
						rrnaf[_j++] = amb[ran%3];
					}
					break;
	
					case 'D':
					{
							      /*A  G  T*/
						char amb[3] = { 0, 2, 3 };
						rrnaf[_j++] = amb[ran%3];
					}
					break;

					case 'H':
					{
							      /*A  C  T*/
						char amb[3] = { 0, 1, 3 };
						rrnaf[_j++] = amb[ran%3];
					}
					break;
			
					case 'K':
					{
							      /*G  T*/
						char amb[2] = { 2, 3 };
						rrnaf[_j++] = amb[ran%2];
					}	
					break;

					case 'M':
					{
							      /*A  C*/
						char amb[2] = { 0, 1 };
						rrnaf[_j++] = amb[ran%2];
					}
					break;

					case 'N':
					{
							      /*A  C  G  T*/
						char amb[4] = { 0, 1, 2, 3 };
						rrnaf[_j++] = amb[ran%4];
					}
					break;

					case 'R':
					{
							      /*A  G*/
						char amb[2] = { 0, 2 };
						rrnaf[_j++] = amb[ran%2];
					}
					break;

					case 'S':
					{
							      /*C  G*/
						char amb[2] = { 1, 2 };
						rrnaf[_j++] = amb[ran%2];
					}
					break;

					case 'V':
					{
							      /*A  C  G*/
						char amb[3] = { 0, 1, 2 };
						rrnaf[_j++] = amb[ran%3];
					}
					break;

					case 'W':
					{
							      /*A  T*/
						char amb[2] = { 0, 3 };
						rrnaf[_j++] = amb[ran%2];
					}
					break;

					case 'Y':
					{
							      /*C  T*/
						char amb[2] = { 1, 3 };
						rrnaf[_j++] = amb[ran%2];
					}
					break;

					default: 
					{
						cout << (char)*mynt << " is not a legal character.\n";
						exit(1); 
					}
					break;
				}//~switch(*nt)
			}//~assign iupac code

	    		rrnalen++;
	  	}

		/* for each minlenread-nt window, make sure the number of s-mers does not decrease
                   below the threshold */

		/* begin at the 1st letter ( forward s-mer ) */
		char* kmer_win_f = rrnaf;
		/* begin at the s'th letter ( reverse s-mer ) */
		char* kmer_win_r = rrnaf+LNWIN-1;

		/* the numerical values of s-mers */
		unsigned int lookupf = 0;
		unsigned int lookupr = 0;
	
		/* there are q = (readlen-LNWIN+1)*2 s-mer windows in a minlenread-nt window, 
		   forward + reverse; this array is sorted and the minimum number of s-mers 
		   existing in the lookup table must be ceil(RATIO*q) for the minlenread-nt 
		   window  */
		unsigned int kmers_occurrences[number_of_kmers];

		/*********** compute the first minoccur value *********************************/

		/* create the first forward and reverse s-mers of the rrna strand */
		for ( int j = 0; j < PARTIALWIN; j++ )	
		{	
			(lookupf <<= 2) |= (int)*kmer_win_f++;
			(lookupr <<= 2) |= (int)*kmer_win_r--;
		}
	
		/* reset the reverse s-mer pointer to the next new character in the rRNA */
		kmer_win_r = rrnaf+LNWIN;

		/* a queue to keep the list of s-mers in a minlenread-nt window */
		queue<unsigned int> kmers_queue;
		
		/* add the first forward & reverse s-mer to the queue */
		kmers_queue.push( lookupf );
		kmers_queue.push( lookupr );

		/* complete the kmer_occurrences array */
		kmers_occurrences[0] = kmerf[lookupf].count;
		kmers_occurrences[1] = kmerr[lookupr].count;

		/* pointer to the next free spot in the s-mer occurrences array */
		unsigned int* ptr_kmer = &kmers_occurrences[2];

		/* number of minlenread-nt windows on the rRNA strand */	
		int testwindows = rrnalen - minlenread + 1;

		/* the number of s-mer windows in a minlenread-mer window */
		int kmers_in_read_len = minlenread - LNWIN + 1;

		/* loop through the 134 s-mers and add them to the queue & kmers_occurrences */
		for ( int p = 1; p < kmers_in_read_len; p++ )
		{
			(( lookupf <<= 2 ) &= mask32 ) |= (int)*kmer_win_f++;
			(lookupr >>= 2) |= ((int)*(kmer_win_r)++ << (LNWIN-2));

			kmers_queue.push( lookupf );
			kmers_queue.push( lookupr );

			*ptr_kmer++ = kmerf[lookupf].count;
			*ptr_kmer++ = kmerr[lookupr].count;
		} //~for every forward & reverse s-mer in the minlenread-nt window

		/* determine the `minoccur' value and build the `minoccurrence_tree' */

		/* sort the kmers_occurrences array */
		qsort( kmers_occurrences, number_of_kmers, sizeof(unsigned int), _qsort );
			
		/* (re)calculate minoccur */
		if ( kmers_occurrences[minoccur_win] <= minoccur ) minoccur = kmers_occurrences[minoccur_win];
		
		/* the number of active windows in the new minlenread-nt window is higher than the min_num_win */
		else
		{
			int reduce = minoccur_win;
			while ( (reduce > 0) && (kmers_occurrences[--reduce] >= minoccur) ) num_win++; 
		}

		/* put all kmer_occurrences up to minoccur into the minoccurence_tree set */
		int roof = number_of_kmers-num_win;
		for ( int l = 0; l < roof; l++ )
			minoccurrence_tree.insert( kmers_occurrences[l] );	

		/* loop through all of the minlenread-nt windows on the rRNA strand, updating the
                   kmers_queue and kmers_occurrences lists */
		for ( int q = 1; q < testwindows; q++ )
		{
			/* find the number of occurrences of the first s-mer */
			unsigned int out_1 = kmerf[kmers_queue.front()].count;

			/* remove the first s-mer from the queue */
			kmers_queue.pop();

			unsigned int out_2 = kmerr[kmers_queue.front()].count;
			kmers_queue.pop();

			/* find the hashid of the newly added s-mers */
			(( lookupf <<= 2 ) &= mask32 ) |= (int)*kmer_win_f++;
			(lookupr >>= 2) |= ((int)*(kmer_win_r)++ << (LNWIN-2));

			/* add them to the end of the queue */
			kmers_queue.push( lookupf );
			kmers_queue.push( lookupr );

			/* find the number of occurrences of the two inserted s\2-mers */
			unsigned int in_1 = kmerf[lookupf].count;
			unsigned int in_2 = kmerr[lookupr].count;

			/* there are 9 unique possibilities of where two s\2-mers are 
                           removed and where two are added,

			   (1)    iioo		
				|---------|---------|	num_win = num_win

			   (2)	   iio	      o
				|---------|---------|   num_win - 1

			   (3)     ii	      oo	
				|---------|---------|   num_win - 2

			   (4)	   i	     ioo
				|---------|---------|   num_win - 1

			   (5)              iioo
				|---------|---------|   num_win = num_win
	
			   (6)     ooi        i
 				|---------|---------|   num_win + 1

			   (7)     oo         ii
				|---------|---------|   num_win + 2

			   (s)	   o	     oii
				|---------|---------|   num_win + 1

			   (9)     oi	      io
				|---------|---------|   num_win = num_win

		         */
			 
			 /* the first 3 cases */
			 if ( (in_1 < minoccur) && (in_2 < minoccur) )
			 {
				minoccurrence_tree.insert( in_1 );
				minoccurrence_tree.insert( in_2 );

				/* case (1) */
				if ( (out_1 < minoccur) && (out_2 < minoccur) )
				{
					/* remove the out_1 and out_2 occurrences from mininum occurrence array
					   (only 1 copy) */
					minoccurrence_tree.erase( minoccurrence_tree.find( out_1 ) );
					minoccurrence_tree.erase( minoccurrence_tree.find( out_2 ) );
				}//~case (1)
				
				/* case (2) */
				else if ( (out_1 < minoccur) ^ (out_2 < minoccur) )
				{
					/* remove the `out' s\2-mer < minoccur from occurrences array */
					minoccurrence_tree.erase( minoccurrence_tree.find ( out_1 > out_2 ? out_2 : out_1 ) );

					num_win--;			

					/* the number of active s-mers in a minlenread-window is
					   less than the minimum amount, edit the minoccur value */
					if ( num_win < min_num_win )
					{
						/* minoccur is set to the next highest occurrence in
                                                   the minlenread-nt window */
						minoccur = *minoccurrence_tree.rbegin();
						minoccurrence_tree.erase( minoccurrence_tree.find(minoccur) );
						num_win++;
					}
				}//~case (2)

				/* case (3) */
				else if ( (out_1 >= minoccur) && (out_2 >= minoccur) )
				{
					num_win-=2;

					if ( num_win < min_num_win )
					{
						/* minoccur is set to the second highest occurrence in
                                                   the minlenread-nt window, we want to erase only one copy */
						minoccur = *minoccurrence_tree.rbegin();
						minoccurrence_tree.erase( minoccurrence_tree.find(minoccur) );
						minoccur = *minoccurrence_tree.rbegin();
						minoccurrence_tree.erase( minoccurrence_tree.find(minoccur) );
						num_win+=2;
					}					

				}//~case (3)
			 }//~ the first 3 cases

			/* cases (6) and (7) */
			else if ( (out_1 < minoccur) && (out_2 < minoccur) )
			{
				minoccurrence_tree.erase( minoccurrence_tree.find ( out_1 ) );
				minoccurrence_tree.erase( minoccurrence_tree.find ( out_2 ) );

				/* case (6) */
				if ( (in_1 < minoccur) ^ (in_2 < minoccur) )
				{
					/* insert the smaller of in_1 and in_2 into the 
                                           minimum occurrence tree */
					minoccurrence_tree.insert( in_1 > in_2 ? in_2 : in_1 ); 
					num_win++;
				}//~case (6)

				/* case (7) */
				else if ( (in_1 >= minoccur) && (in_2 >= minoccur) )
				{
					num_win+=2;
				}//~ case (7)

			 }//~ cases (6) and (7)

			/* cases (5) - nothing happens; and (s) */
			else if ( (in_1 >= minoccur) && (in_2 >= minoccur) )
			{
				/* case (5) */
				if ( (out_1 >= minoccur) && (out_2 >= minoccur) )
				{
					;
				}//~case (5)
				/* case (s) */
				if ( (out_1 < minoccur) ^ (out_2 < minoccur) )
				{
					/* erase the smaller of out_1 and out_2 from the 
                                           minimum occurrence tree */
					minoccurrence_tree.erase( minoccurrence_tree.find( out_1 > out_2 ? out_2 : out_1 ) );
					num_win++;
				}
			 }//~ cases (5) and (s)

			 /* case (4) */
			 else if ( (out_1 >= minoccur) && (out_2 >= minoccur) )
			 {
				if ( (in_1 < minoccur ) ^ (in_2 < minoccur) )
				{
					/* insert the smaller `in' s\2-mer < minoccur into occurrences array */
					minoccurrence_tree.insert( in_1 > in_2 ? in_2 : in_1 );

					if ( --num_win < min_num_win )
					{
						/* minoccur is set to the next highest occurrence in
                                                   the minlenread-nt window */
						minoccur = *minoccurrence_tree.rbegin();
						minoccurrence_tree.erase( minoccurrence_tree.find(minoccur) );
						num_win++;
					}
				}
			 }//~ case(4)

			 /* case (9) */
			 else if ( (out_1 < minoccur) ^ (out_2 < minoccur) )
			 {
				if ( (in_1 < minoccur ) ^ (in_2 < minoccur ) )
				{
					/* erase the smaller `out' s\2-mer < minoccur from occurrences array */
					minoccurrence_tree.erase( minoccurrence_tree.find( out_1 > out_2 ? out_2 : out_1 ) );

					/* insert the smaller `in' s\2-mer < minoccur into occurrences array */
					minoccurrence_tree.insert( in_1 > in_2 ? in_2 : in_1 );
				}				
			 }//~ case (9)
						
		}//~for every minlenread-nt window in the rRNA strand 

		free(rrnaf);
		rrnaf = NULL;

	}//~for every rrna

	/* free memory */
	free(rrnas);

        if ( munmap(rrnatext, filesize ) == -1 )
	{
	     std::cerr << "  Could not munmap rRNAs file!" << std::endl;
	     exit(EXIT_FAILURE);
	}

	/* set the MINCCOUR value */
	MINOCCUR = minoccur;

	return ;

}//~ find_minoccur()


