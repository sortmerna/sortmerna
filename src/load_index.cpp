/**
 * @file preprocess_data.cpp
 * @brief Load data (reads) using mmap.
 * @parblock
 * SortMeRNA - next-generation reads filter for metatranscriptomic or total RNA
 * @copyright 2012-16 Bonsai Bioinformatics Research Group
 * @copyright 2014-16 Knight Lab, Department of Pediatrics, UCSD, La Jolla
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
 * along with SortMeRNA.  If not, see <http://www.gnu.org/licenses/>.
 * @endparblock
 *
 * @contributors Jenya Kopylova, jenya.kopylov@gmail.com
 *               Laurent Noé, laurent.noe@lifl.fr
 *               Pierre Pericard, pierre.pericard@lifl.fr
 *               Daniel McDonald, wasade@gmail.com
 *               Mikaël Salson, mikael.salson@lifl.fr
 *               Hélène Touzet, helene.touzet@lifl.fr
 *               Rob Knight, robknight@ucsd.edu
 */

#include "../include/load_index.hpp"

/*! @brief Map nucleotides to integers.

    Ambiguous letters map to 4. 
    {A/a,C/c,G/g,T/t,U/u} = {0,1,2,3,3} respectively. 
*/
char nt_table[128] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

 /*
 *
 * @function check files: load the @SQ lines for SAM output & input
 *                        all header information to SAM output file
 * @return void
 * @version 2.1 February 7, 2016
 *
 *******************************************************************/
void load_index_stats(vector< pair<string,string> >& myfiles,
                      char** argv,
                      int argc,
                      bool yes_SQ,
                      char* acceptedstrings_sam,
                      long match,
                      long mismatch,
                      long gap_open,
                      long gap_extension,
                      vector<vector<uint32_t> >& skiplengths,
                      vector<uint16_t>& num_index_parts,
                      vector<vector<index_parts_stats> >& index_parts_stats_vec,
                      vector<uint64_t>& full_ref,
                      vector<uint64_t>& full_read,
                      vector<uint32_t>& lnwin,
                      vector<uint32_t>& partialwin,
                      vector<uint32_t>& minimal_score,
                      uint64_t number_total_read,
                      vector<pair<double, double> >& gumbel,
                      vector<uint64_t>& numbvs,
                      vector<uint64_t>& numseq)
{
  ofstream acceptedsam;
  
  if ( samout_gv )
  {
    acceptedsam.open (acceptedstrings_sam);
    if (!acceptedsam.good())
    {
      fprintf(stderr,"  %sERROR%s: could not open SAM output file for writing.\n",startColor,"\033[0m");
      exit(EXIT_FAILURE);
    }
    // @HD header
    else acceptedsam << "@HD\tVN:1.0\tSO:unsorted\n";
  }

  // create and initialize scoring matrix
  long alphabetSize = 4;
  long **scoring_matrix = new long *[alphabetSize];

  for (long i = 0; i < alphabetSize; i++)
  {
    scoring_matrix[i] = new long [alphabetSize];
    for (long j = 0; j < alphabetSize; j++)
    {
      if (i == j)
        scoring_matrix[i][j] = match;
      else
        scoring_matrix[i][j] = mismatch;
    }
  }
    
  // loop through the .stats index files for each database
  for (uint16_t index_num = 0; index_num < (uint16_t)myfiles.size(); index_num++)
  {
    ifstream stats( (char*)(myfiles[index_num].second + ".stats").c_str(), ios::in | ios::binary );
    if ( !stats.good() )
    {
      fprintf(stderr,"\n  %sERROR%s: The index '%s' does not exist.\n",startColor,"\033[0m",(char*)(myfiles[index_num].second + ".stats").c_str());
      fprintf(stderr,"  Make sure you have constructed your index using the command `indexdb'. See `indexdb -h' for help.\n\n");
      exit(EXIT_FAILURE);
    }     
    // read the file size for file used to build the index
    size_t filesize = 0;
    stats.read(reinterpret_cast<char*>(&filesize), sizeof(size_t)); 
    // read the fasta file name used to build the index
    uint32_t fastafile_len = 0;
    stats.read(reinterpret_cast<char*>(&fastafile_len), sizeof(uint32_t));
    char fastafile_name[2000];
    stats.read(reinterpret_cast<char*>(fastafile_name), sizeof(char)*fastafile_len);
    // compute reference database file size for this index
    FILE *fastafile = fopen ((char*)(myfiles[index_num].first).c_str(),"r");
    if ( fastafile == NULL )
    {
      fprintf(stderr,"    %sERROR%s: could not open FASTA reference file: %s .\n",startColor,"\033[0m",(char*)(myfiles[index_num].first).c_str());
      exit(EXIT_FAILURE);
    }
    fseek(fastafile,0L,SEEK_END);
    size_t sz = ftell(fastafile);
    fclose(fastafile);
    if ( sz != filesize )
    {
      fprintf(stderr,"    %sERROR%s: Based on file size, the FASTA file (%s) passed to --ref <FASTA file, index name>\n",startColor,"\033[0m",(char*)(myfiles[index_num].first).c_str());
      fprintf(stderr,"    does not appear to be the same FASTA file (%s) used to build the index %s.\n",fastafile_name,(char*)(myfiles[index_num].second).c_str());
      fprintf(stderr,"    Check your --ref list of files and corresponding indexes.\n\n");
      exit(EXIT_FAILURE);
    }
    // A,C,G,T background frequencies to compute the Gumbel parameters lambda and K
    double background_freq_gv[4] = {0};
    // A/C/G/T distribution frequencies
    stats.read(reinterpret_cast<char*>(&background_freq_gv), sizeof(double)*4);
    // total length of sequences in the complete database
    stats.read(reinterpret_cast<char*>(&full_ref[index_num]), sizeof(uint64_t));
    // sliding window length lnwin & initialize
    stats.read(reinterpret_cast<char*>(&lnwin[index_num]), sizeof(uint32_t));
    // total number of reference sequences in one complete reference database
    stats.read(reinterpret_cast<char*>(&numseq[index_num]), sizeof(uint64_t));
    partialwin[index_num] = lnwin[index_num]/2;
    // number of bitvectors at depth > 0 in [w_1] reverse or [w_2] forward
    numbvs[index_num] = 4*(partialwin[index_num]-3);
    // set the window shift for different seed lengths (if not set by user, or one of the lengths is <= 0)
    if ( (skiplengths[index_num][0] == 0) || (skiplengths[index_num][1] == 0) || (skiplengths[index_num][2] == 0) )
    {
      skiplengths[index_num][0] = lnwin[index_num];
      skiplengths[index_num][1] = partialwin[index_num];
      skiplengths[index_num][2] = 3;
    }    
    // number of index parts
    stats.read(reinterpret_cast<char*>(&num_index_parts[index_num]), sizeof(uint16_t));
    vector<index_parts_stats> hold;
    // information on the location and size of sequences used to build each index part
    for (uint16_t j = 0; j < num_index_parts[index_num]; j++)
    {
      index_parts_stats stats_hold;
      stats.read(reinterpret_cast<char*>(&stats_hold), sizeof(index_parts_stats));
      hold.push_back(stats_hold);
    }  
    index_parts_stats_vec.push_back(hold);
    // Gumbel parameters
    long **substitutionScoreMatrix = scoring_matrix;    
    long gapOpen1 = gap_open;
    long gapOpen2 = gap_open;
    long gapEpen1 = gap_extension;
    long gapEpen2 = gap_extension;
    bool insertions_after_deletions = false;
    double max_time = -1; // required if radomization parameters are set
    double max_mem = 500;
    double eps_lambda = 0.001;
    double eps_K = 0.005;
    long randomSeed = 182345345;
    double *letterFreqs1 = new double[alphabetSize];
    double *letterFreqs2 = new double[alphabetSize];
    long number_of_samples = 14112;
    long number_of_samples_for_preliminary_stages = 39;
    for (long i = 0; i < alphabetSize; i++ )
    {
      // background probabilities for ACGT based on reference file
      letterFreqs1[i] = background_freq_gv[i];
      letterFreqs2[i] = background_freq_gv[i];
    }
    // create an object to store the Gumbel parameters
    AlignmentEvaluer GumbelCalculator_obj;
    // set the randomization parameters
    // (will yield the same Lamba and K values on subsequent runs with the same input files)
    GumbelCalculator_obj.set_gapped_computation_parameters_simplified(max_time,
                                                                      number_of_samples,
                                                                      number_of_samples_for_preliminary_stages);
    GumbelCalculator_obj.initGapped(
    alphabetSize,
    substitutionScoreMatrix,
    letterFreqs1,
    letterFreqs2,
    gapOpen1,
    gapEpen1,
    gapOpen2,
    gapEpen2,
    insertions_after_deletions,
    eps_lambda,
    eps_K,
    max_time,
    max_mem,
    randomSeed);
    gumbel[index_num].first = GumbelCalculator_obj.parameters().lambda;
    gumbel[index_num].second = GumbelCalculator_obj.parameters().K;
    delete [] letterFreqs2;
    delete [] letterFreqs1;
    // Shannon's entropy for reference sequence nucleotide distribution
    double entropy_H_gv = -(background_freq_gv[0]*(log(background_freq_gv[0])/log(2)) +
                            background_freq_gv[1]*(log(background_freq_gv[1])/log(2)) +
                            background_freq_gv[2]*(log(background_freq_gv[2])/log(2)) +
                            background_freq_gv[3]*(log(background_freq_gv[3])/log(2)));
    // Length correction for Smith-Waterman alignment score
    uint64_t expect_L = log((gumbel[index_num].second)*full_read[index_num]*full_ref[index_num])/entropy_H_gv;
    // correct the reads & databases sizes for E-value calculation
    if ( full_ref[index_num] > (expect_L*numseq[index_num]) )
      full_ref[index_num]-=(expect_L*numseq[index_num]);
    full_read[index_num]-=(expect_L*number_total_read);
    // minimum score required to reach E-value
    minimal_score[index_num] = (log(evalue/((double)(gumbel[index_num].second)*full_ref[index_num]*full_read[index_num])))/-(gumbel[index_num].first);
    // SAM @SQ data
    if ( samout_gv )
    {
      // number of nucleotide sequences in the reference file
      uint32_t num_sq = 0;
      stats.read(reinterpret_cast<char*>(&num_sq), sizeof(uint32_t)); 
      // loop through each @SQ
      for ( uint32_t j = 0; j < num_sq; j++ )
      {
        // length of the sequence id
        uint32_t len_id = 0;
        stats.read(reinterpret_cast<char*>(&len_id), sizeof(uint32_t));
        // the sequence id string
        char s[len_id+1];
        memset(s,0,len_id+1);
        stats.read(reinterpret_cast<char*>(&s), sizeof(char)*len_id);
        // the length of the sequence itself
        uint32_t len_seq = 0;
        stats.read(reinterpret_cast<char*>(&len_seq), sizeof(uint32_t));
        // @SQ header
        if ( yes_SQ ) acceptedsam << "@SQ\tSN:" << s << "\tLN:" << len_seq << "\n";
      }
    }     
    stats.close();
  }//~for all indexes
  if ( samout_gv )
  {
    // @PG to sam file
    acceptedsam << "@PG\tID:sortmerna\tVN:1.0\tCL:";
    for ( int j = 0; j < argc; j++ ) acceptedsam << argv[j] << " ";
    acceptedsam << endl;
    acceptedsam.close();
  }//~if samout_gv
  // free memory
  for (long i = 0 ; i < alphabetSize; i++)
  {
    delete [] scoring_matrix[i];
  };
  delete [] scoring_matrix;
}//~load_index_stats()


/*
 *
 * @function load index: read from binary file the L/2-mer look-up
 * tables, the 19-mer position tables and the mini-burst tries
 * @return void
 * @version 1.0 Jan 16, 2013
 *
 *******************************************************************/
void
load_index(char* ptr_dbindex,
           string part_str,
           kmer*& lookup_tbl,
           kmer_origin*& positions_tbl,
           uint32_t& number_elements,
           uint32_t lnwin)
{
  // STEP 1: load the kmer 'count' variables (dbname.kmer.dat)
  string ptr_dbindex_str = ptr_dbindex;
  ifstream inkmer( (char*)(ptr_dbindex_str + ".kmer_" + part_str + ".dat").c_str(), ios::in | ios::binary ); 
  if ( !inkmer.good() )
  {
    fprintf(stderr,"\n  ERROR: The index '%s' does not exist.\n", (char*)(ptr_dbindex_str + ".kmer_" + part_str + ".dat").c_str());
    fprintf(stderr,"  Make sure you have constructed your index using the command `indexdb'. See `indexdb -h' for help.\n\n");
    exit(EXIT_FAILURE);
  } 
  lookup_tbl = new kmer[(1<<lnwin)]();
  if ( lookup_tbl == NULL )
  {
    fprintf(stderr,"\n  ERROR: failed to allocate memory for look-up table (paralleltraversal.cpp)\n\n");
    exit(EXIT_FAILURE);
  }
  uint32_t limit = 1<<lnwin; 
  for ( uint32_t i = 0; i < limit; i++ )
    inkmer.read(reinterpret_cast<char*>(&(lookup_tbl[i].count)), sizeof(uint32_t));
  inkmer.close();
  // STEP 2: load the burst tries ( bursttrief.dat, bursttrier.dat )
  ifstream btrie( (char*)(ptr_dbindex_str + ".bursttrie_" + part_str + ".dat").c_str(), ios::in | ios::binary );
  if ( !btrie.good() )
  {
    fprintf(stderr,"\n  ERROR: The index '%s' does not exist.\n", (char*)(ptr_dbindex_str + ".bursttrie_" + part_str + ".dat").c_str());
    fprintf(stderr,"  Make sure you have constructed your index using the command `indexdb'. See `indexdb -h' for help.\n\n");
    exit(EXIT_FAILURE);
  }
  // loop through all 9-mers
  for ( uint32_t i = 0; i < (uint32_t)(1<<lnwin); i++ )
  {
    uint32_t sizeoftries[2] = {0};     
#ifdef see_binary_output
        cout << "9-mer = " << i; //TESTING
#endif
    // ptr to block of memory for two mini-burst tries
    char *dst = NULL;  
    // the size of both mini-burst tries
    for ( int j = 0; j < 2; j++ )
    {
      btrie.read(reinterpret_cast<char*>(&sizeoftries[j]), sizeof(uint32_t));
#ifdef see_binary_output
            if ( j == 0 ) cout << "\tsizeoftrie f = " << sizeoftries[j]; //TESTING
            else cout << "\tsizeoftrie r = " << sizeoftries[j]; //TESTING
#endif 
    }  
#ifdef see_binary_output
        cout << "\tlookup_tbl[i].count = " << lookup_tbl[i].count << endl; //TESTING
#endif
    // allocate contiguous memory for both mini-burst tries if they exist
    if ( lookup_tbl[i].count != 0 )
    {
      dst = new char[(sizeoftries[0]+sizeoftries[1])]();
      if ( dst == NULL )
      {
        fprintf(stderr,"  ERROR: failed to allocate memory for mini-burst tries (paralleltraversal.cpp)\n");
        exit(EXIT_FAILURE);
      }      
      // load 2 burst tries per 9-mer
      for ( int j = 0; j < 2; j++ )
      {
        // mini-burst trie exists
        if ( sizeoftries[j] != 0 )
        {
#ifdef see_binary_output
          if ( j == 0 ) cout << "forward burst-trie \n"; //TESTING
          if ( j == 1 ) cout << "reverse burst-trie \n"; //TESTING
#endif            
          // create a root trie node
          NodeElement newnode[4];             
          // copy the root trie node into the beginning of burst trie array
          memcpy( dst, &newnode[0], sizeof(NodeElement)*4);
          memset( dst, 0, sizeof(NodeElement)*4);
          if ( j == 0 ) lookup_tbl[i].trie_F = (NodeElement*)dst;
          else lookup_tbl[i].trie_R = (NodeElement*)dst;
          // queue to store the trie nodes as we create them
          deque<NodeElement*> nodes;
          nodes.push_back( (NodeElement*)dst );
          ((NodeElement *&) dst)+=4;
          // queue to store the flags of node elements given in the binary file
          deque<char> flags;
          // read the first trie node
          for ( int i = 0; i < 4; i++ )
          {
            char tmp;
            btrie.read(reinterpret_cast<char*>(&tmp), sizeof(char));
            flags.push_back(tmp);
#ifdef see_binary_output
            cout << " " << (int)tmp; //TESTING
#endif
          }
          // build the mini-burst trie
          while ( !nodes.empty() )
          {
            // ptr to traverse each trie node
            NodeElement* node = nodes.front(); 
            // trie node elements
            for ( int i = 0; i < 4; i++ )
            {
              unsigned char flag = flags.front();         
              // what does the node element point to
              switch( flag )
              {
                // set values to 0
                case 0:
                {
                  node->flag = 0;
                  node->size = 0;
                  node->whichnode.trie = NULL;
                }
                  break;
                // trie node
                case 1:
                {
                  // read the trie node
                  for ( int i = 0; i < 4; i++ )
                  {
                    char tmp;
                    btrie.read(reinterpret_cast<char*>(&tmp), sizeof(char));
#ifdef see_binary_output
                    cout << " " << (int)tmp; //TESTING
#endif
                    flags.push_back(tmp);
                  }               
                  node->flag = 1;
                  node->size = 0;
                  NodeElement newnode[4];
                  memcpy( (NodeElement*)dst, &newnode[0], sizeof(NodeElement)*4);
                  nodes.push_back( (NodeElement*)dst );
                  node->whichnode.trie = (NodeElement*)dst;
                  ((NodeElement *&) dst)+=4;            
                }
                  break;
                // bucket
                case 2:
                {
                  uint32_t sizeofbucket = 0;   
                  // read the bucket info
                  btrie.read(reinterpret_cast<char*>(&sizeofbucket), sizeof(uint32_t)); 
#ifdef see_binary_output
                  cout << "\tsizeofbucket = " << sizeofbucket; //TESTING
#endif                      
                  char* bucket = new char[sizeofbucket]();
                  if ( bucket == NULL )
                  {
                    fprintf(stderr,"\n  %sERROR%s: failed to allocate memory for allocate bucket (paralleltraversal.cpp)\n",startColor,"\033[0m");
                    exit(EXIT_FAILURE);
                  }         
                  btrie.read(reinterpret_cast<char*>(bucket), sizeofbucket);       
                  // copy the bucket into the burst trie array
                  memcpy( (void*)dst, (void*)bucket, sizeofbucket);     
                  delete [] bucket;
                  bucket = NULL;    
                  // assign pointers from trie node to the bucket
                  node->flag = flag;
                  node->whichnode.bucket = dst;
                  node->size = sizeofbucket;
                  dst = ((char *)dst)+sizeofbucket;        
                }
                  break;
                // ?
                default:
                {
                  fprintf(stderr, "\n  %sERROR%s: flag is set to %d (load_index)\n",startColor,"\033[0m",flag);
                  exit(EXIT_FAILURE);
                }
                  break;
              }
              flags.pop_front();
              node++;               
            }//~loop through 4 node elements in a trie node 
#ifdef see_binary_output
            cout << "\n"; //TESTING
#endif            
            nodes.pop_front(); 
          }//~while !nodes.empty()
        }//~if mini-burst trie exists
        else
        {
          if ( j == 0 ) lookup_tbl[i].trie_F = NULL;
          else lookup_tbl[i].trie_R = NULL;
        }
      }//~for both mini-burst tries
    }//~if ( sizeoftries != 0 )
    else
    {
      lookup_tbl[i].trie_F = NULL;
      lookup_tbl[i].trie_R = NULL;
    }
  }//~for all 9-mers in the look-up table
  btrie.close();
  // STEP 3: load the position reference tables (pos.dat)
  ifstream inreff( (char*)(ptr_dbindex_str + ".pos_" + part_str + ".dat").c_str(), ios::in | ios::binary );
  if ( !inreff.good() )
  {
    fprintf(stderr,"\n  ERROR: The database name '%s' does not exist.\n\n", (char*)(ptr_dbindex_str + ".pos_" + part_str + ".dat").c_str());
    exit(EXIT_FAILURE);
  }
  uint32_t size = 0; 
  inreff.read(reinterpret_cast<char*>(&number_elements), sizeof(uint32_t));
  positions_tbl = new kmer_origin[number_elements]();
  if ( positions_tbl == NULL )
  {
    fprintf(stderr,"  ERROR: could not allocate memory for positions_tbl (main(), paralleltraversal.cpp)\n");
    exit(EXIT_FAILURE);
  }
  for ( uint32_t i = 0; i < number_elements; i++ )
  {
    /* the number of positions */
    inreff.read(reinterpret_cast<char*>(&size), sizeof(uint32_t));
    positions_tbl[i].size = size;      
    /* the sequence seq_pos array */
    positions_tbl[i].arr = new seq_pos[size]();
    if ( positions_tbl[i].arr == NULL )
    {
      fprintf(stderr, "  ERROR: could not allocate memory for positions_tbl (paralleltraversal.cpp)\n");
      exit(EXIT_FAILURE);
    }
    inreff.read(reinterpret_cast<char*>(positions_tbl[i].arr), sizeof(seq_pos)*size);
  }
  inreff.close(); 
  return ;
}//~load_index()


/*
 *
 * @function load_ref: load fasta reference sequences into memory
 * for SW alignment
 * @return void
 * @version 1.0 June 12, 2013
 *
 *******************************************************************/
void
load_ref(char* ptr_dbfile,
         char* buffer,
         char** reference_seq,
         uint64_t* reference_seq_len,
         uint64_t seq_part_size,
         uint64_t numseq_part,
         uint64_t start_part,
         bool load_for_search)
{
  FILE *fp = fopen(ptr_dbfile,"r");
  if ( fp == NULL )
  {
    fprintf(stderr,"  %sERROR%s: [Line %d: %s] could not open file %s\n",
                   startColor,"\033[0m", __LINE__, __FILE__, ptr_dbfile);
    exit(EXIT_FAILURE);
  }
  // set the file pointer to the first sequence added to the index for this index file section
  if ( fseek(fp,start_part,SEEK_SET) != 0 )
    {
      fprintf(stderr,"  %sERROR%s: [Line %d: %s] could not locate the sequences used to construct the index.\n",
                     startColor,"\033[0m", __LINE__, __FILE__);
      fprintf(stderr,"  Check that your --ref <FASTA file, index name> correspond correctly for the FASTA file: %s.\n",
              ptr_dbfile);
    }
  // load references sequences into memory, skipping the new lines & spaces in the fasta format
  uint64_t num_seq_read = 0;
  char *s = buffer;
  int i = 0;
  int j = 0;
  char c = fgetc(fp);
  if ( load_for_search )
  {
    do
    {
      // the tag
      reference_seq[i++] = s;
      while ( c != '\n' )
      {
        *s++ = c;
        c = fgetc(fp);
      }
      // new line
      *s++ = c;
      if ( *s == '\n' )
      {
        fprintf(stderr,"  %sERROR%s: [Line %d: %s] your reference sequences are not in FASTA format "
                       "(there is an extra new line).",startColor,"\033[0m", __LINE__, __FILE__);
        exit(EXIT_FAILURE);
      }
      // the sequence
      reference_seq[i++] = s;
      c = fgetc(fp);
      do
      {
        if ( c != '\n' && c != ' ' )
        {
          // keep record of ambiguous character for alignment
          *s++ = nt_table[(int)c];
          // record the sequence length as we read it
          reference_seq_len[j]++;
        }
        c = fgetc(fp);
      } while ( (c != '>') && (c != EOF) );
      *s++ = '\n';
      j++;
      num_seq_read++;
    } while ( (num_seq_read != numseq_part) && (c != EOF) );
  }
  else
  {
    do
    {
      // the tag
      reference_seq[i++] = s;
      while ( c != '\n' )
      {
        *s++ = c;
        c = fgetc(fp);
      }
      // new line
      *s++ = c;
      // the sequence
      reference_seq[i++] = s;
      c = fgetc(fp);
      do
      {
        if ( c != '\n' && c != ' ' )
        {
          // keep record of ambiguous character for alignment
          *s++ = nt_table[(int)c];
        }
        c = fgetc(fp);
      } while ( (c != '>') && (c != EOF) );
      *s++ = '\n';
      j++;
      num_seq_read++;
    } while ( (num_seq_read != numseq_part) && (c != EOF) );
  }
    
  fclose(fp);
}//~load_ref