#ifndef INCLUDED_SLS_ALP_DATA
#define INCLUDED_SLS_ALP_DATA

/* $Id: $
* ===========================================================================
*
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's offical duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*
*  Although all reasonable efforts have been taken to ensure the accuracy
*  and reliability of the software and data, the NLM and the U.S.
*  Government do not and cannot warrant the performance or results that
*  may be obtained by using this software or data. The NLM and the U.S.
*  Government disclaim all warranties, express or implied, including
*  warranties of performance, merchantability or fitness for any particular
*  purpose.
*
*  Please cite the author in any work or product based on this material.
*
* ===========================================================================*/

/*****************************************************************************

File name: sls_alp_data.hpp

Author: Sergey Sheetlin

Contents: Contains input data

******************************************************************************/


//if _MSDOS_ is defined then it is Windows program, if not - then it is UNIX program
//#ifndef _MSDOS_
//#define _MSDOS_
//#endif

#include <complex>
#include <iostream>
#include <map>
#include <vector>
#include <fstream>
#include <float.h>
#include <ctime>
#include <stdlib.h>
#include <limits>
#include <climits>
#include <cstring>
#include <cstdlib>


#ifndef _MSDOS_ //UNIX program
#include <sys/time.h>
#else
#include <sys/timeb.h>

#define _CRTDBG_MAP_ALLOC
#include <crtdbg.h>

#endif

#include "gumbel_params.hpp"

#include "sls_alp_regression.hpp"

#include "njn_random.hpp"
#include "njn_uniform.hpp"



const double mb_bytes=1048576.0;

namespace Sls {

	static long int small_long=(long int)((double)LONG_MIN/2.0);
	static double dbl_max_log=log(DBL_MAX);

	struct struct_for_randomization
	{
		long int d_random_factor;
		vector<long int> d_first_stage_preliminary_realizations_numbers_ALP;
		vector<long int> d_preliminary_realizations_numbers_ALP;
		vector<long int> d_preliminary_realizations_numbers_killing;
		long int d_total_realizations_number_with_ALP;
		long int d_total_realizations_number_with_killing;
	};



	struct error//struct to handle exceptions
	{
		std::string st;
		error(std::string st_,long int error_code_){st=st_;error_code=error_code_;};
		long int error_code;
		//if=0 Result was calculated and returned

		//if =1 //Computation stoped because time
			// or memory requirements exceeded
			// user-specified thresholds

		//if =2 //Computation stopped due to different
			// reasons than time or memory
			// repeating computation with the same
			// input parameters may be successful
		//if =3 //Result can not be computed for current
			// input parameters

		//if =4 //Other cases
		//if =41 //memory allocation error
	};

	struct error_for_single_realization//struct to handle exceptions during calclation of single realization
	{
		std::string st;
		error_for_single_realization(){};
	};


	struct data_for_lambda_equation//struct for lambda_equation
	{
		long int d_number_of_AA;//number of AA
		long int** d_smatr;//scoring matrix
		double *d_RR1;//AA probabilities
		double *d_RR2;//AA probabilities
	};


	class alp_data;

	template<typename T> class array_positive{
	public:
		array_positive(alp_data *alp_data_)// constructor
		{ 
			d_elem=NULL;
			d_alp_data=alp_data_; 
			if(!d_alp_data)
			{
				throw error("Unexpected error",4);
			};
			d_dim=-1;
			d_step=10;
		};

		~array_positive();


		void increment_array();
		

		inline void set_elem(
			long int ind_,
			T elem_)
		{
			while(ind_>d_dim)
			{
				increment_array();
			};

			d_elem[ind_]=elem_;
		};

		inline void increase_elem_by_1(
			long int ind_)
		{
			while(ind_>d_dim)
			{
				increment_array();
			};

			d_elem[ind_]++;
		};

		inline void increase_elem_by_x(
			long int ind_,
			T x_)
		{
			while(ind_>d_dim)
			{
				increment_array();
			};

			d_elem[ind_]+=x_;
		};



	public:
		
		long int d_step;
		long int d_dim;//dimension of the array is d_dim+1
		T * d_elem;
		alp_data *d_alp_data;//initial data
	};


	template<typename T> class array{
	public:
		array(alp_data *alp_data_)// constructor
		{ 
			d_elem=NULL;
			d_alp_data=alp_data_; 
			d_dim=-1;
			d_ind0=0;
			d_step=10;
			d_dim_plus_d_ind0=d_dim+d_ind0;
		};   

		~array();

		void increment_array_on_the_rigth();

		void increment_array_on_the_left();


		inline void set_elem(
			long int ind_,
			T elem_)
		{

			while(ind_>d_dim_plus_d_ind0)
			{
				increment_array_on_the_rigth();
			};

			while(ind_<d_ind0)
			{
				increment_array_on_the_left();
			};

			d_elem[ind_-d_ind0]=elem_;
		};

		inline void increase_elem_by_1(
			long int ind_)
		{
			while(ind_>d_dim_plus_d_ind0)
			{
				increment_array_on_the_rigth();
			};

			while(ind_<d_ind0)
			{
				increment_array_on_the_left();
			};

			d_elem[ind_-d_ind0]++;
		};

		
	public:
		
		long int d_step;
		long int d_dim;//dimension of the array is d_dim+1
		long int d_ind0;//the leftmost index of the array
		long int d_dim_plus_d_ind0;
		T * d_elem;
		alp_data *d_alp_data;//initial data
	};


	struct q_elem
	{
		long int d_a;
		long int d_b;
	};

	class importance_sampling{

	public:
		importance_sampling(
		alp_data *alp_data_,
		long int open_,
		long int open1_,
		long int open2_,

		long int epen_,
		long int epen1_,
		long int epen2_,

		long int number_of_AA_,
		long int **smatr_,
		double *RR1_,
		double *RR2_);

		double d_mu;
		double d_nu;
		double d_eta;
		double d_mu_SI;
		double d_mu_DS;
		double d_mu_ID;
		double d_mu_IS;
		double d_mu_SD;
		q_elem * d_elements;
		double * d_elements_values;


		double d_for_D[3];
		double d_for_I[2];
		double d_for_S[3];

		char d_for_D_states[3];
		char d_for_I_states[2];
		char d_for_S_states[3];

		double **d_exp_s;
		double d_lambda;
		double d_ungap_lambda;



		~importance_sampling();

		static double lambda_equation(double x_,void* func_number_);


		long int d_is_number_of_AA;
		alp_data *d_alp_data;//initial data

	};



	class alp_data{

	

	public:

		alp_data(//constructor
			long int rand_,//randomization number
			string randout_,//if defined, then the program outputs complete randomization information into a file

			long int open_,//gap opening penalty
			long int open1_,//gap opening penalty for a gap in the sequence #1
			long int open2_,//gap opening penalty for a gap in the sequence #2

			long int epen_,//gap extension penalty
			long int epen1_,//gap extension penalty for a gap in the sequence #1
			long int epen2_,//gap extension penalty for a gap in the sequence #2

			string smatr_file_name_,//scoring matrix file name
			string RR1_file_name_,//probabilities1 file name
			string RR2_file_name_,//probabilities2 file name
			double max_time_,//maximum allowed calculation time in seconds
			double max_mem_,//maximum allowed memory usage in MB
			double eps_lambda_,//relative error for lambda calculation
			double eps_K_,//relative error for K calculation
			bool insertions_after_deletions_,//if true, then insertions after deletions are allowed
			long int match,//NEW
			long int mismatch,//NEW
			double A_,//NEW
			double C_,//NEW
			double G_,//NEW
			double T_//NEW
			);


		~alp_data();//destructor

		inline double ran2()//generates the next random value
		{
			return Njn::Uniform::variate <double> (0,1);

			//double rand_C=(double)((double)rand()/(double)RAND_MAX);
			//return rand_C;	
		};

		static void read_smatr(
			string smatr_file_name_,
			long int **&smatr_,
			long int &number_of_AA_smatr_);

		void check_out_file(
			string out_file_name_);



		static double round(//returns nearest integer to x_
			const double &x_);

		static string long_to_string(//convert interer ot string
			long int number_);

		static char digit_to_string(//convert interer ot string
			long int digit_);

		static void get_current_time(
			double &seconds_);

				static bool the_value_is_double(
				string str_,
				double &val_);

				static bool the_value_is_long(
				string str_,
				long int &val_);





		static void read_RR(
			string RR_file_name_,
			double *&RR_,
			double *&RR_sum_,
			long int *&RR_sum_elements_,
			long int &number_of_AA_RR_);

		double get_allocated_memory_in_MB();

		static void assert_mem(void *pointer_);

	

		template<typename T>
		static void get_memory_for_matrix(
		long int dim_,
		T ** &matr_)
		{
			matr_=NULL;
			bool ee_error_flag=false;
			error ee_error("",0);

			try
			{
			try
			{

				long int i;
				matr_=new T *[dim_];
				assert_mem(matr_);

				for(i=0;i<dim_;i++)
				{
					matr_[i]=NULL;
				};

				for(i=0;i<dim_;i++)
				{
					matr_[i]=new T [dim_];
					assert_mem(matr_[i]);
				};
				//d_memory_size_in_MB+=(double)sizeof(T)*(double)dim_*(double)dim_/mb_bytes;

			}
			catch (error er)
			{
				ee_error_flag=true;
				ee_error=er;		
			};
			}
			catch (...)
			{ 
				ee_error_flag=true;
				ee_error=error("Internal error in the program\n",4);
			};

			//memory release

			if(ee_error_flag)
			{

				if(matr_)
				{
					long int i;
					for(i=0;i<dim_;i++)
					{
						if(matr_[i])
						{
							delete[]matr_[i];matr_[i]=NULL;
						};
					};

					delete[]matr_;matr_=NULL;
				};

				throw error(ee_error.st,ee_error.error_code);
			};

		};

		template<typename T>
		static void delete_memory_for_matrix(
		long int dim_,
		T ** &matr_)
		{
			long int i;
			if(matr_)
			{
				for(i=0;i<dim_;i++)
				{
					delete []matr_[i];matr_[i]=NULL;
				};
				delete []matr_;matr_=NULL;
			};

			//d_memory_size_in_MB-=(double)sizeof(T)*(double)dim_*(double)dim_/mb_bytes;
		};

	static long int random_long(
	double value_,
	long int dim_);


	template<typename T>
	static T random_long(
	double value_,
	long int dim_,
	double *sum_distr_,
	T* elements_)//sum_distr_[dim_-1] must be equal to 1
	{
		if(value_<0||value_>1)	
		{
			throw error("Unexpected error in alp_data::random_long\n",4);
		};

		long int v1=0;
		long int v2=dim_;

		while(v2-v1>1)
		{
			long int v3=(long int)(round(double(v2+v1)/2.0));
			/*
			if(sum_distr_[v3-1]==value_)
			{
				v1=v3-1;
				v2=v3;
				break;
			};
			*/

			if(sum_distr_[v3-1]>=value_)
			{
				v2=v3;
			}
			else
			{
				v1=v3;
			};
		};

		return elements_[v2-1];

	};


	template<class T>
	static inline T Tmax(T i_, T j_)
	{
		if(i_>j_)
		{
			return i_;
		};
		return j_;
	};

	template<class T>
	static inline T Tmin(T i_, T j_)
	{
		if(i_<j_)
		{
			return i_;
		};
		return j_;
	};


	template<class T>
	static inline T Tmax(T x_,T y_,T z_)
	{
		return Tmax(Tmax(x_,y_),z_);
	};

	template<class T>
	static inline T Tmin(T x_,T y_,T z_)
	{
		return Tmin(Tmin(x_,y_),z_);
	};

	template<class T>
	static inline T Tmax(T x_,T y_,T z_,T w_)
	{
		return Tmax(Tmax(x_,y_),Tmax(z_,w_));
	};

	template<class T>
	static inline T Tmin(T x_,T y_,T z_,T w_)
	{
		return Tmin(Tmin(x_,y_),Tmin(z_,w_));
	};


	



	public:


	
	//input parameters
	long int d_open;//gap opening penalty
	long int d_open1;//gap opening penalty for a gap in the sequence #1
	long int d_open2;//gap opening penalty for a gap in the sequence #2


	long int d_epen;//gap extension penalty
	long int d_epen1;//gap extension penalty for a gap in the sequence #1
	long int d_epen2;//gap extension penalty for a gap in the sequence #2


	double d_max_time;//maximum allowed calculation time in seconds
	double d_max_mem;//maximum allowed memory usage in MB
	double d_eps_lambda;//relative error for lambda calculation
	double d_eps_K;//relative error for K calculation

	bool d_insertions_after_deletions;//if true, then insertions after deletions are allowed

	long int d_match;//NEW
	long int d_mismatch;//NEW
	double d_A_;//NEW
	double d_C_;//NEW
	double d_G_;//NEW
	double d_T_;//NEW

	
	//additional parameters

	bool d_smatr_symmetric_flag;//true if the scoring matrix is symmetric

	long int d_number_of_AA;//number of AA
	long int d_number_of_AA_smatr;

	long int** d_smatr;//scoring matrix

	double *d_RR1;//AA probabilities
	double *d_RR1_sum;//probability distribution function for d_RR
	long int *d_RR1_sum_elements;//numbers of AA corresponded to d_RR

	double *d_RR2;//AA probabilities
	double *d_RR2_sum;//probability distribution function for d_RR
	long int *d_RR2_sum_elements;//numbers of AA corresponded to d_RR

	long int d_random_factor;
	string d_randout;//if defined, then the program outputs complete randomization information into a file


	double d_memory_size_in_MB;//approximate current allocated memory size

	importance_sampling *d_is;//data for the importance sampling

	double *d_r_i_dot;
	double *d_r_dot_j;

	long int d_minimum_realizations_number;

	bool d_sentinels_flag;

	

	//for debugging
	long int d_dim1_tmp;
	long int d_dim2_tmp;

	long int d_realizations_number2;

	double d_time_before1;

	struct_for_randomization *d_rand_all;
	bool d_rand_flag;
	



private:

	#ifndef _MSDOS_ //UNIX program

	#else
		_CrtMemState d_s1, d_s2, d_s3;
	#endif

	


	};

	//array_positive functions
	template<class T>
	array_positive<T>::~array_positive()
	{
		delete[]d_elem;d_elem=NULL;
		if(d_alp_data)
		{
			d_alp_data->d_memory_size_in_MB-=(double)sizeof(T)*(double)(d_dim+1)/mb_bytes;
		};

	};


	template<class T>
	void array_positive<T>::increment_array()
	{
		bool ee_error_flag=false;
		error ee_error("",0);
		T *d_elem_new=NULL;

		try
		{
		try
		{

			d_dim+=d_step;

			d_elem_new=new T[d_dim+1];
			alp_data::assert_mem(d_elem_new);

			long int i;
			for(i=0;i<d_dim+1-d_step;i++)
			{
				d_elem_new[i]=d_elem[i];
			};

			for(i=d_dim+1-d_step;i<d_dim+1;i++)
			{
				d_elem_new[i]=0;
			};


			delete[]d_elem;d_elem=NULL;
			if(d_alp_data)
			{
				d_alp_data->d_memory_size_in_MB+=(double)sizeof(T)*(double)d_step/mb_bytes;
			};

			d_elem=d_elem_new;d_elem_new=NULL;

		}
		catch (error er)
		{
			ee_error_flag=true;
			ee_error=er;
		};
		}
		catch (...)
		{ 
			ee_error_flag=true;
			ee_error=error("Internal error in the program\n",4);
		};

		//memory release

		if(ee_error_flag)
		{
			delete[]d_elem_new;d_elem_new=NULL;
			throw error(ee_error.st,ee_error.error_code);
		};
		
	};

	//array functions

	template<class T>
	array<T>::~array()
	{
		delete[]d_elem;d_elem=NULL;
		if(d_alp_data)
		{
			d_alp_data->d_memory_size_in_MB-=(double)sizeof(T)*(double)(d_dim+1)/mb_bytes;
		};

	};

	template<class T>
	void array<T>::increment_array_on_the_rigth()
	{
		bool ee_error_flag=false;
		error ee_error("",0);
		T *d_elem_new=NULL;

		try
		{
		try
		{


			d_dim+=d_step;

			d_elem_new=new T[d_dim+1];
			alp_data::assert_mem(d_elem_new);

			long int i;
			for(i=0;i<d_dim+1-d_step;i++)
			{
				d_elem_new[i]=d_elem[i];
			};

			for(i=d_dim+1-d_step;i<d_dim+1;i++)
			{
				d_elem_new[i]=0;
			};

			d_dim_plus_d_ind0=d_dim+d_ind0;

			if(d_alp_data)
			{
				d_alp_data->d_memory_size_in_MB+=(double)sizeof(T)*(double)d_step/mb_bytes;
			};


			delete[]d_elem;d_elem=NULL;
			d_elem=d_elem_new;d_elem_new=NULL;


		}
		catch (error er)
		{
			ee_error_flag=true;
			ee_error=er;		
		};
		}
		catch (...)
		{ 
			ee_error_flag=true;
			ee_error=error("Internal error in the program\n",4);
		};

		//memory release

		if(ee_error_flag)
		{
			delete[]d_elem_new;d_elem_new=NULL;
			throw error(ee_error.st,ee_error.error_code);
		};

	};

	template<class T>
		void array<T>::increment_array_on_the_left()
	{
		bool ee_error_flag=false;
		error ee_error("",0);
		T *d_elem_new=NULL;

		try
		{
		try
		{
			d_dim+=d_step;
			d_ind0-=d_step;

			d_elem_new=new T[d_dim+1];
			alp_data::assert_mem(d_elem_new);

			long int i;

			for(i=0;i<d_step;i++)
			{
				d_elem_new[i]=0;
			};

			for(i=0;i<d_dim+1-d_step;i++)
			{
				d_elem_new[i+d_step]=d_elem[i];
			};

			if(d_alp_data)
			{
				d_alp_data->d_memory_size_in_MB+=(double)sizeof(T)*(double)d_step/mb_bytes;
			};

			delete[]d_elem;d_elem=NULL;
			d_elem=d_elem_new;d_elem_new=NULL;


		}
		catch (error er)
		{
			ee_error_flag=true;
			ee_error=er;		
		};
		}
		catch (...)
		{
			ee_error_flag=true;
			ee_error=error("Internal error in the program\n",4);
		};

		//memory release

		if(ee_error_flag)
		{
			delete[]d_elem_new;d_elem_new=NULL;
			throw error(ee_error.st,ee_error.error_code);
		};

	};


	}


#endif
