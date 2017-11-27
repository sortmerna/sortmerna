/**
* FILE: output.cpp
* Created: Nov 26, 2017 Sun
*/
#include "output.hpp"



void Output::init(Readstats & readstats)
{
	// attach pid to output files
	char pidStr[4000];
	if (pid_gv)
	{
		int32_t pid = getpid();
		sprintf(pidStr, "%d", pid);
	}

	// associate the streams with reference sequence file names
	if (opts.ptr_filetype_ar != NULL)
	{
		if (fastxout_gv)
		{
			// fasta/fastq output
			acceptedstrings.assign(opts.ptr_filetype_ar);
			if (pid_gv)
			{
				acceptedstrings.append("_");
				acceptedstrings.append(pidStr);
			}
			acceptedstrings.append(".");
			acceptedstrings.append(readstats.suffix.c_str());

			acceptedreads.open(acceptedstrings);
			acceptedreads.close();
		}

		if (opts.samout)
		{
			// sam output
			acceptedstrings_sam.assign(opts.ptr_filetype_ar);
			if (pid_gv)
			{
				acceptedstrings_sam.append("_");
				acceptedstrings_sam.append(pidStr);
			}
			acceptedstrings_sam.append(".sam");
			acceptedsam.open(acceptedstrings_sam);
			acceptedsam.close();
		}

		if (opts.blastout)
		{
			// blast output
			acceptedstrings_blast.assign(opts.ptr_filetype_ar);
			if (pid_gv)
			{
				acceptedstrings_blast.append("_");
				acceptedstrings_blast.append(pidStr);
			}
			acceptedstrings_blast.append(".blast");
			acceptedblast.open(acceptedstrings_blast);
			acceptedblast.close();
		}

		if (logout_gv)
		{
			// statistics file output
			ofstream logstream;
			logoutfile.assign(opts.ptr_filetype_ar);
			if (pid_gv)
			{
				logoutfile.append("_");
				logoutfile.append(pidStr);
			}
			logoutfile.append(".log");

			logstream.open(logoutfile);
			logstream.close();
		}

		if (otumapout_gv)
		{
			// OTU map output file
			ofstream otumap;
			acceptedotumap_file.assign(opts.ptr_filetype_ar);
			if (pid_gv)
			{
				acceptedotumap_file.append("_");
				acceptedotumap_file.append(pidStr);
			}
			acceptedotumap_file.append("_otus.txt");
			otumap.open(acceptedotumap_file);
			otumap.close();
		}

		if (de_novo_otu_gv)
		{
			ofstream denovo_otu;
			denovo_otus_file.assign(opts.ptr_filetype_ar);
			if (pid_gv)
			{
				denovo_otus_file.append("_");
				denovo_otus_file.append(pidStr);
			}
			denovo_otus_file.append("_denovo.");
			denovo_otus_file.append(readstats.suffix.c_str());

			denovo_otu.open(denovo_otus_file);
			denovo_otu.close();
		}
	}//~if ( ptr_filetype_ar != NULL ) 

	if (opts.ptr_filetype_or != NULL)
	{
		if (fastxout_gv)
		{
			// output stream for other reads
			ofstream otherreads;
			// add suffix database name to accepted reads file
			if (pid_gv)
			{
				strcat(opts.ptr_filetype_or, "_");
				strcat(opts.ptr_filetype_or, pidStr);
			}
			strcat(opts.ptr_filetype_or, ".");
			strcat(opts.ptr_filetype_or, readstats.suffix.c_str());
			// create the other reads file
			otherreads.open(opts.ptr_filetype_or);
			otherreads.close();
		}
	}
} // ~Output::init