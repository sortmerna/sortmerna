<!-- dx-header -->
# App for running Sortmerna (DNAnexus Platform App)

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://wiki.dnanexus.com/.
<!-- /dx-header -->

This application is for running SortmeRNA 3.
It can be launched either from DNANexus GUI or on command line using 'dx' utility.

Required parameters:

**REFS** a set of file with Reference sequences
**READS** or **READS_GZ** the file with Reads to be aligned.

# Command Line example

```
dx run sortmerna-3.run (*launch the App*)
Entering interactive mode for input selection.

Input:   REFS (REFS)
Class:   array:file

Enter file values, one at a time (^D or <ENTER> to finish, <TAB> twice for compatible files in
current directory, '?' for more options)
REFS[0]: file-FJZz7300V8Pb4KZgBGkXbZBq
REFS[1]: file-FJZz7Fj0V8PX7G23K0Yf2Z2J
REFS[2]: file-FJZz78j0V8Pp0Y6Z07zK7K86
REFS[3]:

Select an optional parameter to set by its # (^D or <ENTER> to finish):

 [0] READS (READS)
 [1] READS_GZ (READS_GZ)
 [2] SAM (SAM) [default=false]
 [3] FASTX (FASTX) [default=true]
 [4] BLAST (BLAST) [default="1 cigar qcov qstrand"]
 [5] advanced (advanced) [default="--paired_out -v"]
 [6] (-h) Sortmerna help (help_smr) [default=false]

Optional param #: 0

Input:   READS (READS)
Class:   file

Enter file ID or path (<TAB> twice for compatible files in current directory, '?' for
more options)
READS: file-FJbP2xj0V8Pqp9VgJjV0jbyq

Select an optional parameter to set by its # (^D or <ENTER> to finish):

 [0] READS (READS) [={"$dnanexus_link": "file-FJbP2xj0V8Pqp9VgJjV0jbyq"}]
 [1] READS_GZ (READS_GZ)
 [2] SAM (SAM) [default=false]
 [3] FASTX (FASTX) [default=true]
 [4] BLAST (BLAST) [default="1 cigar qcov qstrand"]
 [5] advanced (advanced) [default="--paired_out -v"]
 [6] (-h) Sortmerna help (help_smr) [default=false]

Optional param #:

Using input JSON:
{
    "REFS": [
        {
            "$dnanexus_link": "file-FJZz7300V8Pb4KZgBGkXbZBq"
        },
        {
            "$dnanexus_link": "file-FJZz7Fj0V8PX7G23K0Yf2Z2J"
        },
        {
            "$dnanexus_link": "file-FJZz78j0V8Pp0Y6Z07zK7K86"
        }
    ],
    "READS": {
        "$dnanexus_link": "file-FJbP2xj0V8Pqp9VgJjV0jbyq"
    }
}

Confirm running the executable with this input [Y/n]:
Calling app-FKF6Jq00Y21JkJBvKG0xG3Bb with output destination project-FGyB2f00g35X2Pvv891Gf7JK:/

Job ID: job-FKF6q180g35y63fF4kZy0F8Q
Watch launched job now? [Y/n]

Job Log
-------
Watching job job-FKF6q180g35y63fF4kZy0F8Q. Press Ctrl+C to stop.
```