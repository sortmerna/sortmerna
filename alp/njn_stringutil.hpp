#ifndef INCLUDED_NJN_STRINGUTIL
#define INCLUDED_NJN_STRINGUTIL

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

File name: njn_stringutil.hpp

Author: John Spouge

Contents: 

******************************************************************************/


#include <vector>
#include <sstream>
#include <cstring>
using namespace std;

namespace Njn {
	namespace StringUtil {


        bool isAlpha (const char *symbols_); // Is every symbol in symbols_ alphabetic ?
        bool isNoWhiteSpace (const char *symbols_); // Is there no white space ?
        bool isLower (const char *symbols_); // Is every symbol in symbols_ lower-case ?
        bool isUpper (const char *symbols_); // Is every symbol in symbols_ upper-case ?

        void toLower (char *symbols_);
        void toUpper (char *symbols_);

        void eraseInitialWhiteSpace (char *symbols_); // removes initial white-space
        void eraseFinalWhiteSpace (char *symbols_); // removes final white-space
        void eraseWhiteSpace (char *symbols_, bool eraseBlankOnly_ = false); // removes white-space
        void eraseCarriageReturn (char *symbols_); // removes anomalous UNIX carriage return

        void whiteSpace2UnderScore (char *symbols_); // puts single underscore between words

        void eraseInitialChar (char *symbols_, const char *c_); // removes initial char's c_
        void eraseFinalChar (char *symbols_, const char *c_); // removes final char's c_
        void eraseChar (char *symbols_, const char *c_); // removes char's c_
        void substituteChar (char *symbols_, const char cOut_, const char cIn_); // replaces all occurrences of cIn_ with cOut_

        void split (std::vector <std::string> *strVec, const std::string &str_, const std::string &split_); 
        size_t splitCount (const std::string &str_, const std::string &split_); 
        // splits str_ at every occurrence of split_, throwing away split_ and putting the pieces into strVec_

        inline bool isAlpha (const std::string &symbols_); // Is every symbol in symbols_ alphabetic ?
        inline bool isNoWhiteSpace (const std::string &symbols_); // Is there no white space ?
        bool isAllWhiteSpace (const std::string &symbols_); // Is everything white space ?
        inline bool isLower (const std::string &symbols_); // Is every symbol in symbols_ lower-case ?
        inline bool isUpper (const std::string &symbols_); // Is every symbol in symbols_ upper-case ?

        inline void toLower (std::string &symbols_);
        inline void toUpper (std::string &symbols_);

        bool isReplicate ( // Are some of the symbols_ replicated ?
        const std::string &symbols_,
        bool upperEqualsLower_ = false); // Does case matter?, i.e., a == A ?

        void eraseInitialWhiteSpace (std::string &symbols_); // removes initial white-space
        void eraseFinalWhiteSpace (std::string &symbols_); // removes final white-space
        void eraseWhiteSpace (std::string &symbols_, bool eraseBlankOnly_ = false); // removes white-space
        void eraseCarriageReturn (std::string &symbols_); // removes anomalous UNIX carriage return

        void eraseInitialChar (std::string &symbols_, const std::string &c_); // removes initial char's c_
        void eraseFinalChar (std::string &symbols_, const std::string &c_); // removes final char's c_
        void eraseChar (std::string &symbols_, const std::string &c_); // removes char's c_
        void substituteChar (std::string &symbols_, const char cOut_, const char cIn_); // replaces all occurrences of cIn_ with cOut_

        void whiteSpace2UnderScore (std::string &symbols_); // puts single underscore between words

		}
	}

//
// There are no more declarations beyond this point.
//

namespace Njn {
	namespace StringUtil {


        bool isAlpha (const std::string &symbols_) // Is every symbol in symbols_ alphabetic ?
        {
            return isAlpha (symbols_.c_str ());
        }

        bool isNoWhiteSpace (const std::string &symbols_) // Is there no white space ?
        {
            return isNoWhiteSpace (symbols_.c_str ());
        }

        bool isLower (const std::string &symbols_) // Is every symbol in symbols_ lower-case ?
        {
            return isLower (symbols_.c_str ());
        }

        bool isUpper (const std::string &symbols_) // Is every symbol in symbols_ upper-case ?
        {
            return isUpper (symbols_.c_str ());
        }

        void toLower (std::string &symbols_)
        {
            for (std::string::iterator i = symbols_.begin (); i != symbols_.end (); i++) *i = tolower (*i);
        }

        void toUpper (std::string &symbols_)
        {
            for (std::string::iterator i = symbols_.begin (); i != symbols_.end (); i++) *i = toupper (*i);
        }


		}
	}

#endif //! INCLUDED
