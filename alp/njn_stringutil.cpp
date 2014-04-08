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

File name: njn_stringutil.cpp

Author: John Spouge

Contents: 

******************************************************************************/

#include "njn_stringutil.hpp"
#include <assert.h>
#include <algorithm>




using namespace Njn;


bool StringUtil::isAlpha (const char *symbols_) // Is every symbol in symbols_ alphabetic ?
{
    const char *c = '\0';
    for (c = symbols_; *c; c++) {
        if (! isalpha (*c)) return false;
    }
    assert (! *c); // should be at the end of the line

    return true; 
}

bool StringUtil::isNoWhiteSpace (const char *symbols_) // Is there no white space ?
{
    const char *c = '\0';

    for (c = symbols_; *c; c++) { // find next whitespace
        if (isspace (*c)) return false;
    }
    assert (! *c); // should be at the end of the line

    return true;
}

bool StringUtil::isAllWhiteSpace (const std::string &symbols_) // Is everything white space ?
{
    for (string::const_iterator c = symbols_.begin (); c != symbols_.end (); c++) { // find next whitespace
        if (! isspace (*c)) return false;
    }

    return true;
}

bool StringUtil::isLower (const char *symbols_) // Is every symbol in symbols_ lower-case ?
{
    const char *c = '\0';
    for (c = symbols_; *c; c++) {
        if (! islower (*c)) return false;
    }
    assert (! *c); // should be at the end of the line

    return true; 
}

bool StringUtil::isUpper (const char *symbols_) // Is every symbol in symbols_ upper-case ?
{
    const char *c = '\0';
    for (c = symbols_; *c; c++) {
        if (! isupper (*c)) return false;
    }
    assert (! *c); // should be at the end of the line

    return true; 
}

void StringUtil::toLower (char *symbols_)
{
    for (char *i = symbols_; *i; i++) *i = tolower (*i);
}

void StringUtil::toUpper (char *symbols_)
{
    for (char *i = symbols_; *i; i++) *i = toupper (*i);
}

bool StringUtil::isReplicate ( // Are some of the symbols_ replicated ?
const std::string &symbols_,
bool upperEqualsLower_) // Does case matter?, i.e., a == A ?
{
    string symbols (symbols_);
    if (upperEqualsLower_) {
        for (string::iterator i = symbols.begin (); i != symbols.end (); i++) {
            *i = toupper (*i);
        }
    }
    sort <string::iterator> (symbols.begin (), symbols.end ());
    return unique <string::iterator> (symbols.begin (), symbols.end ()) != symbols.end (); 
}

void StringUtil::eraseInitialWhiteSpace (char *symbols_)
// returns a char * symbols_ without initial white-space
{ 
    char *i = '\0';
    char *j = '\0';


    for (i = symbols_; *i && isspace (*i); i++) ;
  
    if (i == symbols_) return;

    for (j = symbols_; *i; i++, j++) 
    {
        *j = *i; // copy and skip white space
    }

    *j = '\0';
}

void StringUtil::eraseInitialWhiteSpace (std::string &symbols_) // erases initial white-space
{
    char *str = new char [symbols_.size () + 1];

    strcpy (str, symbols_.c_str ());
    eraseInitialWhiteSpace (str);
    symbols_ = str;

    delete [] str; str = 0;
}

void StringUtil::eraseFinalWhiteSpace (char *symbols_)
// returns a char * symbols_ without initial white-space
{
    reverse (symbols_, symbols_ + strlen (symbols_));
    eraseInitialWhiteSpace (symbols_);
    reverse (symbols_, symbols_ + strlen (symbols_));
}

void StringUtil::eraseFinalWhiteSpace (std::string &symbols_) // erases initial white-space
{
    char *str = new char [symbols_.size () + 1];

    strcpy (str, symbols_.c_str ());
    eraseFinalWhiteSpace (str);
    symbols_ = str;

    delete [] str; str = 0;
}

void StringUtil::eraseWhiteSpace (char *symbols_, bool eraseBlankOnly_)
// returns a char * symbols_ without white-space
{ 
    char *i = '\0';
    char *j = '\0';

    for (i = j = symbols_; *i; i++) {
        if (! isspace (*i) || *i != ' ' && eraseBlankOnly_) { // copy unless white space
            *j = *i; 
            j++;
        }
    }
    *j = '\0';
}

void StringUtil::eraseWhiteSpace (std::string &symbols_, bool eraseBlankOnly_) // erases white-space
{
    char *str = new char [symbols_.size () + 1];

    strcpy (str, symbols_.c_str ());
    eraseWhiteSpace (str, eraseBlankOnly_);
    symbols_ = str;

    delete [] str; str = 0;
}

void StringUtil::eraseInitialChar (char *symbols_, const char *c_) // removes initial char's c_
{ 
    char *i = '\0';
    char *j = '\0';
    const char *s = '\0';

    for (i = symbols_; *i; i++) 
    {
        bool found = false;

        for (s = c_; *s; s++)
        {
            if (*i == *s)
            {
                found = true;
                break;
            }
        }

        if (! found) break;
    }
  
    if (i == symbols_) return;

    for (j = symbols_; *i; i++, j++) 
    {
        *j = *i; // copy and skip white space
    }

    *j = '\0';
}

void StringUtil::eraseInitialChar (std::string &symbols_, const std::string &c_) // removes initial char's c_
{
    char *str = new char [symbols_.size () + 1];

    strcpy (str, symbols_.c_str ());
    eraseInitialChar (str, c_.c_str ());
    symbols_ = str;

    delete [] str; str = 0;
}

void StringUtil::eraseFinalChar (char *symbols_, const char *c_) // removes final char's c_
{
    reverse (symbols_, symbols_ + strlen (symbols_));
    eraseInitialChar (symbols_, c_);
    reverse (symbols_, symbols_ + strlen (symbols_));
}

void StringUtil::eraseFinalChar (std::string &symbols_, const std::string &c_) // removes final char's c_
{
    char *str = new char [symbols_.size () + 1];

    strcpy (str, symbols_.c_str ());
    eraseFinalChar (str, c_.c_str ());
    symbols_ = str;

    delete [] str; str = 0;
}

void StringUtil::eraseChar (char *symbols_, const char *c_) // removes char's c_
{ 
    char *i = '\0';
    char *j = '\0';
    const char *s = '\0';

    for (i = j = symbols_; *i; i++) 
    {
        bool found = false;

        for (s = c_; *s; s++)
        {
            if (*i == *s)
            {
                found = true;
                break;
            }
        }

        if (! found) 
        {
            *j = *i; // copy 
            j++;
        }
    }
    *j = '\0';
}

void StringUtil::eraseChar (std::string &symbols_, const std::string &c_) // removes char's c_
{
    char *str = new char [symbols_.size () + 1];

    strcpy (str, symbols_.c_str ());
    eraseChar (str, c_.c_str ());
    symbols_ = str;

    delete [] str; str = 0;
}

void StringUtil::substituteChar (char *symbols_, char cOut_, const char cIn_) // replaces all occurrences of cIn_ with cOut_
{ 
    char *i = '\0';
    const char *s = '\0';

    for (i = symbols_; *i; i++) 
    {
        if (*i == cIn_) *i = cOut_; 
    }
}

void StringUtil::substituteChar (std::string &symbols_, char cOut_, char cIn_) // replaces all occurrences of cIn_ with cOut_
{
    char *str = new char [symbols_.size () + 1];

    strcpy (str, symbols_.c_str ());
    substituteChar (str, cOut_, cIn_);
    symbols_ = str;

    delete [] str; str = 0;
}

void StringUtil::eraseCarriageReturn (char *symbols_) // erases white-space
{
	// Erase the first carriage return character '\r',
	//     to avoid problems when transferring files between DOS and UNIX.
    if (! symbols_ || ! *symbols_) return;
    if (symbols_ [strlen (symbols_) - 1] == '\r') symbols_ [strlen (symbols_) - 1] = '\0';
}

void StringUtil::eraseCarriageReturn (std::string &symbols_) // erases white-space
{
	// Erase the first carriage return character '\r',
	//     to avoid problems when transferring files between DOS and UNIX.
    if (symbols_.empty ()) return;
    if (*symbols_.rbegin () == '\r') symbols_.erase (symbols_.size () - 1);
}

void StringUtil::whiteSpace2UnderScore (char *symbols_) // puts single underscore between words
{
    assert (strlen (symbols_) != 0);

    stringstream sstr; 

    sstr.str (symbols_);
    sstr.clear ();

    sstr >> skipws;

    string word;

    sstr >> word;

    if (sstr.fail ()) 
    {
        symbols_ [0] = '\0';
        return;
    }

    string str = word;

    while (sstr >> word) 
    {
        str += string ("_") + word;
    }

    strcpy (symbols_, str.c_str ());
}

void StringUtil::whiteSpace2UnderScore (std::string &symbols_) // puts single underscore between words
{
    char *str = new char [symbols_.size () + 1];

    strcpy (str, symbols_.c_str ());
    whiteSpace2UnderScore (str);
    symbols_ = str;

    delete [] str; str = 0;
}

void StringUtil::split (
std::vector <std::string> *strVec_, // vector of pieces
const std::string &str_, // input string (can be part of strVec_)
const std::string &split_) // split points
{
     assert (strVec_);
     assert (split_.length () != 0);

     string str (str_);

     strVec_->clear ();
     string::size_type pos0 = 0;

     for (string::size_type pos = str.find (split_ [0]); pos != string::npos && pos != str.length (); )
     {
          if (str.substr (pos, min (split_.length (), str.length () - pos)) != split_) {

                pos++;
                continue;
          }

          strVec_->push_back (pos == pos0 ? string ("") : str.substr (pos0, pos - pos0));
          pos += split_.length ();
          pos0 = pos;
     }

     strVec_->push_back (str.length () == pos0 ? string ("") : str.substr (pos0, str.length () - pos0));
}

size_t StringUtil::splitCount (
const std::string &str_, // input string (can be part of strVec_)
const std::string &split_) // split points
{
    std::vector <std::string> strVec;

    split (&strVec, str_, split_);

    return strVec.size ();
}
