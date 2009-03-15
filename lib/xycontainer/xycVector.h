/*
Copyright (c) 2009, Universitaet Bremen
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Universitaet Bremen nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*! \author Udo Frese */
/*! \file xycVector.h 
    
    This file contains the template class \c XycVector. 
*/

#ifndef XYCVECTOR
#define XYCVECTOR

#include <stdlib.h>

#ifndef assert
#include <assert.h>
#endif

//! Forward declaration
template<class T> class XycVector;

//! Performs a lexicographical comparison between v1 and v2
/*! The result indicates, wheter \c v1<v2 (-1), \c v1==v2 (0) or
    \c v1>v2 (+1). Comparison is performed lexicographically.

    This means an empty vector is smaller than everything (except
    an empty vector). A vector starting with t1 is smaller than a
    vector starting with t2, iif t1<t2 and vice versa. A vector
    v1 starting with t and a vector v2 starting with t compare the
    same as the corresponding vectors with t removed.
*/
template<class T> int compare (const XycVector<T>& v1, const XycVector<T>& v2);

//! Returns, whether two vectors are equal
/*! Equal means: Same size and the i-th element is == in both
 */
template<class T> bool operator== (const XycVector<T>& v1, const XycVector<T>& v2);

//! Swap two vectors, simply by handing over memory
/*! No copying of data is performed. */
template<class T> int swap (XycVector<T>& v1, XycVector<T>& v2);


//! A container vector mostly compatible to stl::vector
/*! It 
    provides a vector container which is mostly compatible with
    STL::vector<>. However, it uses only a single template level,
    which has the following advantages:

    \li reasonable performance in debug mode
    \li member functions do not call other functions, so you can
        step throught them more easily in the debugger
    \li range checking
    \li new functions \c compact and \c resizeCompactlyWithUndefindedData
        for smaller memory usage when desired

    Function not supported by the STL vector are marked as NONSTD.
 */
template<class T> class XycVector {
protected:
   //! The 0th entry starts here
   T* _begin;

   //! One beyond the last entry
   T* _end;

   //! One beyond the end of the allocated memory
   T* _storageEnd;

 public:
   //! Entries in the container
   typedef T value_type;
   //! Pointer to entries 
   typedef T* pointer;
   //! Reference to entries 
   typedef T& reference;
   //! Const reference to entries
   typedef const T& const_reference;
   //! For compatibility with STL
   typedef unsigned int size_type;
   //! For compatibility with STL
   typedef int difference_type;
   //! Plain pointer used as an iterator
   /*! We use plain pointers as iterators, because this prevents
       the multiple level of templates encountered in STL.
   */
   typedef T* iterator;
   //! Plain pointer used as an const iterator   
   typedef const T* const_iterator;

   //! Empty container without memory
   XycVector() :_begin(NULL), _end(NULL), _storageEnd(NULL){}

   //! Container with n initialized T() entries
   XycVector (int n)
   {
      _begin = new T[n];
      _end = _storageEnd = _begin+n;
   }

   //! Container with n entries, initialized as t
   XycVector (int n, const T& t=T())
   {
      _begin = new T[n];
      _end = _storageEnd = _begin+n;
      for (T* p=_begin; p!=_end; p++) *p = t;
   }

   //! Copy constructor (deep copy)
   XycVector (const XycVector<T>& v2)
   {
      int n = v2.size();
      if (n>0) {
         _begin = new T[n];
         _storageEnd = _end  = _begin+n;
         for (T *p=v2._begin, *p2=_begin; p!=v2._end; p++, p2++) *p2 = *p;
      }
      else _begin = _end = _storageEnd = NULL;
   }

   //! Copy a vector from a range of iterators including \c from, not including \c to
   XycVector (const T* from, const T* to)
   {
      if (_begin!=_end) {
         _begin = new T[to-from];
         for (T *p=from, *p2=_begin; p!=to; p++, p2++) *p2 = *p;
         _end   = _storageEnd = _begin + (to-from);
      }
      else _begin = _end = _storageEnd = NULL;
   }

   //! Assignment operator (deep copy)
   XycVector<T>& operator = (const XycVector<T>& v2)
   {
      resizeWithUndefinedData (v2.size());
      for (T *p=v2._begin, *p2=_begin; p!=v2._end; p++, p2++) *p2 = *p;      
      return *this;
   }

   //! Destructor, frees all elements
   ~XycVector() {if (_begin!=NULL) delete[] _begin;}

   //! Whether there is no entry in the vector
   bool empty() const {return _begin==_end;}

   //! Number of entries
   int size() const {return _end-_begin;}

   //! Capacity reserved
   /*! The vector can grow up to size()==capacity()
       without the need for allocating new memory and
       copying data. */
   int capacity() const {return _storageEnd-_begin;}

   //! NONSTD: return, whether \c idx is a valid index for (*this)[]
   bool idx (int idx) const {return 0<=idx && idx<(int) (_end-_begin);}          

   //! (*this)[i] returns the i-th entry of the vector
   /*! Validity of the index is asserted. */
   T& operator[] (int idx) {
      assert (0<=idx && idx<size());
      return _begin[idx];
   }

   //! (*this)[i] returns the i-th entry of the vector
   /*! Validity of the index is asserted. */
   const T& operator[] (int idx) const {
      assert (0<=idx && idx<size());
      return _begin[idx];
   }

   //! Iterator/pointer to first entry
   T* begin() {return _begin;}

   //! Iterator/pointer to first entry
   const T* begin() const {return _begin;}

   //! Iterator/pointer to one beyond the last entry
   T* end() {return _end;}

   //! Iterator/pointer to one beyond the last entry
   const T* end() const {return _end;}

   //! The last entry
   T& back() {return *(_end-1);}

   //! The last entry
   const T& back() const {return *(_end-1);}

   //! The first entry
   T& front() {return *_begin;}

   //! The first entry
   const T& front() const {return *_begin;}

   //! NONSTD: Erases all entry after \c end (including) 
   /*! Slightly faster than the general \c erase function.*/
   void eraseAfter (T* end)
   {
     _end = end;
   }   


   //! NONSTD: Resizes but does not initinialize
   void resizeWithUndefinedData (int n)
   {
     T* newEnd = _begin+n;     
     if (_storageEnd<newEnd) {
       if (_begin!=NULL) delete[] _begin;
       _begin = new T[n];
       _end   = _storageEnd = _begin + n;       
     }     
     else _end = newEnd;
   }

   //! NONSTD: Resizes and allocates new memory except the current memory exactly fits
   /*! This function allows to avoid having unused memory in a vector. */
   void resizeCompactlyWithUndefindedData (int n)
   {
     T* newEnd = _begin+n;     
     if (_storageEnd!=newEnd) {
       if (_begin!=NULL) delete[] _begin;
       _begin = new T[n];
       _end   = _storageEnd = _begin + n;       
     }     
     else _end = newEnd;     
   }   

   //! Shrink or expand to size \c n initializing new entries with \c t
   /*! With resize the current memory is used, unless \c n is larger than
       capacity. This means, that \c resize never frees unused memory.
       If this is desired, call \c resizeCompactlyWithUndefindedData
   */
   void resize (int n, const T& t = T())
   {
      reserve (n);
      T* nEnd = _begin + n;
      for (T* p=_end; p<nEnd; p++) *p = t;
      _end = nEnd;
   }

   //! Like reserve but does not deallocate the old ptr and returns it instead
   T* internal_reserve (int n)
     {
       if (_storageEnd<_begin+n) { // extend and reallocate
         if (n<8) n = 8;        
         if (n<2*(_storageEnd-_begin)) n = 2*(_storageEnd-_begin);        
         T* newBegin = new T[n];
         for (T* p=_begin, *p2=newBegin; p!=_end; p++, p2++) *p2 = *p;
         T* oldBegin = _begin;         
         _end = newBegin + (_end-_begin);
         _begin = newBegin;
         _storageEnd = _begin+n;
         return oldBegin;         
       }
       else return NULL;       
     }
   
   //! Allocate memory for at least \c n entries
   /*! If \c n is smaller than the current \c capacity()
       nothing is done.
   */
   void reserve (int n)
   {
     T* old = internal_reserve (n);
     if (old!=NULL) delete[] old;     
   }

   //! Make the vector empty
   /*! Does not free any memory */
   void clear () {_end = _begin;}
   
   //! NONSTD: Copy the vector content to new memory of exacttly the right size
   void compact () 
     {
       if (_storageEnd>_end) {
         int n = _end-_begin; 
         T* newBegin;
         if (n>0) {         
           newBegin = new T[n];
           for (T* p=_begin, *p2=newBegin; p!=_end; p++, p2++) *p2 = *p;
         }
         else newBegin = NULL;
         if (_begin!=NULL) delete[] _begin;
         _end = newBegin + n;
         _begin = newBegin;
         _storageEnd = _begin+n;
       }
     }   

   //! Memory used by the vector
   int memory () const
   {
     return capacity()*sizeof(T);
   }   


   //! Append \c t to the end of the vector
   /*! If adding an entry exceeds the vectors capacity, new memory
       is allocated and everything copied.
   */
   void push_back (const T& t)
   {
     if (_storageEnd>_end) {       
       *_end = t;
       _end++;
     }
     else {
       T* old = internal_reserve (_storageEnd-_begin+1);
       *_end = t;
       _end++;
       if (old!=NULL) delete[] old;
     }     
   }

   //! Remove the last entry
   void pop_back () {
      assert (_end!=_begin);
      _end--;
   }

   //! Erase all entries between \c from and \c to (not including)
   void erase (T* from, T* to)
   {
      assert (_begin<=from && from<=to && to<=_end);
      for (T *p=to, *p2=from; p!=_end; p++,p2++) *p2 = *p;
      _end -= (to-from);
   }

   //! Swap \c *this and \c v2 without copying entries
   void swap (XycVector<T>& v2)
   {
   /*
      int n1 = size(), n2 = v2.size();
      if (capacity()<n2) reserve (n2);
      if (v2.capacity()<n1) v2.reserve (n1);
      T* pEnd;
      if (n1>n2) pEnd = _begin + n1;
      else pEnd = _begin + n2;      
      for (T *p=_begin, *p2=v2._begin; p!=pEnd; p++, p2++) swap (*p, *p2);
      _end = _begin + n2;
      v2._end = v2._begin + n1;      
   */
     T* buf = _begin; _begin = v2._begin; v2._begin = buf;   
     buf = _end; _end = v2._end, v2._end = buf;
     buf = _storageEnd; _storageEnd = v2._storageEnd; v2._storageEnd = buf;     
   }

   //! Insert an entry t before \c pos
   /*! Returns a new iterator to the entry
       that was \c pos before.  

       If \c capacity() is too low, new memory is
       allocated and the whole vector copied.
   */
   T* insert (T* pos, const T& t)
   {
      assert (_begin<=pos && pos<=_end);
      if (_storageEnd>_end) {        
        for (T *p = _end, *p2=_end+1; p2!=pos; p--, p2--) *p2 = *p;
        *pos = t;
      }
      else {
        T* old = internal_reserve (size()+1);
        pos = _begin + pos - old;        
        for (T *p = _end, *p2=_end+1; p2!=pos; p--, p2--) *p2 = *p;
        *pos = t;
        if (old!=NULL) delete[] old;        
      }
      
      return pos;        
   }

   //! insert entries \c *from to  *to (exclusive) before \c *pos
   /*! Returns a new iterator to the entry
       that was \c pos before.  

       If \c capacity() is too low, new memory is
       allocated and the whole vector copied.
   */
   void insert (T* pos, const T* from, const T* to)
   {
      assert (_begin<=pos && pos<=_end && from<=to);
      reserve (size()+(to-from));
      for (T *p = _end, *p2=_end+(to-from); p2!=pos; p--, p2--) *p2 = *p;
      for (T *p = from, *p2 = pos; p!=to; p++, p2++) *p2 = *p;
   }


   //! Insert \c n copies of \c t before \c pos
   /*! Returns a new iterator to the entry
       that was \c pos before.  

       If \c capacity() is too low, new memory is
       allocated and the whole vector copied.
   */
   T* insert (T* pos, int n, const T& t)
   {
      assert (_begin<=pos && pos<=_end);
      if (_storageEnd>=_end+n) {
        T *p, *p2, *to=_end+n;      
        for (p = _end, p2=to; p2!=pos; p--, p2--) *p2 = *p;
        for (p2 = pos; p2!=to; p2++) *p2 = t;
      }
      else {
        T* old = internal_reserve (size()+n);
        pos = _begin + pos - old;                
        T *p, *p2, *to=_end+n;      
        for (p = _end, p2=to; p2!=pos; p--, p2--) *p2 = *p;
        for (p2 = pos; p2!=to; p2++) *p2 = t;
        if (old!=NULL) delete[] old;        
      }      
      return pos;      
   }

   //! operator== can directly access internas
   template<class TT> friend bool operator== (const XycVector<TT>& v1, const XycVector<TT>& v2);   
   //! the comparison function can directly access internas
   template<class TT> friend int compare (const XycVector<TT>& v1, const XycVector<TT>& v2);   
   //! the swap function can directly access internas
   template<class TT> friend void swap (XycVector<TT>& v1, XycVector<TT>& v2);   
};


template<class T> bool operator== (const XycVector<T>& v1, const XycVector<T>& v2)
{
   if (v1.size()!=v2.size()) return false;
   for (T *p1=v1._begin, *p2=v2._begin; p1!=v1._end; p1++, p2++)
      if (!(*p1==*p2)) return false;
   return true;
}


template<class T> int compare (const XycVector<T>& v1, const XycVector<T>& v2)
{
   int n = min(v1.size(), v2.size());
   T* pEnd = v1._begin+n;
   for (T* p1=v1._begin, *p2=v2._begin; p1!=pEnd; p1++, p2++)
      if (*p1<*p2) return -1;
      else if (*p2<*p1) return +1;
   if (v1.size()<v2.size()) return -1;
   else if (v1.size()>v2.size()) return +1;
   else return 0;
}


//! Lexicographical comparison. See \c compare
template<class T> bool operator!= (const XycVector<T>& v1, const XycVector<T>& v2)
{
  return !(v1==v2);
}


//! Lexicographical comparison. See \c compare
template<class T> bool operator< (const XycVector<T>& v1, const XycVector<T>& v2)
{
  return compare (v1, v2)<0;
}

//! Lexicographical comparison. See \c compare
template<class T> bool operator<= (const XycVector<T>& v1, const XycVector<T>& v2)
{
  return compare (v1, v2)<=0;
}

//! Lexicographical comparison. See \c compare
template<class T> bool operator> (const XycVector<T>& v1, const XycVector<T>& v2)
{
  return compare (v1, v2)>0;
}

//! Lexicographical comparison. See \c compare
template<class T> bool operator>= (const XycVector<T>& v1, const XycVector<T>& v2)
{
  return compare (v1, v2)>=0;
}

template<class T> void swap (XycVector<T>& v1, XycVector<T>& v2)
{
   T* buf = v1._begin; v1._begin = v2._begin; v2._begin = buf;   
   buf = v1._end; v1._end = v2._end, v2._end = buf;
   buf = v1._storageEnd; v1._storageEnd = v2._storageEnd; v2._storageEnd = buf;   
}

#endif
