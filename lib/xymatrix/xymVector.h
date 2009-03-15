//!\author Udo Frese

/*!\file xymVector.h Contains the class \c XymVector representing a
  vector.
*/

#ifndef XYMVECTOR_H
#define XYMVECTOR_H

#include "xymVectorV.h"

//!Vector owning it's own memory. 
/*! To avoid confusion: This is a matrix for double vectors used in
    linear algebra not for general indexed storage like the STL
    vector.

    In contrast to \c VectorV the stride is always 1, so the access
    operator is quicker. One can reserve extra space, so the vector
    can be enlarged without copying it. (like \c XymMatrixC)
*/
class XymVector : public XymVectorV
{
public:
    //! Construct an empty vector
    XymVector();

    //! Construct a vector with \c n elements.
    /*! If \c initialised is \c true all elementsare initialised to 0
     */
    XymVector(int n, bool initialize=true);

    //! Create a vector with \c n elements and memory reserved for \c nReserve
    //! elements.
    /*! The elements are initialized to 0, if \c initialized==true.
        If \c autoGrow==true an automatic reallocation of the memory
        is performed if a call to \c append() or \c insert(...) exceed
        the reserved memory, otherwise this is an error. Note, that
        this process renders all views invalid.
    */
    XymVector(int n, int nReserve, bool initialize=true, bool autoGrow=false);

    //! Create a copy of 'v'
    XymVector(const XymVectorV& v);

    //! Create a copy of 'v'
    XymVector(const XymVector& v);


    //! Creates a vector from raw data
    /*! The vector has \c n elements, which are copied from \c d,
        where the \c i -th element is \c d[i*stride]. Note: \c stride
        refers to the source data at \c *d. The vector itself has unit
        stride.
    */
    XymVector(double* d, int n, int stride=1);

    //! Destructor
    ~XymVector();
    
    
    //! Returns a view to the subvector starting at entry \c i
    //! with size \c n entries.
    XymVectorV subVector (int i, int n) {
        XYMRANGECHECK (0<=i && i+n<=_n);
        return XymVectorV (_d+i*_stride, n, _stride);
    }

    //! The whole data is copied. If necessary new memory is allocated.
    XymVector& operator = (const XymVector& v);

    //! The whole data is copied. If necessary new memory is allocated.
    XymVector& operator = (const XymVectorV& v);
    
    
    //! stores the vector in the stream
    friend ostream& operator<<(ostream&, const XymVectorV& v);

    //! Reserves memory for being able to store an \c n vector.
    /*! If enough memory is reserved, nothing is done. Otherwise
        new memory is allocated and the matrix is copied.

        \warning This operation renders \em all views generated for
        this vector invalid.  \c autoGrow determines the behaviour,
        when an growing operation exceeds the reserved size. If \c
        autoGrow==true (default for not reserved vectors), new memory
        will be allocated and the vector copied (rendering all matrix
        views invalid).  If \c autoGrow==false an error will be
        thrown.
    */
    void reserve (int n, bool autoGrow=false);

    
    //! Append 'n' entries to the bottom of the vector.
    /*! If \c initialize==true the entries are filled with 0.
     */
    void append (int n, bool initialize=true);
    
    //! Removes the last \c n entries from the vector.
    void deleteLast (int n);

    //! Inserts \c n new entries before entry \c i.
    /*! If \c initialize==true the entries are filled with 0.
     */
    void insert (int i, int n, bool initialize=true);

    //! Deletes \c n entries starting at \c i.
    /*! Note, that since the following entries are moved all views
        refer to the same index, thus to a different entry.  \warning
        This cannot be called 'delete' because 'delete' is a reserved
        word.
    */
    void remove (int i, int n);

    //! Changes the size to \c n
    /*! If the vector is larger than \c n, entries are removed, otherwise
        appended. If \c fill==true, the appended entries are initialized
        with 0. */
    void resize (int n, bool initialize=true);
     

    //! Frees all memory of the vector
    /*! The vector will be in uninitialised state.
     */
    void clear();

    //! Creates a vector with \c n entries and \c reserve reserved entries
    /*! If \c initialize==true the entries are filled with 0.
     */
    void create(int n, bool init=true, int reserve=0);

    //! Set the autoGrow flag
    /*! If \c autoGrow==true an automatic reallocation of the memory
        is performed if a call to \c append() or \c insert(...) exceed
        the reserved memory, otherwise this is an error. Note, that this
        process renders all views invalid.
    */
    void setAutoGrow (bool autoGrow=true);

    //! Returns the meory usage (in bytes) of this vector
    /*! Does not include \c sizeof(*this). */
    int memoryUsage() const;

    //! Sets \c *this(i)=v(idx[i]).
    /*! If the vector is not initialised, it is created with the right size. 
     */
    void extractFrom (const XymVectorV& v, const vector<int>& idx);

    //! Copies data from \c v. If necessary the format is changed.
    void copyFrom (const XymVectorV& v);

    //! Transfers data from \c v. 
    /*! The data is moved from \c v to \c this by moving the pointer
        without copying. As a consequence \c v is invalid afterward
        containing no data. This function is for instane usefull to
        move computation results from a local variable to somewhere
        else, where the loca variable would be destroyed anyway.
     */
    void transferFrom (XymVector& v);
    
    //! loads the vector from the stream
    /*! It allocating memory and changes the vectors size if necessary
     */
    friend istream& operator>>(ifstream&, XymVectorV& v) throw (runtime_error);

protected:
    //! Pointer to the base of the allocated memory
    double *_memBase;

    //! Number of double allocated
    int _nReserved;

    //! Whether to automatically increase the memory size
    bool _autoGrow;
};

//! Load the vector from a stream
istream& operator>>(ifstream&, XymVectorV& v) throw (runtime_error);


#endif  /* XYMVECTOR_H */
