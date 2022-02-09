//------------------------------------------------------------------------------
// LAGraph_Matrix_print:  pretty-print a matrix
//------------------------------------------------------------------------------

// LAGraph, (c) 2021 by The LAGraph Contributors, All Rights Reserved.
// SPDX-License-Identifier: BSD-2-Clause
// Contributed by Tim Davis, Texas A&M University.

//------------------------------------------------------------------------------

// LAGraph_Matrix_print:  pretty-print a matrix.  The type is either derived
// from GxB_Matrix_type (if available) or assumed to be GrB_FP64 otherwise,
// or passed in as a parameter.
// Contributed by Tim Davis, Texas A&M.

#include "LG_internal.h"

#undef  LG_FREE_WORK
#define LG_FREE_WORK                \
{                                   \
    LAGraph_Free ((void **) &I) ;   \
    LAGraph_Free ((void **) &J) ;   \
    LAGraph_Free ((void **) &X) ;   \
}

#undef  LG_FREE_ALL
#define LG_FREE_ALL LG_FREE_WORK

//------------------------------------------------------------------------------
// LG_Matrix_print_TYPE: print with the specified type
//------------------------------------------------------------------------------

#define LG_MATRIX_PRINT(suffix,ctype,gtype,fmt1,fmt2)                       \
int LG_Matrix_print_ ## suffix                                              \
(                                                                           \
    GrB_Matrix A, int pr, FILE *f, char *msg                                \
)                                                                           \
{                                                                           \
    LG_CLEAR_MSG ;                                                          \
    ctype *X = NULL ;                                                       \
    GrB_Index *I = NULL, *J = NULL ;                                        \
    LG_ASSERT (A != NULL && f != NULL, GrB_NULL_POINTER) ;                  \
    if (pr < 0) return (GrB_SUCCESS) ;                                      \
    /* get basic properties */                                              \
    GrB_Index nrows, ncols, nvals ;                                         \
    GrB_TRY (GrB_Matrix_nrows (&nrows, A)) ;                                \
    GrB_TRY (GrB_Matrix_ncols (&ncols, A)) ;                                \
    GrB_TRY (GrB_Matrix_nvals (&nvals, A)) ;                                \
    /* print header line */                                                 \
    FPRINTF (f, "%s matrix: %" PRIu64 "-by-%" PRIu64 " entries: %" PRIu64   \
        "\n", LG_XSTR (gtype), nrows, ncols, nvals) ;                       \
    if (pr <= 1) return (GrB_SUCCESS) ;                                     \
    /* extract tuples */                                                    \
    I = LAGraph_Malloc (nvals, sizeof (GrB_Index)) ;                        \
    J = LAGraph_Malloc (nvals, sizeof (GrB_Index)) ;                        \
    X = LAGraph_Malloc (nvals, sizeof (ctype)) ;                            \
    LG_ASSERT (I != NULL && J != NULL && X != NULL, GrB_OUT_OF_MEMORY) ;    \
    GrB_Info info = GrB_Matrix_extractTuples (I, J, X, &nvals, A) ;         \
    LG_ASSERT_MSG (info != GrB_DOMAIN_MISMATCH,                             \
        GrB_NOT_IMPLEMENTED, "user-defined types not supported") ; /* RETVAL */\
    GrB_TRY (info) ;                                                        \
    /* determine the format */                                              \
    char *format = (pr <= 3) ? fmt1 : fmt2 ;                                \
    bool summary = (pr == 2 || pr == 4) && (nvals > 30) ;                   \
    for (int64_t k = 0 ; k < nvals ; k++)                                   \
    {                                                                       \
        /* print the kth tuple */                                           \
        GrB_Index i = I [k] ;                                               \
        GrB_Index j = J [k] ;                                               \
        ctype     x = X [k] ;                                               \
        FPRINTF (f, "    (%" PRIu64 ", %" PRIu64 ")   ", i, j) ;            \
        FPRINTF (f, format, x) ;                                            \
        FPRINTF (f, "\n") ;                                                 \
        if (summary && k >= 29)                                             \
        {                                                                   \
            /* quit early if a only a summary is requested */               \
            FPRINTF (f, "    ...\n") ;                                      \
            break ;                                                         \
        }                                                                   \
    }                                                                       \
    LG_FREE_WORK ;                                                          \
    return (GrB_SUCCESS) ;                                                  \
}

LG_MATRIX_PRINT (BOOL  , bool    , GrB_BOOL  , "%d"  , "%d"    ) ;
LG_MATRIX_PRINT (INT8  , int8_t  , GrB_INT8  , "%d"  , "%d"    ) ;
LG_MATRIX_PRINT (INT16 , int16_t , GrB_INT16 , "%d"  , "%d"    ) ;
LG_MATRIX_PRINT (INT32 , int32_t , GrB_INT32 , "%" PRId32, "%" PRId32  ) ;
LG_MATRIX_PRINT (INT64 , int64_t , GrB_INT64 , "%" PRId64, "%" PRId64  ) ;
LG_MATRIX_PRINT (UINT8 , uint8_t , GrB_UINT8 , "%d"  , "%d"    ) ;
LG_MATRIX_PRINT (UINT16, uint16_t, GrB_UINT16, "%d"  , "%d"    ) ;
LG_MATRIX_PRINT (UINT32, uint32_t, GrB_UINT32, "%" PRIu32, "%" PRIu32  ) ;
LG_MATRIX_PRINT (UINT64, uint64_t, GrB_UINT64, "%" PRIu64, "%" PRIu64  ) ;
LG_MATRIX_PRINT (FP32  , float   , GrB_FP32  , "%g"  , "%0.7g" ) ;
LG_MATRIX_PRINT (FP64  , double  , GrB_FP64  , "%g"  , "%0.15g") ;
#if 0
LG_MATRIX_PRINT (FC32  , GxB_FC32_t, GxB_FC32, ...) ;
LG_MATRIX_PRINT (FC64  , GxB_FC64_t, GxB_FC64, ...) ;
#endif

#undef  LG_FREE_WORK
#define LG_FREE_WORK ;
#undef  LG_FREE_ALL
#define LG_FREE_ALL ;

//------------------------------------------------------------------------------
// LAGraph_Matrix_print: automatically determine the type
//------------------------------------------------------------------------------

int LAGraph_Matrix_print    // TODO rename LAGraph_Matrix_Print
(
    // input:
    GrB_Matrix A,       // matrix to pretty-print to the file
    // TODO: use an enum for pr
    int pr,             // print level: -1 nothing, 0: one line, 1: terse,
                        //      2: summary, 3: all,
                        //      4: as 2 but with %0.15g for float/double
                        //      5: as 3 but with %0.15g for float/double
    FILE *f,            // file to write it to, must be already open; use
                        // stdout or stderr to print to those locations.
    char *msg
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    LG_CLEAR_MSG ;
    LG_ASSERT (A != NULL && f != NULL, GrB_NULL_POINTER) ;

    //--------------------------------------------------------------------------
    // determine the type
    //--------------------------------------------------------------------------

    GrB_Type type ;
    #if LG_SUITESPARSE
        // SuiteSparse:GraphBLAS: query the type and print accordingly
        GrB_TRY (GxB_Matrix_type (&type, A)) ;
    #else
        // no way to determine the type with pure GrB*; print as if FP64
        type = GrB_FP64 ;
    #endif

    //--------------------------------------------------------------------------
    // print the matrix
    //--------------------------------------------------------------------------

    if (type == GrB_BOOL)
    {
        return (LG_Matrix_print_BOOL (A, pr, f, msg)) ;
    }
    else if (type == GrB_INT8) 
    {
        return (LG_Matrix_print_INT8 (A, pr, f, msg)) ;
    }
    else if (type == GrB_INT16) 
    {
        return (LG_Matrix_print_INT16 (A, pr, f, msg)) ;
    }
    else if (type == GrB_INT32) 
    {
        return (LG_Matrix_print_INT32 (A, pr, f, msg)) ;
    }
    else if (type == GrB_INT64) 
    {
        return (LG_Matrix_print_INT64 (A, pr, f, msg)) ;
    }
    else if (type == GrB_UINT8) 
    {
        return (LG_Matrix_print_UINT8 (A, pr, f, msg)) ;
    }
    else if (type == GrB_UINT16) 
    {
        return (LG_Matrix_print_UINT16 (A, pr, f, msg)) ;
    }
    else if (type == GrB_UINT32) 
    {
        return (LG_Matrix_print_UINT32 (A, pr, f, msg)) ;
    }
    else if (type == GrB_UINT64) 
    {
        return (LG_Matrix_print_UINT64 (A, pr, f, msg)) ;
    }
    else if (type == GrB_FP32) 
    {
        return (LG_Matrix_print_FP32 (A, pr, f, msg)) ;
    }
    else if (type == GrB_FP64) 
    {
        return (LG_Matrix_print_FP64 (A, pr, f, msg)) ;
    }
    #if 0
    else if (type == GxB_FC32)
    {
        return (LG_Matrix_print_FC32 (A, pr, f, msg)) ;
    }
    else if (type == GxB_FC32)
    {
        return (LG_Matrix_print_FC64 (A, pr, f, msg)) ;
    }
    #endif
    else
    {
        LG_ASSERT_MSG (false,
            GrB_NOT_IMPLEMENTED, "user-defined types not supported") ; // RETVAL
        return (GrB_SUCCESS) ;
    }
}

