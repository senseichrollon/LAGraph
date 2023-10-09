#define LG_FREE_WORK                \
{                                   \
    GrB_free (&h_old) ;                 \
    GrB_free (&a_old) ;                 \
}

#define LG_FREE_ALL                 \
{                                   \
    LG_FREE_WORK ; \
    GrB_free(&h)                  \
    GrB_free(&a)  \
}

#include "LG_internal.h"
#include "LAGraphX.h"
#include <math.h>

int LaGr_HITS(
    GrB_Vector * hubs,
    GrB_Vector* authorities,
    int * iters,
    const LAGraph_Graph G,
    float tol,
    int itermax,
    char *msg 
) {


    LG_ASSERT (hubs != NULL, GrB_NULL_POINTER) ;
    LG_ASSERT (authorities != NULL, GrB_NULL_POINTER) ;
    LG_TRY (LAGraph_CheckGraph (G, msg)) ;
    GrB_Matrix AT ;
    if (G->kind == LAGraph_ADJACENCY_UNDIRECTED ||
        G->is_symmetric_structure == LAGraph_TRUE)
    {
        // A and A' have the same structure
        AT = G->A ;
    }
    else
    {
        // A and A' differ
        AT = G->AT ;
        LG_ASSERT_MSG (AT != NULL,
            LAGRAPH_NOT_CACHED, "G->AT is required") ;
    }

    //initializations
    GrB_Vector h = NULL, a = NULL, h_old = NULL, a_old=NULL, h_square=NULL, a_square=NULL, t=NULL;
    GrB_Index n;
    (*hubs) = NULL;
    (*authorities) = NULL;
    GRB_TRY(GrB_Matrix_nrows(&n, AT))
    float rdiff = 1;


    GRB_TRY (GrB_Vector_new (&h_old, GrB_FP32, n)) ;
    GRB_TRY (GrB_Vector_new (&a_old, GrB_FP32, n)) ;

    GRB_TRY (GrB_Vector_new (&h_square, GrB_FP32, n)) ;
    GRB_TRY (GrB_Vector_new (&a_square, GrB_FP32, n)) ;

    GRB_TRY (GrB_Vector_new (&h, GrB_FP32, n));
    GRB_TRY (GrB_Vector_new (&a, GrB_FP32, n)) ;

   // GRB_TRY (GrB_Vector_new (&t, GrB_FP32, n)) ;



    for((*iters) = 0; (*iters) < itermax && rdiff > tol; (*iters)++) {
        //save old values of h and a
        GrB_Vector temp = h_old ; h_old = h ; h = temp ;
        temp = a_old ; a_old = a ; a = temp ;

        //a = AT . h
        GRB_TRY(GrB_mxv(&a, NULL,NULL, LAGraph_plus_second_fp32, AT, h_old, NULL));
        //h = A . a
        GRB_TRY(GrB_mxv(&h, NULL,NULL, LAGraph_plus_second_fp32, G->A, a_old, NULL));
       
        //scale a

        //a_square = a*a
        GRB_TRY(GrB_eWiseMult(&a_square,NULL, NULL, GrB_TIMES_FP32, a, a, NULL));
        float norm;
        // norm = sum(a_square)
        GRB_TRY (GrB_reduce(&norm, NULL, GrB_PLUS_MONOID_FP32, a_square NULL)) ;
        norm = sqrt(norm);
        
        //a = a/norm
        GRB_TRY(GrB_Vector_assign(a,NULL, GrB_DIV_FP32, &norm, NULL, NULL, NULL));

        //scale h

        // h_square = h*h
        GRB_TRY(GrB_eWiseMult(&h_square,NULL, NULL, GrB_TIMES_FP32, h, h, NULL));
        //norm = sum(h_square)
        GRB_TRY (GrB_reduce(&norm, NULL, GrB_PLUS_MONOID_FP32, h_square NULL)) ;
        norm = sqrt(norm);

        // h = h/norm
        GRB_TRY(GrB_assign(h,NULL, GrB_DIV_FP32, &norm, NULL, NULL, NULL));


        //deal with tolerance

        //a_old -= a
        GRB_TRY (GrB_assign (a_old, NULL, GrB_MINUS_FP32, a, GrB_ALL, n, NULL));
        //a_old = abs(a_old)
        GRB_TRY (GrB_apply (a_old, NULL, NULL, GrB_ABS_FP32, a_old, NULL));
        //rdiff = sum(a_old)
        GRB_TRY (GrB_reduce (&rdiff, NULL, GrB_PLUS_MONOID_FP32, t, NULL));

        //h_old -= h
        GRB_TRY (GrB_assign (h_old, NULL, GrB_MINUS_FP32, h, GrB_ALL, n, NULL));
        // h_old = abs(h_old)
        GRB_TRY (GrB_apply (h_old, NULL, NULL, GrB_ABS_FP32, h_old, NULL));
        //rdiff += sum(h_old)
        GRB_TRY (GrB_reduce (&rdiff, GrB_PLUS_FP32, GrB_PLUS_MONOID_FP32, t, NULL)) ;


        //rdiff = rdiff/2
        rdiff /= 2;
    }

    (*hubs) = h;
    (*authorities) = a;
    LG_FREE_WORK;
    return (GrB_SUCCESS);
}