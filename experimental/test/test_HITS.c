//------------------------------------------------------------------------------
// test_HITS.c: Test suite for HITS algorithm using GraphBLAS and LAGraph
//------------------------------------------------------------------------------

#include <stdio.h>
#include <acutest.h>
#include "LAGraph_test.h"

#define LEN 512
char msg[LAGRAPH_MSG_LEN];
char filename[LEN + 1];
LAGraph_Graph G = NULL;

// Utility function to compare vectors with expected results
float difference(GrB_Vector vector, double *expected_result, GrB_Index n) {
    GrB_Vector diff = NULL;
    float err = 0.0;
    OK(GrB_Vector_new(&diff, GrB_FP32, n));
    for (GrB_Index i = 0; i < n; i++) {
        OK(GrB_Vector_setElement_FP64(diff, expected_result[i], i));
    }
    OK(GrB_eWiseAdd(diff, NULL, NULL, GrB_MINUS_FP32, vector, diff, NULL));
    OK(GrB_apply(diff, NULL, NULL, GrB_ABS_FP32, diff, NULL));
    OK(GrB_reduce(&err, NULL, GrB_MAX_MONOID_FP32, diff, NULL));
    OK(GrB_free(&diff));
    return err;
}

// Test function for a specific graph
void test_HITS_on_graph(const char *graph_file, double *expected_hubs, double *expected_authorities, GrB_Index n) {
    GrB_Vector hubs = NULL, authorities = NULL;
    GrB_Matrix A = NULL;
    int iters = 0;
    float tol = 1e-4;
    int itermax = 100;

    // Load the graph
    snprintf(filename, LEN, LG_DATA_DIR "%s", graph_file);
    FILE *f = fopen(filename, "r");
    TEST_CHECK(f != NULL);
    OK(LAGraph_MMRead(&A, f, msg));
    OK(fclose(f));
    OK (LAGraph_New (&G, &A, LAGraph_ADJACENCY_DIRECTED, msg)) ;
    TEST_CHECK (A == NULL) ;    // A has been moved into G->A
    OK (LAGraph_Cached_OutDegree (G, msg)) ;
    OK (LAGraph_Cached_InDegree (G, msg)) ;
    OK (LAGraph_Cached_AT(G, msg)) ;
    // Run HITS algorithm
    OK(LAGr_HITS(&hubs, &authorities, &iters, G, tol, itermax, msg));

    // Compare results with expected values
    float err_hubs = difference(hubs, expected_hubs, n);
    float err_auth = difference(authorities, expected_authorities, n);
    TEST_CHECK(err_hubs < 1e-4 && err_auth < 1e-4);

    // Clean up
    OK(GrB_free(&hubs));
    OK(GrB_free(&authorities));
    OK(LAGraph_Delete(&G, msg));
}

// Test cases for different graphs
void test_HITS(void) {
    LAGraph_Init(msg);

    // Test case for Graph 1
    // Define expected hubs and authorities values for the first graph
    double structure_hubs[] ={
        0.13865498151112815,
        0.13865498151112823,
        -7.52772710372611e-17,
        0.19608775534362824,
        -7.52772710372611e-17,
        0.15423823730232433,
        0.3723640443317912
    };

    double structure_authorities[] = {
    0.0884024498169609,
    0.06250997173907646,
    0.3258111125561342,
    0.230383247074376,
    0.23038324707437605,
    -1.5901485815586466e-16,
    0.06250997173907653};
    test_HITS_on_graph("structure.mtx", structure_hubs, structure_authorities, sizeof(structure_hubs)/sizeof(structure_hubs[0]));

    // Additional test cases can be added here
    double karate_authorities[] = {
        0.07141272880825196,
        0.05342723123552997,
        0.06371906455637479,
        0.042422737124709016,
        0.01526095970620749,
        0.01596691350305964,
        0.015966913503059642,
        0.03434316721905367,
        0.0456819251197503,
        0.020625667749388638,
        0.015260959706207484,
        0.010617891511071214,
        0.01692545079230685,
        0.045494864068056355,
        0.020370345825614276,
        0.02037034582561427,
        0.004748031847301578,
        0.018561637037432098,
        0.020370345825614273,
        0.02971333388643479,
        0.02037034582561427,
        0.018561637037432088,
        0.02037034582561428,
        0.030156497509356398,
        0.011460952230971697,
        0.011893664396281381,
        0.015182734330338388,
        0.026813494117104757,
        0.0263315057779539,
        0.027111539628217676,
        0.035106237976714395,
        0.038375741862956045,
        0.06200184647383098,
        0.0750029421565755
    };
    double karate_hubs[] = {
            0.07141272880825197,
            0.053427231235529976,
            0.06371906455637479,
            0.04242273712470899,
            0.015260959706207477,
            0.01596691350305964,
            0.015966913503059642,
            0.034343167219053644,
            0.045681925119750305,
            0.020625667749388638,
            0.01526095970620748,
            0.010617891511071217,
            0.01692545079230686,
            0.045494864068056355,
            0.020370345825614276,
            0.020370345825614276,
            0.004748031847301572,
            0.018561637037432077,
            0.020370345825614276,
            0.02971333388643479,
            0.020370345825614276,
            0.018561637037432077,
            0.020370345825614276,
            0.030156497509356388,
            0.011460952230971698,
            0.011893664396281383,
            0.015182734330338387,
            0.026813494117104767,
            0.02633150577795389,
            0.02711153962821767,
            0.03510623797671438,
            0.038375741862956045,
            0.062001846473830974,
            0.07500294215657549
        };
    test_HITS_on_graph("karate.mtx", karate_hubs, karate_authorities, sizeof(karate_hubs)/sizeof(karate_hubs[0]));
    LAGraph_Finalize(msg);
}

// List of tests to run
TEST_LIST = {
    {"test_HITS", test_HITS},
    {NULL, NULL}
};

