
/*
 * This test compares the performance of multiply add "as used in convolution reverbs"
 * with 2 different implementations : contiguous memory / non contiguous memory
 * (variable parameters are block size / number of blocks)
 *
 * Tests for sizes < and > page size
 */
TEST(PerfMultiplyAdd, measure) {
    imajuscule::profiling::measure_madd::run_multiplyadd_test();
}
