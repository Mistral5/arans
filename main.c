#define ARANS_STATIC

//includes
#include <stdio.h>
#include <string.h>

#include "platform.h"

#include "arans_8.h"
//#include "arans_4x4.h"
//#include "arans_3x5.h"
//#include "arans_2x6.h"
//#include "arans_2x3x3.h"
//#include "arans_2x2x4.h"
//#include "arans_2x2x2x2.h"
//#include "arans_3x5_clear.h"
//#include "arans_SIMD.h"

//main function
int main(int argc, char *argv[]) {
    if (argc != 4) {
        printf("Missing argument!\n");
        return 0;
    }

    // 0 - mode arg is wrong, 1 - encoding, 2 - decoding
    uint8_t mode = 0;

    if (strcmp(argv[1], "enc") == 0)
        mode = 1;

    if (strcmp(argv[1], "dec") == 0)
        mode = 2;

    if (!mode) {
        printf("Check the selected mode!\n");
        return 0;
    }

    FILE *in_file = fopen(argv[2], "rb");
    if (!in_file) {
        printf("File not found!\n");
        return 0;
    }

    FILE *out_file = fopen(argv[3], "wb");
    if (!out_file) {
        fclose(in_file);
        printf("Unable to write!\n");
        return 0;
    }

    fseek(in_file, 0, SEEK_END);
    size_t in_size = ftell(in_file);
    fseek(in_file, 0, SEEK_SET);
    size_t out_size = 0;

    unsigned char *in = (unsigned char *) malloc(in_size);
    unsigned char *out = NULL;

    if (!in) {
        fclose(in_file);
        fclose(out_file);
        printf("Allocate failed!\n");
        return 0;
    }

    size_t res = fread(in, in_size, 1, in_file);

    double start_execution_time;
    uint64_t start_clocks;
    uint64_t clocks;
    double execution_time;

    if (mode == 1) {
        out_size = in_size + 128;
        out = (unsigned char *) malloc(out_size);

        struct Arans arans;
        aransInit(&arans);

        start_execution_time = timer();
        start_clocks = __rdtsc();

        //do encoding
        out_size = aransEncode(&arans, out, out_size, in, in_size);

        clocks = __rdtsc() - start_clocks;
        execution_time = timer() - start_execution_time;
    }

    if (mode == 2) {
        out_size = aransGetOutFileSize(in);
        out = (unsigned char *) malloc(out_size);

        struct Arans arans;
        aransInit(&arans);

        start_execution_time = timer();
        start_clocks = __rdtsc();

        //do decoding
        aransDecode(&arans, out, out_size, in, in_size);

        clocks = __rdtsc() - start_clocks;
        execution_time = timer() - start_execution_time;
    }

    //write to output file
    fwrite(out, out_size, 1, out_file);

    //print decompression ratio
    printf("%zu to %zu (%.1f%%)\n", in_size, out_size, 100.0 * (double) out_size / (double) in_size);

    printf("%" PRIu64" clocks, %.1f clocks/symbol (%5.1fMiB/s)\n", clocks,
           (double) clocks / (double) in_size,
           (double) in_size / (execution_time * 1048576.0));

    if (mode == 1) {
        FILE *fstat = fopen("stat.txt", "a");
        fprintf(fstat, "%zu\t%zu\n", in_size, out_size);
        fclose(fstat);
    }

    FILE *f_speed_stat = fopen("speed_stat.txt", "a");
    fprintf(f_speed_stat, "%.1f\n", (double) clocks / (double) in_size);
    fclose(f_speed_stat);

    //free buffers
    free(in);
    free(out);

    //close files
    fclose(in_file);
    fclose(out_file);

    return 0;
}