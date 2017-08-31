/*
c_multiply.c

simple C function that alters data passed in via a pointer

    used to see how we can do this with Cython/numpy

*/

#include "huffman.h"
#include "DrsCalib.h"
#include "factfits.h"

extern "C"{
    void c_multiply (double* array, double multiplier, int m, int n) {

        int i, j ;
        int index = 0 ;

        for (i = 0; i < m; i++) {
            for (j = 0; j < n; j++) {
                array[index] = array[index]  * multiplier ;
                index ++ ;
                }
            }
        return ;
    }

    void Decode_dom(const uint8_t *bufin, size_t bufinlen, int16_t *bufout, size_t bufoutlen) {
        int64_t i = 0;
        const Huffman::Decoder decoder(bufin, i);
        uint16_t * mybuf = (uint16_t *)bufout;
        uint16_t * mybuf_end = (uint16_t *)(bufout+bufoutlen);
        decoder.Decode(bufin+i, bufin+bufinlen, mybuf, mybuf_end);
        return ;
    }

    void remove_spikes_4_dom(
        float* calib_data,
        size_t number_of_pixel,
        uint32_t roi
    )
    {
        for (size_t ch=0; ch<number_of_pixel; ch++)
        {
            DrsCalibrate::RemoveSpikes4(
                calib_data+ch*roi,
                roi);
        }
    }
}
