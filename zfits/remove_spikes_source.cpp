#include "huffman.h"
#include "DrsCalib.h"
#include "factfits.h"

extern "C"{
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
