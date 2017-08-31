/*
 * factfits.h
 *
 *  Created on: May 26, 2013
 *      Author: lyard
 */

#ifndef MARS_FACTFITS
#define MARS_FACTFITS

#include "zfits.h"

class factfits : public zfits
{
public:
    // Default constructor
    factfits(const std::string& fname, const std::string& tableName="", bool force=false) :
        zfits(fname, tableName, force),
        fOffsetCalibration(0),
        fOffsetStartCellData(0),
        fOffsetData(0),
        fNumRoi(0)
    {
        if (init())
            readDrsCalib(fname);
    }

    // Alternative constructor
    factfits(const std::string& fname, const std::string& fout, const std::string& tableName, bool force=false) :
        zfits(fname, fout, tableName, force),
        fOffsetCalibration(0),
        fOffsetStartCellData(0),
        fOffsetData(0),
        fNumRoi(0)
    {
        if (init())
            readDrsCalib(fname);
    }

        const std::vector<int16_t> &GetOffsetCalibration() const { return fOffsetCalibration; }

private:

    void StageRow(size_t row, char* dest)
    {
        zfits::StageRow(row, dest);

        // This file does not contain fact data or no calibration to be applied
        if (fOffsetCalibration.empty())
            return;

        //re-get the pointer to the data to access the offsets
        const uint8_t offset = (row*fTable.bytes_per_row)%4;

        int16_t *startCell = reinterpret_cast<int16_t*>(fBufferRow.data() + offset + fOffsetStartCellData);
        int16_t *data      = reinterpret_cast<int16_t*>(fBufferRow.data() + offset + fOffsetData);

         /*
         for (uint32_t i=0; i<1440*1024; i+=1024, startCell++)
         {
             if (*startCell < 0)
             {
                 data += fNumRoi;
                 continue;
             }
             for (uint32_t j=0; j<fNumRoi; j++, data++)
                 *data += fOffsetCalibration[i + (*startCell+j)%1024];
         }
         */

        // This version is faster because the compilers optimization
        // is not biased by the evaluation of %1024
        for (int ch=0; ch<1440; ch++)
        {
            if (startCell[ch]<0)
            {
                data += fNumRoi;
                continue;
            }

            const int16_t modStart = startCell[ch] % 1024;
            const int16_t *off = fOffsetCalibration.data() + ch*1024;

            const int16_t *cal = off+modStart;
            const int16_t *end_stride = data+fNumRoi;

            if (modStart+fNumRoi>1024)
            {
                while (cal<off+1024)
                    *data++ += *cal++;

                cal = off;
            }
            while (data<end_stride)
                *data++ += *cal++;
        }

    }

    bool init()
    {
        if (!HasKey("NPIX") || !HasKey("NROI"))
            return false;

        if (Get<uint16_t>("NPIX")!=1440)
            return false;

        fNumRoi = Get<uint16_t>("NROI");
        if (fNumRoi>1024)
            return false;

        // check column details for Data
        const Table::Columns::const_iterator it = fTable.cols.find("Data");
        if (it==fTable.cols.end() || it->second.num!=1440*fNumRoi || it->second.type!='I')
            return false;

        // check column details for StartCellData
        const Table::Columns::const_iterator is = fTable.cols.find("StartCellData");
        if (is==fTable.cols.end() || is->second.num!=1440 || is->second.type!='I')
            return false;

        fOffsetStartCellData = is->second.offset;
        fOffsetData          = it->second.offset;

        return true;
    }

    //  Read the Drs calibration data
    void readDrsCalib(const std::string& fileName)
    {
        //should not be mandatory, but improves the perfs a lot when reading not compressed, gzipped files
        if (!IsCompressedFITS())
            return;

        zfits calib(fileName, "ZDrsCellOffsets");

        if (calib.bad())
        {
            clear(rdstate()|std::ios::badbit);
            return;
        }

        if (calib.eof())
            return;

        // Check correct size and format of table
        if (calib.GetNumRows() != 1)
        {
            clear(rdstate()|std::ios::badbit);
#ifdef __EXCEPTIONS
            throw std::runtime_error("Table 'ZDrsCellOffsets' found, but not with one row as expected");
#else
            gLog << ___err___ << "ERROR - Table 'ZDrsCellOffsets' found, but not with one row as expected" << std::endl;
            return;
#endif
        }
        if (calib.GetStr("TTYPE1") != "OffsetCalibration")
        {
            clear(rdstate()|std::ios::badbit);
#ifdef __EXCEPTIONS
            throw std::runtime_error("Table 'ZDrsCellOffsets' found, but first column is not the one expected");
#else
            gLog << ___err___ << "ERROR - Table 'ZDrsCellOffsets' found, but first column is not the one expected" << std::endl;
            return;
#endif
        }
        bool isColumnPresent = false;
        if (calib.HasKey("TFORM1") && calib.GetStr("TFORM1") == "1474560I") isColumnPresent = true;
        if (calib.HasKey("ZFORM1") && calib.GetStr("ZFORM1") == "1474560I") isColumnPresent = true;
        if (!isColumnPresent)  // 1024*1440
        {
            clear(rdstate()|std::ios::badbit);
#ifdef __EXCEPTIONS
            throw std::runtime_error("Table 'ZDrsCellOffsets' has wrong column format (TFROM1)");
#else
            gLog << ___err___ << "ERROR - Table 'ZDrsCellOffsets' has wrong column format (TFORM1)" << std::endl;
            return;
#endif
        }

        fOffsetCalibration.resize(1024*1440);

        calib.SetPtrAddress("OffsetCalibration", fOffsetCalibration.data());
        if (calib.GetNextRow())
            return;

        clear(rdstate()|std::ios::badbit);

#ifdef __EXCEPTIONS
        throw std::runtime_error("Reading column 'OffsetCalibration' failed.");
#else
        gLog << ___err___ << "ERROR - Reading column 'OffsetCalibration' failed." << std::endl;
#endif

    }

    std::vector<int16_t> fOffsetCalibration; ///< integer values of the drs calibration used for compression

    size_t fOffsetStartCellData;
    size_t fOffsetData;

    uint16_t fNumRoi;

#warning Time marker channels currently unhandled

}; //class factfits

#endif
