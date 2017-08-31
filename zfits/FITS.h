/*
 *
 * FITS.h
 *
 * Global fits header
 *
 * Author: lyard
 *
 */

#ifndef MARS_FITS
#define MARS_FITS

#include <stdint.h>
#include <string.h>

#include <vector>
#include <string>

#ifndef __CINT__
#include <unordered_set>
#endif

namespace FITS
{
    static inline bool IsReservedKeyWord(const std::string &key)
    {
#ifndef __CINT__
        static const std::unordered_set<std::string> keys =
        {
            "DATASUM",  "END",      "EXTNAME",  "PCOUNT",   "NAXIS",
            "NAXIS1",   "NAXIS2",   "RAWSUM",   "SIMPLE",   "TFIELDS",
            "THEAP",    "XTENSION", "ZHEAPPTR", "ZNAXIS1",  "ZNAXIS2",
            "ZPCOUNT",  "ZRATIO",   "ZSHRINK",  "ZTABLE",   "ZTILELEN",
        };

        static const std::unordered_set<std::string> short_keys =
        {
            "TFORM", "TUNIT", "TTYPE", "ZCTYP", "ZFORM",
        };

        if (keys.find(key)!=keys.end())
            return true;

        const std::string five = key.substr(0, 5);
        return short_keys.find(five)!=short_keys.end();
#endif
    }

    static inline std::string CommentFromType(char type)
    {
        std::string comment;

        switch (type)
        {
        case 'L': comment = "[1-byte BOOL]";  break;
        case 'A': comment = "[1-byte CHAR]";  break;
        case 'B': comment = "[1-byte BOOL]";  break;
        case 'I': comment = "[2-byte INT]";   break;
        case 'J': comment = "[4-byte INT]";   break;
        case 'K': comment = "[8-byte INT]";   break;
        case 'E': comment = "[4-byte FLOAT]"; break;
        case 'D': comment = "[8-byte FLOAT]"; break;
        case 'Q': comment = "[var. Length]"; break;
        }

        return comment;
    }

    static inline uint32_t SizeFromType(char type)
    {
        size_t size = 0;

        switch (type)
        {
        case 'L': 
        case 'A': 
        case 'B': size =  1; break;
        case 'I': size =  2; break;
        case 'J': 
        case 'E': size =  4; break;
        case 'K': 
        case 'D': size =  8; break;
        case 'Q': size = 16; break;
        }

        return size;
    }

    //Identifier of the compression schemes processes
    enum CompressionProcess_t
    {
        kFactRaw       = 0x0,
        kFactSmoothing = 0x1,
        kFactHuffman16 = 0x2
    };

    //ordering of the columns / rows
    enum RowOrdering_t
    {
        kOrderByCol = 'C',
        kOrderByRow = 'R'
    };

#ifdef __CINT__
    // CINT doesn't like the packed attribute...
    // Therefore we give another hint of the size of the structure
    struct TileHeader  { char dummy[16]; };
    struct BlockHeader { char dummy[10]; };
#else
    //Structure helper for tiles headers
    struct TileHeader
    {
      char     id[4];
      uint32_t numRows;
      uint64_t size;

      TileHeader() {}

      TileHeader(uint32_t nRows,
                 uint64_t s) : id{'T', 'I', 'L', 'E'},
                                 numRows(nRows),
                                 size(s)
      { };
    } __attribute__((__packed__));

    //Structure helper for blocks headers and compresion schemes
    struct BlockHeader
    {
        uint64_t      size;
        char          ordering;
        unsigned char numProcs;
        uint16_t      processings[];

        BlockHeader(uint64_t      s=0,
                    char          o=kOrderByRow,
                    unsigned char n=1) : size(s),
                                         ordering(o),
                                         numProcs(n)
        {
        }
    } __attribute__((__packed__));
#endif

    //Helper structure to simplify the initialization and handling of compressed blocks headers
    struct Compression
    {
        std::vector<uint16_t> sequence;
        BlockHeader header;

        Compression(const std::vector<uint16_t> &seq, const RowOrdering_t &order=kOrderByCol)
            : sequence(seq), header(0, order, seq.size())
        {
        }

        Compression(const CompressionProcess_t &compression=kFactRaw, const RowOrdering_t &order=kOrderByCol)
            : sequence(1), header(0, order, 1)
        {
            sequence[0] = compression;
        }

#ifdef __MARS__ // needed for CINT
        Compression(const int &compression)
            : sequence(1), header(0, kOrderByCol, 1)
        {
            sequence[0] = compression;
        }
#endif

        RowOrdering_t getOrdering() const { return RowOrdering_t(header.ordering); }
        uint32_t getSizeOnDisk() const { return sizeof(BlockHeader) + sizeof(uint16_t)*header.numProcs; }
        CompressionProcess_t getProc(uint32_t i) const { return CompressionProcess_t(sequence[i]); }
        uint16_t getNumProcs() const { return header.numProcs; }

        void SetBlockSize(uint64_t size) { header.size = size; }
        void Memcpy(char *dest) const
        {
            memcpy(dest, &header, sizeof(BlockHeader));
            memcpy(dest+sizeof(BlockHeader), sequence.data(), header.numProcs*sizeof(uint16_t));
        }
    };
};

#endif //_FITS_H_

