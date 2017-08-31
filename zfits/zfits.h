/*
 * zfits.h
 *
 *  Created on: May 16, 2013
 *      Author: lyard
 */

#ifndef MARS_zfits
#define MARS_zfits

#include "fits.h"
#include "huffman.h"

#include "FITS.h"

class zfits : public fits
{
public:

    // Basic constructor
    zfits(const std::string& fname, const std::string& tableName="", bool force=false)
        : fCatalogInitialized(false), fNumTiles(0), fNumRowsPerTile(0), fCurrentRow(-1), fHeapOff(0), fTileSize(0)
    {
        open(fname.c_str());
        Constructor(fname, "", tableName, force);
//        InitCompressionReading();
    }

    // Alternative constructor
    zfits(const std::string& fname, const std::string& fout, const std::string& tableName, bool force=false)
        : fCatalogInitialized(false), fNumTiles(0), fNumRowsPerTile(0), fCurrentRow(-1), fHeapOff(0), fTileSize(0)
    {
        open(fname.c_str());
        Constructor(fname, fout, tableName, force);
//        InitCompressionReading();
    }

    //  Skip the next row
    bool SkipNextRow()
    {
        if (!fTable.is_compressed)
            return fits::SkipNextRow();

        fRow++;
        return true;
    }

    virtual bool IsFileOk() const
    {
        if (!HasKey("RAWSUM"))
            return fits::IsFileOk();

        const bool rawsum = GetStr("RAWSUM") == std::to_string((long long int)fRawsum.val());
        return fits::IsFileOk() && rawsum;
    };

    size_t GetNumRows() const
    {
        return fTable.Get<size_t>(fTable.is_compressed ? "ZNAXIS2" : "NAXIS2");
    }

    size_t GetBytesPerRow() const
    {
        return fTable.Get<size_t>(fTable.is_compressed ? "ZNAXIS1" : "NAXIS1");
    }

protected:

    //  Stage the requested row to internal buffer
    //  Does NOT return data to users
    virtual void StageRow(size_t row, char* dest)
    {
        if (!fTable.is_compressed)
        {
            fits::StageRow(row, dest);
            return;
        }
        ReadBinaryRow(row, dest);
    }

private:

    // Do what it takes to initialize the compressed structured
    void InitCompressionReading()
    {
        fCatalogInitialized = true;

        if (!fTable.is_compressed)
            return;

        //The constructor may have failed
        if (!good())
            return;

        if (fTable.is_compressed)
            for (auto it=fTable.sorted_cols.cbegin(); it!= fTable.sorted_cols.cend(); it++)
            {
                if (it->comp == kCompFACT)
                    continue;

                clear(rdstate()|std::ios::badbit);
#ifdef __EXCEPTIONS
                throw std::runtime_error("Only the FACT compression scheme is handled by this reader.");
#else
                gLog << ___err___ << "ERROR - Only the FACT compression scheme is handled by this reader." << std::endl;
                return;
#endif
            }

        fColumnOrdering.resize(fTable.sorted_cols.size(), FITS::kOrderByRow);

        //Get compressed specific keywords
        fNumTiles       = fTable.is_compressed ? GetInt("NAXIS2") : 0;
        fNumRowsPerTile = fTable.is_compressed ? GetInt("ZTILELEN") : 0;

        //read the file's catalog
        ReadCatalog();

        //give it some space for uncompressing
        AllocateBuffers();
    }

    // Copy decompressed data to location requested by user
    void MoveColumnDataToUserSpace(char* dest, const char* src, const Table::Column& c)
    {
        if (!fTable.is_compressed)
        {
            fits::MoveColumnDataToUserSpace(dest, src, c);
            return;
        }

        memcpy(dest, src, c.num*c.size);
    }

    bool  fCatalogInitialized;

    std::vector<char> fBuffer;           ///<store the uncompressed rows
    std::vector<char> fTransposedBuffer; ///<intermediate buffer to transpose the rows
    std::vector<char> fCompressedBuffer; ///<compressed rows
    std::vector<char> fColumnOrdering;   ///< ordering of the column's rows. Can change from tile to tile.

    size_t fNumTiles;       ///< Total number of tiles
    size_t fNumRowsPerTile; ///< Number of rows per compressed tile
    int64_t fCurrentRow;    ///< current row in memory signed because we need -1
    size_t fShrinkFactor;   ///< shrink factor

    streamoff fHeapOff;           ///< offset from the beginning of the file of the binary data
    streamoff fHeapFromDataStart; ///< offset from the beginning of the data table

    std::vector<std::vector<std::pair<int64_t, int64_t>>> fCatalog;     ///< Catalog, i.e. the main table that points to the compressed data.
    std::vector<size_t>                                   fTileSize;    ///< size in bytes of each compressed tile
    std::vector<std::vector<size_t>>                      fTileOffsets; ///< offset from start of tile of a given compressed column

    Checksum fRawsum;   ///< Checksum of the uncompressed, raw data

    // Get buffer space
    void AllocateBuffers()
    {
        uint32_t buffer_size = fTable.bytes_per_row*fNumRowsPerTile;
        uint32_t compressed_buffer_size = fTable.bytes_per_row*fNumRowsPerTile +
            //use a bit more memory for block headers. 256 char coding the compression sequence max.
            fTable.num_cols*(sizeof(FITS::BlockHeader)+256) +
            //a bit more for the tile headers
            sizeof(FITS::TileHeader) +
            //and a bit more for checksuming
            8;

        if (buffer_size % 4 != 0)
            buffer_size += 4 - (buffer_size%4);

        if (compressed_buffer_size % 4 != 0)
            compressed_buffer_size += 4 - (compressed_buffer_size%4);

        fBuffer.resize(buffer_size);

        fTransposedBuffer.resize(buffer_size);
        fCompressedBuffer.resize(compressed_buffer_size);
    }

    // Read catalog data. I.e. the address of the compressed data inside the heap
    void ReadCatalog()
    {
        std::vector<char> readBuf(16);
        fCatalog.resize(fNumTiles);

        const streampos catalogStart = tellg();

        fChkData.reset();

        //do the actual reading
        for (uint32_t i=0;i<fNumTiles;i++)
            for (uint32_t j=0;j<fTable.num_cols;j++)
            {
                read(readBuf.data(), 2*sizeof(int64_t));
                fChkData.add(readBuf);
                //swap the bytes
                int64_t tempValues[2] = {0,0};
                revcpy<8>(reinterpret_cast<char*>(tempValues), readBuf.data(), 2);
                if (tempValues[0] < 0 || tempValues[1] < 0)
                {
                    clear(rdstate()|std::ios::badbit);
#ifdef __EXCEPTIONS
                    throw std::runtime_error("Negative value in the catalog");
#else
                    gLog << ___err___ << "ERROR - negative value in the catalog" << std::endl;
                    return;
#endif
                }
                //add catalog entry
                fCatalog[i].emplace_back(tempValues[0], tempValues[1]);
            }

        //see if there is a gap before heap data
        fHeapOff = tellg()+fTable.GetHeapShift();
        fHeapFromDataStart = fNumTiles*fTable.num_cols*2*sizeof(int64_t) + fTable.GetHeapShift();

        //check if the catalog has been shrinked
        fShrinkFactor = HasKey("ZSHRINK") ? GetUInt("ZSHRINK") : 1;

        if (fNumRowsPerTile%fShrinkFactor)
        {
            clear(rdstate()|std::ios::badbit);
#ifdef __EXCEPTIONS
            throw std::runtime_error("Rows per tile and shrink factor do not match");
#else
            gLog << ___err___ << "ERROR - Rows per tile and shrink factor do not match" << std::endl;
            return;
#endif
        }

        if (fShrinkFactor>0)
            fNumRowsPerTile /= fShrinkFactor;

        //compute the total size of each compressed tile
        fTileSize.resize(fNumTiles);
        fTileOffsets.resize(fNumTiles);
        for (uint32_t i=0;i<fNumTiles;i++)
        {
            fTileSize[i] = 0;
            for (uint32_t j=0;j<fTable.num_cols;j++)
            {
                fTileSize[i] += fCatalog[i][j].first;
                fTileOffsets[i].emplace_back(fCatalog[i][j].second - fCatalog[i][0].second);
            }
        }

        if (!fCopy.is_open())
            return;

        //write catalog and heap gap to target file
        seekg(catalogStart);

        const size_t catSize = fTable.GetHeapShift() + fTable.total_bytes;

        std::vector<char> buf(catSize);
        read(buf.data(), catSize);

        fCopy.write(buf.data(), catSize);
        if (!fCopy)
            clear(rdstate()|std::ios::badbit);
    }

    //overrides fits.h method with empty one
    //work is done in ReadBinaryRow because it requires volatile data from ReadBinaryRow
    virtual void WriteRowToCopyFile(size_t row)
    {
        if (row == fRow+1)
            fRawsum.add(fBufferRow);
    }

    // Compressed version of the read row, even files with shrunk catalogs
    // can be read fully sequentially so that streaming, e.g. through
    // stdout/stdin, is possible.
    bool ReadBinaryRow(const size_t &rowNum, char *bufferToRead)
    {
        if (rowNum >= GetNumRows())
            return false;

        if (!fCatalogInitialized)
            InitCompressionReading();

        // Book keeping, where are we?
        const int64_t requestedTile      = rowNum        / fNumRowsPerTile;
        const int64_t currentTile        = fCurrentRow   / fNumRowsPerTile;

        const int64_t requestedSuperTile = requestedTile / fShrinkFactor;
        //const int64_t currentSuperTile   = currentTile   / fShrinkFactor;

        const int64_t requestedSubTile   = requestedTile % fShrinkFactor;
        //const int64_t currentSubTile     = currentTile   % fShrinkFactor;

        // Is this the first tile we read at all?
        const bool isFirstTile = fCurrentRow<0;

        // Is this just the next tile in the sequence?
        const bool isNextTile = requestedTile==currentTile+1 || isFirstTile;

        fCurrentRow = rowNum;

        // Do we have to read a new tile from disk?
        if (requestedTile!=currentTile || isFirstTile)
        {
            //skip to the beginning of the tile
            const int64_t superTileStart = fCatalog[requestedSuperTile][0].second - sizeof(FITS::TileHeader);

            std::vector<size_t> offsets = fTileOffsets[requestedSuperTile];

            // If this is a sub tile we might have to step forward a bit and
            // seek for the sub tile. If we were just reading the previous one
            // we can skip that.
            if (!isNextTile || isFirstTile)
            {
                // step to the beginnig of the super tile
                seekg(fHeapOff+superTileStart);

                // If there are sub tiles we might have to seek through the super tile
                for (uint32_t k=0; k<requestedSubTile; k++)
                {
                    // Read header
                    FITS::TileHeader header;
                    read((char*)&header, sizeof(FITS::TileHeader));

                    // Skip to the next header
                    seekg(header.size-sizeof(FITS::TileHeader), cur);
                }
            }

            // this is now the beginning of the sub-tile we want to read
            const int64_t subTileStart = tellg() - fHeapOff;
            // calculate the 32 bits offset of the current tile.
            const uint32_t offset = (subTileStart + fHeapFromDataStart)%4;

            // start of destination buffer (padding comes later)
            char *destBuffer = fCompressedBuffer.data()+offset;

            // Store the current tile size once known
            size_t currentTileSize = 0;

            // If this is a request for a sub tile which is not cataloged
            // recalculate the offsets from the buffer, once read
            if (requestedSubTile>0)
            {
                // Read header
                read(destBuffer, sizeof(FITS::TileHeader));

                // Get size of tile
                currentTileSize = reinterpret_cast<FITS::TileHeader*>(destBuffer)->size;

                // now read the remaining bytes of this tile
                read(destBuffer+sizeof(FITS::TileHeader), currentTileSize-sizeof(FITS::TileHeader));

                // Calculate the offsets recursively
                offsets[0] = 0;

                //skip through the columns
                for (size_t i=0; i<fTable.num_cols-1; i++)
                {
                    //zero sized column do not have headers. Skip it
                    if (fTable.sorted_cols[i].num == 0)
                    {
                        offsets[i+1] = offsets[i];
                        continue;
                    }

                    const char *pos = destBuffer + offsets[i] + sizeof(FITS::TileHeader);
                    offsets[i+1] = offsets[i] + reinterpret_cast<const FITS::BlockHeader*>(pos)->size;
                }
            }
            else
            {
                // If we are reading the first tile of a super tile, all information
                // is already available.
                currentTileSize = fTileSize[requestedSuperTile] + sizeof(FITS::TileHeader);
                read(destBuffer, currentTileSize);
            }


            // If we are reading sequentially, calcualte checksum
            if (isNextTile)
            {
                // Padding for checksum calculation
                memset(fCompressedBuffer.data(),   0, offset);
                memset(destBuffer+currentTileSize, 0, fCompressedBuffer.size()-currentTileSize-offset);
                fChkData.add(fCompressedBuffer);
            }

            // Check if we are writing a copy of the file
            if (isNextTile && fCopy.is_open() && fCopy.good())
            {
                fCopy.write(fCompressedBuffer.data()+offset, currentTileSize);
                if (!fCopy)
                    clear(rdstate()|std::ios::badbit);
            }
            else
                if (fCopy.is_open())
                    clear(rdstate()|std::ios::badbit);


            // uncompress  the buffer
            const uint32_t thisRoundNumRows = (GetNumRows()<fCurrentRow + fNumRowsPerTile) ? GetNumRows()%fNumRowsPerTile : fNumRowsPerTile;
            if (!UncompressBuffer(offsets, thisRoundNumRows, offset+sizeof(FITS::TileHeader)))
                return false;

            // pointer to column (source buffer)
            const char *src = fTransposedBuffer.data();

            uint32_t i=0;
            for (auto it=fTable.sorted_cols.cbegin(); it!=fTable.sorted_cols.cend(); it++, i++)
            {
                char *buffer = fBuffer.data() + it->offset; // pointer to column (destination buffer)

                switch (fColumnOrdering[i])
                {
                case FITS::kOrderByRow:
                    // regular, "semi-transposed" copy
                    for (char *dest=buffer; dest<buffer+thisRoundNumRows*fTable.bytes_per_row; dest+=fTable.bytes_per_row) // row-by-row
                    {
                        memcpy(dest, src, it->bytes);
                        src += it->bytes;  // next column
                    }
                    break;

                case FITS::kOrderByCol:
                    // transposed copy
                    for (char *elem=buffer; elem<buffer+it->bytes; elem+=it->size) // element-by-element (arrays)
                    {
                        for (char *dest=elem; dest<elem+thisRoundNumRows*fTable.bytes_per_row; dest+=fTable.bytes_per_row) // row-by-row
                        {
                                memcpy(dest, src, it->size);
                                src += it->size; // next element
                        }
                    }
                    break;

                default:
                    clear(rdstate()|std::ios::badbit);

                    std::ostringstream str;
                    str << "Unkown column ordering scheme found (i=" << i << ", " << fColumnOrdering[i] << ")";
#ifdef __EXCEPTIONS
                    throw std::runtime_error(str.str());
#else
                    gLog << ___err___ << "ERROR - " << str.str() << std::endl;
                    return false;
#endif
                };
            }
        }

        //Data loaded and uncompressed. Copy it to destination
        memcpy(bufferToRead, fBuffer.data()+fTable.bytes_per_row*(fCurrentRow%fNumRowsPerTile), fTable.bytes_per_row);
        return good();
    }

    // Read a bunch of uncompressed data
    uint32_t UncompressUNCOMPRESSED(char*       dest,
                                    const char* src,
                                    uint32_t    numElems,
                                    uint32_t    sizeOfElems)
    {
        memcpy(dest, src, numElems*sizeOfElems);
        return numElems*sizeOfElems;
    }

    // Read a bunch of data compressed with the Huffman algorithm
    uint32_t UncompressHUFFMAN16(char*       dest,
                                 const char* src,
                                 uint32_t    numChunks)
    {
        std::vector<uint16_t> uncompressed;

        //read compressed sizes (one per row)
        const uint32_t* compressedSizes = reinterpret_cast<const uint32_t*>(src);
        src += sizeof(uint32_t)*numChunks;

        //uncompress the rows, one by one
        uint32_t sizeWritten = 0;
        for (uint32_t j=0;j<numChunks;j++)
        {
            Huffman::Decode(reinterpret_cast<const unsigned char*>(src), compressedSizes[j], uncompressed);

            memcpy(dest, uncompressed.data(), uncompressed.size()*sizeof(uint16_t));

            sizeWritten += uncompressed.size()*sizeof(uint16_t);
            dest        += uncompressed.size()*sizeof(uint16_t);
            src         += compressedSizes[j];
        }
        return sizeWritten;
    }

    // Apply the inverse transform of the integer smoothing
    uint32_t UnApplySMOOTHING(int16_t*   data,
                              uint32_t   numElems)
    {
        //un-do the integer smoothing
        for (uint32_t j=2;j<numElems;j++)
            data[j] = data[j] + (data[j-1]+data[j-2])/2;

        return numElems*sizeof(uint16_t);
    }

    // Data has been read from disk. Uncompress it !
    bool UncompressBuffer(const std::vector<size_t> &offsets,
                          const uint32_t &thisRoundNumRows,
                          const uint32_t offset)
    {
        char *dest = fTransposedBuffer.data();

        //uncompress column by column
        for (uint32_t i=0; i<fTable.sorted_cols.size(); i++)
        {
            const fits::Table::Column &col = fTable.sorted_cols[i];
            if (col.num == 0)
                continue;

            //get the compression flag
            const int64_t compressedOffset = offsets[i]+offset;

            const FITS::BlockHeader* head = reinterpret_cast<FITS::BlockHeader*>(&fCompressedBuffer[compressedOffset]);

            fColumnOrdering[i] = head->ordering;

            const uint32_t numRows = (head->ordering==FITS::kOrderByRow) ? thisRoundNumRows : col.num;
            const uint32_t numCols = (head->ordering==FITS::kOrderByCol) ? thisRoundNumRows : col.num;

            const char *src = fCompressedBuffer.data()+compressedOffset+sizeof(FITS::BlockHeader)+sizeof(uint16_t)*head->numProcs;

            for (int32_t j=head->numProcs-1;j >= 0; j--)
            {
                uint32_t sizeWritten=0;

                switch (head->processings[j])
                {
                case FITS::kFactRaw:
                    sizeWritten = UncompressUNCOMPRESSED(dest, src, numRows*numCols, col.size);
                    break;

                case FITS::kFactSmoothing:
                    sizeWritten = UnApplySMOOTHING(reinterpret_cast<int16_t*>(dest), numRows*numCols);
                    break;

                case FITS::kFactHuffman16:
                    sizeWritten = UncompressHUFFMAN16(dest, src, numRows);
                    break;

                default:
                    clear(rdstate()|std::ios::badbit);

                    std::ostringstream str;
                    str << "Unknown processing applied to data (col=" << i << ", proc=" << j << "/" << (int)head->numProcs;
#ifdef __EXCEPTIONS
                    throw std::runtime_error(str.str());
#else
                    gLog << ___err___ << "ERROR - " << str.str() << std::endl;
                    return false;
#endif
                }
                //increment destination counter only when processing done.
                if (j==0)
                    dest+= sizeWritten;
            }
        }

        return true;
    }

    void CheckIfFileIsConsistent(bool update_catalog=false)
    {
        //goto start of heap
        const streamoff whereAreWe = tellg();
        seekg(fHeapOff);

        //init number of rows to zero
        uint64_t numRows = 0;

        //get number of columns from header
        const size_t numCols = fTable.num_cols;

        std::vector<std::vector<std::pair<int64_t, int64_t> > > catalog;

        FITS::TileHeader tileHead;
        FITS::BlockHeader columnHead;

        streamoff offsetInHeap = 0;
        //skip through the heap
        while (true)
        {
            read((char*)(&tileHead), sizeof(FITS::TileHeader));
            //end of file
            if (!good())
                break;

            //padding or corrupt data
            if (memcmp(tileHead.id, "TILE", 4))
            {
                clear(rdstate()|std::ios::badbit);
                break;
            }

            //a new tile begins here
            catalog.emplace_back();
            offsetInHeap += sizeof(FITS::TileHeader);

            //skip through the columns
            for (size_t i=0;i<numCols;i++)
            {
                //zero sized column do not have headers. Skip it
                if (fTable.sorted_cols[i].num == 0)
                {
                    catalog.back().emplace_back(0,0);
                    continue;
                }

                //read column header
                read((char*)(&columnHead), sizeof(FITS::BlockHeader));

                //corrupted tile
                if (!good())
                    break;

                catalog.back().emplace_back((int64_t)(columnHead.size),offsetInHeap);
                offsetInHeap += columnHead.size;
                seekg(fHeapOff+offsetInHeap);
            }

            //if we ain't good, this means that something went wrong inside the current tile.
            if (!good())
            {
                catalog.pop_back();
                break;
            }
            //current tile is complete. Add rows
            numRows += tileHead.numRows;
        }

        if (numRows != fTable.num_rows)
        {
            clear(rdstate()|std::ios::badbit);
            std::ostringstream str;
            str << "Heap data does not agree with header: " << numRows << " calculated vs " << fTable.num_rows << " from header.";
#ifdef __EXCEPTIONS
                    throw std::runtime_error(str.str());
#else
                    gLog << ___err___ << "ERROR - " << str.str() << std::endl;
                    return;
#endif
        }

        if (update_catalog)
        {
            fCatalog = catalog;
            //clear the bad bit before seeking back (we hit eof)
            clear();
            seekg(whereAreWe);
            return;
        }

        if (catalog.size() != fCatalog.size())
        {
                    clear(rdstate()|std::ios::badbit);
#ifdef __EXCEPTIONS
                    throw std::runtime_error("Heap data does not agree with header.");
#else
                    gLog << ___err___ << "ERROR - Heap data does not agree with header." << std::endl;
                    return;
#endif
        }

        for (uint32_t i=0;i<catalog.size(); i++)
            for (uint32_t j=0;j<numCols;j++)
            {
                if (catalog[i][j].first  != fCatalog[i][j].first ||
                    catalog[i][j].second != fCatalog[i][j].second)
                {
                    clear(rdstate()|std::ios::badbit);
#ifdef __EXCEPTIONS
                    throw std::runtime_error("Heap data does not agree with header.");
#else
                    gLog << ___err___ << "ERROR - Heap data does not agree with header." << std::endl;
                    return;
#endif
                }
            }
        //go back to start of heap
        //clear the bad bit before seeking back (we hit eof)
        clear();
        seekg(whereAreWe);
    }

};//class zfits

#endif 
