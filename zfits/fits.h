#ifndef MARS_fits
#define MARS_fits

#include <stdint.h>

#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "FITS.h"
#include "checksum.h"

class fits : public ifstream
{
public:
    //I know I know, you're going to yiell that this does not belong here.
    //It will belong in the global scope eventually, and it makes the coding of zfits much simpler this way.
    enum Compression_t
    {
        kCompUnknown,
        kCompFACT
    };

    enum fitsstate
    {
        throwbit = 8
    };

    struct Entry
    {
        char        type;
        std::string value;
        std::string comment;
        std::string fitsString;

        template<typename T>
            T Get() const
        {
            T t;

            std::istringstream str(value);
            str >> t;

            return t;
        }
    };

    struct Table
    {
        off_t offset;

        bool is_compressed;

        std::string name;
        size_t bytes_per_row;
        size_t num_rows;
        size_t num_cols;
        size_t total_bytes; // NAXIS1*NAXIS2

        struct Column
        {
            size_t offset;
            size_t num;
            size_t size;
            size_t bytes;  // num*size
            char   type;
            std::string unit;
            Compression_t comp;
        };

        typedef std::map<std::string, Entry>  Keys;
        typedef std::map<std::string, Column> Columns;
        typedef std::vector<Column> SortedColumns;

        Columns cols;
        SortedColumns sorted_cols;
        Keys    keys;

        int64_t datasum;

        std::string Trim(const std::string &str, char c=' ') const
        {
            // Trim Both leading and trailing spaces
            const size_t pstart = str.find_first_not_of(c); // Find the first character position after excluding leading blank spaces
            const size_t pend   = str.find_last_not_of(c);  // Find the first character position from reverse af

            // if all spaces or empty return an empty string
            if (std::string::npos==pstart || std::string::npos==pend)
                return std::string();

            return str.substr(pstart, pend-pstart+1);
        }

        bool Check(const std::string &key, char type, const std::string &value="") const
        {
            const Keys::const_iterator it = keys.find(key);
            if (it==keys.end())
            {
                std::ostringstream str;
                str << "Key '" << key << "' not found.";
                throw std::runtime_error(str.str());
            }

            if (it->second.type!=type)
            {
                std::ostringstream str;
                str << "Wrong type for key '" << key << "': expected " << type << ", found " << it->second.type << ".";
                throw std::runtime_error(str.str());
            }

            if (!value.empty() && it->second.value!=value)
            {
                std::ostringstream str;
                str << "Wrong value for key '" << key << "': expected " << value << ", found " << it->second.value << ".";
                throw std::runtime_error(str.str());
            }

            return true;
        }

        Keys ParseBlock(const std::vector<std::string> &vec) const
        {
            Keys rc;

            for (unsigned int i=0; i<vec.size(); i++)
            {
                const std::string key = Trim(vec[i].substr(0,8));
                // Keywords without a value, like COMMENT / HISTORY
                if (vec[i].substr(8,2)!="= ")
                    continue;

                char type = 0;

                std::string com;
                std::string val = Trim(vec[i].substr(10));

                if (val[0]=='\'')
                {
                    // First skip all '' in the string
                    size_t p = 1;
                    while (1)
                    {
                        const size_t pp = val.find_first_of('\'', p);
                        if (pp==std::string::npos)
                            break;

                        p = val[pp+1]=='\'' ? pp+2 : pp+1;
                    }

                    // Now find the comment
                    const size_t ppp = val.find_first_of('/', p);

                    // Set value, comment and type
                    // comments could be just spaces. take care of this.
                    if (ppp!=std::string::npos && val.size()!=ppp+1)
                        com = Trim(val.substr(ppp+1));

                    val  = Trim(val.substr(1, p-2));
                    type = 'T';
                }
                else
                {
                    const size_t p = val.find_first_of('/');

                    if (p!=std::string::npos && val.size()!=p+1)
                        com = Trim(val.substr(p+2));

                    val = Trim(val.substr(0, p));

                    if (val.empty() || val.find_first_of('T')!=std::string::npos || val.find_first_of('F')!=std::string::npos)
                        type = 'B';
                    else
                        type = val.find_last_of('.')==std::string::npos ? 'I' : 'F';
                }

                const Entry e = { type, val, com, vec[i] };

                rc[key] = e;
            }

            return rc;
        }

        Table() : offset(0), is_compressed(false) { }
        Table(const std::vector<std::string> &vec, off_t off) : offset(off),
            keys(ParseBlock(vec))
        {
            is_compressed = HasKey("ZTABLE") && Check("ZTABLE", 'B', "T");

            if (!Check("XTENSION", 'T', "BINTABLE") ||
                !Check("NAXIS",    'I', "2")        ||
                !Check("BITPIX",   'I', "8")        ||
                !Check("GCOUNT",   'I', "1")        ||
                !Check("EXTNAME",  'T')             ||
                !Check("NAXIS1",   'I')             ||
                !Check("NAXIS2",   'I')             ||
                !Check("TFIELDS",  'I'))
                return;

            if (is_compressed)
            {
                if (!Check("ZNAXIS1", 'I') ||
                    !Check("ZNAXIS2", 'I') ||
                    !Check("ZPCOUNT", 'I', "0"))
                    return;
            }
            else
            {
                if (!Check("PCOUNT", 'I', "0"))
                    return;
            }

            total_bytes   = Get<size_t>("NAXIS1")*Get<size_t>("NAXIS2");
            bytes_per_row = is_compressed ? Get<size_t>("ZNAXIS1") : Get<size_t>("NAXIS1");
            num_rows      = is_compressed ? Get<size_t>("ZNAXIS2") : Get<size_t>("NAXIS2");
            num_cols      = Get<size_t>("TFIELDS");
            datasum       = Get<int64_t>("DATASUM", -1);
            size_t bytes = 0;

            const std::string tFormName = is_compressed ? "ZFORM" : "TFORM";
            for (long long unsigned int i=1; i<=num_cols; i++)
            {
                const std::string num(std::to_string(i));

                if (!Check("TTYPE"+num, 'T') ||
                    !Check(tFormName+num, 'T'))
                    return;

                const std::string id   = Get<std::string>("TTYPE"+num);
                const std::string fmt  = Get<std::string>(tFormName+num);
                const std::string unit = Get<std::string>("TUNIT"+num, "");
                const std::string comp = Get<std::string>("ZCTYP"+num, "");

                Compression_t compress = kCompUnknown;
                if (comp == "FACT")
                    compress = kCompFACT;

                std::istringstream sin(fmt);
                int n = 0;
                sin >> n;
                if (!sin)
                    n = 1;

                const char type = fmt[fmt.length()-1];

                size_t size = 0;
                switch (type)
                {
                    // We could use negative values to mark floats
                    // otheriwse we could just cast them to int64_t?
                case 'L':                  // logical
                case 'A':                  // char
                case 'B': size = 1; break; // byte
                case 'I': size = 2; break; // short
                case 'J': size = 4; break; // int
                case 'K': size = 8; break; // long long
                case 'E': size = 4; break; // float
                case 'D': size = 8; break; // double
                // case 'X': size =  n; break; // bits (n=number of bytes needed to contain all bits)
                // case 'C': size =  8; break; // complex float
                // case 'M': size = 16; break; // complex double
                // case 'P': size =  8; break; // array descriptor (32bit)
                // case 'Q': size = 16; break; // array descriptor (64bit)
                default:
                    {
                        std::ostringstream str;
                        str << "FITS format TFORM='" << fmt << "' not yet supported.";
                        throw std::runtime_error(str.str());
                    }
                }

                const Table::Column col = { bytes, size_t(n), size, n*size, type, unit, compress};

                cols[id] = col;
                sorted_cols.emplace_back(col);
                bytes  += n*size;
            }

            if (bytes!=bytes_per_row)
            {
                std::ostringstream str;
                str << "Sum of bytes in columns [" << bytes << "] does not match (Z)NAXIS2 [" << bytes_per_row << "].";

                throw std::runtime_error(str.str());
            }

            name = Get<std::string>("EXTNAME");
        }

        void PrintKeys(bool display_all=false) const
        {
        }

        void PrintColumns() const
        {
        }

        operator bool() const { return !name.empty(); }

        bool HasKey(const std::string &key) const
        {
            return keys.find(key)!=keys.end();
        }

        bool HasColumn(const std::string& col) const
        {
            return cols.find(col)!=cols.end();
        }

        const Columns &GetColumns() const
        {
            return cols;
        }

        const Keys &GetKeys() const
        {
            return keys;
        }

        // Values of keys are always signed
        template<typename T>
            T Get(const std::string &key) const
        {
            const Keys::const_iterator it = keys.find(key);
            if (it==keys.end())
            {
                std::ostringstream str;
                str << "Key '" << key << "' not found.";
                throw std::runtime_error(str.str());
            }
            return it->second.Get<T>();
        }

        // Values of keys are always signed
        template<typename T>
            T Get(const std::string &key, const T &deflt) const
        {
            const Keys::const_iterator it = keys.find(key);
            return it==keys.end() ? deflt :it->second.Get<T>();
        }

        size_t GetN(const std::string &key) const
        {
            const Columns::const_iterator it = cols.find(key);
            return it==cols.end() ? 0 : it->second.num;
        }



        // There may be a gap between the main table and the start of the heap:
        // this computes the offset
        streamoff GetHeapShift() const
        {
            if (!HasKey("ZHEAPPTR"))
                return 0;

            const size_t shift = Get<size_t>("ZHEAPPTR");
            return shift <= total_bytes ? 0 : shift - total_bytes;
        }

        // return total number of bytes 'all inclusive'
        streamoff GetTotalBytes() const
        {
            //get offset of special data area from start of main table
            const streamoff shift = GetHeapShift();

            //and special data area size
            const streamoff size  = HasKey("PCOUNT") ? Get<streamoff>("PCOUNT") : 0;

            // Get the total size
            const streamoff total = total_bytes + size + shift;

            // check for padding
            if (total%2880==0)
                return total;

            // padding necessary
            return total + (2880 - (total%2880));
        }
    };

    void Exception(const std::string &txt)
    {
        if (exceptions()&throwbit)
            throw std::runtime_error(txt);
    }

    // Public for the root dictionary
    typedef std::pair<void*, Table::Column> Address;
    typedef std::vector<Address> Addresses;
    typedef std::unordered_map<std::string, void*> Pointers;

protected:
    std::ofstream fCopy;
    std::vector<std::string> fListOfTables; // List of skipped tables. Last table is open table

    Table fTable;

    //map<void*, Table::Column> fAddresses;
    Addresses fAddresses;

    Pointers fPointers;

    std::vector<std::vector<char>> fGarbage;

    std::vector<char> fBufferRow;
    std::vector<char> fBufferDat;

    size_t fRow;

    Checksum fChkHeader;
    Checksum fChkData;

    bool ReadBlock(std::vector<std::string> &vec)
    {
        int endtag = 0;
        for (int i=0; i<36; i++)
        {
            char c[81];
            c[80] = 0;
            read(c, 80);
            if (!good())
                break;

            fChkHeader.add(c, 80);

            std::string str(c);

            if (endtag==2 || str=="END                                                                             ")
            {
                endtag = 2; // valid END tag found
                continue;
            }

            if (endtag==1 || str=="                                                                                ")
            {
                endtag = 1; // end tag not found, but expected to be there
                continue;
            }

            vec.emplace_back(str);
        }

        // Make sure that no empty vector is returned
        if (endtag && vec.size()%36==0)
            vec.emplace_back("END     = '' / ");

        return endtag==2;
    }

    std::string Compile(const std::string &key, int16_t i=-1) const
    {
        return i<0 ? key : key+std::to_string(i);
    }

    void Constructor(const std::string &fname, std::string fout="", const std::string& tableName="", bool force=false)
    {
        char simple[10];
        read(simple, 10);
        if (!good())
            return;

        EnableAddressExceptions();

        if (memcmp(simple, "SIMPLE  = ", 10))
        {
            clear(rdstate()|std::ios::badbit);
            throw std::runtime_error("File is not a FITS file.");
        }

        seekg(0);

        while (good())
        {
            std::vector<std::string> block;
            while (1)
            {
                // If we search for a table, we implicitly assume that
                // not finding the table is not an error. The user
                // can easily check that by eof() && !bad()
                peek();
                if (eof() && !bad() && !tableName.empty())
                {
                    break;
                }
                // FIXME: Set limit on memory consumption
                const int rc = ReadBlock(block);
                if (!good())
                {
                    clear(rdstate()|std::ios::badbit);
                    throw std::runtime_error("FITS file corrupted.");
                }

                if (block.size()%36)
                {
                    if (!rc && !force)
                    {
                        clear(rdstate()|std::ios::badbit);
                        throw std::runtime_error("END keyword missing in FITS header.");
                    }
                    break;
                }
            }

            if (block.empty())
                break;

            if (block[0].substr(0, 9)=="SIMPLE  =")
            {
                fChkHeader.reset();
                continue;
            }

            if (block[0].substr(0, 9)=="XTENSION=")
            {
                fTable = Table(block, tellg());
                fRow   = (size_t)-1;

                if (!fTable)
                {
                    clear(rdstate()|std::ios::badbit);
                    return;
                }

                const std::string &tname = fTable.Get<std::string>("EXTNAME");

                fListOfTables.emplace_back(tname);

                // Check for table name. Skip until eof or requested table are found.
                // skip the current table?
                if ((!tableName.empty() && tableName!=tname) || (tableName.empty() && "ZDrsCellOffsets"==tname))
                {
                    const streamoff skip = fTable.GetTotalBytes();
                    seekg(skip, std::ios_base::cur);

                    fChkHeader.reset();

                    continue;
                }

                fBufferRow.resize(fTable.bytes_per_row + 8-fTable.bytes_per_row%4);
                fBufferDat.resize(fTable.bytes_per_row);

                break;
            }
        }

        if (fout.empty())
            return;

        if (*fout.rbegin()=='/')
        {
            const size_t p = fname.find_last_of('/');
            fout.append(fname.substr(p+1));
        }

        fCopy.open(fout);
        if (!fCopy)
        {
            clear(rdstate()|std::ios::badbit);
            throw std::runtime_error("Could not open output file.");
        }

        const streampos p = tellg();
        seekg(0);

        std::vector<char> buf(p);
        read(buf.data(), p);

        fCopy.write(buf.data(), p);
        if (!fCopy)
            clear(rdstate()|std::ios::badbit);
    }

public:
    fits(const std::string &fname, const std::string& tableName="", bool force=false) : ifstream(fname.c_str())
    {
        Constructor(fname, "", tableName, force);
        if ((fTable.is_compressed ||fTable.name=="ZDrsCellOffsets") && !force)
        {
            throw std::runtime_error("Trying to read a compressed fits with the base fits class. Use factfits instead.");
            clear(rdstate()|std::ios::badbit);
        }
    }

    fits(const std::string &fname, const std::string &fout, const std::string& tableName, bool force=false) : ifstream(fname.c_str())
    {
        Constructor(fname, fout, tableName, force);
        if ((fTable.is_compressed || fTable.name=="ZDrsCellOffsets") && !force)
        {
            throw std::runtime_error("Trying to read a compressed fits with the base fits class. Use factfits instead.");
            clear(rdstate()|std::ios::badbit);
        }
    }

    fits() : ifstream()
    {

    }

    ~fits()
    {
        std::copy(std::istreambuf_iterator<char>(*this),
                  std::istreambuf_iterator<char>(),
                  std::ostreambuf_iterator<char>(fCopy));
    }

    virtual void StageRow(size_t row, char* dest)
    {
        // if (row!=fRow+1) // Fast seeking is ensured by ifstream
        seekg(fTable.offset+row*fTable.bytes_per_row);
        read(dest, fTable.bytes_per_row);
        //fin.clear(fin.rdstate()&~ios::eofbit);
    }

    virtual void WriteRowToCopyFile(size_t row)
    {
        if (row==fRow+1)
        {
            const uint8_t offset = (row*fTable.bytes_per_row)%4;

            fChkData.add(fBufferRow);
            if (fCopy.is_open() && fCopy.good())
                fCopy.write(fBufferRow.data()+offset, fTable.bytes_per_row);
            if (!fCopy)
                clear(rdstate()|std::ios::badbit);
        }
        else
            if (fCopy.is_open())
                clear(rdstate()|std::ios::badbit);
    }

    void ZeroBufferForChecksum(std::vector<char>& vec, const uint64_t extraZeros=0)
    {
        auto ib = vec.begin();
        auto ie = vec.end();

        *ib++ = 0;
        *ib++ = 0;
        *ib++ = 0;
        *ib   = 0;

        for (uint64_t i=0;i<extraZeros+8;i++)
            *--ie = 0;
    }

    uint8_t ReadRow(size_t row)
    {
        // For the checksum we need everything to be correctly aligned
        const uint8_t offset = (row*fTable.bytes_per_row)%4;

        ZeroBufferForChecksum(fBufferRow);

        StageRow(row, fBufferRow.data()+offset);

        WriteRowToCopyFile(row);

        fRow = row;

        return offset;
    }

    template<size_t N>
        void revcpy(char *dest, const char *src, const int &num)
    {
        const char *pend = src + num*N;
        for (const char *ptr = src; ptr<pend; ptr+=N, dest+=N)
            std::reverse_copy(ptr, ptr+N, dest);
    }

    virtual void MoveColumnDataToUserSpace(char *dest, const char *src, const Table::Column& c)
    {
        // Let the compiler do some optimization by
        // knowing that we only have 1, 2, 4 and 8
        switch (c.size)
        {
        case 1: memcpy   (dest, src, c.bytes); break;
        case 2: revcpy<2>(dest, src, c.num);   break;
        case 4: revcpy<4>(dest, src, c.num);   break;
        case 8: revcpy<8>(dest, src, c.num);   break;
        }
    }

    virtual bool GetRow(size_t row, bool check=true)
    {
        if (check && row>=fTable.num_rows)
            return false;

        const uint8_t offset = ReadRow(row);
        if (!good())
            return good();

        const char *ptr = fBufferRow.data() + offset;

        for (Addresses::const_iterator it=fAddresses.cbegin(); it!=fAddresses.cend(); it++)
        {
            const Table::Column &c = it->second;

            const char *src = ptr + c.offset;
            char *dest = reinterpret_cast<char*>(it->first);

            MoveColumnDataToUserSpace(dest, src, c);
        }

        return good();
    }

    bool GetNextRow(bool check=true)
    {
        return GetRow(fRow+1, check);
    }

    virtual bool SkipNextRow()
    {
        seekg(fTable.offset+(++fRow)*fTable.bytes_per_row);
        return good();
    }

    static bool Compare(const Address &p1, const Address &p2)
    {
        return p1.first>p2.first;
    }

    template<class T, class S>
    const T &GetAs(const std::string &name)
    {
        return *reinterpret_cast<S*>(fPointers[name]);
    }

    void EnableAddressExceptions(bool b=true)
    {
        if (b)
            exceptions(iostate(throwbit));
        else
            exceptions(iostate(exceptions()&~throwbit));
    }

    void DisableAddressExceptions()
    {
        EnableAddressExceptions(false);
    }

    void *SetPtrAddress(const std::string &name)
    {
        if (fTable.cols.count(name)==0)
        {
            std::ostringstream str;
            str << "SetPtrAddress('" << name << "') - Column not found.";
            Exception(str.str());
            return NULL;
        }

        Pointers::const_iterator it = fPointers.find(name);
        if (it!=fPointers.end())
            return it->second;

        fGarbage.emplace_back(fTable.cols[name].bytes);

        void *ptr = fGarbage.back().data();

        fPointers[name] = ptr;
        fAddresses.emplace_back(ptr, fTable.cols[name]);
        sort(fAddresses.begin(), fAddresses.end(), Compare);
        return ptr;
    }

    template<typename T>
    bool SetPtrAddress(const std::string &name, T *ptr, size_t cnt)
    {
        if (fTable.cols.count(name)==0)
        {
            std::ostringstream str;
            str << "SetPtrAddress('" << name << "') - Column not found.";
            Exception(str.str());
            return false;
        }

        if (sizeof(T)!=fTable.cols[name].size)
        {
            std::ostringstream str;
            str << "SetPtrAddress('" << name << "') - Element size mismatch: expected "
                << fTable.cols[name].size << " from header, got " << sizeof(T);
            Exception(str.str());
            return false;
        }

        if (cnt!=fTable.cols[name].num)
        {
            std::ostringstream str;
            str << "SetPtrAddress('" << name << "') - Element count mismatch: expected "
                << fTable.cols[name].num << " from header, got " << cnt;
            Exception(str.str());
            return false;
        }

        //fAddresses[ptr] = fTable.cols[name];
        fPointers[name] = ptr;
        fAddresses.emplace_back(ptr, fTable.cols[name]);
        sort(fAddresses.begin(), fAddresses.end(), Compare);
        return true;
    }

    template<class T>
    bool SetRefAddress(const std::string &name, T &ptr)
    {
        return SetPtrAddress(name, &ptr, sizeof(ptr)/sizeof(T));
    }

    template<typename T>
    bool SetVecAddress(const std::string &name, std::vector<T> &vec)
    {
        return SetPtrAddress(name, vec.data(), vec.size());
    }

    template<typename T>
        T Get(const std::string &key) const
    {
        return fTable.Get<T>(key);
    }

    template<typename T>
        T Get(const std::string &key, const std::string &deflt) const
    {
        return fTable.Get<T>(key, deflt);
    }

    bool SetPtrAddress(const std::string &name, void *ptr, size_t cnt=0)
    {
        if (fTable.cols.count(name)==0)
        {
            std::ostringstream str;
            str <<"SetPtrAddress('" << name << "') - Column not found.";
            Exception(str.str());
            return false;
        }

        if (cnt && cnt!=fTable.cols[name].num)
        {
            std::ostringstream str;
            str << "SetPtrAddress('" << name << "') - Element count mismatch: expected "
                << fTable.cols[name].num << " from header, got " << cnt;
            Exception(str.str());
            return false;
        }

        fPointers[name] = ptr;
        fAddresses.emplace_back(ptr, fTable.cols[name]);
        sort(fAddresses.begin(), fAddresses.end(), Compare);
        return true;
    }

    bool     HasKey(const std::string &key) const { return fTable.HasKey(key); }
    bool     HasColumn(const std::string& col) const { return fTable.HasColumn(col);}
    const Table::Columns &GetColumns() const { return fTable.GetColumns();}
    const Table::SortedColumns& GetSortedColumns() const { return fTable.sorted_cols;}
    const Table::Keys &GetKeys() const { return fTable.GetKeys();}

    int64_t     GetInt(const std::string &key) const { return fTable.Get<int64_t>(key); }
    uint64_t    GetUInt(const std::string &key) const { return fTable.Get<uint64_t>(key); }
    double      GetFloat(const std::string &key) const { return fTable.Get<double>(key); }
    std::string GetStr(const std::string &key) const { return fTable.Get<std::string>(key); }

    size_t GetN(const std::string &key) const
    {
        return fTable.GetN(key);
    }

    size_t GetRow() const { return fRow==(size_t)-1 ? 0 : fRow; }

    operator bool() const { return fTable && fTable.offset!=0; }

    void PrintKeys(bool all_keys=false) const { fTable.PrintKeys(all_keys); }
    void PrintColumns() const { fTable.PrintColumns(); }

    bool IsHeaderOk() const { return fTable.datasum<0?false:(fChkHeader+Checksum(fTable.datasum)).valid(); }
    virtual bool IsFileOk() const { return (fChkHeader+fChkData).valid(); }

    bool IsCompressedFITS() const { return fTable.is_compressed;}

    virtual size_t GetNumRows() const
    {
        return fTable.Get<size_t>("NAXIS2");
    }

    virtual size_t GetBytesPerRow() const
    {
        return fTable.Get<size_t>("NAXIS1");
    }

    const std::vector<std::string> &GetTables() const
    {
        return fListOfTables;
    }
};

#endif
