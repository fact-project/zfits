#ifndef FACT_huffman
#define FACT_huffman

#include <string.h>
#include <stdint.h>
#include <stdexcept>

#include <set>
#include <string>
#include <vector>
#include <stdio.h>
#define MAX_SYMBOLS (1<<16)

// ================================================================

namespace Huffman
{
    static unsigned long numbytes_from_numbits(unsigned long numbits)
    {
        return numbits / 8 + (numbits % 8 ? 1 : 0);
    }

    struct TreeNode
    {
        TreeNode *parent;
        union
        {
            struct
            {
                TreeNode *zero, *one;
            };
            uint16_t symbol;
        };

        size_t count;
        bool isLeaf;

        TreeNode(uint16_t sym, size_t cnt=0) : parent(0), isLeaf(true)
        {
            symbol = sym;
            count  = cnt;
        }

        TreeNode(TreeNode *n0=0, TreeNode *n1=0) : parent(0), isLeaf(false)
        {
            count = n0 && n1 ?  n0->count + n1->count : 0;
            zero  = n0 && n1 ? (n0->count > n1->count ? n0 : n1) : NULL;
            one   = n0 && n1 ? (n0->count > n1->count ? n1 : n0) : NULL;

            if (n0)
                n0->parent = this;

            if (n1)
                n1->parent = this;
        }

        ~TreeNode()
        {
            if (isLeaf)
                return;

            if (zero)
                delete zero;
            if (one)
                delete one;
        }

        bool operator() (const TreeNode *hn1, const TreeNode *hn2) const
        {
            return hn1->count < hn2->count;
        }
    };


    struct Encoder
    {
        struct Code
        {
            size_t bits;
            uint8_t numbits;

            Code() : numbits(0) { }
        };

        size_t count;
        Code lut[1<<16];

        bool CreateEncoder(const TreeNode *n, size_t bits=0, uint8_t nbits=0)
        {
            if (n->isLeaf)
            {
                if (nbits>sizeof(size_t)*8)
                    throw std::runtime_error("Too many different symbols - this should not happen!");
                lut[n->symbol].bits    = bits;
                lut[n->symbol].numbits = nbits==0 ? 1 : nbits;
                count++;
                return true;
            }

            return
                CreateEncoder(n->zero, bits,              nbits+1) &&
                CreateEncoder(n->one,  bits | (1<<nbits), nbits+1);

        }

        void WriteCodeTable(std::string &out) const
        {
            out.append((char*)&count, sizeof(size_t));

            for (uint32_t i=0; i<MAX_SYMBOLS; i++)
            {
                const Code &n = lut[i];
                if (n.numbits==0)
                    continue;

                // Write the 2 byte symbol.
                out.append((char*)&i, sizeof(uint16_t));
                if (count==1)
                    return;

                // Write the 1 byte code bit length.
                out.append((char*)&n.numbits, sizeof(uint8_t));

                // Write the code bytes.
                uint32_t numbytes = numbytes_from_numbits(n.numbits);
                out.append((char*)&n.bits, numbytes);
            }
        }

        void Encode(std::string &out, const uint16_t *bufin, size_t bufinlen) const
        {
            if (count==1)
                return;

            uint8_t curbyte = 0;
            uint8_t curbit  = 0;

            for (uint32_t i=0; i<bufinlen; ++i)
            {
                const uint16_t &symbol = bufin[i];

                const Code *code = lut+symbol;

                uint8_t nbits = code->numbits;
                const uint8_t *bits = (uint8_t*)&code->bits;

                while (nbits>0)
                {
                    // Number of bits available in the current byte
                    const uint8_t free_bits = 8 - curbit;

                    // Write bits to current byte
                    curbyte |= *bits<<curbit;

                    // If the byte has been filled, put it into the output buffer
                    // If the bits exceed the current byte step to the next byte
                    // and fill it properly
                    if (nbits>=free_bits)
                    {
                        out += curbyte;
                        curbyte = *bits>>free_bits;

                        bits++;
                    }

                    // Adapt the number of available bits, the number of consumed bits
                    // and the bit-pointer accordingly
                    const uint8_t consumed = nbits>8 ? 8 : nbits;
                    nbits  -= consumed;
                    curbit += consumed;
                    curbit %= 8;
                }
            }

            // If the buffer-byte is half-full, also add it to the output buffer
            if (curbit>0)
                out += curbyte;
        }

        Encoder(const uint16_t *bufin, size_t bufinlen) : count(0)
        {
            uint64_t counts[MAX_SYMBOLS];
            memset(counts, 0, sizeof(uint64_t)*MAX_SYMBOLS);

            // Count occurances
            for (const uint16_t *p=bufin; p<bufin+bufinlen; p++)
                counts[*p]++;

            // Copy all occuring symbols into a sorted list
            std::multiset<TreeNode*, TreeNode> set;
            for (int i=0; i<MAX_SYMBOLS; i++)
                if (counts[i])
                    set.insert(new TreeNode(i, counts[i]));

            // Create the tree bottom-up
            while (set.size()>1)
            {
                auto it = set.begin();

                auto it1 = it++;
                auto it2 = it;

                TreeNode *nn = new TreeNode(*it1, *it2);

                set.erase(it1, ++it2);

                set.insert(nn);
            }

            // get the root of the tree
            const TreeNode *root = *set.begin();

            CreateEncoder(root);

            // This will delete the whole tree
            delete root;
        }

    };



    struct Decoder
    {
        uint16_t symbol;
        uint8_t nbits;
        bool isLeaf;

        Decoder *lut;

        Decoder() : isLeaf(false), lut(NULL)
        {
        }

        ~Decoder()
        {
            if (lut)
                delete [] lut;
        }

        void Set(uint16_t sym, uint8_t n=0, size_t bits=0)
        {
            if (!lut)
                lut = new Decoder[256];

            if (n>8)
            {
                lut[bits&0xff].Set(sym, n-8, bits>>8);
                return;
            }

            const int nn = 1<<(8-n);

            for (int i=0; i<nn; i++)
            {
                const uint8_t key = bits | (i<<n);

                lut[key].symbol = sym;
                lut[key].isLeaf = true;
                lut[key].nbits  = n;
            }
        }

        void Build(const TreeNode &p, uint64_t bits=0, uint8_t n=0)
        {
            if (p.isLeaf)
            {
                Set(p.symbol, n, bits);
                return;
            }

            Build(*p.zero, bits,          n+1);
            Build(*p.one,  bits | (1<<n), n+1);
        }

        Decoder(const TreeNode &p) : symbol(0), isLeaf(false), lut(NULL)
        {
            Build(p);
        }

        const uint8_t *Decode(const uint8_t *in_ptr, const uint8_t *in_end,
                              uint16_t *out_ptr, const uint16_t *out_end) const
        {
            Decoder const *p = this;

            if (in_ptr==in_end)
            {
                while (out_ptr < out_end)
                    *out_ptr++ = p->lut->symbol;
                return in_ptr;
            }

            uint8_t curbit = 0;
            while (in_ptr<in_end && out_ptr<out_end)
            {
                const uint16_t *two = (uint16_t*)in_ptr;

                const uint8_t curbyte = (*two >> curbit);

                if (!p->lut)
                    throw std::runtime_error("Unknown bitcode in stream!");

                p = p->lut + curbyte;
                if (!p->isLeaf)
                {
                    in_ptr++;
                    continue;
                }

                *out_ptr++ = p->symbol;
                curbit += p->nbits;

                p = this;

                if (curbit>=8)
                {
                    curbit %= 8;
                    in_ptr++;
                }

            }

            return curbit ? in_ptr+1 : in_ptr;
        }

        Decoder(const uint8_t* bufin, int64_t &pindex) : isLeaf(false), lut(NULL)
        {
            // FIXME: Sanity check for size missing....

            // Read the number of entries.
            size_t count=0;
            memcpy(&count, bufin + pindex, sizeof(count));
            pindex += sizeof(count);
            // Read the entries.
            for (size_t i=0; i<count; i++)
            {
                uint16_t sym;
                memcpy(&sym, bufin + pindex, sizeof(uint16_t));
                pindex += sizeof(uint16_t);

                if (count==1)
                {
                    Set(sym);
                    break;
                }

                uint8_t numbits;
                memcpy(&numbits, bufin + pindex, sizeof(uint8_t));
                pindex += sizeof(uint8_t);

                const uint8_t numbytes = numbytes_from_numbits(numbits);

                if (numbytes>sizeof(size_t))
                    throw std::runtime_error("Number of bytes for a single symbol exceeds maximum.");
                size_t bits=0;
                memcpy(&bits, bufin+pindex, numbytes);
                pindex += numbytes;

                Set(sym, numbits, bits);
            }
        }
    };

    inline bool Encode(std::string &bufout, const uint16_t *bufin, size_t bufinlen)
    {
        const Encoder encoder(bufin, bufinlen);

        bufout.append((char*)&bufinlen, sizeof(size_t));
        encoder.WriteCodeTable(bufout);
        encoder.Encode(bufout, bufin, bufinlen);

        return true;
    }

    inline int64_t Decode(const uint8_t *bufin,
                          size_t         bufinlen,
                          std::vector<uint16_t> &pbufout)
    {
        int64_t i = 0;

        // Read the number of data bytes this encoding represents.
        size_t data_count = 0;
        memcpy(&data_count, bufin, sizeof(size_t));
        i += sizeof(size_t);


        pbufout.resize(data_count);

        const Decoder decoder(bufin, i);

        const uint8_t *in_ptr =
            decoder.Decode(bufin+i, bufin+bufinlen,
                           pbufout.data(), pbufout.data()+data_count);

        return in_ptr-bufin;
    }

};

#endif
