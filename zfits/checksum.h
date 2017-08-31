#ifndef MARS_checksum
#define MARS_checksum

#ifndef __CINT__
#include <arpa/inet.h>
#endif

#include <stdint.h>

class Checksum
{
public:
    uint64_t buffer;

    void reset() { buffer = 0; }
    Checksum() : buffer(0) { }
    Checksum(const Checksum &sum) : buffer(sum.buffer) { }
    Checksum(uint64_t v) : buffer(((v>>16)&0xffff) | ((v&0xffff)<<32)) { }

    uint32_t val() const { return (((buffer&0xffff)<<16) | ((buffer>>32)&0xffff)); }

    bool valid() const { return buffer==0xffff0000ffff; }

    void HandleCarryBits()
    {
        while (1)
        {
            const uint64_t carry = ((buffer>>48)&0xffff) | ((buffer&0xffff0000)<<16);
            if (!carry)
                break;

            buffer = (buffer&0xffff0000ffff) + carry;
        }
    }

    Checksum &operator+=(const Checksum &sum)
    {
        buffer += sum.buffer;
        HandleCarryBits();
        return *this;
    }

    Checksum operator+(Checksum sum) const
    {
        return (sum += *this);
    }


    bool add(const char *buf, size_t len, bool big_endian = true)
    {
        // Avoid overflows in carry bits
        if (len>262140) // 2^18-4
        {
            add(buf, 262140);
            return add(buf+262140, len-262140);
        }

        if (len%4>0)
        {
            std::ostringstream sout;
            sout << "Length " << len << " not dividable by 4";

#ifdef __EXCEPTIONS
            throw std::runtime_error(sout.str());
#else
            gLog << ___err___ << "ERROR - " << sout.str() << std::endl;
            return false;
#endif
        }

        const uint16_t *sbuf = reinterpret_cast<const uint16_t *>(buf);

        uint32_t *hilo  = reinterpret_cast<uint32_t*>(&buffer);


        const uint16_t *end = sbuf + len/2;

        if (big_endian)
            addLoopSwapping(sbuf, end, hilo);
        else
            addLoop(sbuf, end, hilo);
        /*const uint16_t *end = sbuf + len/2;
        while (1)
        {
            if (sbuf==end)
                break;

            hilo[0] += ntohs(*sbuf++);

            if (sbuf==end)
                break;

            hilo[1] += ntohs(*sbuf++);
        }*/

        HandleCarryBits();

        return true;
    }

    void addLoopSwapping(const uint16_t *sbuf, const uint16_t *end, uint32_t* hilo)
    {
        /*
        for (size_t i = 0; i < len/2; i++)
        {
            //swap the bytes of the 32 bits value. but...
            //the hi and lo values are stored in fits-like order. do not swap them
            hilo[i%2] += ntohs(sbuf[i]); //(sbuf[i]&0xff00)>>8 | (sbuf[i]&0x00ff)<<8;
        }*/

        // This is about as twice as fast as the loop above
        // ntohs is CPU optimized, i%2 doesn't need to be computed
        while (1)
        {
            if (sbuf==end)
                break;

            hilo[0] += ntohs(*sbuf++);

            if (sbuf==end)
                break;

            hilo[1] += ntohs(*sbuf++);
        }
    }

    void addLoop(const uint16_t *sbuf, const uint16_t *end, uint32_t* hilo)
    {
        while (1)
        {
            if (sbuf==end)
                break;

            hilo[0] += ntohs(*sbuf++);

            if (sbuf==end)
                break;

            hilo[1] += ntohs(*sbuf++);
        }
    }

    bool add(const std::vector<char> &v, bool big_endian = true)
    {
        return add(v.data(), v.size(), big_endian);
    }

    std::string str(bool complm=true) const
    {
        std::string rc(16,0);

        const uint8_t exclude[13] =
        {
            0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f, 0x40,
            0x5b, 0x5c, 0x5d, 0x5e, 0x5f, 0x60
        };

        const uint32_t value = complm ? ~val() : val();   // complement each bit of the value

        for (int ii = 0; ii < 4; ii++)
        {
            uint8_t byte = (value >> (24 - (8 * ii)));

            const uint8_t quotient  = byte / 4 + '0';
            const uint8_t remainder = byte % 4;

            uint32_t ch[4] = { uint32_t(quotient+remainder), quotient, quotient, quotient };

            // avoid ASCII  punctuation
            while (1)
            {
                bool check = false;
                for (int kk = 0; kk < 13; kk++)
                {
                    for (int jj = 0; jj < 4; jj += 2)
                    {
                        if (ch[jj] != exclude[kk] &&  ch[jj+1] != exclude[kk])
                            continue;

                        ch[jj]++;
                        ch[jj+1]--;
                        check=true;
                    }
                }

                if (!check)
                    break;
            }

            for (int jj = 0; jj < 4; jj++)        // assign the bytes
                rc[4*jj+ii] = ch[jj];
        }

        const char lastChar = rc[15];
        for (int i=15;i>0;i--)
            rc[i] = rc[i-1];
        rc[0] = lastChar;
        return rc;
/*
        uint8_t *p = reinterpret_cast<uint8_t*>(&value);
        //swap the bytes of the value
        uint8_t temp;
        temp = p[0];
        p[0] = p[3];
        p[3] = temp;
        temp = p[1];
        p[1] = p[2];
        p[2] = temp;

        for (int i=0; i<4; i++)
        {
            rc[i+ 0] = '0' + p[i]/4 + p[i]%4;
            rc[i+ 4] = '0' + p[i]/4;
            rc[i+ 8] = '0' + p[i]/4;
            rc[i+12] = '0' + p[i]/4;
        }

        while(1)
        {
            bool ok  = true;
            for (int i=0; i<16-4; i++)
            {
                for (int j=0; j<13; j++)
                    if (rc[i]==exclude[j] || rc[(i+4)%16]==exclude[j])
                    {
                        rc[i]++;
                        rc[(i+4)%16]--;
                        ok = false;
                    }
            }

            if (ok)
                break;
        }


        return rc;
 */
   }
};

#endif
