#ifndef MARS_DrsCalib
#define MARS_DrsCalib

#include <math.h>   // fabs
#include <errno.h>  // errno
#include <algorithm>  // sort

class DrsCalibrate
{
protected:
    int64_t fNumEntries;

    size_t fNumSamples;
    size_t fNumChannels;

    std::vector<int64_t> fSum;
    std::vector<int64_t> fSum2;

public:
    DrsCalibrate() : fNumEntries(0), fNumSamples(0), fNumChannels(0)
    {
        fSum.reserve(1024*1440);
        fSum2.reserve(1024*1440);
    }

    void Reset()
    {
        fNumEntries  = 0;
        fNumSamples  = 0;
        fNumChannels = 0;

        fSum.clear();
        fSum2.clear();
    }

    void InitSize(uint16_t channels, uint16_t samples)
    {
        fNumChannels = channels;
        fNumSamples  = samples;

        fSum.assign(samples*channels, 0);
        fSum2.assign(samples*channels, 0);
    }

    void AddRel(const int16_t *val, const int16_t *start)
    {
        /*
        for (size_t ch=0; ch<fNumChannels; ch++)
        {
            const int16_t &spos = start[ch];
            if (spos<0)
                continue;

            const size_t pos = ch*1024;
            for (size_t i=0; i<1024; i++)
            {
                // Value is relative to trigger
                // Abs is corresponding index relative to DRS pipeline
                const size_t rel = pos +  i;
                const size_t abs = pos + (spos+i)%1024;

                const int64_t v = val[rel];

                fSum[abs]  += v;
                fSum2[abs] += v*v;
            }
        }*/

        // This version is 2.5 times faster because the compilers optimization
        // is not biased by the evaluation of %1024
        for (size_t ch=0; ch<fNumChannels; ch++)
        {
            const int16_t &spos = start[ch];
            if (spos<0)
                continue;

            const size_t pos = ch*1024;

            const int16_t *beg_val  = val          + pos;
            int64_t       *beg_sum  = fSum.data()  + pos;
            int64_t       *beg_sum2 = fSum2.data() + pos;

            const int16_t *pval  = beg_val;          // val[rel]
            int64_t       *psum  = beg_sum  + spos;  // fSum[abs]
            int64_t       *psum2 = beg_sum2 + spos;  // fSum2[abs]

            while (psum<beg_sum+1024)
            {
                const int64_t v = *pval++;

                *psum++  += v;
                *psum2++ += v*v;
            }

            psum  = beg_sum;
            psum2 = beg_sum2;

            while (pval<beg_val+1024)
            {
                const int64_t v = *pval++;

                *psum++  += v;
                *psum2++ += v*v;
            }
        }

        fNumEntries++;
    }

    void AddRel(const int16_t *val,    const int16_t *start,
                const int32_t *offset, const int64_t scale)
    {
        /*
        for (size_t ch=0; ch<fNumChannels; ch++)
        {
            const int16_t spos = start[ch];
            if (spos<0)
                continue;

            const size_t pos = ch*1024;

            for (size_t i=0; i<fNumSamples; i++)
            {
                // Value is relative to trigger
                // Offset is relative to DRS pipeline
                // Abs is corresponding index relative to DRS pipeline
                const size_t rel = pos +  i;
                const size_t abs = pos + (spos+i)%1024;

                const int64_t v = int64_t(val[rel])*scale-offset[abs];

                fSum[abs]  += v;
                fSum2[abs] += v*v;
            }
        }*/

        // This version is 2.5 times faster because the compilers optimization
        // is not biased by the evaluation of %1024
        for (size_t ch=0; ch<fNumChannels; ch++)
        {
            const int16_t &spos = start[ch];
            if (spos<0)
                continue;

            const size_t pos = ch*1024;

            const int16_t *beg_val    = val          + pos;
            const int32_t *beg_offset = offset       + pos;
            int64_t       *beg_sum    = fSum.data()  + pos;
            int64_t       *beg_sum2   = fSum2.data() + pos;


            const int16_t *pval    = beg_val;            // val[rel]
            const int32_t *poffset = beg_offset + spos;  // offset[abs]
            int64_t       *psum    = beg_sum    + spos;  // fSum[abs]
            int64_t       *psum2   = beg_sum2   + spos;  // fSum2[abs]

            while (psum<beg_sum+1024)
            {
                const int64_t v = int64_t(*pval++)*scale - *poffset++;

                *psum++  += v;
                *psum2++ += v*v;
            }

            psum    = beg_sum;
            psum2   = beg_sum2;
            poffset = beg_offset;

            while (pval<beg_val+1024)
            {
                const int64_t v = int64_t(*pval++)*scale - *poffset++;

                *psum++  += v;
                *psum2++ += v*v;
            }
        }

        fNumEntries++;
    }

    void AddAbs(const int16_t *val,    const int16_t *start,
                const int32_t *offset, const int64_t scale)
    {
        /*
        // 1440 without tm, 1600 with tm
        for (size_t ch=0; ch<fNumChannels; ch++)
        {
            const int16_t spos = start[ch];
            if (spos<0)
                continue;

            const size_t pos = ch*fNumSamples;
            const size_t drs = ch>1439 ? ((ch-1440)*9+8)*1024 : ch*1024;

            for (size_t i=0; i<fNumSamples; i++)
            {
                // Value is relative to trigger
                // Offset is relative to DRS pipeline
                // Abs is corresponding index relative to DRS pipeline
                const size_t rel = pos +  i;
                const size_t abs = drs + (spos+i)%1024;

                const int64_t v = int64_t(val[rel])*scale-offset[abs];

                fSum[rel]  += v;
                fSum2[rel] += v*v;
            }
            }*/

        // This version is 1.5 times faster because the compilers optimization
        // is not biased by the evaluation of %1024
        for (size_t ch=0; ch<fNumChannels; ch++)
        {
            const int16_t &spos = start[ch];
            if (spos<0)
                continue;

            const size_t pos = ch*fNumSamples;

            const int32_t *beg_offset = offset       + ch*1024;
            const int16_t *beg_val    = val          + pos;
            int64_t *beg_sum          = fSum.data()  + pos;
            int64_t *beg_sum2         = fSum2.data() + pos;


            const int16_t *pval    = beg_val;             // val[rel]
            const int32_t *poffset = beg_offset + spos;   // offset[abs]
            int64_t *psum          = beg_sum;             // fSum[rel]
            int64_t *psum2         = beg_sum2;            // fSum2[rel]

            if (spos+fNumSamples>1024)
            {
                while (poffset<beg_offset+1024)
                {
                    const int64_t v = int64_t(*pval++)*scale - *poffset++;

                    *psum++  += v;
                    *psum2++ += v*v;
                }

                poffset = beg_offset;
            }

            while (psum<beg_sum+fNumSamples)
            {
                const int64_t v = int64_t(*pval++)*scale - *poffset++;

                *psum++  += v;
                *psum2++ += v*v;
            }
        }

        fNumEntries++;
    }


    static void ApplyCh(
        float *vec,
        const int16_t *val,
        int16_t start,
        uint32_t roi,
        const int32_t *offset,
        const int64_t scaleabs,
        const int64_t *gain,
        const int64_t scalegain
    )
    {
        if (start<0)
        {
            memset(vec, 0, roi);
            return;
        }
        /*
        for (size_t i=0; i<roi; i++)
        {
            // Value is relative to trigger
            // Offset is relative to DRS pipeline
            // Abs is corresponding index relative to DRS pipeline
            const size_t abs = (start+i)%1024;

            const int64_t v =
                + int64_t(val[i])*scaleabs-offset[abs]
                ;

            const int64_t div = gain[abs];
            vec[i] = div==0 ? 0 : double(v)*scalegain/div;
        }
        */

        // This version is faster because the compilers optimization
        // is not biased by the evaluation of %1024
        // (Here we are dominated by numerics... improvement ~10%)
        const int32_t *poffset = offset + start; // offset[abs]
        const int64_t *pgain   = gain   + start; // gain[abs]
        const int16_t *pval    = val;            // val[rel]
        float         *pvec    = vec;            // vec[rel]

        if (start+roi>1024)
        {
            while (poffset<offset+1024)
            {
                const int64_t v =
                    + int64_t(*pval++)*scaleabs - *poffset++
                    ;

                *pvec++ = *pgain==0 ? 0 : double(v)*scalegain / *pgain;

                pgain++;
            }

            poffset = offset;
            pgain   = gain;
        }

        while (pvec<vec+roi)
        {
            const int64_t v =
                + int64_t(*pval++)*scaleabs - *poffset++
                ;

            *pvec++ = *pgain==0 ? 0 : double(v)*scalegain / *pgain;

            pgain++;
        }
    }

    static void ApplyCh(float *vec, const int16_t *val, int16_t start, uint32_t roi,
                        const int32_t *offset, const int64_t scaleabs,
                        const int64_t *gain,   const int64_t scalegain,
                        const int64_t *trgoff, const int64_t scalerel)
    {
        if (start<0)
        {
            memset(vec, 0, roi);
            return;
        }
        /*
        for (size_t i=0; i<roi; i++)
        {
            // Value is relative to trigger
            // Offset is relative to DRS pipeline
            // Abs is corresponding index relative to DRS pipeline
            const size_t abs = (start+i)%1024;

            const int64_t v =
                + (int64_t(val[i])*scaleabs-offset[abs])*scalerel
                - trgoff[i]
                ;

            const int64_t div = gain[abs]*scalerel;
            vec[i] = div==0 ? 0 : double(v)*scalegain/div;
        }
        */
        // (Here we are dominated by numerics... improvement ~10%)
        const int32_t *poffset = offset + start; // offset[abs]
        const int64_t *pgain   = gain   + start; // gain[abs]
        const int16_t *pval    = val;            // val[rel]
        const int64_t *ptrgoff = trgoff;         // trgoff[rel]
        float         *pvec    = vec;            // vec[rel]

        if (start+roi>1024)
        {
            while (poffset<offset+1024)
            {
                const int64_t v =
                    + (int64_t(*pval++)*scaleabs - *poffset++)*scalerel
                    - *ptrgoff++;
                ;

                const int64_t div = *pgain * scalerel;
                *pvec++ = div==0 ? 0 : double(v)*scalegain / div;

                pgain++;
            }

            poffset = offset;
            pgain   = gain;
        }

        while (pvec<vec+roi)
        {
            const int64_t v =
                + (int64_t(*pval++)*scaleabs - *poffset++)*scalerel
                - *ptrgoff++;
            ;

            const int64_t div = *pgain * scalerel;
            *pvec++ = div==0 ? 0 : double(v)*scalegain / div;

            pgain++;
        }
    }

    static double FindStep(const size_t ch0, const float *vec, int16_t roi, const int16_t pos, const uint16_t *map=NULL)
    {
        // We have about 1% of all cases which are not ahndled here,
        // because the baseline jumps up just before the readout window
        // and down just after it. In this cases we could determine the jump
        // from the board time difference or throw the event away.
        if (pos==0 || pos>=roi)
            return 0;

        double step = 0; // before
        double rms  = 0; // before
        int    cnt  = 0;

        // Exclude TM channel
        for (int p=0; p<8; p++)
        {
            const size_t hw = ch0+p;
            const size_t sw = (map?map[hw]:hw)*roi + pos;

            const double diff = vec[sw]-vec[sw-1];

            step += diff;
            rms  += (vec[sw]-vec[sw-1])*(vec[sw]-vec[sw-1]);

            cnt++;
        }

        return cnt==0 ? 0 : step/cnt;
    }

    static void SubtractStep(const size_t ch0, const double avg, float *vec, int16_t roi, int32_t pos, const uint16_t *map=NULL)
    {
        if (pos==0 || pos>=roi)
            return;

        const int begin = avg>0 ? pos :   0;
        const int end   = avg>0 ? roi : pos;

        const double sub = fabs(avg);

        for (int p=0; p<9; p++)
        {
            for (int j=begin; j<end; j++)
            {
                const size_t hw = ch0+p;
                const size_t sw = (map?map[hw]:hw)*roi + j;

                vec[sw] -= sub;
            }
        }
    }

    struct Step
    {
        Step() : avg(0), rms(0), pos(0), cnt(0) { }
        double   avg;
        double   rms;
        double   pos;
        uint16_t cnt;

        static bool sort(const Step &s, const Step &r) { return s.avg<r.avg; }
    };

    static Step AverageSteps(const std::vector<Step>::iterator beg, const std::vector<Step>::iterator end)
    {
        Step rc;
        for (auto it=beg; it!=end; it++)
        {
            rc.pos += it->pos;
            rc.avg += it->avg;
            rc.rms += it->avg*it->avg;
        }

        rc.cnt = end-beg;

        rc.pos /= rc.cnt;
        rc.avg /= rc.cnt;
        rc.rms /= rc.cnt;
        rc.rms -= rc.avg*rc.avg;
        rc.rms  = rc.rms<0 ? 0 : sqrt(rc.rms);

        return rc;
    }


    static Step CorrectStep(float *vec, uint16_t nch, uint16_t roi,
                            const int16_t *prev, const int16_t *start,
                            const int16_t offset, const uint16_t *map=NULL)
    {

        std::vector<Step> list;
        list.reserve(nch);

        // Fill steps into array
        // Exclude broken pixels?
        // Remove maximum and minimum patches (4max and 4min)?
        for (size_t ch=0; ch<nch; ch += 9)
        {
            if (prev[ch]<0 || start[ch]<0)
                continue;

            const int16_t dist = (prev[ch]-start[ch]+1024+offset)%1024;
            const double  step = FindStep(ch, vec, roi, dist, map);
            if (step==0)
                continue;

            Step rc;
            rc.pos = dist;
            rc.avg = step;
            list.push_back(rc);
        }

        if (list.empty())
            return Step();

        Step rc = AverageSteps(list.begin(), list.begin()+list.size());;

        if (rc.avg==0)
            return Step();

        // std::cout << "   A0=" << rc.avg << "    rms=" << rc.rms << std::endl;
        if (rc.rms>5)
        {
            sort(list.begin(), list.end(), Step::sort);

            //for (auto it=list.begin(); it!=list.end(); it++)
            //    std::cout << "     " << it->avg << std::endl;

            const size_t skip = list.size()/10;
            rc = AverageSteps(list.begin()+skip, list.begin()+list.size()-skip);

            // std::cout << "   A1=" << rc.avg << "    rms=" << rc.rms << std::endl;
        }

        for (size_t ch=0; ch<nch; ch += 9)
        {
            const int16_t dist = (prev[ch]-start[ch]+1024+offset)%1024;
            SubtractStep(ch, rc.avg, vec, roi, dist, map);
        }

        return rc;
    }

    static void RemoveSpikes(float *p, uint32_t roi)
    {
        if (roi<4)
            return;

        for (size_t i=1; i<roi-2; i++)
        {
            if (p[i]-p[i-1]>25 && p[i]-p[i+1]>25)
            {
                p[i] = (p[i-1]+p[i+1])/2;
            }

            if (p[i]-p[i-1]>22 && fabs(p[i]-p[i+1])<4 && p[i+1]-p[i+2]>22)
            {
                p[i] = (p[i-1]+p[i+2])/2;
                p[i+1] = p[i];
            }
        }
    }

    static void RemoveSpikes2(float *p, uint32_t roi)//from Werner
    {
        if (roi<4)
            return;

        std::vector<float> Ameas(p, p+roi);

        std::vector<float> diff(roi);
        for (size_t i=1; i<roi-1; i++)
            diff[i] = (p[i-1] + p[i+1])/2 - p[i];

        //std::vector<float> N1mean(roi);
        //for (size_t i=1; i<roi-1; i++)
        //    N1mean[i] = (p[i-1] + p[i+1])/2;

        const float fract = 0.8;

        for (size_t i=0; i<roi-3; i++)
        {
            if (diff[i]<5)
                continue;

            if (Ameas[i+2] - (Ameas[i] + Ameas[i+3])/2 > 10)
            {
                p[i+1]=   (Ameas[i+3] - Ameas[i])/3 + Ameas[i];
                p[i+2]= 2*(Ameas[i+3] - Ameas[i])/3 + Ameas[i];

                i += 3;

                continue;
            }

            if ( (diff[i+1]<-diff[i]*fract*2) && (diff[i+2]>10) )
            {
                p[i+1]    = (Ameas[i]+Ameas[i+2])/2;
                diff[i+2] = (p[i+1] + Ameas[i+3])/2 - Ameas[i+2];

                i += 2;
            }

            // const float x = Ameas[i] - N1mean[i];
            // if (x > -5.)
            //     continue;

            // if (Ameas[i+2] - (Ameas[i] + Ameas[i+3])/2. > 10.)
            // {
            //     p[i+1]=   (Ameas[i+3] - Ameas[i])/3 + Ameas[i];
            //     p[i+2]= 2*(Ameas[i+3] - Ameas[i])/3 + Ameas[i];
            //     i += 3;
            //     continue;
            // }

            // const float xp  = Ameas[i+1] - N1mean[i+1];
            // const float xpp = Ameas[i+2] - N1mean[i+2];

            // if ( (xp > -2.*x*fract) && (xpp < -10.) )
            // {
            //     p[i+1] = N1mean[i+1];
            //     N1mean[i+2] = Ameas[i+1] - Ameas[i+3]/2;
            //
            //     i += 2;
            // }
        }
    }

    static void RemoveSpikes3(float *vec, uint32_t roi)//from Werner
    {
        if (roi<4)
            return;

        const float SingleCandidateTHR = -10.;
        const float DoubleCandidateTHR =  -5.;

        const std::vector<float> src(vec, vec+roi);

        std::vector<float> diff(roi);
        for (size_t i=1; i<roi-1; i++)
            diff[i] = src[i] - (src[i-1] + src[i+1])/2;

        // find the spike and replace it by mean value of neighbors
        for (unsigned int i=1; i<roi-3; i++)
        {
            // Speed up (no leading edge)
            if (diff[i]>=DoubleCandidateTHR)
                continue;

            //bool checkDouble = false;

            // a single spike candidate
            if (diff[i]<SingleCandidateTHR)
            {
                // check consistency with a single channel spike
                if (diff[i+1] > -1.6*diff[i])
                {
                    vec[i+1] = (src[i] + src[i+2]) / 2;

                    i += 2;

                    /*** NEW ***/
                    continue;
                    /*** NEW ***/
                }
                /*
                else
                {
                    // do nothing - not really a single spike,
                    // but check if it is a double
                    checkDouble = true;
                }*/
            }

            // a double spike candidate
            //if (diff[i]>DoubleCandidateTHR || checkDouble == 1)
            {
                // check the consistency with a double spike
                if ((diff[i+1] > -DoubleCandidateTHR) &&
                    (diff[i+2] > -DoubleCandidateTHR))
                {
                    vec[i+1] =   (src[i+3] - src[i])/3 + src[i];
                    vec[i+2] = 2*(src[i+3] - src[i])/3 + src[i];

                    //vec[i]   = (src[i-1] + src[i+2]) / 2.;
                    //vec[i+1] = (src[i-1] + src[i+2]) / 2.;

                    //do not care about the next sample it was the spike
                    i += 3;
                }
            }
        }
    }

    static void RemoveSpikes4(float *ptr, uint32_t roi)
    {
        if (roi<7)
            return;

        for (uint32_t i=0; i<roi-6; i++)
        {
            double d10, d21, d32, d43, d54;

            // ============================================
            d43 = ptr[i+4]-ptr[i+3];
            d54 = ptr[i+5]-ptr[i+4];

            if ((d43>35 && -d54>35) || (d43<-35 && -d54<-35))
            {
                ptr[i+4] = (ptr[i+3]+ptr[i+5])/2;
            }

            // ============================================
            d32 = ptr[i+3]-ptr[i+2];
            d54 = ptr[i+5]-ptr[i+4];

            if ((d32>9   && -d54>13  && d32-d54>31)/* || (d32<-13 && -d54<-13 && d32+d54<-63)*/)
            {
                double avg0 = (ptr[i+2]+ptr[i+5])/2;
                double avg1 = (ptr[i+3]+ptr[i+4])/2;

                ptr[i+3] = ptr[i+3] - avg1+avg0;
                ptr[i+4] = ptr[i+4] - avg1+avg0;
            }

            // ============================================
            d21 = ptr[i+2]-ptr[i+1];
            d54 = ptr[i+5]-ptr[i+4];

            if (d21>15 && -d54>17)
            {
                double avg0 = (ptr[i+1]+ptr[i+5])/2;
                double avg1 = (ptr[i+2]+ptr[i+3]+ptr[i+4])/3;

                ptr[i+2] = ptr[i+2] - avg1+avg0;
                ptr[i+3] = ptr[i+3] - avg1+avg0;
                ptr[i+4] = ptr[i+4] - avg1+avg0;
            }

            // ============================================
            d10 = ptr[i+1]-ptr[i];
            d54 = ptr[i+5]-ptr[i+4];

            if (d10>18 && -d54>20)
            {
                double avg0 = (ptr[i]+ptr[i+5])/2;
                double avg1 = (ptr[i+1]+ptr[i+2]+ptr[i+3]+ptr[i+4])/4;

                ptr[i+1] = ptr[i+1] - avg1+avg0;
                ptr[i+2] = ptr[i+2] - avg1+avg0;
                ptr[i+3] = ptr[i+3] - avg1+avg0;
                ptr[i+4] = ptr[i+4] - avg1+avg0;
            }
        }
    }

    static void SlidingAverage(float *const vec, const uint32_t roi, const uint16_t w)
    {
        if (w==0 || w>roi)
            return;

        for (float *pix=vec; pix<vec+1440*roi; pix += roi)
        {
            for (float *ptr=pix; ptr<pix+roi-w; ptr++)
            {
                for (float *p=ptr+1; p<ptr+w; p++)
                    *ptr += *p;
                *ptr /= w;
            }
        }
    }

    std::pair<std::vector<double>,std::vector<double> > GetSampleStats() const
    {
        if (fNumEntries==0)
            return make_pair(std::vector<double>(),std::vector<double>());

        std::vector<double> mean(fSum.size());
        std::vector<double> error(fSum.size());

        std::vector<int64_t>::const_iterator it = fSum.begin();
        std::vector<int64_t>::const_iterator i2 = fSum2.begin();
        std::vector<double>::iterator        im = mean.begin();
        std::vector<double>::iterator        ie = error.begin();

        while (it!=fSum.end())
        {
            *im = /*cnt<fResult.size() ? fResult[cnt] :*/ double(*it)/fNumEntries;
            *ie = sqrt(double(*i2*int64_t(fNumEntries) - *it * *it))/fNumEntries;

            im++;
            ie++;
            it++;
            i2++;
        }


        /*
         valarray<double> ...

         mean /= fNumEntries;
         error = sqrt(error/fNumEntries - mean*mean);
         */

        return make_pair(mean, error);
    }

    void GetSampleStats(float *ptr, float scale) const
    {
        const size_t sz = fNumSamples*fNumChannels;

        if (fNumEntries==0)
        {
            memset(ptr, 0, sizeof(float)*sz*2);
            return;
        }

        std::vector<int64_t>::const_iterator it = fSum.begin();
        std::vector<int64_t>::const_iterator i2 = fSum2.begin();

        while (it!=fSum.end())
        {
            *ptr      = scale*double(*it)/fNumEntries;
            *(ptr+sz) = scale*sqrt(double(*i2*fNumEntries - *it * *it))/fNumEntries;

            ptr++;
            it++;
            i2++;
        }
    }

    static double GetPixelStats(float *ptr, const float *data, uint16_t roi, uint16_t begskip=0, uint16_t endskip=0)
    {
        if (roi==0)
            return -1;

        // Skip first 10 samples
        const uint beg = roi>begskip ? begskip : 0;
        const uint end = roi-beg>endskip ? roi-endskip : roi;
        const uint len = end-beg;

        double max = 0;
        double patch = 0;
        for (uint i=0; i<1440; i++)
        {
            const float *vec = data+i*roi;

            uint   pos  = beg;
            double sum  = vec[beg];
            double sum2 = vec[beg]*vec[beg];

            for (uint j=beg+1; j<end; j++)
            {
                sum  += vec[j];
                sum2 += vec[j]*vec[j];

                if (vec[j]>vec[pos])
                    pos = j;
            }
            sum  /= len;
            sum2 /= len;
            sum2 -= sum*sum;

            if (i%9!=8)
                patch += vec[pos];
            else
            {
                if (patch > max)
                    max = patch;
                patch = 0;
            }

            *(ptr+0*1440+i) = sum;
            *(ptr+1*1440+i) = sum2<0 ? 0 : sqrt(sum2);
            *(ptr+2*1440+i) = vec[pos];
            *(ptr+3*1440+i) = pos;
        }

        return max/8;
    }

    static void GetPixelMax(float *max, const float *data, uint16_t roi, int32_t first, int32_t last)
    {
        if (roi==0 || first<0 || last<0 || first>=roi || last>=roi || last<first)
            return;

        for (int i=0; i<1440; i++)
        {
            const float *beg = data+i*roi+first;
            const float *end = data+i*roi+last;

            const float *pmax = beg;

            for (const float *ptr=beg+1; ptr<=end; ptr++)
                if (*ptr>*pmax)
                    pmax = ptr;

            max[i] = *pmax;
        }
    }

    const std::vector<int64_t> &GetSum() const { return fSum; }

    int64_t GetNumEntries() const { return fNumEntries; }
};

class DrsCalibrateTime
{
public:
    int64_t fNumEntries;

    size_t fNumSamples;
    size_t fNumChannels;

    std::vector<std::pair<double, double>> fStat;

public:
    DrsCalibrateTime() : fNumEntries(0), fNumSamples(0), fNumChannels(0)
    {
        InitSize(160, 1024);
    }

    DrsCalibrateTime(const DrsCalibrateTime &p) : fNumEntries(p.fNumEntries), fNumSamples(p.fNumSamples), fNumChannels(p.fNumChannels), fStat(p.fStat)
    {
    }
    virtual ~DrsCalibrateTime()
    {
    }

    double Sum(uint32_t i) const { return fStat[i].first; }
    double W(uint32_t i) const { return fStat[i].second; }

    virtual void InitSize(uint16_t channels, uint16_t samples)
    {
        fNumChannels = channels;
        fNumSamples  = samples;

        fNumEntries  = 0;

        fStat.clear();

        fStat.resize(samples*channels);
    }

    void Reset()
    {
        for (auto it=fStat.begin(); it!=fStat.end(); it++)
        {
            it->first = 0;
            it->second = 0;
        }
    }

    void AddT(const float *val, const int16_t *start, signed char edge=0)
    {
        if (fNumSamples!=1024 || fNumChannels!=160)
            return;

        // Rising or falling edge detection has the advantage that
        // we are much less sensitive to baseline shifts

        for (size_t ch=0; ch<160; ch++)
        {
            const size_t tm = ch*9+8;

            const int16_t spos = start[tm];
            if (spos<0)
                continue;

            const size_t pos = ch*1024;

            double  p_prev =  0;
            int32_t i_prev = -1;

            for (size_t i=0; i<1024-1; i++)
            {
                const size_t rel = tm*1024 + i;

                const float &v0 = val[rel];  //-avg;
                const float &v1 = val[rel+1];//-avg;

                // If edge is positive ignore all falling edges
                if (edge>0 && v0>0)
                    continue;

                // If edge is negative ignore all falling edges
                if (edge<0 && v0<0)
                    continue;

                // Check if there is a zero crossing
                if ((v0<0 && v1<0) || (v0>0 && v1>0))
                    continue;

                // Calculate the position p of the zero-crossing
                // within the interval [rel, rel+1] relative to rel
                // by linear interpolation.
                const double p = v0==v1 ? 0.5 : v0/(v0-v1);

                // If this was at least the second zero-crossing detected
                if (i_prev>=0)
                {
                    // Calculate the distance l between the
                    // current and the last zero-crossing
                    const double l = i+p - (i_prev+p_prev);

                    // By summation, the average length of each
                    // cell is calculated. For the first and last
                    // fraction of a cell, the fraction is applied
                    // as a weight.
                    const double w0 = 1-p_prev;
                    fStat[pos+(spos+i_prev)%1024].first  += w0*l;
                    fStat[pos+(spos+i_prev)%1024].second += w0;

                    for (size_t k=i_prev+1; k<i; k++)
                    {
                        fStat[pos+(spos+k)%1024].first  += l;
                        fStat[pos+(spos+k)%1024].second += 1;
                    }

                    const double w1 = p;
                    fStat[pos+(spos+i)%1024].first  += w1*l;
                    fStat[pos+(spos+i)%1024].second += w1;
                }

                // Remember this zero-crossing position
                p_prev = p;
                i_prev = i;
            }
        }
        fNumEntries++;
    }

    void FillEmptyBins()
    {
        for (int ch=0; ch<160; ch++)
        {
            const auto beg = fStat.begin() + ch*1024;
            const auto end = beg + 1024;

            double   avg = 0;
            uint32_t num = 0;
            for (auto it=beg; it!=end; it++)
            {
                if (it->second<fNumEntries-0.5)
                    continue;

                avg += it->first / it->second;
                num++;
            }
            avg /= num;

            for (auto it=beg; it!=end; it++)
            {
                if (it->second>=fNumEntries-0.5)
                    continue;

                // {
                //     result[i+1].first  = *is2;
                //     result[i+1].second = *iw2;
                // }
                // else
                // {
                it->first  = avg*fNumEntries;
                it->second = fNumEntries;
                // }
            }
        }
    }

    DrsCalibrateTime GetComplete() const
    {
        DrsCalibrateTime rc(*this);
        rc.FillEmptyBins();
        return rc;
    }

    void CalcResult()
    {
        for (int ch=0; ch<160; ch++)
        {
            const auto beg = fStat.begin() + ch*1024;
            const auto end = beg + 1024;

            // First calculate the average length s of a single
            // zero-crossing interval in the whole range [0;1023]
            // (which is identical to the/ wavelength of the
            // calibration signal)
            double s = 0;
            double w = 0;
            for (auto it=beg; it!=end; it++)
            {
                s += it->first;
                w += it->second;
            }
            s /= w;

            // Dividing the average length s of the zero-crossing
            // interval in the range [0;1023] by the average length
            // in the interval [0;n] yields the relative size of
            // the interval in the range [0;n].
            //
            // Example:
            // Average [0;1023]: 10.00  (global interval size in samples)
            // Average [0;512]:   8.00  (local interval size in samples)
            //
            // Globally, on average one interval is sampled by 10 samples.
            // In the sub-range [0;512] one interval is sampled on average
            // by 8 samples.
            // That means that the interval contains 64 periods, while
            // in the ideal case (each sample has the same length), it
            // should contain 51.2 periods.
            // So, the sampling position 512 corresponds to a time 640,
            // while in the ideal case with equally spaces samples,
            // it would correspond to a time 512.
            //
            // The offset (defined as 'ideal - real') is then calculated
            // as 512*(1-10/8) = -128, so that the time is calculated as
            // 'sampling position minus offset'
            //
            double sumw = 0;
            double sumv = 0;
            int n = 0;

            // Sums about many values are numerically less stable than
            // just sums over less. So we do the exercise from both sides.
            // First from the left
            for (auto it=beg; it!=end-512; it++, n++)
            {
                const double valv = it->first;
                const double valw = it->second;

                it->first  = sumv>0 ? n*(1-s*sumw/sumv) : 0;

                sumv += valv;
                sumw += valw;
            }

            sumw = 0;
            sumv = 0;
            n = 1;

            // Second from the right
            for (auto it=end-1; it!=beg-1+512; it--, n++)
            {
                const double valv = it->first;
                const double valw = it->second;

                sumv += valv;
                sumw += valw;

                it->first  = sumv>0 ? n*(s*sumw/sumv-1) : 0;
            }

            // A crosscheck has shown, that the values from the left
            // and right perfectly agree over the whole range. This means
            // the a calculation from just one side would be enough, but
            // doing it from both sides might still make the numerics
            // a bit more stable.
        }
    }

    DrsCalibrateTime GetResult() const
    {
        DrsCalibrateTime rc(*this);
        rc.CalcResult();
        return rc;
    }

    double Offset(uint32_t ch, double pos) const
    {
        const auto p = fStat.begin() + ch*1024;

        const uint32_t f = floor(pos);

        const double v0 = p[f].first;
        const double v1 = p[(f+1)%1024].first;

        return v0 + fmod(pos, 1)*(v1-v0);
    }

    double Calib(uint32_t ch, double pos) const
    {
        return pos-Offset(ch, pos);
    }


};

struct DrsCalibration
{
    std::vector<int32_t> fOffset;
    std::vector<int64_t> fGain;
    std::vector<int64_t> fTrgOff;

    int64_t fNumOffset;
    int64_t fNumGain;
    int64_t fNumTrgOff;

    uint32_t fStep;
    uint16_t fRoi;   // Region of interest for trgoff
    uint16_t fNumTm; // Number of time marker channels in trgoff

    std::string fDateObs;
    std::string fDateRunBeg[3];
    std::string fDateRunEnd[3];
    std::string fDateEnd;

//    uint16_t fDAC[8];

    DrsCalibration() :
        fOffset(1440*1024, 0),
        fGain(1440*1024, 4096),
        fTrgOff (1600*1024, 0),
        fNumOffset(1),
        fNumGain(2000),
        fNumTrgOff(1),
        fStep(0),
        fDateObs("1970-01-01T00:00:00"),
        fDateEnd("1970-01-01T00:00:00")
    {
        for (int i=0; i<3; i++)
        {
            fDateRunBeg[i] = "1970-01-01T00:00:00";
            fDateRunEnd[i] = "1970-01-01T00:00:00";
        }
    }

    DrsCalibration(const DrsCalibration &cpy) :
        fOffset(cpy.fOffset),
        fGain(cpy.fGain),
        fTrgOff(cpy.fTrgOff),
        fNumOffset(cpy.fNumOffset),
        fNumGain(cpy.fNumGain),
        fNumTrgOff(cpy.fNumTrgOff),
        fStep(cpy.fStep),
        fRoi(cpy.fRoi),
        fNumTm(cpy.fNumTm),
        fDateObs(cpy.fDateObs),
        fDateEnd(cpy.fDateEnd)
    {
        for (int i=0; i<3; i++)
        {
            fDateRunBeg[i] = cpy.fDateRunBeg[i];
            fDateRunEnd[i] = cpy.fDateRunEnd[i];
        }
    }

    void Clear()
    {
        // Default gain:
        // 0.575*[45590]*2.5V / 2^16 = 0.99999 V
        fOffset.assign(1440*1024, 0);
        fGain.assign  (1440*1024, 4096);
        fTrgOff.assign(1600*1024, 0);

        fNumOffset = 1;
        fNumGain   = 2000;
        fNumTrgOff = 1;

        fStep = 0;

        fDateObs = "1970-01-01T00:00:00";
        fDateEnd = "1970-01-01T00:00:00";

        for (int i=0; i<3; i++)
        {
            fDateRunBeg[i] = "1970-01-01T00:00:00";
            fDateRunEnd[i] = "1970-01-01T00:00:00";
        }
    }


};

#endif
