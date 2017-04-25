//#ifndef SAMPLE_H
//#define SAMPLE_H
#pragma once

class Sample {
    private:
        int kt;
    public:
        double dJ;
        Sample(const int _kt, const double _dJ);

        int getKT();
        //inline bool operator<(const Sample& rhs);
        inline bool operator<(const Sample& rhs ){
            if( dJ < rhs.dJ ){
                return true;
            }else{
                return false;
            }
        }
        inline bool operator>(const Sample& rhs ){
            if( dJ > rhs.dJ ){
                return true;
            }else{
                return false;
            }
        }
        //inline bool operator<(const Sample& rhs);
};

//#endif

