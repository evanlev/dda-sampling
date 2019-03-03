#ifndef SAMPLE_H
#define SAMPLE_H

class Sample {
    private:
        int kt;
    public:
        double dJ;
        Sample() {
            // empty
        }
        Sample(const int _kt, const double _dJ);

        int getIndex() const;

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
};

#endif