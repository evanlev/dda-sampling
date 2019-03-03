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

        int getKTIndex() const;
};

#endif