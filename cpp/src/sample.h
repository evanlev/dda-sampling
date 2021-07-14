#ifndef SAMPLE_H
#define SAMPLE_H

class Sample
{
private:
    int kt;

public:
    double dJ;
    Sample()
    {
        // empty
    }
    Sample(const int kt, const double dJ);

    int getKTIndex() const;
};

#endif