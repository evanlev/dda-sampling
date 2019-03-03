#include "sample.h"

// Sample class    
Sample::Sample(const int _kt, const double _dJ) : kt(_kt), dJ(_dJ){
    // Empty
}
int Sample::getIndex() const{
    return kt;
}



