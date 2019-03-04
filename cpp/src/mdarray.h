#ifndef MDARRAY_H
#define MDARRAY_H 1

#include <string>
#include <iostream>
#include <fstream>      // std::ifstream
#include <sstream>      // std::stringstream
#include <vector>
#include "debug.h"
#include "multind.h"

// TODO
// array(sub) = ***
// sub is a long *

template<int Dims, typename T>
class MDArray 
{
    private:
        long len;
    public:
        const static int kDims = Dims;
        T *data;
        long strs[kDims];
        long dims[kDims];

    void init(const long _dims[]){
        debug_printf(DP_DEBUG3, "Copying dims\n");
        memcpy(dims, _dims, kDims*sizeof(long));
        debug_print_dims(DP_DEBUG3, kDims, dims);
        md_calc_strides(kDims, strs, dims, 1);
        len = md_calc_size(kDims, dims);
        debug_printf(DP_DEBUG3, "Allocating %d\n", len);
        data = new T[len];
        debug_printf(DP_DEBUG3, "Done with MDArray constructor\n");
    }

    long length() const {
        return len;
    }

    // Constructor
    MDArray(){
        data = nullptr;
    }
    MDArray(const long _dims[]){
        init(_dims);
    }

    // Read weighting w. First 4 doubles are size, remaining are the values
    MDArray<kDims, T>(const std::string &fname){
        debug_printf(DP_DEBUG3, "Reading %s...\n", fname.c_str());
        int i = 0;
        using namespace std;

        std::string fname_dat = fname + ".dat";
        std::string fname_hdr = fname + ".hdr";

        try{
            // Read header file
            ifstream file_hdr (fname_hdr.c_str() );
            if (file_hdr.is_open())
            {
                long sized;
                vector<long> dims1;
                while( file_hdr >> sized ){
                    dims1.push_back(sized);
                }
                // Read dimensions, pad with ones
                long dims[kDims];
                for( unsigned int i = 0 ; i < kDims ; i++ ){
                    dims[i] = 1u;
                }
                for( unsigned int i = 0 ; i < dims1.size() ; i++ )
                    dims[i] = dims1[i];

                this->init(dims);
            }else{
                throw 1;
            }
                
            std::ifstream file_dat (fname_dat.c_str(), ios::in|ios::binary|ios::ate);
            if (file_dat.is_open())
            {
                // Read image data
                streampos size = file_dat.tellg();
                file_dat.seekg (0, ios::beg);
                assert(size % sizeof(T) == 0);
                if(size / sizeof(T) != this->length()){
                    cerr << "Size: " << size << " / " << sizeof(T) << " = " << size / sizeof(T) << ", len: " << this->length() << endl;
                    assert(0);
                }
                file_dat.read ((char *) this->data, size);
                file_dat.close();
            }else{
                throw 2;
            }
        }catch( int e ){
            if( e == 1 ){
                cerr << "Error reading header: " << fname << endl;
            }else if(e == 2){
                cerr << "Error reading data: " << fname << endl;
            }else{
                cerr << "Error reading ????: " << fname << endl;
            }
            exit(0);
        }catch( ... ){
            assert(0);
        }
        debug_printf(DP_DEBUG3, "Done!\n");
    }

    // Destructor
    ~MDArray(){
        if( data ){
            delete [] data;
        }
    }

    void Clear() {
        memset(data, 0, sizeof(T) * len);
    }

    // set data
    void set_data(const T *ptr){
        memcpy((void *) data, (void *) ptr, sizeof(T) * len);
    }

    T sum(){
        T res = 0;
        for( long i = 0 ; i < len ; i++ ){
            res += data[i];
        }
        return res;
    }

    std::string toString(){
        std::stringstream ss;
        if( len > 300 ){
            //cout << "len > 300 " << endl;
            ss << "array of size [";
            for( long i = 0 ; i < kDims ; i++ )
                ss << dims[i] << " ";
            ss << "]";
        }else{
            ss << "[";
            for( long i = 0 ; i < len ; i++ ){
                ss << data[i] << " ";
            }
            ss << "];\n";
        }

        return ss.str();
    }

    /*
    T operator()(const long *subs){
        long ind = sub2ind<long>(D, strs, subs);
        return data[ind];
    }
    */

    // For [] operator
#if 0
    class Proxy
    {
        MDArray<T> &a;
        int idx;
        public:
            Proxy(MDArray<T> &a, int idx) : a(a), idx(idx) {}
            long& operator= (long x) { a.subs[idx] = x; return a.subs[idx]; }
            // TODO: need to overload << to print a Proxy 
    };
    Proxy operator[] (int index) { return Proxy(*this, index); }
#endif

    const T& operator[] (const int index) const {
        assert( index < len && index >= 0 );
        return this->data[index]; 
    }

    T& operator[] (const int index) { 
        assert( index < len && index >= 0 );
        return this->data[index];
    }

    void Write(const std::string &filename) const {
        debug_printf(DP_DEBUG1, "Writing %s...", filename.c_str());
        using namespace std;

        string filename_dat = filename + ".dat";
        string filename_hdr = filename + ".hdr";
        
        ofstream myfile;
        myfile.open (filename_hdr.c_str());
        for( int i = 0 ; i < kDims ; i++ ){
            myfile << this->dims[i];
            myfile << " ";
        }
        myfile.close();

        FILE *fp = fopen(filename_dat.c_str(), "wb");
        if( fp == NULL ){
            cerr << "Failed to open fwrite cfl: " << filename_dat << endl;
        }
        size_t ret_code = fwrite(this->data, 1, this->len*sizeof(T), fp);
        if( ret_code != (size_t) this->len * sizeof(T) ){
            cerr << "Failed to write: " << filename_dat << endl;
        }
        fclose(fp);
        debug_printf(DP_DEBUG1, "Done!\n");
    }

    
};

template <uint32_t Dims, typename T>
std::ostream& operator<<(std::ostream &os, const MDArray<Dims, T> & dt){
    os << dt.toString();
    return os;
}

#endif // MDARRAY_H

