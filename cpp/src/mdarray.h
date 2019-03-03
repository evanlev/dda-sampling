#ifndef MDARRAY_H
#define MDARRAY_H 1

#include "debug.h"

// TODO
// array(sub) = ***
// sub is a long *

template<typename T>
class MDArray 
{
    public:
        T *data;
        long D;
        long strs[DIMS];
        long dims[DIMS];
        size_t el_size;
        long len;

    long length(){
        return len;
    }

    // Constructor
    MDArray(){
        data = NULL;
    }
    MDArray(const long _D, const long *_dims){
        D = _D;
        for( unsigned int i = 0 ; i < DIMS ; i++ ){
            dims[i] = 1;
        }
        debug_printf(DP_DEBUG3, "Copying dims\n");
        memcpy(dims, _dims, D*sizeof(long));
        debug_print_dims(DP_DEBUG3, DIMS, dims);
        md_calc_strides(DIMS, strs, dims, 1);
        el_size = sizeof(T);
        len = md_calc_size(DIMS, dims);
        debug_printf(DP_DEBUG3, "Allocating %d\n", len);
        data = new T[len];
        debug_printf(DP_DEBUG3, "Done with MDArray constructor\n");
    }
    // Destructor
    ~MDArray(){
        if( data ){
            delete [] data;
        }
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

    string toString(){
        std::stringstream ss;
        if( len > 300 ){
            //cout << "len > 300 " << endl;
            ss << "array of size [";
            for( long i = 0 ; i < D ; i++ )
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

    T operator[] (int index) { return this->data[index]; }

    static void write_array(string filename, MDArray<T> *array){
        debug_printf(DP_DEBUG1, "Writing %s...", filename.c_str());

        string filename_dat = filename + ".dat";
        string filename_hdr = filename + ".hdr";
        
        ofstream myfile;
        myfile.open (filename_hdr.c_str());
        for( int i = 0 ; i < array->D ; i++ ){
            myfile << array->dims[i];
            myfile << " ";
        }
        myfile.close();

        FILE *fp = fopen(filename_dat.c_str(), "wb");
        if( fp == NULL ){
            cerr << "Failed to open fwrite cfl: " << filename_dat << endl;
        }
        size_t ret_code = fwrite(array->data, 1, array->len*sizeof(T), fp);
        if( ret_code != (size_t) array->len * sizeof(T) ){
            cerr << "Failed to write: " << filename_dat << endl;
        }
        fclose(fp);
        debug_printf(DP_DEBUG1, "Done!\n");
    }

    // Read weighting w. First 4 doubles are size, remaining are the values
    static MDArray<T>* read_array(const string fname){
        debug_printf(DP_DEBUG3, "Reading %s...\n", fname.c_str());
        int i = 0;

        string fname_dat = fname + ".dat";
        string fname_hdr = fname + ".hdr";

        MDArray<T> *res = NULL;

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
                // Read dimensions
                long ndims;
                long dims[dims1.size()];
                for( unsigned int i = 0 ; i < dims1.size() ; i++ )
                    dims[i] = dims1[i];

                res = new MDArray<T>(dims1.size(), dims);
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
                if(size / sizeof(T) != res->length()){
                    cerr << "Size: " << size << " / " << sizeof(T) << " = " << size / sizeof(T) << ", len: " << res->length() << endl;
                    assert(0);
                }
                file_dat.read ((char *) res->data, size);
                file_dat.close();

                /*
                cout << "the entire file content is in memory \n";
                for(i=0; i<=10; i++){
                    T value = wblock [i];
                    cout << "value ("<<i<<")=" << value << "\n";
                }
                */
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
        return res;
    }
};

template <typename T>
ostream& operator<<(ostream &os, const MDArray<T> & dt){
    os << dt.toString();
    return os;
}

#endif // MDARRAY_H

