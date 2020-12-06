#ifndef MDARRAY_H
#define MDARRAY_H 1

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>      // std::ifstream
#include <memory>
#include <sstream>      // std::stringstream
#include <vector>
#include "debug.h"
#include "multind.h"

template<unsigned int D, typename T>
class MDArray 
{
private:
    long m_len;
    std::array<long, D> m_dims; // TODO private
public:
    static constexpr int kDims = D;
    std::vector<T> m_data;
    long m_strs[D]; // TODO private

    // Constructors
    MDArray() = default;

    MDArray(const std::array<long, D>& dims)
    {
        Resize(dims);
    }

    void Resize(const std::array<long, kDims>& dims)
    {
        m_dims = dims; // memcpy(m_dims, dims, kDims * sizeof(long));
        md_calc_strides(kDims, m_strs, m_dims.data(), 1);
        m_len = md_calc_size(kDims, m_dims.data());
        m_data.resize(m_len);
    }

    T kthLargest(long k) const
    {
        assert(k < m_data.size());
        assert(k >= 0);
        std::vector<T> dataCopy = m_data;
        std::nth_element(dataCopy.begin(), dataCopy.begin() + dataCopy.size() - k, dataCopy.end());
        return dataCopy[dataCopy.size() - k - 1];
    }

    const std::array<long, kDims>& Dims() const
    {
        return m_dims;
    }

    long Length() const
    {
        return m_data.size();
    }

    void ReadHeader(const std::string &fname_hdr, std::array<long, kDims> dims)
    {
        // Read header file
        std::ifstream file_hdr(fname_hdr.c_str());
        if (file_hdr.is_open())
        {
            long sized;
            std::vector<long> dims1;
            while (file_hdr >> sized)
            {
                dims1.push_back(sized);
            }
            // Read dimensions, pad with ones
            for (unsigned int i = 0; i < kDims; i++)
            {
                dims[i] = 1u;
            }
            for (unsigned int i = 0; i < dims1.size(); i++)
            {
                dims[i] = dims1[i];
            }
        }
        else
        {
            throw 1;
        }
    }

    void ReadData(const std::string &fname_dat)
    {
        std::ifstream file_dat(fname_dat.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
        if (file_dat.is_open())
        {
            // Read image data
            std::streampos size = file_dat.tellg();
            file_dat.seekg(0, std::ios::beg);
            assert(size % sizeof(T) == 0);
            if (size / sizeof(T) != Length())
            {
                std::cerr << "Size: " << size << " / " << sizeof(T) << " = " << size / sizeof(T) << ", len: " << Length() << "\n";
                assert(0);
            }
            file_dat.read((char *)m_data.data(), size);
            file_dat.close();
        }
        else
        {
            throw 2;
        }
    }

    // Read weighting w. First 4 doubles are size, remaining are the values
    MDArray(const std::string &fname)
    {
        debug_printf(DP_DEBUG3, "Reading %s...\n", fname.c_str());
        const std::string fname_dat = fname + ".dat";
        const std::string fname_hdr = fname + ".hdr";

        try
        {
            std::array<long, kDims> dims;
            ReadHeader(fname_hdr, dims);
            Resize(dims);
            ReadData(fname_dat);
        }
        catch (int e)
        {
            if (e == 1)
            {
                std::cerr << "Error reading header: " << fname << std::endl;
            }
            else if (e == 2)
            {
                std::cerr << "Error reading data: " << fname << std::endl;
            }
            else
            {
                std::cerr << "Error reading ????: " << fname << std::endl;
            }
            exit(0);
        }
        catch (...)
        {
            assert(0);
        }
        debug_printf(DP_DEBUG3, "Done!\n");
    }

    void Clear()
    {
        // TODO
        memset(m_data.data(), 0, sizeof(T) * m_data.size());
    }

    // set data
    void set_data(const T *ptr)
    {
        memcpy(m_data.data(), (void *)ptr, sizeof(T) * m_data.size());
    }

    T sum()
    {
        // TODO std::accumulate
        T res = 0;
        for (long i = 0; i < m_data.size(); i++)
        {
            res += m_data[i];
        }
        return res;
    }

    std::string toString(){
        std::stringstream ss;
        if (m_len > 300)
        {
            ss << "array of size [";
            for (long i = 0; i < kDims; i++)
                ss << m_dims[i] << " ";
            ss << "]";
        }
        else
        {
            ss << "[";
            for (long i = 0; i < m_len; i++)
            {
                ss << m_data[i] << " ";
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

    const T &operator[](const int index) const
    {
        assert(index < m_len && index >= 0);
        return m_data[index];
    }

    T &operator[](const int index)
    {
        assert(index < m_len && index >= 0);
        return m_data[index];
    }

    void Write(const std::string &filename) const
    {
        debug_printf(DP_DEBUG1, "Writing %s...", filename.c_str());

        const std::string filename_dat = filename + ".dat";
        const std::string filename_hdr = filename + ".hdr";

        std::ofstream myfile;
        myfile.open(filename_hdr.c_str());
        for (int i = 0; i < kDims; i++)
        {
            myfile << m_dims[i];
            myfile << " ";
        }
        myfile.close();

        FILE *fp = fopen(filename_dat.c_str(), "wb");
        if (fp == NULL)
        {
            std::cerr << "Failed to open fwrite cfl: " << filename_dat << std::endl;
        }
        size_t ret_code = fwrite(m_data.data(), 1, m_len * sizeof(T), fp);
        if (ret_code != (size_t)m_len * sizeof(T))
        {
            std::cerr << "Failed to write: " << filename_dat << std::endl;
        }
        fclose(fp);
        debug_printf(DP_DEBUG1, "Done!\n");
    }
};

template <unsigned int Dims, typename T>
std::ostream& operator<<(std::ostream &os, const MDArray<Dims, T> & dt){
    os << dt.toString();
    return os;
}

#endif // MDARRAY_H

