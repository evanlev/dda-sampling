#ifndef MDARRAY_H
#define MDARRAY_H 1

#include <algorithm>
#include <numeric>
#include <string>
#include <iostream>
#include <fstream>      // std::ifstream
#include <memory>
#include <sstream>      // std::stringstream
#include <vector>
#include "debug.h"
#include "multind.h"

template<size_t D, typename T>
class MDArray 
{
public:
    static constexpr size_t kDims = D;

    // Constructors
    MDArray() = default;

    explicit MDArray(const std::array<long, D>& dims)
    {
        Resize(dims);
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

    T* Data()
    {
        return m_data.data();
    }

    void Resize(const std::array<long, kDims>& dims)
    {
        m_dims = dims;
        m_strs = md_calc_strides(m_dims, 1);
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

    const std::array<long, D>& Dims() const
    {
        return m_dims;
    }

    const std::array<long, D>& Strides() const
    {
        return m_strs;
    }

    long Length() const
    {
        return m_data.size();
    }

    void Clear()
    {
        std::fill(m_data.begin(), m_data.end(), static_cast<T>(0));
    }

    T sum() const
    {
        return std::accumulate(m_data.begin(), m_data.end(), 0);
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

    const T &operator[](const int index) const
    {
        return const_cast<MDArray<D, T>&>(*this)[index];
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
private:
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

    long m_len;
    std::array<long, D> m_dims;
    std::array<long, D> m_strs;
    std::vector<T> m_data;
};

template <unsigned int Dims, typename T>
std::ostream& operator<<(std::ostream &os, const MDArray<Dims, T> & dt){
    os << dt.toString();
    return os;
}

#endif // MDARRAY_H

