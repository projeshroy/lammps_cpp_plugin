#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

template <class T, int D>
class EigenMultiDimArray
{
      public:
            std::vector<int> arraydims;
	    std::vector<int> arraydims_offset;
	    int total_size = 1;
            Eigen::Matrix<T, 1, Eigen::Dynamic> DummyMD;
            
	    int getMDSize()
            {
	     total_size = 1;
	     for(int i = 0; i < D; i++)
	        total_size *= arraydims[i];
             return total_size;
            }
            
            void resize_DummyMD(std::vector<int> idarray)
            {
	       arraydims = idarray;
	       arraydims_offset = idarray;

	       total_size = getMDSize();
	       arraydims_offset[0] = 1;

               if(D > 1){
	         for(int i = 1; i < D; i++)
	             arraydims_offset[i] = arraydims[i-1] * arraydims_offset[i-1];
                }

               DummyMD.resize(total_size);
            }

	    void resize_DummyMD(std::vector<int> idarray, T value)
            {
	       arraydims = idarray;
	       arraydims_offset = idarray;

	       total_size = getMDSize();
	       arraydims_offset[0] = 1;

               if(D > 1){
	         for(int i = 1; i < D; i++)
	             arraydims_offset[i] = arraydims[i-1] * arraydims_offset[i-1];
                }

               DummyMD.resize(total_size);

               #pragma omp parallel for
	       for(int i = 0; i < total_size; i++)
	           DummyMD[i] = value;
            }

//basics operations ............................................
	    int gid( std::vector<int> idarray)
	    {
                int id = 0;
                for(int i = 0; i < D; i++)
                   id += idarray[i] * arraydims_offset[i];
                return id;
	    }

            void s(std::vector<int> idarray, T value)
            {
              int id = gid(idarray);
	      DummyMD[id] = value;
	    }

	    void s(T value)
	    {
	      #pragma omp parallel for
	      for(int i = 0; i < total_size; i++)
	          DummyMD[i] = value;
	    }

	    T g(std::vector<int> idarray)
	    {
	      int id = gid(idarray);
              return DummyMD[id];
	    }

            T* p(std::vector<int> idarray)
	    {
              int id = gid(idarray);
              return DummyMD[id];
            }
};

