#pragma once
#include <Eigen/Core>
#include <Eigen/Dense>

template <class T>
class EigenDummyTensor
{
      public:
            unsigned int l1,l2,l3;

            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> DummyTensor;
            EigenDummyTensor(unsigned int size1, unsigned int size2, unsigned int size3)
            {
                l1 = size1; l2 = size2; l3 = size3;
                DummyTensor.resize(l1, l2*l3);
            }
            EigenDummyTensor(){}
            ~EigenDummyTensor(){}



            Eigen::Vector3i
            getTensorSize()
            {
                Eigen::Vector3i size;
                size << l1, l2 , l3;
                return size;
            }

            void resize_DummyTensor(unsigned int i, unsigned int j, unsigned int k)
            {
                l1 = i; l2 = j; l3 = k;
                DummyTensor.resize(l1, l2*l3);
		DummyTensor.setZero();
            }
            
	    void increment( unsigned int i, unsigned int j, unsigned int k, T value)
	    {
	       DummyTensor(i, j+k*l2) += value;
	    }

            void check(unsigned int i, unsigned int j, unsigned int k)
            {
            if((i >= l1) || (j >= l2) ||(k >= l3)){
               std::cout << " Tensor index greater than size " << std::endl;
               throw  std::exception(); }
            }

            void setTensor(unsigned int i, unsigned int j, unsigned int k, T value)
            {
                check(i, j, k);
                DummyTensor(i, j+k*l2) = value;
            }

            T getTensor(unsigned int i, unsigned int j, unsigned int k)
            {
                check(i, j, k);
                T return_value = DummyTensor(i, j+k*l2);
                return return_value;
            }

            T sumTensor()
            {
            T sum_DummyTensor = DummyTensor.sum();
            return sum_DummyTensor;
            }

            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
            getTensorBlock(unsigned int iini, unsigned int ifin,
                           unsigned int jini, unsigned int jfin,
                           unsigned int kini, unsigned int kfin)
            {
               Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> return_matrix;
               if(((ifin - iini+1) <= 0) ||
                  ((jfin - jini+1) <= 0) ||
                  ((kfin - kini+1) <= 0)){
                   std::cout << " Improper block operation ! "  << std::endl;
                   throw  std::exception(); }
               return_matrix.resize((ifin - iini+1), (jfin - jini+1) * (kfin - kini+1));

            #pragma omp parallel for schedule(static)
            for (int i = 0; i <ifin-iini+1; i++){
	        int count = -1;
                for (int j = 0; j <jfin-jini+1; j++){
                    for (int k = 0; k <kfin-kini+1; k++){
		        count++;
                        return_matrix(i, count) = DummyTensor(i+iini, j+jini+(k+kini)*l2);
                        }
                }
            }

            return return_matrix;
            }

            void setTensorBlock(unsigned int iini, unsigned int ifin,
                           unsigned int jini, unsigned int jfin,
                           unsigned int kini, unsigned int kfin,
                           Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tensor_block)
            {
                if(((ifin - iini+1) <= 0) ||
                   ((jfin - jini+1) <= 0) ||
                   ((kfin - kini+1) <= 0)){
                    std::cout << " Improper block operation ! "  << std::endl;
                    throw  std::exception(); }
                else if ((tensor_block.rows() != (ifin - iini+1)) ||
                         (tensor_block.cols()!= (jfin - jini+1)*(kfin - kini+1))){
                            std::cout << " Given block does not have proper size ! "  << std::endl;
                            throw  std::exception(); }

                #pragma omp parallel for schedule(static) 
                for (int i = 0; i <ifin - iini+1; i++){
		    int count = -1;
                    for (int j = 0; j <jfin - jini+1; j++){
                        for (int k = 0; k <kfin - kini+1; k++){
			    count++;
                            DummyTensor(i+iini, j+jini+(k+kini)*l2) = tensor_block(i, count);
                        }
                    }
                }
            }
};

