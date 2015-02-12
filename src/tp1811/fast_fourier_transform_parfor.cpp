#include "fourier_transform.hpp"
#include "tbb/parallel_for.h"
#include <cmath>
#include <cassert>

namespace hpce{
	namespace tp1811{

class fast_fourier_transform_parfor
	: public fourier_transform
{
protected:
	/* Standard radix-2 FFT only supports binary power lengths */
	virtual size_t calc_padded_size(size_t n) const
	{
		assert(n!=0);
		
		size_t ret=1;
		while(ret<n){
			ret<<=1;
		}
		
		return ret;
	}
	

	virtual void forwards_impl(
		size_t n,	const std::complex<double> &wn,
		const std::complex<double> *pIn, size_t sIn,
		std::complex<double> *pOut, size_t sOut
	) const 
	{
		assert(n>0);
		const int DEFAULT_LOOP_K = 32;//
		char *chunk_size=getenv("HPCE_FFT_LOOP_K");//
		size_t chunk_size_int = (chunk_size==NULL) ? DEFAULT_LOOP_K : atoi(chunk_size);//
		
		if (n == 1){
			pOut[0] = pIn[0];
		}else if (n == 2){
			pOut[0] = pIn[0]+pIn[sIn];
			pOut[sOut] = pIn[0]-pIn[sIn];
		}else{
			size_t m = n/2;
			forwards_impl(m,wn*wn,pIn,2*sIn,pOut,sOut);
			forwards_impl(m,wn*wn,pIn+sIn,2*sIn,pOut+sOut*m,sOut);
			 
			std::complex<double> w=std::complex<double>(1.0, 0.0);
			
			if(m<=chunk_size_int){
				for (size_t j=0;j<m;j++){
				std::complex<double> t1 = w*pOut[m+j];
				std::complex<double> t2 = pOut[j]-t1;
				pOut[j] = pOut[j]+t1;                 /*  pOut[j] = pOut[j] + w^i pOut[m+j] */
				pOut[j+m] = t2;                          /*  pOut[j] = pOut[j] - w^i pOut[m+j] */
				w = w*wn;
				}
			}
			
			else{
				typedef tbb::blocked_range<size_t> my_range_t;
				my_range_t range(0, m, chunk_size_int);
				auto f=[=](const my_range_t &chunk_size_int){
					std::complex<double> w2 = pow(wn, chunk_size_int.begin());
					for(size_t i=chunk_size_int.begin(); i!=chunk_size_int.end(); i++){
						std::complex<double> t1 = w2*pOut[m+i];
						std::complex<double> t2 = pOut[i]-t1;
						pOut[i] = pOut[i]+t1;                 /*  pOut[j] = pOut[j] + w^i pOut[m+j] */
						pOut[i+m] = t2;                          /*  pOut[j] = pOut[j] - w^i pOut[m+j] */
						w2 = w2*wn;
					}
				};
				tbb::parallel_for(range, f, tbb::simple_partitioner());
			}
		}
	}
	
	virtual void backwards_impl(
		size_t n,	const std::complex<double> &wn,
		const std::complex<double> *pIn, size_t sIn,
		std::complex<double> *pOut, size_t sOut
	) const 
	{
		complex_t reverse_wn=1.0/wn;
		forwards_impl(n, reverse_wn, pIn, sIn, pOut, sOut);
		
		double scale=1.0/n;
		for(size_t i=0;i<n;i++){
			pOut[i]=pOut[i]*scale;
		}
	}
	
public:
	virtual std::string name() const
	{ return "hpce.fast_fourier_transform"; }
	
	virtual bool is_quadratic() const
	{ return false; }
};

std::shared_ptr<fourier_transform> Create_fast_fourier_transform_parfor()
{
	return std::make_shared<fast_fourier_transform_parfor>();
}
	}
}; // namespace hpce
