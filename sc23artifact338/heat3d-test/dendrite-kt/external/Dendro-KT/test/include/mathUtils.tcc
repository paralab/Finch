/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains useful math helper routines
*
*/


template <typename T>
T normLInfty(T *vec1, T *vec2, unsigned int n)
{
    assert(n>0);
    if((std::isnan(vec1[0])) ||(std::isnan(vec2[0]))) return NAN;

    T l1=fabs(vec1[0]-vec2[0]);
    for(unsigned int i=1;i<n;i++)
    {
        if((std::isnan(vec1[i])) ||(std::isnan(vec2[i])))
            return NAN;

        if(l1<(fabs(vec1[i]-vec2[i])))
            l1=fabs(vec1[i]-vec2[i]);

    }


    return l1;

}


template <typename T>
T normL2(T * vec1,T* vec2, unsigned int n)
{

    T l2=0.0;
    for(unsigned int i=0;i<n;i++)
    {
        if((std::isnan(vec1[i])) ||(std::isnan(vec2[i])))
            return NAN;

        l2+=pow((vec1[i]-vec2[i]),2);
    }


    return sqrt(l2);

}

template <typename T>
T normL2(T * vec,unsigned int n)
{
    T l2=0;
    for(unsigned int i=0;i<n;i++)
    {
        if((std::isnan(vec[i])))
            return NAN;

        l2+=pow((vec[i]),2);
    }



    return sqrt(l2);

}

template <typename T>
T normLInfty(T * vec,unsigned int n)
{
    if((std::isnan(vec[0]))) return NAN;
    T linf=fabs(vec[0]);
    for(unsigned int i=1;i<n;i++)
    {
        if((std::isnan(vec[i]))) return NAN;

        if(linf<fabs(vec[i]))linf=fabs(vec[i]);
    }


    return linf;
}



template <typename T>
T vecMin(T * vec,unsigned int n)
{
    if((std::isnan(vec[0]))) return NAN;
    T min=fabs(vec[0]);
    for(unsigned int i=1;i<n;i++)
    {
        if((std::isnan(vec[i]))) return NAN;

        if(min>fabs(vec[i]))min=fabs(vec[i]);
    }


    return min;
}


template <typename T>
T vecMax(T * vec,unsigned int n)
{
    if((std::isnan(vec[0]))) return NAN;
    T max=fabs(vec[0]);
    for(unsigned int i=1;i<n;i++)
    {
        if((std::isnan(vec[i]))) return NAN;

        if(max<fabs(vec[i]))max=fabs(vec[i]);
    }


    return max;
}

template <typename T>
T normLInfty(T * vec,unsigned int n,MPI_Comm comm)
{
    T linf=normLInfty(vec,n);
    T linf_g=0;

    par::Mpi_Reduce(&linf,&linf_g,1,MPI_MAX,0,comm);
    return linf_g;
}


template <typename T>
T normL2(T * vec1,T* vec2, unsigned int n,MPI_Comm comm)
{

    int rank,npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    T l2=normL2(vec1,vec2,n);
    l2=l2*l2;

    T l2_sum=0;
    par::Mpi_Reduce(&l2,&l2_sum,1,MPI_SUM,0,comm);

    return sqrt(l2_sum);

}

template <typename T>
T normLInfty(T *vec1, T *vec2, unsigned int n, MPI_Comm comm)
{
    int rank,npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    assert(n>0);

    T l1=normLInfty(vec1,vec2,n);
    T l1_max=0;
    par::Mpi_Reduce(&l1,&l1_max,1,MPI_MAX,0,comm);

    return (l1_max);
}


template <typename T>
T normL2(T * vec,unsigned int n,MPI_Comm comm)
{
    T l2=normL2(vec,n);
    l2=l2*l2;

    T l2_sum=0;
    par::Mpi_Reduce(&l2,&l2_sum,1,MPI_SUM,0,comm);

    return sqrt(l2_sum);

}



template <typename T>
T vecMin(T * vec,unsigned int n,MPI_Comm comm)
{
   T min=vecMin(vec,n);
   T min_g=0;
   par::Mpi_Reduce(&min,&min_g,1,MPI_MIN,0,comm);

   return min_g;

}



template <typename T>
T vecMax(T * vec,unsigned int n,MPI_Comm comm)
{
    T max=vecMax(vec,n);
    T max_g=0;
    par::Mpi_Reduce(&max,&max_g,1,MPI_MAX,0,comm);

    return max_g;
}


template <typename T>
T dot(const T* v1, const T*v2,const unsigned int n)
{
    T ans=0;
    for(unsigned int i=0;i<n;i++)
        ans+=(v1[i]*v2[i]);

    return ans;
}

template <typename T>
T dot(const T* v1, const T*v2,const unsigned int n,MPI_Comm comm)
{
    T ans=dot(v1,v2,n);
    T ans_g;
    par::Mpi_Reduce(&ans,&ans_g,1,MPI_SUM,0,comm);
    return ans_g;
}


template<typename T>
void mul(const T alpha, const T* v, const unsigned int n, T* out)
{
    for(unsigned int i=0;i<n;i++)
        out[i]=alpha*v[i];
}

template <typename T>
T add(const T* v1, const T*v2,const unsigned int n, T* out)
{
    for(unsigned int i=0;i<n;i++)
        out[i]=v1[i]+v2[i];
}

template <typename T>
T subt(const T* v1, const T*v2,const unsigned int n, T* out)
{
    for(unsigned int i=0;i<n;i++)
        out[i]=v1[i]-v2[i];
}

template <typename T>
constexpr T intPow(T b, unsigned p, T A)
{
  return (!p ? A : intPow<T>(b, p-1, b*A));
}

template <typename T>
constexpr T intFactorial(T f, T A)
{
  return (f <= 1 ? A : intFactorial<T>(f-1, f*A));
}


// Note: No special handling for integer overflow, we hope n and k are small.
template <typename T>
constexpr T intCombination(T n, T k, T p, T q)
{
  return (k <= 0 ? p / q : intCombination<T>(n, k-1, p * (n-k+1), q * k));
}


template <typename T>
T intCombinationSum(T n, T k, T &combo)
{
  T accum = 0;
  combo = 1;
  for (T kk = 1; kk <= k; kk++)
  {
    accum += combo;
    combo = combo * (n-kk+1) / kk;
  }
  return accum;
}


template <typename T, unsigned int L>
void incrementBaseB(std::array<T,L> &digits, T B, unsigned int lstart)
{
  unsigned int l = lstart;
  do
  {
    digits[l]++;
    if (digits[l] == B)
      digits[l] = 0;
  }
  while (digits[l] == 0 && (l++, l < L));
}


template <typename T, unsigned int L>
void incrementFor(std::array<T,L> &digits,
                  const std::array<T,L> &limits,
                  unsigned int d)
{
  do
  {
    digits[d]++;
    if (digits[d] == limits[d])
      digits[d] = 0;
  }
  while (digits[d] == 0 && (d++, d < L));
}

