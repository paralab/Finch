/**
 * @file templateUtils.h
 * @author Masado Ishii  -- UofU SoC,
 * @date 2019-10-25
 * @brief Some compact template metaprogramming utilities.
 */

#ifndef DENDRO_KT_TEMPLATE_UTILS_H
#define DENDRO_KT_TEMPLATE_UTILS_H

#include <tuple>

namespace litb// adapted from (Johannes Schaub - litb) https://stackoverflow.com/a/7858971
{
  //
  // TupleExtractor / ConstTupleExtractor
  //
  // Most compact pre- c++14 solution I have found to apply std::get<>().
  // The problem is, we want to stuff tuple elements into function parameters.
  // Which would be trivial, except when using a parameter pack.
  //
  // Usage:
  //
  //   std::tuple<double, int> arg_tuple;
  //
  //   struct : public litb::TupleExtractor<void, double, int>
  //   {
  //     virtual void user_defined(double &f, int &i)
  //     {
  //       f = 2.0; i = 5; return;
  //     }
  //   } extractor;
  //
  //   struct : public litb::ConstTupleExtractor<double, double, int>
  //   {
  //     virtual double user_defined(const double &f, const int &i)
  //     {
  //       return pow(f, i);
  //     }
  //   } cextractor;
  //
  //   extractor.applied_to(arg_tuple);                   // assigns (2.0, 5).
  //   double result = cextractor.applied_to(arg_tuple);  // pow(2.0, 5).
  //

  // Type that will be used to identify an integer pack.
  template <int ...Seq>
  struct IntSeq { };

  // Pushes N-1 and decrements N.
  template <int N, int ...IntTail>
  struct GenSeq : GenSeq<N-1, N-1, IntTail...> { };

  // Base case specialization: When N is 0, all ints have been pushed to tail.
  template <int ...IntTail>
  struct GenSeq<0, IntTail...> {
    typedef IntSeq<IntTail...> type;
  };

  // TupleExtractor
  template <typename RetType, typename ...Types>
  struct TupleExtractor
  {
    public:
      RetType applied_to(std::tuple<Types...> &tup) {
        return extract(tup, typename GenSeq<sizeof...(Types)>::type());
      }
    private:
      template <int ...Seq>
      RetType extract(std::tuple<Types...> &tup, IntSeq<Seq...>) {
        return user_defined(std::get<Seq>(tup)...);
      }
    protected:
      virtual RetType user_defined(Types& ...args) = 0;
  };

  // ConstTupleExtractor
  template <typename RetType, typename ...Types>
  struct ConstTupleExtractor
  {
    public:
      RetType applied_to(const std::tuple<Types...> &tup) {
        return extract(tup, typename GenSeq<sizeof...(Types)>::type());
      }
    private:
      template <int ...Seq>
      RetType extract(const std::tuple<Types...> &tup, IntSeq<Seq...>) {
        return user_defined(std::get<Seq>(tup)...);
      }
    protected:
      virtual RetType user_defined(const Types& ...args) = 0;
  };

  template <typename ...>
  struct Typelist { };

  // BiTupleExtractor
  template <typename RetType, typename ...Types>
  struct BiTupleExtractor
  {
    public:
      RetType applied_to(std::tuple<Types...> &tup1,
                         std::tuple<Types...> &tup2) {
        return extract(tup1, tup2, typename GenSeq<sizeof...(Types)>::type());
      }
    private:
      template <int ...Seq>
      RetType extract(std::tuple<Types...> &tup1,
                      std::tuple<Types...> &tup2, IntSeq<Seq...>) {
        return user_defined(std::get<Seq>(tup1)...,
                            std::get<Seq>(tup2)...);
      }
    protected:
      virtual RetType user_defined(Types& ...args1, Types& ...args2) = 0;
  };


}// litb ( https://stackoverflow.com/a/7858971)

#endif//DENDRO_KT_TEMPLATE_UTILS_H
