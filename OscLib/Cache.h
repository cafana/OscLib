
#ifndef OSCLIBCACHE_H
#define OSCLIBCACHE_H

#include <unordered_map>
#include <Eigen/Eigen>

// We want to put ArrayXd into an unordered_map, so define hash and equality
namespace std
{
  template<> struct hash<Eigen::ArrayXd>
  {
    size_t operator()(const Eigen::ArrayXd& x) const
    {
      // Adapted from the boost `hash_combine` function
      // http://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine

      size_t seed = 0;
      for(int i = 0; i < x.size(); ++i) {
        seed ^= std::hash<double>()(x[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2));
      }
      return seed;
    }
  };

  template<> struct equal_to<Eigen::ArrayXd>
  {
    bool operator()(const Eigen::ArrayXd& a, const Eigen::ArrayXd& b) const
    {
      return (a == b).all();
    }
  };
}

namespace osc {
namespace analytic {
  
  template<class T> class Probs
  {
  public:
    Probs(T ee, T me, T em, T mm)
      : Pee(ee), Pme(me), Pem(em), Pmm(mm)
    {
    }

    inline __attribute__((always_inline)) T P(int from, int to) const {
      // convert flavours to indices into matrix
      const int i0 = (from-12)/2;
      const int i1 = (to-12)/2;

      // Exploit unitarity
      switch(i0*3+i1){
      case 0: return Pee;
      case 1: return Pme;
      case 2: return 1-Pee-Pme; // Pte
      case 3: return Pem;
      case 4: return Pmm;
      case 5: return 1-Pem-Pmm; // Ptm
      case 6: return 1-Pee-Pem; // Pet
      case 7: return 1-Pme-Pmm; // Pmt
      case 8: return Pee+Pem+Pme+Pmm-1; // Ptt
      default: abort();
      }
    }

  protected:
    T Pee, Pme, Pem, Pmm;
  };
  
  template<class KT, class VT> class ProbCache : public std::unordered_map<KT, Probs<VT>> {};
}
}

#endif

