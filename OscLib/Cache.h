
#ifndef OSCLIBCACHE_H
#define OSCLIBCACHE_H

#include <unordered_map>
#include <Eigen/Eigen>
#include <iostream> 

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
      // exploit 3f unitarity
      Pte = 1-Pee-Pme;
      Ptm = 1-Pem-Pmm;
      Pet = 1-Pee-Pem;
      Pmt = 1-Pme-Pmm;
      Ptt = Pee+Pem+Pme+Pmm-1;
      Pes = 0;
      Pms = 0;
      Pts = 0;
    }

    Probs(T ee, T me, T te,
          T em, T mm, T tm,
          T et, T mt, T tt)
      : Pee(ee), Pme(me), Pte(te), Pem(em), Pmm(mm), Ptm(tm), Pet(et), Pmt(mt), Ptt(tt)
    {
      Pes = Pee+Pem+Pet;
      Pms = Pme+Pmm+Pmt;
      Pts = Pte+Ptm+Ptt;
    }

    inline __attribute__((always_inline)) T P(int from, int to) const {
      switch(from*100+to){
      case 1212: return Pee;
      case 1412: return Pme;
      case 1612: return Pte;
      case 1214: return Pem;
      case 1414: return Pmm;
      case 1614: return Ptm;
      case 1216: return Pet;
      case 1416: return Pmt;
      case 1616: return Ptt;
      case 1200: return Pes;
      case 1400: return Pms;
      case 1600: return Pts;
      default: abort();
      }
    }

  protected:
    T Pee, Pme, Pte, Pem, Pmm, Ptm, Pet, Pmt, Ptt;
    T Pes, Pms, Pts; // osc probability to active states

  };

  template<class KT, class VT> class ProbCache : public std::unordered_map<KT, Probs<VT>> {};


}
}

#endif

