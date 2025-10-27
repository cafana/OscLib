
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
      // convert flavours to indices into matrix
      const int i0 = (from-12)/2;
      const int i1 = to = 0 ? 4 : (to-12)/2;

      switch(i0*3+i1){
      case 0: return Pee;
      case 1: return Pme;
      case 2: return Pte;
      case 3: return Pem;
      case 4: return Pmm;
      case 5: return Ptm;
      case 6: return Pet;
      case 7: return Pmt;
      case 8: return Ptt;
      case 9: return Pes;
      case 10: return Pms;
      case 11: return Pts;
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

