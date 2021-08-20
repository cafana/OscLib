#include <cassert>
#include <vector>

namespace osc
{
  struct Layer
  {
    Layer(double _L, double _rho) : L(_L), rho(_rho)
    {
    }

    double L, rho;
  };

  struct Bin
  {
    Bin(double _low, double _center, double _high)
      : low(_low), center(_center), high(_high)
    {
    }

    double low, center, high;
  };

  // TODO figure out the best name for this concept
  class Path
  {
  public:
    double GetLOverE() const
    {
      assert(fLOverEs.size() == 1 || (fEs.size() == 1 && fLayers.size() == 1));
      if(fLOverEs.size() == 1) return fLOverEs[0].center;
      return fLayers[0].L / fEs[0].center;
    }

    void GetLOverERange(double& lo, double& hi) const
    {
      assert(fLOverEs.size() == 1 || (fEs.size() == 1 && fLayers.size() == 1));
      if(fLOverEs.size() == 1){
        lo = fLOverEs[0].low;
        hi = fLOverEs[0].high;
      }
      else{
        lo = fLayers[0].L / fEs[0].high;
        hi = fLayers[0].L / fEs[0].low;
      }
      assert(lo >= 0 && hi >= 0);
    }

    double GetE() const
    {
      assert(fEs.size() == 1);
      const double E = fEs[0].center;
      assert(E >= 0);
      return E;
    }

    void GetERange(double& lo, double& hi) const
    {
      assert(fEs.size() == 1);
      lo = fEs[0].low;
      hi = fEs[0].high;
      assert(lo >= 0 && hi >= 0);
    }

    double GetL() const
    {
      assert(fLayers.size() == 1);
      return fLayers[0].L;
    }

    double GetRho() const
    {
      assert(fLayers.size() == 1);
      return fLayers[0].rho;
    }

    std::vector<double> GetEs() const
    {
      assert(!fEs.empty());
      std::vector<double> Es;
      Es.reserve(fEs.size());
      for(const Bin& E: fEs){
        assert(E.center >= 0);
        Es.push_back(E.center);
      }
      return Es;
    }

    const std::vector<Layer>& GetLayers() const
    {
      assert(!fLayers.empty());
      return fLayers;
    }

  protected:
    Path(){}

    std::vector<Bin> fEs;
    std::vector<Layer> fLayers;
    std::vector<Bin> fLOverEs;
  };

  // Glorified constructors

  struct Energy: public Path
  {
    Energy(double E)
    {
      fEs.emplace_back(-1, E, -1);
    }

    Energy(const std::vector<double>& Es)
    {
      fEs.reserve(Es.size());
      for(double E: Es) fEs.emplace_back(-1, E, -1);
    }
  };

  struct EnergyEdges: public Path
  {
    EnergyEdges(double Elo, double Ehi)
    {
      fEs.emplace_back(Elo, -1, Ehi);
    }

    EnergyEdges(const std::vector<double>& edges)
    {
      fEs.reserve(edges.size()-1);
      for(unsigned int i = 0; i+1 < edges.size(); ++i){
        fEs.emplace_back(edges[i], -1, edges[i+1]);
      }
    }
  };

  struct LOverE: public Path
  {
    LOverE(double LoE)
    {
      fLOverEs.emplace_back(-1, LoE, -1);
    }
    
    LOverE(const std::vector<double>& LoEs)
    {
      fLOverEs.reserve(LoEs.size());
      for(double LoE: LoEs) fLOverEs.emplace_back(-1, LoE, -1);
    }
  };

  struct LOverEEdges: public Path
  {
    LOverEEdges(double LoElo, double LoEhi)
    {
      fLOverEs.emplace_back(LoElo, -1, LoEhi);
    }

    LOverEEdges(const std::vector<double>& edges)
    {
      fLOverEs.reserve(edges.size()-1);
      for(unsigned int i = 0; i+1 < edges.size(); ++i){
        fLOverEs.emplace_back(edges[i], -1, edges[i+1]);
      }
    }
  };

  // TODO version that allows setting low+center+high for OscCurve interface

  // TODO full Layers interface for atmospherics
}
