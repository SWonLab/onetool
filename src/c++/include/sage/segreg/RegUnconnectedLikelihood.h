#ifndef REG_UNCONNECTED_LIKELIHOOD
#define REG_UNCONNECTED_LIKELIHOOD

#include "sage/segreg/RegPenetranceCalculator.h"

namespace SAGE {
namespace SEGREG {

class RegUnconnectedLikelihood
{
  public:
    typedef RegPenetranceCalculator::penetrance_info   penetrance_info;

    RegUnconnectedLikelihood(const RegPenetranceCalculator& pen,
                             const model& modp);

    // -------------------------------------------------
    // unconnected likelihood for unconnected individual
    // -------------------------------------------------

    log_double unconnected_likelihood(int                      genotype, 
				      FPED::MemberConstPointer memit)
    {
        log_double p(1.0);

        double psi = mod.freq_sub_model.prob(genotype);

        p = psi * pen.get_penetrance(penetrance_info(*memit,genotype));

        return p;
    }

  protected:  

    const RegPenetranceCalculator& pen;
    const model &                  mod;
};

}
}

#include "sage/segreg/RegUnconnectedLikelihood.ipp"

#endif

