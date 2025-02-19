#ifndef ABLATELIBRARY_ESSENTIALGHOST_HPP
#define ABLATELIBRARY_ESSENTIALGHOST_HPP

#include <mathFunctions/fieldFunction.hpp>
#include "ghost.hpp"

namespace ablate::flow::boundaryConditions {

class EssentialGhost : public Ghost {
   private:
    static PetscErrorCode EssentialGhostUpdate(PetscReal time, const PetscReal* c, const PetscReal* n, const PetscScalar* a_xI, PetscScalar* a_xG, void* ctx);

    const std::shared_ptr<ablate::mathFunctions::FieldFunction> boundaryFunction;

   public:
    EssentialGhost(std::string boundaryName, std::vector<int> labelId, std::shared_ptr<ablate::mathFunctions::FieldFunction> boundaryFunction, std::string labelName = {});
};
}  // namespace ablate::flow::boundaryConditions
#endif  // ABLATELIBRARY_ESSENTIALGHOST_HPP
