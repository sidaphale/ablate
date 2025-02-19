#include "densityMassFractions.hpp"
#include <mathFunctions/functionPointer.hpp>

ablate::flow::fieldFunctions::DensityMassFractions::DensityMassFractions(std::shared_ptr<ablate::flow::fieldFunctions::CompressibleFlowState> flowStateIn)
    : ablate::mathFunctions::FieldFunction("densityYi",
                                           std::make_shared<ablate::mathFunctions::FunctionPointer>(ablate::flow::fieldFunctions::CompressibleFlowState::ComputeDensityYiFromState, flowStateIn.get())),
      flowState(flowStateIn) {}

#include "parser/registrar.hpp"
REGISTER(ablate::mathFunctions::FieldFunction, ablate::flow::fieldFunctions::DensityMassFractions, "initializes the densityYi conserved field variables based upon a CompressibleFlowState",
         ARG(ablate::flow::fieldFunctions::CompressibleFlowState, "state", "The CompressibleFlowState used to initialize"));