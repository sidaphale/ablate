#include "dmViewFromOptions.hpp"
#include <utilities/petscError.hpp>
#include <utilities/petscOptions.hpp>
#include "flow/flow.hpp"
ablate::monitors::DmViewFromOptions::DmViewFromOptions(Scope scope, std::string options, std::string optionNameIn)
    : petscOptions(nullptr), optionName(optionNameIn.empty() ? "-CallDmViewFromOptions" : optionNameIn), scope(scope) {
    // Set the options
    if (!options.empty()) {
        PetscOptionsCreate(&petscOptions) >> checkError;

        // build the string
        std::string optionString = optionName + " " + options;
        PetscOptionsInsertString(petscOptions, optionString.c_str());
    }
}

ablate::monitors::DmViewFromOptions::~DmViewFromOptions() {
    if (petscOptions) {
        ablate::utilities::PetscOptionsDestroyAndCheck("DmViewFromOptions", &petscOptions);
    }
}
void ablate::monitors::DmViewFromOptions::Register(std::shared_ptr<Monitorable> monitorableObject) {
    if (scope == Scope::INITIAL) {
        // if the scope is initial, dm plex only once during register
        // this probe will only work with fV flow with a single mpi rank for now.  It should be replaced with DMInterpolationEvaluate
        auto flow = std::dynamic_pointer_cast<ablate::flow::Flow>(monitorableObject);
        if (!flow) {
            throw std::invalid_argument("The DmViewFromOptions monitor can only be used with ablate::flow::Flow");
        }

        DMViewFromOptions(flow->GetDM()) >> checkError;
    }
}

PetscErrorCode ablate::monitors::DmViewFromOptions::DMViewFromOptions(DM dm) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    PetscViewer viewer;
    PetscBool flg;
    PetscViewerFormat format;

    ierr = PetscOptionsGetViewer(PetscObjectComm((PetscObject)dm), petscOptions, NULL, optionName.c_str(), &viewer, &format, &flg);
    CHKERRQ(ierr);

    if (flg) {
        ierr = PetscViewerPushFormat(viewer, format);
        CHKERRQ(ierr);
        ierr = PetscObjectView((PetscObject)dm, viewer);
        CHKERRQ(ierr);
        ierr = PetscViewerFlush(viewer);
        CHKERRQ(ierr);
        ierr = PetscViewerPopFormat(viewer);
        CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&viewer);
        CHKERRQ(ierr);
    }

    PetscFunctionReturn(0);
}

PetscErrorCode ablate::monitors::DmViewFromOptions::CallDmViewFromOptions(TS ts, PetscInt steps, PetscReal time, Vec u, void* mctx) {
    PetscFunctionBeginUser;
    PetscErrorCode ierr;
    DM dm;
    ierr = TSGetDM(ts, &dm);
    CHKERRQ(ierr);

    DmViewFromOptions* monitor = (DmViewFromOptions*)mctx;
    if (monitor->scope == Scope::MONITOR) {
        ierr = monitor->DMViewFromOptions(dm);
    }
    CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

std::ostream& ablate::monitors::operator<<(std::ostream& os, const ablate::monitors::DmViewFromOptions::Scope& v) {
    switch (v) {
        case DmViewFromOptions::Scope::INITIAL:
            return os << "initial";
        case DmViewFromOptions::Scope::MONITOR:
            return os << "monitor";
        default:
            return os;
    }
}

std::istream& ablate::monitors::operator>>(std::istream& is, ablate::monitors::DmViewFromOptions::Scope& v) {
    std::string enumString;
    is >> enumString;

    if (enumString == "initial") {
        v = DmViewFromOptions::Scope::INITIAL;
    } else if (enumString == "monitor") {
        v = DmViewFromOptions::Scope::MONITOR;
    } else {
        throw std::invalid_argument("Unknown norm type " + enumString);
    }
    return is;
}

#include "parser/registrar.hpp"
REGISTER(ablate::monitors::Monitor, ablate::monitors::DmViewFromOptions, "replicates the DMViewFromOptions function in PETSC",
         ENUM(ablate::monitors::DmViewFromOptions::Scope, "scope", "determines if DMViewFromOptions is called initially (initial) or every time step (monitor)"),
         OPT(std::string, "options", "if provided these options are used for the DMView call, otherwise global options is used"),
         OPT(std::string, "optionName", "if provided the optionsName is used for DMViewFromOptions.  Needed if using global options."));
