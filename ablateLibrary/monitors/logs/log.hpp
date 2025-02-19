#ifndef ABLATELIBRARY_LOG_HPP
#define ABLATELIBRARY_LOG_HPP

#include <petsc.h>
#include <vector>

namespace ablate::monitors::logs {
class Log {
   public:
    virtual ~Log() = default;

    // each log must support the print command
    virtual void Printf(const char*, ...) = 0;
    virtual void Print(const char* value) { Printf(value); }
    virtual void Print(const char* name, std::size_t num, const double*, const char* format = nullptr);
    virtual void Initialize(MPI_Comm comm = MPI_COMM_SELF) = 0;

    // built in support calls
    void Print(const char* name, const std::vector<double>& values, const char* format = nullptr);
};
}  // namespace ablate::monitors::logs

#endif  // ABLATELIBRARY_LOG_HPP
