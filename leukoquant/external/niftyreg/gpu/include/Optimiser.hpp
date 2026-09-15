/** @file Optimiser.hpp
 * @author Marc Modat
 * @date 20/07/2012
 */

#pragma once

#include "_reg_tools.h"

/* *************************************************************** */
namespace NiftyReg {
/* *************************************************************** */
/** @brief Interface between the registration class and the optimiser
 */
class InterfaceOptimiser {
public:
    // [leukoquant patch] Pure-virtual interface with no virtual destructor --
    // the same missing-virtual-destructor pattern fixed in reg_aladin (see
    // that class's own patch comment for the real crash this caused).
    // reg_base (and therefore reg_f3d/reg_f3d2) derive from this, so without
    // a virtual destructor here, ANY caller deleting one of those through a
    // reg_base<T>*/InterfaceOptimiser* base pointer would hit the identical
    // heap-corrupting bug reg_aladin_sym did -- not actively triggered
    // today (GIF's own reg_f3d usage happens to use matching declared/
    // allocated types), but the same latent defect. Fixed at the root
    // interface so the whole class hierarchy is safe for polymorphic
    // deletion, not just the one call site that happened to crash.
    virtual ~InterfaceOptimiser() {}
    /// @brief Returns the registration current objective function value
    virtual double GetObjectiveFunctionValue() = 0;
    /// @brief The transformation parameters are optimised
    virtual void UpdateParameters(float) = 0;
    /// @brief The best objective function values are stored
    virtual void UpdateBestObjFunctionValue() = 0;
};
/* *************************************************************** */
/** @class Optimiser
 * @brief Standard gradient ascent optimisation
 */
template <class T>
class Optimiser {
protected:
    bool isSymmetric;
    size_t dofNumber;
    size_t dofNumberBw;
    size_t ndim;
    T *currentDof; // pointer to the cpp nifti image array
    T *currentDofBw; // pointer to the cpp nifti image array (backwards)
    T *bestDof;
    T *bestDofBw;
    T *gradient;
    T *gradientBw;
    bool optimiseX;
    bool optimiseY;
    bool optimiseZ;
    size_t maxIterationNumber;
    size_t currentIterationNumber;
    double bestObjFunctionValue;
    double currentObjFunctionValue;
    InterfaceOptimiser *intOpt;

#ifdef NR_TESTING
public:
#endif
    /// @brief Update the gradient array
    virtual void UpdateGradientValues() {}

public:
    Optimiser();
    virtual ~Optimiser();
    virtual void StoreCurrentDof();
    virtual void RestoreBestDof();
    virtual size_t GetDofNumber() {
        return this->dofNumber;
    }
    virtual size_t GetDofNumberBw() {
        return this->dofNumberBw;
    }
    virtual size_t GetNDim() {
        return this->ndim;
    }
    virtual size_t GetVoxNumber() {
        return this->dofNumber / this->ndim;
    }
    virtual size_t GetVoxNumberBw() {
        return this->dofNumberBw / this->ndim;
    }
    virtual T* GetBestDof() {
        return this->bestDof;
    }
    virtual T* GetBestDofBw() {
        return this->bestDofBw;
    }
    virtual T* GetCurrentDof() {
        return this->currentDof;
    }
    virtual T* GetCurrentDofBw() {
        return this->currentDofBw;
    }
    virtual T* GetGradient() {
        return this->gradient;
    }
    virtual T* GetGradientBw() {
        return this->gradientBw;
    }
    virtual bool GetOptimiseX() {
        return this->optimiseX;
    }
    virtual bool GetOptimiseY() {
        return this->optimiseY;
    }
    virtual bool GetOptimiseZ() {
        return this->optimiseZ;
    }
    virtual size_t GetMaxIterationNumber() {
        return this->maxIterationNumber;
    }
    virtual size_t GetCurrentIterationNumber() {
        return this->currentIterationNumber;
    }
    virtual size_t ResetCurrentIterationNumber() {
        return this->currentIterationNumber = 0;
    }
    virtual double GetBestObjFunctionValue() {
        return this->bestObjFunctionValue;
    }
    virtual void SetBestObjFunctionValue(double i) {
        this->bestObjFunctionValue = i;
    }
    virtual double GetCurrentObjFunctionValue() {
        return this->currentObjFunctionValue;
    }
    virtual void IncrementCurrentIterationNumber() {
        this->currentIterationNumber++;
    }
    virtual void Initialise(size_t nvox,
                            int ndim,
                            bool optX,
                            bool optY,
                            bool optZ,
                            size_t maxIt,
                            size_t startIt,
                            InterfaceOptimiser *intOpt,
                            T *cppData,
                            T *gradData,
                            size_t nvoxBw,
                            T *cppDataBw,
                            T *gradDataBw);
    virtual void Optimise(T maxLength,
                          T smallLength,
                          T& startLength);
    virtual void Perturbation(float length);
    // Reset any accumulated search-direction history. No-op for memoryless optimisers.
    virtual void RestartOptimisation() {}
};
/* *************************************************************** */
/** @class ConjugateGradient
 * @brief Conjugate gradient ascent optimisation
 */
template <class T>
class ConjugateGradient: public Optimiser<T> {
protected:
    T *array1 = nullptr;
    T *array1Bw = nullptr;
    T *array2 = nullptr;
    T *array2Bw = nullptr;
    bool firstCall = true;

#ifdef NR_TESTING
public:
#endif
    virtual void UpdateGradientValues() override;

public:
    ConjugateGradient() { NR_FUNC_CALLED(); }
    virtual ~ConjugateGradient();
    virtual void Initialise(size_t nvox,
                            int ndim,
                            bool optX,
                            bool optY,
                            bool optZ,
                            size_t maxIt,
                            size_t startIt,
                            InterfaceOptimiser *intOpt,
                            T *cppData,
                            T *gradData,
                            size_t nvoxBw,
                            T *cppDataBw,
                            T *gradDataBw) override;
    virtual void Optimise(T maxLength,
                          T smallLength,
                          T& startLength) override;
    virtual void Perturbation(float length) override;
    virtual void RestartOptimisation() override;
};
/* *************************************************************** */
} // namespace NiftyReg
/* *************************************************************** */
