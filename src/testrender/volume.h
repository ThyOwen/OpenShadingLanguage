// Copyright Contributors to the Open Shading Language project.
// SPDX-License-Identifier: BSD-3-Clause
// https://github.com/AcademySoftwareFoundation/OpenShadingLanguage

#pragma once

#include <OSL/dual_vec.h>
#include <BSDL/static_virtual.h>


#include "sampling.h"

OSL_NAMESPACE_BEGIN

// StaticVirtual generates a switch/case dispatch method for us given
// a list of possible subtypes. We just need to forward declare them. 
// See shading.h for the same simplmentation 

struct HenyeyGreenstein;
struct IsotropicPhase;
using AbstractPhaseFunction = bsdl::StaticVirtual<HenyeyGreenstein, IsotropicPhase>;

struct PhaseFunction : public AbstractPhaseFunction { 
    struct Sample {

        OSL_HOSTDEVICE Sample() : wi(0.0f), weight(0.0f), pdf(0.0f) 
        {
        }
        OSL_HOSTDEVICE Sample(const Sample& o)
            : wi(o.wi), weight(o.weight), pdf(o.pdf)
        {
        }
        OSL_HOSTDEVICE Sample(Vec3 wi, Color3 w, float pdf)
            : wi(wi), weight(w), pdf(pdf)
        {
        }

        Vec3 wi;
        Color3 weight;
        float pdf;      
    };
        
    template<typename LOBE> OSL_HOSTDEVICE PhaseFunction(LOBE* lobe) : AbstractPhaseFunction(lobe)
    {
    }
    // Default implementations, to be overriden by subclasses

    OSL_HOSTDEVICE float eval(const Vec3& wo, const Vec3& wi) const
    {
        return {};
    }

    OSL_HOSTDEVICE float pdf(const Vec3& wo, const Vec3& wi) const
    {
        return {};
    }

    OSL_HOSTDEVICE  Sample sample(const Vec3& wo, float rx, float ry) const
    {
        return {};
    }


    // And the "virtual" versions of the above. They are implemented via
    // dispatch with a lambda, but it has to be written after subclasses
    // with their inline methods have been defined. See shading.cpp
    OSL_HOSTDEVICE float eval_vrtl(const Vec3& wo, const Vec3& wi) const;
    OSL_HOSTDEVICE float pdf_vrtl(const Vec3& wo, const Vec3& wi) const;
    OSL_HOSTDEVICE Sample sample_vrtl(const Vec3& wo, float rx, float ry) const;
#ifdef __CUDACC__
    // Copied from BSDF in shading.h
    // TODO: This is a total hack to avoid a misaligned address error
    // that sometimes occurs with the EnergyCompensatedOrenNayar BSDF.
    // It's not clear what the issue is or why this fixes it, but that
    // will take a bit of digging.
    int pad;
#endif
};



struct IsotropicPhase : public PhaseFunction {

    OSL_HOSTDEVICE IsotropicPhase()
        : PhaseFunction(this)
    {
    }

    OSL_HOSTDEVICE float eval(const Vec3& /*wo*/, const Vec3& /*wi*/) const 
    {
        return 0.25f * float(M_1_PI);  // 1 / (4π)
    }

    OSL_HOSTDEVICE PhaseFunction::Sample sample(const Vec3& /*wo*/, float rx, float ry) const 
    {

        float z   = 1.0f - 2.0f * rx;
        float r   = sqrtf(OIIO::clamp(1.0f - z * z, 0.0f, 1.0f));
        float phi = 2.0f * float(M_PI) * ry;
        Vec3 wi   = Vec3(r * cosf(phi), r * sinf(phi), z);

        float pdf     = 0.25f * float(M_1_PI);
        Color3 weight = Color3(1.0f);
        return PhaseFunction::Sample(wi, weight, pdf);
    }

    OSL_HOSTDEVICE float pdf(const Vec3& /*wo*/, const Vec3& /*wi*/) const 
    {
        return 0.25f * float(M_1_PI);
    }
};

struct HenyeyGreenstein : public PhaseFunction {
    const float g;
    OSL_HOSTDEVICE HenyeyGreenstein(float g) 
        : PhaseFunction(this),
        g(g) 
    {
    }

    static OSL_HOSTDEVICE float PhaseHG(float cos_theta, float g) {
        float denom = 1 + g * g + 2 * g * cos_theta;
        return (1 - g * g) / (4 * M_PI * denom * sqrtf(denom));
    }

    OSL_HOSTDEVICE float eval(const Vec3& wo, const Vec3& wi) const 
    {
        return PhaseHG(dot(wo, wi), g);
    }

    OSL_HOSTDEVICE PhaseFunction::Sample sample(const Vec3& wo, float rx, float ry) const 
    {
        TangentFrame frame = TangentFrame::from_normal(wo);

        float cos_theta;
        if (abs(g) < 1e-3f) {
            cos_theta = 1.0f - 2.0f * rx;
        } else {
            float sqr_term = (1 - g * g) / (1 - g + 2 * g * rx);
            cos_theta = (1 + g * g - sqr_term * sqr_term) / (2 * g);
            cos_theta = OIIO::clamp(cos_theta, -1.0f, 1.0f);
        }

        float sin_theta =  sqrtf(OIIO::clamp(1.0f - cos_theta * cos_theta, 0.0f, 1.0f));
        float phi = 2 * M_PI * ry;
        Vec3 local_wi = Vec3(sin_theta * cosf(phi), sin_theta * sinf(phi), cos_theta);

        Vec3 wi = frame.toworld(local_wi);
        float pdf_val = PhaseHG(cos_theta, g);

        Color3 weight = Color3(1.0f);
        return PhaseFunction::Sample(wi, weight, pdf_val);
    }

    OSL_HOSTDEVICE float pdf(const Vec3& wo, const Vec3& wi) const 
    {
        return eval(wo, wi);
    }
};



struct MediumProperties {
    //https://graphics.pixar.com/library/ProductionVolumeRendering/paper.pdf
    Color3 sigma_t       = Color3(0.0f); // extinction coefficient
    Color3 sigma_s       = Color3(0.0f); //scattering 
    float medium_g       = 0.0f;  // volumetric anisotropy
    float refraction_ior = 1.0f;
    int priority         = 0;

    OSL_HOSTDEVICE MediumProperties() :
        sigma_t(0),
        sigma_s(0),
        medium_g(0),
        refraction_ior(1),
        priority(0)
    {
    }

    OSL_HOSTDEVICE MediumProperties(Color3 sigma_t, Color3 sigma_s, 
                                    float medium_g, float refraction_ior, int priority) :
        sigma_t(sigma_t),
        sigma_s(sigma_s),
        medium_g(medium_g),
        refraction_ior(refraction_ior),
        priority(priority)
    {
    }


    OSL_HOSTDEVICE float sample_distance(float rnd) const
    {
        float sigma_t = max_sigma_t();
        if (sigma_t <= 0.0f)
            return std::numeric_limits<float>::infinity();

        // Exponential sampling: t = -log(1 - ξ) / σ_t
        return -logf(1.0f - rnd) / sigma_t;
    }

    //transmittance using Beer-Lambert law
    OSL_HOSTDEVICE Color3 transmittance(float distance) const
    {
        return Color3(expf(-sigma_t.x * distance),
                      expf(-sigma_t.y * distance),
                      expf(-sigma_t.z * distance));
    }

    OSL_HOSTDEVICE bool is_vacuum() const
    {
        return sigma_t.x <= 0.0f && sigma_t.y <= 0.0f && sigma_t.z <= 0.0f;
    }

    OSL_HOSTDEVICE float max_sigma_t() const
    {
        return std::max(std::max(sigma_t.x, sigma_t.y), sigma_t.z);
    }
};

struct MediumStack {
    OSL_HOSTDEVICE MediumStack() : depth(0) {}

    OSL_HOSTDEVICE void push(const MediumProperties& vol)
    {
        if (depth < MaxEntries) {
            stack[depth] = vol;
            depth++;
        }
    }

    OSL_HOSTDEVICE void pop()
    {
        if (depth > 0)
            depth--;
    }

    OSL_HOSTDEVICE const MediumProperties* current() const
    {
        return depth > 0 ? &stack[depth - 1] : nullptr;
    }

    OSL_HOSTDEVICE bool in_medium() const { return depth > 0; }

    OSL_HOSTDEVICE int size() const { return depth; }

//    private:

    enum { MaxEntries = 8 };

    MediumProperties stack[MaxEntries];
    int depth;
};

OSL_NAMESPACE_END
