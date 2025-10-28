// Copyright Contributors to the Open Shading Language project.
// SPDX-License-Identifier: BSD-3-Clause
// https://github.com/AcademySoftwareFoundation/OpenShadingLanguage

#pragma once

#include <OSL/dual_vec.h>

#include "bsdl_config.h"
#include <BSDL/static_virtual.h>

#include "raytracer.h"
#include "shading.h"

OSL_NAMESPACE_BEGIN

// StaticVirtual generates a switch/case dispatch method for us given
// a list of possible subtypes. We just need to forward declare them.
// See shading.h for the same simplmentation
struct HomogeneousVolume;
using AbstractMedium = bsdl::StaticVirtual<HomogeneousVolume>;


struct Medium : public AbstractMedium {
    struct Sample {
        OSL_HOSTDEVICE Sample()
            : t(0.0f), transmittance(0.0f), pdf(0.0f)
        {
        }
        OSL_HOSTDEVICE Sample(const Sample& o)
            : t(o.t), scatter(o.scatter), transmittance(o.transmittance)
        {
        }
        OSL_HOSTDEVICE Sample(float t, bool scatter, Color3 transmittance, float pdf)
            : t(t), scatter(scatter), transmittance(transmittance)
        {
        }
        float t;
        bool scatter; 
        Color3 transmittance;
    };
    
    template<typename LOBE>
    OSL_HOSTDEVICE Medium(LOBE* lobe) : AbstractMedium(lobe)
    {
    }

    OSL_HOSTDEVICE Sample sample(const Ray &ray, Sampler &sampler) const;
    OSL_HOSTDEVICE Color3 transmittance(float distance) const;

    int priority = -1;
};



struct HomogeneousVolume final : public Medium, VolumeParams {

    OSL_HOSTDEVICE HomogeneousVolume(const VolumeParams& params)
        : Medium(this), VolumeParams(params)
    {
    }
    OSL_HOSTDEVICE 

    OSL_HOSTDEVICE Medium::Sample sample(const Ray &ray, Sampler &sampler, Intersection& hit) const
    {
        // sample distance
        Vec3 rand_vol = sampler.get();

        float sigma_t_max = max_sigma_t();
        if (sigma_t_max <= 0.0f) {
            // No extinction, treat as vacuum
            return Medium::Sample { hit.t, false, Color3(1.0f) };
        }

        // Exponential sampling: t = -log(1 - ξ) / σ_t
        float t_volume = -logf(1.0f - rand_vol.x) / sigma_t_max;

        // Determine if volume scattering occurs before surface hit
        bool volume_scatter = (t_volume < hit.t);
        float t_event       = volume_scatter ? t_volume : hit.t;

        // Apply transmittance up to the event (always needed)
        Color3 tr = transmittance(t_event);

        return Medium::Sample { t_volume, volume_scatter, tr };
    }

    OSL_HOSTDEVICE Color3 transmittance(float distance) const
    { //Beer-Lambert law
        return Color3(expf(-sigma_t.x * distance),
                      expf(-sigma_t.y * distance),
                      expf(-sigma_t.z * distance));
    }
};

/*
struct MediumStack {

    OSL_HOSTDEVICE MediumStack() : depth(0) {}

    OSL_HOSTDEVICE void push(const VolumeParams& vol)
    {
        if (depth >= MaxEntries)
            return;

        int insert_pos = depth;

        for (int i = 0; i < depth; ++i) {
            if (vol.priority < stack[i].priority) {
                insert_pos = i;
                break;
            }
        }

        for (int j = depth; j > insert_pos; --j) {
            stack[j] = stack[j - 1];
        }

        stack[insert_pos] = vol;
        depth++;
    }

    OSL_HOSTDEVICE void pop()
    {
        if (depth > 0)
            depth--;
    }

    OSL_HOSTDEVICE const VolumeParams* current() const
    {
        return depth > 0 ? &stack[depth - 1] : nullptr;
    }

    enum { MaxEntries = 16 };

    VolumeParams stack[MaxEntries];
    int depth;
};
*/

struct MediumStack {

    OSL_HOSTDEVICE MediumStack() : depth(0), num_bytes(0) {}

    OSL_HOSTDEVICE const Medium* current() const {
        return depth > 0 ? mediums[depth - 1] : nullptr;
    }

    OSL_HOSTDEVICE bool in_medium() const { return depth > 0; }

    OSL_HOSTDEVICE int size() const { return depth; }

    template<typename Medium_Type, typename... Medium_Args>
    OSL_HOSTDEVICE bool push_medium(Medium_Args&&... args) {

        static_assert(sizeof(decltype(Medium_Type::priority)) > 0,
                      "Medium_Type must have a 'priority' field");


        Medium_Type* new_medium = new (pool + num_bytes)
            Medium_Type(std::forward<Medium_Args>(args)...);

        if (depth >= MaxEntries)
            return false;

        if (num_bytes + sizeof(Medium_Type) > MaxSize)
            return false;

        int insert_pos = depth;

        for (int i = 0; i < depth; ++i) {
            if (new_medium->priority < mediums[i]->priority) {
                insert_pos = i;
                break;
            }
        }

        for (int j = depth; j > insert_pos; --j) {
            mediums[j] = mediums[j - 1];
        }

        mediums[insert_pos] = new_medium;
        depth++;
        num_bytes += sizeof(Medium_Type);
        return true;
    }

    OSL_HOSTDEVICE void pop_medium()
    {
        if (depth > 0)
            depth--;
    }

private:
    /// Never try to copy this struct because it would invalidate the bsdf pointers
    OSL_HOSTDEVICE MediumStack(const MediumStack& c);
    OSL_HOSTDEVICE MediumStack& operator=(const MediumStack& c);

    enum { MaxEntries = 8 };
    enum { MaxSize = 256 * sizeof(float) };

    Medium* mediums[MaxEntries];
    char pool[MaxSize];
    int depth, num_bytes;
};


OSL_NAMESPACE_END
