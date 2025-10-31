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
struct EmptyVolume;
using AbstractMedium = bsdl::StaticVirtual<HomogeneousVolume, EmptyVolume>;

template <typename, typename = void>
struct has_equal : std::false_type {};

template <typename T>
struct has_equal<T, std::void_t<decltype(std::declval<const T&>() == std::declval<const T&>())>> : std::true_type {};

struct Medium : public AbstractMedium {
    struct Sample {
        OSL_HOSTDEVICE Sample()
            : scatter(false), t(0.0f), transmittance(0.0f), weight(0.0f)
        {
        }
        OSL_HOSTDEVICE Sample(const Sample& o)
            : scatter(o.scatter), t(o.t), transmittance(o.transmittance), weight(o.weight)
        {
        }
        OSL_HOSTDEVICE Sample(bool scatter, float t, Color3 transmittance, Color3 weight)
            : scatter(scatter), t(t), transmittance(transmittance), weight(weight)
        {
        }
        bool scatter; 
        float t;
        Color3 transmittance;
        Color3 weight;
    };
    
    template<typename LOBE>
    OSL_HOSTDEVICE Medium(LOBE* lobe) : AbstractMedium(lobe)
    {
    }

    OSL_HOSTDEVICE Sample sample(const Ray &ray, Sampler &sampler, Intersection& hit) const {
        return {};
    }

    OSL_HOSTDEVICE Color3 transmittance(float distance) const {
        return {};
    }

    OSL_HOSTDEVICE Sample sample_vrtl(const Ray &ray, Sampler &sampler, Intersection& hit) const;
    OSL_HOSTDEVICE Color3 transmittance_vrtl(float distance) const;

    template<typename Param_Type>
    OSL_HOSTDEVICE bool has_params_as(const Param_Type &params) const {
        if constexpr (std::is_base_of_v<Param_Type, std::decay_t<decltype(*this)>>) {
            auto lhsParam = static_cast<const Param_Type*>(this);
            if (auto rhsParam = dynamic_cast<const Param_Type*>(&params)) {
                if constexpr (has_equal<Param_Type>::value) {
                    return *lhsParam == *rhsParam;
                }
            }
        }
        return false;
    }

    OSL_HOSTDEVICE bool operator==(const Medium& rhs) const {
        using Self = std::decay_t<decltype(*this)>;

        if (typeid(Self) == typeid(rhs)) { // only compare if both objects have the same conrete type
            const Self& lhsRef = static_cast<const Self&>(*this);
            const Self& rhsRef = static_cast<const Self&>(rhs);

            if constexpr (has_equal<Self>::value) {
                return lhsRef.operator==(rhsRef); // call the subclass's own operator== (not Medium's)
            }
        }
        return false;
    }

    int priority = -1;
};

struct HomogeneousVolume final : public Medium, VolumeParams {

    OSL_HOSTDEVICE HomogeneousVolume(const VolumeParams& params)
        : Medium(this), VolumeParams(params)
    {
    }

    OSL_HOSTDEVICE Medium::Sample sample(const Ray &ray, Sampler &sampler, Intersection& hit) const
    {
        Vec3 rand_vol = sampler.get();
        float max_sigma_t = std::max(std::max(sigma_t.x, sigma_t.y), sigma_t.z);
    /*
        if (max_sigma_t <= 0.0f) {
            // vacuum: no attenuation, pdf irrelevant
            return Medium::Sample { hit.t, false, Color3(1.0f) };
        }
    */
        float t_volume = -logf(1.0f - rand_vol.x) / max_sigma_t;
        bool volume_scatter = (t_volume < hit.t);

        Color3 weight;
        Color3 tr;
        
        if (volume_scatter) {
            tr = transmittance(t_volume);

            Color3 albedo = Color3(
                sigma_s.x / sigma_t.x,
                sigma_s.y / sigma_t.y,
                sigma_s.z / sigma_t.z
            );

            weight = albedo / tr;
        } else {
            tr = transmittance(hit.t);
            weight = Color3(
                1.0 / tr.x,
                1.0 / tr.y,
                1.0 / tr.z
            );
        }
    /*
        if (pdf <= 0.0f) {
            // ensure Tr is black in this degenerate case if we expected absorption
            return Medium::Sample { hit.t, false, tr };
        }
    */
        return Medium::Sample { volume_scatter, t_volume, tr, weight };
    }

    OSL_HOSTDEVICE Color3 transmittance(float distance) const
    { //Beer-Lambert law
        return Color3(expf(-sigma_t.x * distance),
                      expf(-sigma_t.y * distance),
                      expf(-sigma_t.z * distance));
    }
};

struct EmptyVolume final : public Medium {
    OSL_HOSTDEVICE EmptyVolume() :
        Medium(this)
    {
    }
};

struct MediumStack {

    OSL_HOSTDEVICE MediumStack() : depth(0), num_bytes(0) {}

    OSL_HOSTDEVICE const Medium* current() const {
        return depth > 0 ? mediums[depth - 1] : nullptr;
    }

    OSL_HOSTDEVICE bool in_medium() const { return depth > 0; }

    OSL_HOSTDEVICE int size() const { return depth; }

    template<typename Medium_Type, typename... Medium_Args>
    OSL_HOSTDEVICE bool push_medium(Medium_Args&&... args) {

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
