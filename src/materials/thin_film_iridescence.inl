#include "../microfacet.h"
#include <math.h>
#define SPECTRUM_SAMPLES 3
Spectrum atan(const Spectrum& y, const Spectrum& x) {
    Spectrum value;
    for (int i = 0; i < SPECTRUM_SAMPLES; i++)
        value[i] = std::atan2(y[i], x[i]);
    return value;
}
Spectrum sqr(const Spectrum& s) {
    Spectrum value;
    for (int i = 0; i < SPECTRUM_SAMPLES; i++)
        value[i] = s[i] * s[i];
    return value;
}
void fresnelConductorExact(float cosThetaI,
    float eta, float k,
    double& Rp2, double& Rs2) {
    /* Modified from "Optics" by K.D. Moeller, University Science Books, 1988 */

    float cosThetaI2 = cosThetaI * cosThetaI,
        sinThetaI2 = 1 - cosThetaI2,
        sinThetaI4 = sinThetaI2 * sinThetaI2;

    float temp1 = eta * eta - k * k - sinThetaI2,
        a2pb2 = sqrt(temp1 * temp1 + 4 * k * k * eta * eta),
        a = sqrt(0.5f * (a2pb2 + temp1));

    float term1 = a2pb2 + cosThetaI2,
        term2 = 2 * a * cosThetaI;

    Rs2 = (term1 - term2) / (term1 + term2);

    float term3 = a2pb2 * cosThetaI2 + sinThetaI4,
        term4 = term2 * sinThetaI2;

    Rp2 = Rs2 * (term3 - term4) / (term3 + term4);
}

inline void fresnelPhaseExact(const Spectrum& cost,
    const Spectrum& eta1,
    const Spectrum& eta2, const Spectrum& kappa2,
    Spectrum& phiP, Spectrum& phiS) {
    const Spectrum sinThetaSqr = make_const_spectrum(1) - sqr(cost);
    const Spectrum A = sqr(eta2) * (make_const_spectrum(1) - sqr(kappa2)) - sqr(eta1) * sinThetaSqr;
    const Spectrum B = sqrt(sqr(A) + sqr(2 * sqr(eta2) * kappa2));
    const Spectrum U = sqrt((A + B) / 2.0);
    const Spectrum V = sqrt((B - A) / 2.0);

    phiS = atan(2 * eta1 * V * cost, sqr(U) + sqr(V) - sqr(eta1 * cost));
    phiP = atan(2 * eta1 * sqr(eta2) * cost * (2 * kappa2 * U - (make_const_spectrum(1) - sqr(kappa2)) * V),
        sqr(sqr(eta2) * (make_const_spectrum(1) + sqr(kappa2)) * cost) - sqr(eta1) * (sqr(U) + sqr(V)));
}


Spectrum IrisdescenceTerm(Real eta1, Real eta2, Real eta3, Real kappa3, Real filmheight, bool spectralAntialiasing, bool useGaussianFit, Real h_dot_in, const Real& wavelengths)
{
    Spectrum R12p, T121p, R23p, R12s, T121s, R23s, ct2;

    for (int i = 0; i < SPECTRUM_SAMPLES; ++i)
    {
        const float scale = eta1 / eta2; 
        const float cosThetaTSqr = 1 - (1 - pow(h_dot_in,2)) * pow(scale,2);

        if (cosThetaTSqr <= 0)
        {
            R12s[i] = 1.0;
            R12p[i] = 1.0;

            // Compute the transmission coefficients
            T121p[i] = 0.0;
            T121s[i] = 0.0;
        }
        else
        {
            ct2[i] = std::sqrt(cosThetaTSqr);
            fresnelConductorExact(h_dot_in, eta2 / eta1, 0.0, R12p[i], R12s[i]);

            // Reflected part by the base
            fresnelConductorExact(ct2[i], eta3 / eta2, kappa3 / eta2, R23p[i], R23s[i]);

            // Compute the transmission coefficients
            T121p[i] = 1.0 - R12p[i];
            T121s[i] = 1.0 - R12s[i];

        }
    }

    const Spectrum D = 2.0 * eta2 * filmheight * ct2;
    const Spectrum Dphi = 2.0 * c_PI * D / wavelengths;
    Spectrum phi21p= make_zero_spectrum(), phi21s= make_zero_spectrum(), phi23p= make_zero_spectrum(), phi23s= make_zero_spectrum(), I = make_zero_spectrum();
    Spectrum r123s, r123p, Rs, cosP, irid;
    return(make_const_spectrum(1));
}
Spectrum eval_op::operator()(const ThinFilmIridescence& bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
        dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // If we are going into the surface, then we use normal eta
    // (internal/external), otherwise we use external/internal.
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Real wavelengths[3] ={580.0, 550.0, 450.0};
    Real eta1 = bsdf.eta1;
    Real eta2 = bsdf.eta2;
    Real eta3 = bsdf.eta3;
    Real kappa = bsdf.kappa3;

    Real filmheight= bsdf.filmheight;
    bool spectralAntialiasing = bsdf.spectralAntialiasing;
    bool useGaussianFit = bsdf.useGaussianFit;

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    }
    else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // Compute F / D / G
    // Note that we use the incoming direction
    // for evaluating the Fresnel reflection amount.
    // We can also use outgoing direction -- then we would need to
    // use 1/bsdf.eta and we will get the same result.
    // However, using the incoming direction allows
    // us to use F to decide whether to reflect or refract during sampling.
    Real h_dot_in = dot(half_vector, dir_in);
    Spectrum F = IrisdescenceTerm(eta1, eta2, eta3, kappa, filmheight ,spectralAntialiasing, useGaussianFit , h_dot_in ,*wavelengths);
    Real D = GTR2(dot(frame.n, half_vector), roughness);
    Real G = smith_masking_gtr2(to_local(frame, dir_in), roughness) *
        smith_masking_gtr2(to_local(frame, dir_out), roughness);

    return((F * D * G) / (4 * fabs(dot(frame.n, dir_in))));


}

Real pdf_sample_bsdf_op::operator()(const ThinFilmIridescence& bsdf) const {
       bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1 - anisotropic * Real(0.9));
    constexpr Real min_alpha = Real(0.0001);
    Real alpha_x = max(min_alpha, roughness * roughness / aspect);
    Real alpha_y = max(min_alpha, roughness * roughness * aspect);
    return 0;
}

std::optional<BSDFSampleRecord>
sample_bsdf_op::operator()(const ThinFilmIridescence& bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }

    // Homework 1: implement this!

    // Fetch the texture values for later use
    Real specular_transmission =
        eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = Real(0.25) * eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);

    Real aspect = sqrt(1 - anisotropic * Real(0.9));
    constexpr Real min_alpha = Real(0.0001);
    Real alpha_x = max(min_alpha, roughness * roughness / aspect);
    Real alpha_y = max(min_alpha, roughness * roughness * aspect);
    Real alpha_c = (1 - clearcoat_gloss) * Real(0.1) + clearcoat_gloss * Real(0.001);

    Real diffuse_weight = (1 - metallic) * (1 - specular_transmission);
    Real metallic_weight = (1 - specular_transmission * (1 - metallic));
    Real glass_weight = (1 - metallic) * specular_transmission;
    Real clearcoat_weight = clearcoat;

    // Two cases: 1) if we are coming from "outside" the surface, 
    // sample all lobes
    if (dot(vertex.geometric_normal, dir_in) >= 0) {
        Real total_weight = diffuse_weight + metallic_weight + glass_weight + clearcoat_weight;
        Real diffuse_prob = diffuse_weight / total_weight;
        Real metallic_prob = metallic_weight / total_weight;
        Real glass_prob = glass_weight / total_weight;
        // Real clearcoat_prob = clearcoat_weight / total_weight;
        if (rnd_param_w <= diffuse_prob) {
            return BSDFSampleRecord{
                to_world(vertex.shading_frame, sample_cos_hemisphere(rnd_param_uv)),
                Real(0) /* eta */, Real(1) /* roughness */ };
        }
        else if (rnd_param_w <= (diffuse_prob + metallic_prob)) { // metallic
            // Visible normal sampling

            // Convert the incoming direction to local coordinates
            Vector3 local_dir_in = to_local(vertex.shading_frame, dir_in);
            Vector3 local_micro_normal =
                sample_visible_normals(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

            // Transform the micro normal to world space
            Vector3 half_vector = to_world(vertex.shading_frame, local_micro_normal);
            // Reflect over the world space normal
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            return BSDFSampleRecord{
                reflected, Real(0) /* eta */, roughness /* roughness */ };
        }
        else if (rnd_param_w <= (diffuse_prob + metallic_prob + glass_prob)) { // glass
            if (glass_prob <= 0) {
                // Just to be safe numerically.
                return {};
            }
            // Visible normal sampling

            // Convert the incoming direction to local coordinates
            Vector3 local_dir_in = to_local(vertex.shading_frame, dir_in);
            Vector3 local_micro_normal =
                sample_visible_normals(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

            // Transform the micro normal to world space
            Vector3 half_vector = to_world(vertex.shading_frame, local_micro_normal);
            // Flip half-vector if it's below surface
            if (dot(half_vector, frame.n) < 0) {
                half_vector = -half_vector;
            }

            // Now we need to decide whether to reflect or refract.
            // We do this using the Fresnel term.
            Real h_dot_in = dot(half_vector, dir_in);
            Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
            Real F = fresnel_dielectric(h_dot_in, eta);
            // rescale rnd_param_w from
            // (diffuse_prob + metallic_prob, diffuse_prob + metallic_prob + glass_prob]
            // to
            // (0, 1]
            Real u = (rnd_param_w - (diffuse_prob + metallic_prob)) / glass_prob;
            if (u <= F) {
                // Reflect over the world space normal
                Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
                return BSDFSampleRecord{
                    reflected, Real(0) /* eta */, roughness };
            }
            else {
                // Refraction
                Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
                if (h_dot_out_sq <= 0) {
                    return {};
                }
                // flip half_vector if needed
                if (h_dot_in < 0) {
                    half_vector = -half_vector;
                }
                Real h_dot_out = sqrt(h_dot_out_sq);
                Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
                return BSDFSampleRecord{ refracted, eta, roughness };
            }
        }
        else { // clearcoat
            // Only importance sampling D

            // Appendix B.2 Burley's note
            Real alpha2 = alpha_c * alpha_c;
            // Equation 5
            Real cos_h_elevation =
                sqrt(max(Real(0), (1 - pow(alpha2, 1 - rnd_param_uv[0])) / (1 - alpha2)));
            Real sin_h_elevation = sqrt(max(1 - cos_h_elevation * cos_h_elevation, Real(0)));
            Real h_azimuth = 2 * c_PI * rnd_param_uv[1];
            Vector3 local_micro_normal{
                sin_h_elevation * cos(h_azimuth),
                sin_h_elevation * sin(h_azimuth),
                cos_h_elevation
            };
            // Transform the micro normal to world space
            Vector3 half_vector = to_world(frame, local_micro_normal);

            // Reflect over the world space normal
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            return BSDFSampleRecord{
                reflected, Real(0) /* eta */, sqrt(alpha_c) /* roughness */
            };
        }
    }
    else {
        // 2) otherwise, only consider the glass lobes.

        // Convert the incoming direction to local coordinates
        Vector3 local_dir_in = to_local(vertex.shading_frame, dir_in);
        Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha_x, alpha_y, rnd_param_uv);

        // Transform the micro normal to world space
        Vector3 half_vector = to_world(vertex.shading_frame, local_micro_normal);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }

        // Now we need to decide whether to reflect or refract.
        // We do this using the Fresnel term.
        Real h_dot_in = dot(half_vector, dir_in);
        Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
        Real F = fresnel_dielectric(h_dot_in, eta);
        Real u = rnd_param_w;
        if (u <= F) {
            // Reflect over the world space normal
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            return BSDFSampleRecord{
                reflected, Real(0) /* eta */, roughness /* roughness */ };
        }
        else {
            // Refraction
            Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
            if (h_dot_out_sq <= 0) {
                return {};
            }
            // flip half_vector if needed
            if (h_dot_in < 0) {
                half_vector = -half_vector;
            }
            Real h_dot_out = sqrt(h_dot_out_sq);
            Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
            return BSDFSampleRecord{ refracted, eta, roughness };
        }
    }
}

TextureSpectrum get_texture_op::operator()(const ThinFilmIridescence& bsdf) const {
    return bsdf.base_color;
}
