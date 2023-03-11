#include "../microfacet.h"

Spectrum eval_op::operator()(const ThinFilmIridescence& bsdf) const {
 //
    return(make_zero_spectrum());
}

Real pdf_sample_bsdf_op::operator()(const ThinFilmIridescence& bsdf) const {
   

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
