//
// Created by root on 1/29/23.
//

#ifndef CU_CAMERA_H
#define CU_CAMERA_H

class camera {
private:
    point3 origin;
    point3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 u, v, w;
    float lens_radius;
public:
    __device__ camera(point3 lookfrom = {0, 0, 0}, point3 lookat = {0, 0, -1},
                      vec3 vup = {0, 1, 0}, float vfov = 90.f,
                      float aspect_ratio = 16.f / 9.f, float aperture = 0., float focus_dist = 1.f) {
        auto theta = degrees_to_radians(vfov);
        auto h = tan(theta / 2.f);
        auto viewport_height = 2.f * h;
        auto viewport_width = aspect_ratio * viewport_height;

        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        origin = lookfrom;
        horizontal = focus_dist * viewport_width * u;
        vertical = focus_dist * viewport_height * v;
        lower_left_corner = origin - horizontal / 2 - vertical / 2 - focus_dist * w;
        lens_radius = aperture / 2;
    }

    __device__ ray get_ray(float s, float t, curandState *st) {
        vec3 rd = lens_radius * random_in_unit_disk(st);
        vec3 offset = u * rd.x() + v * rd.y();
        return ray(origin + offset,
                   lower_left_corner + s * horizontal + t * vertical - origin - offset);
    }
};

#endif //CU_CAMERA_H
