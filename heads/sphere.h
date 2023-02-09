//
// Created by root on 1/29/23.
//

#ifndef CU_SPHERE_H
#define CU_SPHERE_H

#include "vec3.h"
#include "hittable.h"

class sphere : public hittable{
public:
    point3 center;
    float radius;
    material* mat_ptr;
public:
    __device__ sphere(){}
    __device__ sphere(point3 cen,float r,material* mat)
    :center{cen},radius{r},mat_ptr{mat}{}

    __device__ bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const override{
        vec3 oc = r.origin() - center;
        auto a = r.direction().length_squared();
        auto half_b = dot(oc,r.direction());
        auto c = oc.length_squared() - radius * radius;
        auto discriminant = half_b * half_b - a * c;
        if (discriminant < 0) {
            return false;
        }
        auto sqrtd = sqrt(discriminant);

        // Find the nearest root that lies in the acceptable range.
        auto root = (-half_b - sqrtd) /a;
        if (root < t_min || t_max < root){
            root = (-half_b + sqrtd) /a;
            if (root < t_min || t_max < root){
                return false;
            }
        }
        rec.t = root;
        rec.p = r.at(root);
        vec3 outward_normal = (rec.p - center) / radius;
        rec.set_face_normal(r,outward_normal);
        rec.mat_ptr = mat_ptr;

        return true;
    }
};

#endif //CU_SPHERE_H
