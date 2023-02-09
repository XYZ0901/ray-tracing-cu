//
// Created by root on 1/29/23.
//

#ifndef CU_HITTABLE_LIST_H
#define CU_HITTABLE_LIST_H

#include "hittable.h"

class hittable_list : public hittable {
public:
    hittable** objects;
    int list_size{};
public:
    __device__ hittable_list() {};

    __device__ hittable_list(hittable**l,int n) { objects = l;list_size = n; }

    __device__ bool hit(const ray &r, float t_min, float t_max, hit_record &rec) const override {
        hit_record temp_rec;
        bool hit_anything = false;
        auto closest_so_far = t_max;

        for (int i =0;i<list_size;++i ){
//            if(objects[i] == nullptr){
//                continue;
//            }
            if (objects[i]->hit(r,t_min,closest_so_far,temp_rec)){
                hit_anything = true;
                closest_so_far = temp_rec.t;
                rec = temp_rec;
            }
        }
        return hit_anything;
    }
};

#endif //CU_HITTABLE_LIST_H
