//
// Created by root on 1/28/23.
//

#include <iostream>
#include <cfloat>

#include <curand_kernel.h>

#include "heads/rtweekend.h"

#include "heads/global.h"
#include "heads/vec3.h"
#include "heads/color.h"
#include "heads/ray.h"
#include "heads/hittable.h"
#include "heads/sphere.h"
#include "heads/hittable_list.h"
#include "heads/camera.h"
#include "heads/material.h"

#define sample_per_pixel_float 320.f
#define sample_per_pixel_int 320
#define _aspect_ratio (3.f/2.f)
#define _R cos(pi/4.f)

__device__ color ray_color(const ray &r, hittable **world, curandState *st) {
    ray cur_ray = r;
    vec3 cur_attenuation{1., 1., 1.};

    for (int i = 0; i < 50; i++) {
        hit_record rec;
        if ((*world)->hit(cur_ray, 0.1f, FLT_MAX, rec)) {
            ray scattered;
            vec3 attenuation;
            if (rec.mat_ptr->scatter(cur_ray, rec, attenuation, scattered, st)) {
                cur_attenuation = cur_attenuation * attenuation;
                cur_ray = scattered;
            } else {
                return {0, 0, 0};
            }
        } else {
            vec3 unit_direction = unit_vector(cur_ray.direction());
            auto t = 0.5f * (unit_direction.y() + 1.f);
            vec3 c = (1.f - t) * color(1.f, 1.f, 1.f) + t * color(.5f, .7f, 1.f);
            return cur_attenuation * c;
        }
    }
    return {0, 0, 0};
}

__global__ void render(vec3 *fb, camera **cam, hittable **world) {
    int i = blockIdx.x;
    int j = blockIdx.y;
    int t = threadIdx.x;

    __shared__ float cache[sample_per_pixel_int * 3];

    int max_x = gridDim.x;
    int max_y = gridDim.y;

    int pixel_index = (j * max_x + i);

    curandState st;
    curand_init(pixel_index * sample_per_pixel_float + t, 0, 0, &st);

    float u = float(i + curand_uniform(&st)) / float(max_x - 1);
    float v = float(j + curand_uniform(&st)) / float(max_y - 1);

    ray r = (*cam)->get_ray(u, v, &st);
    color pixel_color = ray_color(r, world, &st);

    cache[t * 3 + 0] = pixel_color.x();
    cache[t * 3 + 1] = pixel_color.y();
    cache[t * 3 + 2] = pixel_color.z();

    __syncthreads();

    int idx = blockDim.x >> 1;
    while (idx != 0) {
        if (t < idx) {
            cache[t * 3 + 0] += cache[(t + idx) * 3 + 0];
            cache[t * 3 + 1] += cache[(t + idx) * 3 + 1];
            cache[t * 3 + 2] += cache[(t + idx) * 3 + 2];
        }
        __syncthreads();
        idx >>= 1;
    }

    if (t == 0) {
        color pixel = color{cache[0],
                            cache[1],
                            cache[2]} / sample_per_pixel_float;
        fb[pixel_index] = pixel;
    }
}

#define RND (curand_uniform(&st))

//__global__ void create_world(hittable **d_list, hittable **d_world, camera **cam) {
//    if (threadIdx.x == 0 && blockIdx.x == 0) {
//        curandState st;
//        curand_init(1234, 0, 0, &st);
//        auto ground_material = new lambertian({.5, .5, .5});
//        d_list[0] = new sphere({0, -1000, 0}, 1000, ground_material);
////        int i = 1;
//        for (int a = -11; a < 11; a++) {
//            for (int b = -11; b < 11; b++) {
//                auto choose_mat = RND;
//                point3 center(a + .9 * RND, .2, b + .9 * RND);
//                int idx = ((a + 11) * 22 + b + 11) + 1;
//                d_list[idx] = nullptr;
//                if ((center - point3(4, .2, 0)).length() > .9) {
//                    material *sphere_material;
//                    if (choose_mat < .8) {
//                        auto albedo = _RANDVEC3 * _RANDVEC3;
//                        sphere_material = new lambertian(albedo);
//                    } else if (choose_mat < .95) {
//                        auto albedo = _RANDVEC3 / 2 + color{.5, .5, .5};
//                        auto fuzz = curand_uniform(&st) / 2;
//                        sphere_material = new metal(albedo, fuzz);
//                    } else {
//                        sphere_material = new dielectric(1.5);
//                    }
//                    d_list[idx] = new sphere(center, .2, sphere_material);
//                }
//            }
//        }
//
//        d_list[484 + 1] = new sphere({0, 1, 0}, 1, new dielectric(1.5));
//        d_list[484 + 2] = new sphere({-4, 1, 0}, 1, new lambertian({.4, .2, .1}));
//        d_list[484 + 3] = new sphere({4, 1, 0}, 1, new metal({.7, .6, .5}, 0.));
//
//        *d_world = new hittable_list(d_list, (484 + 1 + 3));
//
//        point3 lookfrom(13, 2, 3);
//        point3 lookat(0, 0, 0);
//        vec3 vup(0, 1, 0);
//        auto dist_to_focus = 10.f;
//        auto aperture = .1f;
//
//        *cam = new camera(lookfrom, lookat,
//                          vup, 20, _aspect_ratio, aperture, dist_to_focus);
//    }
//}

__global__ void create_world(hittable **d_list, hittable **d_world, camera **cam) {

    if (threadIdx.x == 0 && blockIdx.x == 0) {
        curandState st;
        curand_init(0, 0, 0, &st);

        d_list[0] = new sphere(vec3(0, -1000.0, -1), 1000,
                               new lambertian(vec3(0.5, 0.5, 0.5)));

        int i = 1;
        for (int a = -11; a < 11; a++) {
            for (int b = -11; b < 11; b++) {
                float choose_mat = RND;
                vec3 center(a + RND, 0.2, b + RND);
                if (choose_mat < 0.8f) {
                    d_list[i++] = new sphere(center, 0.2,
                                             new lambertian(vec3(RND * RND, RND * RND, RND * RND)));
                } else if (choose_mat < 0.95f) {
                    d_list[i++] = new sphere(center, 0.2,
                                             new metal(vec3(0.5f * (1.0f + RND), 0.5f * (1.0f + RND),
                                                            0.5f * (1.0f + RND)), 0.5f * RND));
                } else {
                    d_list[i++] = new sphere(center, 0.2, new dielectric(1.5));
                }
            }
        }

        d_list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
        d_list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
        d_list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

        *d_world = new hittable_list(d_list, 22 * 22 + 1 + 3);

        vec3 lookfrom(13, 2, 3);
        vec3 lookat(0, 0, 0);
        vec3 vup{0, 1, 0};
        float dist_to_focus = 10.0;
        float aperture = 0.;
        *cam = new camera(lookfrom,
                          lookat,
                          vup,
                          20.0,
                          _aspect_ratio,
                          aperture,
                          dist_to_focus);
    }

}

__global__ void free_world(hittable **d_list, hittable **d_world, camera **cam) {
    for (int i = 0; i < 22*22+1+3; i++) {
//        if (i != 3)
//        if (d_list[i] == nullptr) {
//            continue;
//        }
        delete ((sphere *) d_list[i])->mat_ptr;
        delete d_list[i];
    }

    delete *d_world;
    delete *cam;
}

__global__ void render_init(int max_x, int max_y, curandState *rand_state) {
    int i = blockIdx.x;
    int j = blockIdx.y;

    int pixel_idx = j * max_x + i;

    curand_init(1999, pixel_idx, 0, &rand_state[pixel_idx]);
}

int main(int argc, char **argv) {
    // Image
    const auto aspect_ratio = _aspect_ratio;
    const int image_width = 1200;
    const int image_height = static_cast<int>(image_width / aspect_ratio);

    int num_pixels = image_height * image_width;
    size_t fb_size = num_pixels * sizeof(vec3);

    // allocate FB
    vec3 *fb;
    checkCudaErrors(cudaMalloc((void **) &fb, fb_size));

    // World
    hittable **d_list;
    int num_hitables = 22 * 22 + 1 + 3;
    checkCudaErrors(cudaMalloc((void **) &d_list, num_hitables * sizeof(hittable *)));
    hittable **d_world;
    checkCudaErrors(cudaMalloc((void **) &d_world, sizeof(hittable *)));
    // Camera
    camera **d_cam;
    checkCudaErrors(cudaMalloc((void **) &d_cam, sizeof(camera *)));
    create_world<<<1, 1>>>(d_list, d_world, d_cam);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());

    dim3 block_size(image_width, image_height);
    dim3 thread_size(sample_per_pixel_int);

//    curandState *d_rand_state;
//    checkCudaErrors(cudaMalloc((void **) &d_rand_state, num_pixels * sizeof(curandState)));
//    render_init<<<block_size, 1>>>(image_width, image_height, d_rand_state);
//    checkCudaErrors(cudaGetLastError());
//    checkCudaErrors(cudaDeviceSynchronize());

    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));

    checkCudaErrors(cudaEventRecord(start));

    render<<<block_size, thread_size>>>(fb, d_cam, d_world);

    checkCudaErrors(cudaEventRecord(stop));
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaEventSynchronize(stop));
    float elapsed_time;
    checkCudaErrors(cudaEventElapsedTime(&elapsed_time, start, stop));
    printf("Time = %g ms.\n", elapsed_time);
    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));

    free_world<<<1, 1>>>(d_list, d_world, d_cam);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaFree(d_list));
    checkCudaErrors(cudaFree(d_world));
    checkCudaErrors(cudaFree(d_cam));

    auto rb = (vec3 *) calloc(num_pixels, sizeof(vec3));
    checkCudaErrors(cudaMemcpy(rb, fb, num_pixels * sizeof(vec3), cudaMemcpyDeviceToHost));

    checkCudaErrors(cudaFree(fb));

    // Render
//    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";

    FILE *fp = fopen("./image.ppm", "wb");
    (void) fprintf(fp, "P6\n%d %d\n255\n", image_width, image_height);
    for (int j = image_height - 1; j >= 0; --j) {
//        std::cerr << "\rScanlines remaining: " << j << ' ';
//        std::cerr.flush();
        for (int i = 0; i < image_width; ++i) {
            size_t pixel_idx = (j * image_width + i);
            color pixel_color = rb[pixel_idx];
            pixel_color.e[0] = sqrt(pixel_color.x());
            pixel_color.e[1] = sqrt(pixel_color.y());
            pixel_color.e[2] = sqrt(pixel_color.z());
            write_color(fp, pixel_color);
//            write_color(std::cout, pixel_color);
        }
        UpdateProgress(j / (float) image_height);
    }
    fclose(fp);

    return 0;
}
