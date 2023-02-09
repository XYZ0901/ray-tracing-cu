//
// Created by root on 1/28/23.
//

#ifndef CU_COLOR_H
#define CU_COLOR_H

#include "vec3.h"

#include <iostream>

void write_color(FILE *fp, color pixel_color) {
    static unsigned char colors[3];
    colors[0] = static_cast<char>(255.999 * pixel_color.x());
    colors[1] = static_cast<char>(255.999 * pixel_color.y());
    colors[2] = static_cast<char>(255.999 * pixel_color.z());
    fwrite(colors, 1, 3, fp);
}

void write_color(std::ostream &out, color pixel_color) {
    out << static_cast<int>(255.999 * pixel_color.x()) << ' ' <<
        static_cast<int>(255.999 * pixel_color.y()) << ' ' <<
        static_cast<int>(255.999 * pixel_color.z()) << '\n';
}

#endif //CU_COLOR_H
