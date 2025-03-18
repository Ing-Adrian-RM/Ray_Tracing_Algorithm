///////////////////////////////////////////////////////////////////////////////
//
// Ray Tracing Header
//
///////////////////////////////////////////////////////////////////////////////

//Verification //////////////////////////////////
#ifndef RAY_TRACING_H
#define RAY_TRACING_H


//Includes //////////////////////////////////////
#include "window_module.h"
#include "ray_elements.h"
#include "vector.h"

//Structs ///////////////////////////////////////

//Global Vars ///////////////////////////////////
float epsilon = 0.5;
float t_min, t, t1, t2, a, b, c, disc, k0, k1, k2, F_att, I, E, xw, yw, zw, mag, dot_prod, aux_Kd, aux_Ke, aux_Km;
int shadow_intersection, reflex;

COLOR color_background = (COLOR) {.r=1, .g=1, .b=1};
COLOR color_pixel = (COLOR) {.r=0, .g=0, .b=0};
COLOR color_white = (COLOR) {.r=1, .g=1, .b=1};
COLOR color_black = (COLOR) {.r=0, .g=0, .b=0};
COLOR color_red = (COLOR) {.r=1, .g=0, .b=0};
COLOR color_green = (COLOR) {.r=0, .g=1, .b=0};
COLOR color_blue = (COLOR) {.r=0, .g=0, .b=1};
COLOR color_purple = (COLOR) {.r=1, .g=0, .b=1};
COLOR color_cyan = (COLOR) {.r=0, .g=1, .b=1};
COLOR color_yellow = (COLOR) {.r=1, .g=1, .b=0};
COLOR color_mirror = (COLOR) {.r=1, .g=1, .b=1};
COLOR color_aux = (COLOR) {.r=0, .g=0, .b=0};

VECTOR direction, ray, L, R, V, N, intersection, d, c_min, c_max, aux1, aux2;

//Functions /////////////////////////////////////
void ray_tracing();
COLOR first_intersection(VECTOR eye, VECTOR direction);
void ray_sphere_const (SPHERE_LIST_PTR ptr, VECTOR eye, VECTOR direction);
void ray_cone_const (CONE_LIST_PTR ptr, VECTOR eye, VECTOR direction, VECTOR c_dir);
void ray_cylinder_const (CYLINDER_LIST_PTR ptr, VECTOR eye, VECTOR direction, VECTOR cy_dir);
float ray_t_calculate (float a, float b, float c);
VECTOR ray_intersection (float t, VECTOR origin, VECTOR direction);
COLOR ray_point_color(float t, COLOR figure_color, float Kd, float Ke, float Km, VECTOR N);
float ray_difuse_lighting (VECTOR N, SPOTLIGHT_LIST_PTR ptr);
float ray_specular_lighting (VECTOR N, SPOTLIGHT_LIST_PTR ptr);
COLOR ray_mirror_lighting (VECTOR N, VECTOR direction);
int ray_shadow_intersection (SPOTLIGHT_LIST_PTR ptr);
COLOR ray_calculate_color (COLOR figure_color, COLOR light_color , float I , float E, float Km, COLOR mirror_color);

#endif
///////////////////////////////////////////////////////////////////////////////