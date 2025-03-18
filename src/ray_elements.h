///////////////////////////////////////////////////////////////////////////////
//
// Raytracing Window Header
//
///////////////////////////////////////////////////////////////////////////////

//Verification //////////////////////////////////
#ifndef RAY_ELEMENTS_H
#define RAY_ELEMENTS_H


//Includes //////////////////////////////////////
#include "window_module.h"
#include "vector.h"
#include "ray_tracing.h"

//Structs ///////////////////////////////////////
typedef struct color{
    float r;
    float g;
    float b;
}COLOR, *COLOR_PTR;

typedef struct window{
    VECTOR c_min, c_max;
}WINDOW, *WINDOW_PTR;

typedef struct sphere{
    VECTOR s_center;
    COLOR s_color;
    float R,Kd,Ke,Km;
}SPHERE, *SPHERE_PTR;

typedef struct sphere_list {
    SPHERE sphere;
    struct sphere_list *next;

}SPHERE_LIST, *SPHERE_LIST_PTR;

typedef struct plane{
    float A,B,C,D,Kd,Ke,Km;
    COLOR p_color;
}PLANE, *PLANE_PTR;

typedef struct plane_list {
    PLANE plane;
    struct plane_list *next;

}PLANE_LIST, *PLANE_LIST_PTR;

typedef struct cone{
    VECTOR c_anchor;
    VECTOR c_direction;
    COLOR c_color;
    float C1, C2, c_open,Kd,Ke,Km; // C1 = Axis step ; C2 = Radius step ; C2 = c_open*C1 ; Means that by every C1 units, the radius increse C2 units.
}CONE, *CONE_PTR;

typedef struct cone_list {
    CONE cone;
    struct cone_list *next;

}CONE_LIST, *CONE_LIST_PTR;

typedef struct cylinder{
    VECTOR cy_axis;
    VECTOR cy_direction;
    COLOR cy_color;
    float R,Kd,Ke,Km;
}CYLINDER, *CYLINDER_PTR;

typedef struct cylinder_list {
    CYLINDER cylinder;
    struct cylinder_list *next;

}CYLINDER_LIST, *CYLINDER_LIST_PTR;

typedef struct spotlight{
    VECTOR origin;
    float Is;
    COLOR color;
}SPOTLIGHT, *SPOTLIGHT_PTR;

typedef struct spotlight_list {
    SPOTLIGHT spotlight;
    struct spotlight_list *next;

}SPOTLIGHT_LIST, *SPOTLIGHT_LIST_PTR;

typedef struct illumination{
    SPOTLIGHT_LIST_PTR spots;
    float Ia;
    COLOR color;
}ILLUMINATION, *ILLUMINATION_PTR;

//Global Vars ///////////////////////////////////

SPHERE_LIST_PTR s_list;
PLANE_LIST_PTR p_list;
CONE_LIST_PTR c_list;
CYLINDER_LIST_PTR cy_list;
SPOTLIGHT_LIST_PTR spot_list;
ILLUMINATION_PTR illumination;
WINDOW ray_window;
VECTOR eye;

//Functions /////////////////////////////////////
void window_initialize (VECTOR corner_min, VECTOR corner_max);
void window_pan(float delta_x, float delta_y, float delta_z);
void window_zoom(float scale_factor);
void eye_initialize ();
void eye_pan(float delta_x, float delta_y, float delta_z);
void elements_create_scene(SPHERE_LIST_PTR *s_list, PLANE_LIST_PTR *p_list, CONE_LIST_PTR *c_list, CYLINDER_LIST_PTR *cy_list, SPOTLIGHT_LIST_PTR *spot_list);
PLANE elements_create_plane(VECTOR p0, VECTOR p1, VECTOR p2, COLOR color, float Kd, float Ke, float Km);
SPHERE_LIST_PTR elements_insert_sphere_last(SPHERE_LIST_PTR list, SPHERE sphere);
PLANE_LIST_PTR elements_insert_plane_last(PLANE_LIST_PTR list, PLANE plane);
CONE_LIST_PTR elements_insert_cone_last(CONE_LIST_PTR list, CONE cone);
CYLINDER_LIST_PTR elements_insert_cylinder_last(CYLINDER_LIST_PTR list, CYLINDER cylinder);
SPOTLIGHT_LIST_PTR elements_insert_spot_last(SPOTLIGHT_LIST_PTR list, SPOTLIGHT spot);
VECTOR elements_sphere_normal(SPHERE sphere, VECTOR intersection);
VECTOR elements_cone_normal(CONE_LIST_PTR ptr, VECTOR intersection);
VECTOR elements_cylinder_normal(CYLINDER_LIST_PTR ptr, VECTOR intersection, VECTOR cy_dir);


#endif
///////////////////////////////////////////////////////////////////////////////