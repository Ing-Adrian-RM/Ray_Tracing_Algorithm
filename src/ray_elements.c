///////////////////////////////////////////////////////////////////////////////
//
// Raytracing Window Module
// Author: Adrián Rodríguez Murillo
//
///////////////////////////////////////////////////////////////////////////////

#include "ray_elements.h"

///////////////////////////////////////////////////////////////////////////////
//
// window_initialize -> Create a window for raytracing process.
//
///////////////////////////////////////////////////////////////////////////////
void window_initialize (VECTOR c_min, VECTOR c_max){
    ray_window.c_max.x = c_max.x;
    ray_window.c_max.y = c_max.y;
    ray_window.c_max.z = c_max.z;
    ray_window.c_min.x = c_min.x;
    ray_window.c_min.y = c_min.y;
    ray_window.c_min.z = c_min.z;

}

///////////////////////////////////////////////////////////////////////////////
//
// window_pan -> Traslate de window 
//
///////////////////////////////////////////////////////////////////////////////
void window_pan(float delta_x, float delta_y, float delta_z){
    ray_window.c_min.x = (int)(ray_window.c_min.x + delta_x);
    ray_window.c_min.y = (int)(ray_window.c_min.y + delta_y);
    ray_window.c_min.z = (int)(ray_window.c_min.z + delta_z);

    ray_window.c_max.x = (int)(ray_window.c_max.x + delta_x);
    ray_window.c_max.y = (int)(ray_window.c_max.y + delta_y);
    ray_window.c_max.z = (int)(ray_window.c_max.z + delta_z);
}

///////////////////////////////////////////////////////////////////////////////
//
// window_zoom -> Changes de window`s size
//
///////////////////////////////////////////////////////////////////////////////
void window_zoom(float scale_factor){
    float centro_x = (ray_window.c_max.x + ray_window.c_min.x)/2;
    float centro_y = (ray_window.c_max.y + ray_window.c_min.y)/2;
    //Transfer to origin
    ray_window.c_max.x -= centro_x;
    ray_window.c_min.x -= centro_x;
    ray_window.c_max.y -= centro_y;
    ray_window.c_min.y -= centro_y;

    //Scale
    float x_max_new = ray_window.c_max.x * scale_factor;
    float x_min_new = ray_window.c_min.x * scale_factor;
    float y_max_new = ray_window.c_max.y * scale_factor;
    float y_min_new = ray_window.c_min.y * scale_factor;

    //Move back to original center
    ray_window.c_max.x = (int)(x_max_new + centro_x);
    ray_window.c_min.x = (int)(x_min_new + centro_x);
    ray_window.c_max.y = (int)(y_max_new + centro_y);
    ray_window.c_min.y = (int)(y_min_new + centro_y);
}

///////////////////////////////////////////////////////////////////////////////
//
// eye_initialize -> Create a eye for raytracing process.
//
///////////////////////////////////////////////////////////////////////////////
void eye_initialize (){
    float cMx = (float) ray_window.c_max.x;
    float cMy = (float) ray_window.c_max.y;
    float cmx = (float) ray_window.c_min.x;
    float cmy = (float) ray_window.c_min.y;
    float z = -50;
    eye.x = (cMx + cmx) / 2;
    eye.y = (cMy + cmy) / 2;
    eye.z = z;
}

///////////////////////////////////////////////////////////////////////////////
//
// eye_pan -> Traslate de ray tracing's eye 
//
///////////////////////////////////////////////////////////////////////////////
void eye_pan(float delta_x, float delta_y, float delta_z){
    eye.x = (int)(eye.x + delta_x);
    eye.y = (int)(eye.y + delta_y);
    eye.z = (int)(eye.z + delta_z);
}

///////////////////////////////////////////////////////////////////////////////
//
// element_create_scene-> Initialize the 3D objects in the scene
//
///////////////////////////////////////////////////////////////////////////////
void elements_create_scene(SPHERE_LIST_PTR *s_list, PLANE_LIST_PTR *p_list, CONE_LIST_PTR *c_list, CYLINDER_LIST_PTR *cy_list, SPOTLIGHT_LIST_PTR *spot_list){
    float global_Ke = 1;

    //Initialize the spotlights configuration and add them in the spotligh list
    VECTOR spot_1_origin = (VECTOR) {.x = -50, .y =150, .z = -25};
    SPOTLIGHT spot_1 = (SPOTLIGHT) {.origin = spot_1_origin, .Is = 0.2, .color = color_white};
    //*spot_list = elements_insert_spot_last(*spot_list, spot_1);

    VECTOR spot_2_origin = (VECTOR) {.x = 250, .y = 100, .z = -25};
    SPOTLIGHT spot_2 = (SPOTLIGHT) {.origin = spot_2_origin, .Is = 0.5, .color = color_white};
    *spot_list = elements_insert_spot_last(*spot_list, spot_2);

    //Initialize the illumination configuration
    illumination = (ILLUMINATION_PTR) malloc(sizeof(ILLUMINATION));
    illumination->spots = *spot_list;
    illumination->Ia = 0.2;

    //Initialize the spheres and add them in the spheres list
    SPHERE sphere_1;
    sphere_1.s_center = (VECTOR) {.x = 100, .y = 100, .z = 50};
    sphere_1.s_color = color_red;
    sphere_1.R = 30;
    sphere_1.Kd = 1;
    sphere_1.Ke = global_Ke;
    sphere_1.Km = 0;
    *s_list = elements_insert_sphere_last(*s_list, sphere_1);

    //Initialize the planes and add them in the planes list
    VECTOR p0 = (VECTOR) {.x = 0, .y = 0, .z = 500};
    VECTOR p1 = (VECTOR) {.x = 0, .y = 200, .z = 500};
    VECTOR p2 = (VECTOR) {.x = 200, .y = 0, .z = 500};
    PLANE plane_1 = elements_create_plane(p0, p1, p2, color_red, 1, global_Ke, 0);
    *p_list = elements_insert_plane_last(*p_list, plane_1);

    p0 = (VECTOR) {.x = 0, .y = 0, .z = 0};
    p1 = (VECTOR) {.x = 200, .y = 0, .z = 0};
    p2 = (VECTOR) {.x = 0, .y = 0, .z = 200};
    PLANE plane_2 = elements_create_plane(p0, p1, p2, color_purple, 1, global_Ke, 0);
    *p_list = elements_insert_plane_last(*p_list, plane_2);

    p0 = (VECTOR) {.x = -50, .y = 0, .z = 0};
    p1 = (VECTOR) {.x = -50, .y = 200, .z = 0};
    p2 = (VECTOR) {.x = -50, .y = 0, .z = 200};
    PLANE plane_3 = elements_create_plane(p0, p1, p2, color_cyan, 1, global_Ke, 1);
    *p_list = elements_insert_plane_last(*p_list, plane_3);

    p0 = (VECTOR) {.x = 250, .y = 0, .z = 0};
    p1 = (VECTOR) {.x = 250, .y = 200, .z = 0};
    p2 = (VECTOR) {.x = 250, .y = 0, .z = 200};
    PLANE plane_4 = elements_create_plane(p0, p1, p2, color_yellow, 1, global_Ke, 1);
    *p_list = elements_insert_plane_last(*p_list, plane_4);

    //Initialize the cones and add them in the cones list
    VECTOR c_anchor = (VECTOR) {.x = 100, .y = 100, .z = 50};
    VECTOR c_direction = (VECTOR) {.x = 0, .y = 1, .z = 0};
    float C1 = 1;
    float C2 = 0.3;
    CONE cone_1 = (CONE) {.c_anchor = c_anchor, .c_direction = c_direction, .c_color = color_blue, .C1 = C1, .C2 = C2, .c_open = C2/C1, .Kd = 1, .Ke = global_Ke, 0};
    *c_list = elements_insert_cone_last(*c_list, cone_1);

    //Initialize the cylinders configuration and add them in the cylinders list
    VECTOR cy_axis = (VECTOR) {.x = 0, .y = 0, .z = 150};
    VECTOR cy_direction = (VECTOR) {.x = 0, .y = 1, .z = 0};
    CYLINDER cylinder_1 = (CYLINDER) {.cy_axis = cy_axis, .cy_direction = cy_direction, .cy_color = color_green, .R = 50, .Kd = 1, .Ke = global_Ke, 0};
    *cy_list = elements_insert_cylinder_last(*cy_list, cylinder_1);

    cy_axis = (VECTOR) {.x = 200, .y = 0, .z = 75};
    cy_direction = (VECTOR) {.x = 0, .y = 1, .z = 0};
    CYLINDER cylinder_2 = (CYLINDER) {.cy_axis = cy_axis, .cy_direction = cy_direction, .cy_color = color_green, .R = 20, .Kd = 1, .Ke = global_Ke, 0};
    *cy_list = elements_insert_cylinder_last(*cy_list, cylinder_2);
}

///////////////////////////////////////////////////////////////////////////////
//
// element_create_plane-> Create one plane from two points
//
///////////////////////////////////////////////////////////////////////////////
PLANE elements_create_plane(VECTOR p0, VECTOR p1, VECTOR p2, COLOR color, float Kd, float Ke, float Km){
    PLANE plane;

    VECTOR v1 = (VECTOR){.x = p2.x - p1.x, .y = p2.y - p1.y, .z = p2.z - p1.z};
    VECTOR v2 = (VECTOR){.x = p0.x - p1.x, .y = p0.y - p1.y, .z = p0.z - p1.z};
    VECTOR N = vector_cross_prod(v1, v2);

    if (N.x == 0 && N.y == 0 && N.z == 0) printf("Invalids points for plane. \n");

    plane.A = N.x;
    plane.B = N.y;
    plane.C = N.z;
    plane.D = -1*(plane.A*p0.x + plane.B*p0.y + plane.C*p0.z);
    plane.Kd = Kd;
    plane.Ke = Ke;
    plane.Km = Km;
    plane.p_color = color;    
    return plane;
}

///////////////////////////////////////////////////////////////////////////////
//
// elements_insert_sphere_last-> Insert a sphere in a sphere list
//
///////////////////////////////////////////////////////////////////////////////
SPHERE_LIST_PTR elements_insert_sphere_last(SPHERE_LIST_PTR list, SPHERE sphere){
    SPHERE_LIST_PTR new_node = (SPHERE_LIST_PTR) malloc(sizeof(SPHERE_LIST));
    new_node->sphere = sphere;
    new_node->next = NULL;

    if (list == NULL) {
        return new_node;
    } else {
        for (SPHERE_LIST_PTR ptr = list; ptr != NULL; ptr = ptr->next){
            if (ptr->next == NULL) {
                ptr->next = new_node;
                return list;
            }
        }
    }
    return list;
}

///////////////////////////////////////////////////////////////////////////////
//
// elements_insert_plane_last-> Insert a plane in a plane list
//
///////////////////////////////////////////////////////////////////////////////
PLANE_LIST_PTR elements_insert_plane_last(PLANE_LIST_PTR list, PLANE plane){
    PLANE_LIST_PTR new_node = (PLANE_LIST_PTR) malloc(sizeof(PLANE_LIST));
    new_node->plane = plane;
    new_node->next = NULL;

    if (list == NULL) {
        return new_node;
    } else {
        for (PLANE_LIST_PTR ptr = list; ptr != NULL; ptr = ptr->next){
            if (ptr->next == NULL) {
                ptr->next = new_node;
                return list;
            }
        }
    }
    return list;
}

///////////////////////////////////////////////////////////////////////////////
//
// elements_insert_cone_last-> Insert a cone in a cone list
//
///////////////////////////////////////////////////////////////////////////////
CONE_LIST_PTR elements_insert_cone_last(CONE_LIST_PTR list, CONE cone){
    CONE_LIST_PTR new_node = (CONE_LIST_PTR) malloc(sizeof(CONE_LIST));
    new_node->cone = cone;
    new_node->next = NULL;

    if (list == NULL) {
        return new_node;
    } else {
        for (CONE_LIST_PTR ptr = list; ptr != NULL; ptr = ptr->next){
            if (ptr->next == NULL) {
                ptr->next = new_node;
                return list;
            }
        }
    }
    return list;
}

///////////////////////////////////////////////////////////////////////////////
//
// elements_insert_cylinder_last-> Insert a cylinder in a cylinder list
//
///////////////////////////////////////////////////////////////////////////////
CYLINDER_LIST_PTR elements_insert_cylinder_last(CYLINDER_LIST_PTR list, CYLINDER cylinder){
    CYLINDER_LIST_PTR new_node = (CYLINDER_LIST_PTR) malloc(sizeof(CYLINDER_LIST));
    new_node->cylinder = cylinder;
    new_node->next = NULL;

    if (list == NULL) {
        return new_node;
    } else {
        for (CYLINDER_LIST_PTR ptr = list; ptr != NULL; ptr = ptr->next){
            if (ptr->next == NULL) {
                ptr->next = new_node;
                return list;
            }
        }
    }
    return list;
}

///////////////////////////////////////////////////////////////////////////////
//
// elements_insert_spot_last-> Insert a spot in a spot list
//
///////////////////////////////////////////////////////////////////////////////
SPOTLIGHT_LIST_PTR elements_insert_spot_last(SPOTLIGHT_LIST_PTR list, SPOTLIGHT spot){
    SPOTLIGHT_LIST_PTR new_node = (SPOTLIGHT_LIST_PTR) malloc(sizeof(SPOTLIGHT_LIST));
    new_node->spotlight = spot;
    new_node->next = NULL;

    if (list == NULL) {
        return new_node;
    } else {
        for (SPOTLIGHT_LIST_PTR ptr = list; ptr != NULL; ptr = ptr->next){
            if (ptr->next == NULL) {
                ptr->next = new_node;
                return list;
            }
        }
    }
    return list;
}

///////////////////////////////////////////////////////////////////////////////
//
// element_sphere_normal-> Calculate the normal vector of the sphere
//
///////////////////////////////////////////////////////////////////////////////
VECTOR elements_sphere_normal(SPHERE sphere, VECTOR intersection){
    N.x = (intersection.x - sphere.s_center.x) / sphere.R;
    N.y = (intersection.y - sphere.s_center.y) / sphere.R;
    N.z = (intersection.z - sphere.s_center.z) / sphere.R;
    return N;
}

///////////////////////////////////////////////////////////////////////////////
//
// element_cone_normal-> Calculate the normal vector of the cone
//
///////////////////////////////////////////////////////////////////////////////
VECTOR elements_cone_normal(CONE_LIST_PTR ptr, VECTOR intersection){
    // Vector from vertex to intersection point
    float V_p_x = intersection.x - ptr->cone.c_anchor.x;
    float V_p_y = intersection.y - ptr->cone.c_anchor.y;
    float V_p_z = intersection.z - ptr->cone.c_anchor.z;

    // Projection on the axis of the cone
    float dot_V_p_axis = V_p_x * ptr->cone.c_direction.x + V_p_y * ptr->cone.c_direction.y + V_p_z * ptr->cone.c_direction.z;

    float V_p_parallel_x = dot_V_p_axis * ptr->cone.c_direction.x;
    float V_p_parallel_y = dot_V_p_axis * ptr->cone.c_direction.y;
    float V_p_parallel_z = dot_V_p_axis * ptr->cone.c_direction.z;

    // Paso 3: Perpendicular component
    float V_p_perp_x = V_p_x - V_p_parallel_x;
    float V_p_perp_y = V_p_y - V_p_parallel_y;
    float V_p_perp_z = V_p_z - V_p_parallel_z;

    // Paso 4: Adjustment by cone angle
    N.x = V_p_perp_x - ptr->cone.c_open * V_p_parallel_x;
    N.y = V_p_perp_y - ptr->cone.c_open * V_p_parallel_y;
    N.z = V_p_perp_z - ptr->cone.c_open * V_p_parallel_z;

    // Paso 5: Normal Vector Normalization
    N = vector_norm(N);
    return N;
}

///////////////////////////////////////////////////////////////////////////////
//
// element_cylinder_normal-> Calculate the normal vector of the cylinder
//
///////////////////////////////////////////////////////////////////////////////
VECTOR elements_cylinder_normal(CYLINDER_LIST_PTR ptr, VECTOR intersection, VECTOR cy_dir){
    float aux = ((intersection.x - ptr->cylinder.cy_axis.x) * cy_dir.x + (intersection.y - ptr->cylinder.cy_axis.y) * cy_dir.y + (intersection.z - ptr->cylinder.cy_axis.z) * cy_dir.z);
    N.x = (intersection.x - ptr->cylinder.cy_axis.x) - aux * cy_dir.x;
    N.y = (intersection.y - ptr->cylinder.cy_axis.y) - aux * cy_dir.y;
    N.z = (intersection.z - ptr->cylinder.cy_axis.z) - aux * cy_dir.z;

    // Paso 5: Normal Vector Normalization
    N = vector_norm(N);
    return N;
}
///////////////////////////////////////////////////////////////////////////////