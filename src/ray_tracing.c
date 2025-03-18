///////////////////////////////////////////////////////////////////////////////
//
// Ray Tracing Module
// Author: Adrián Rodríguez Murillo
//
///////////////////////////////////////////////////////////////////////////////

#include "ray_tracing.h"

///////////////////////////////////////////////////////////////////////////////
//
// ray_tracing -> Algorithm for ray tracing technique
//
///////////////////////////////////////////////////////////////////////////////
void ray_tracing(){

    // Initialize the ray tracing window
    c_min = (VECTOR) {.x=0, .y=0, .z=0};
    c_max = (VECTOR) {.x=200, .y=200, .z=0};
    window_initialize (c_min, c_max);

    //Optional operations for the window
    //window_zoom(2);
    //window_pan(100, 100, 0);
    
    //Initialize the ray tracing eye
    eye_initialize();
    
    //Optional operations for the eye
    //eye_pan(0, 50, 0);

    //Reseat the elements lists for prevention
    s_list = NULL;
    p_list = NULL;
    c_list = NULL;
    cy_list = NULL;
    spot_list = NULL;

    //Create the 3D scene
    elements_create_scene(&s_list, &p_list, &c_list, &cy_list, &spot_list);

    //Paint the buffer
    for (int i=0; i < x_res; i++){
        for(int j=0; j < y_res; j++){
            
            //Calculate the coordanates for the intersection between the ray and the window
            xw = ((i + 1/2)*(ray_window.c_max.x - ray_window.c_min.x) + ray_window.c_min.x)/x_res;
            yw = ((j + 1/2)*(ray_window.c_max.y - ray_window.c_min.y) + ray_window.c_min.y)/y_res;
            zw = ray_window.c_max.z;

            //Calculate the direction of the ray
            direction.x = xw - eye.x;
            direction.y = yw - eye.y;
            direction.z = zw - eye.z;
            direction = vector_norm(direction);

            //Calculate the color for every pixel
            reflex = 0;
            buffer[i][j] = first_intersection(eye, direction);
            //printf("color.x: %f, color.y:%f, color.z:%f\n", color_pixel.r, color_pixel.g, color_pixel.b);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//
// first_intersection -> Gets the first intersection for a ray
//
///////////////////////////////////////////////////////////////////////////////
COLOR first_intersection(VECTOR eye, VECTOR direction){
    t_min = INFINITY;

    //Intersection analysis for the sphere's list
    for (SPHERE_LIST_PTR ptr = s_list; ptr != NULL; ptr = ptr->next){
        // Calculate the minor intersection
        t = - INFINITY;
        ray_sphere_const (ptr, eye, direction);
        t = ray_t_calculate(a,b,c);
        if (t >= 0 && t < t_min) {
            t_min = t;
            intersection = ray_intersection(t_min,eye,direction);

            // Calculate the normal vector
            N = elements_sphere_normal(ptr->sphere, intersection);

            //Set the auxiliar variables
            color_aux = ptr->sphere.s_color;
            aux_Kd = ptr->sphere.Kd;
            aux_Ke = ptr->sphere.Ke;
            aux_Km = ptr->sphere.Km;
        }
    }

    //Intersection analysis for the plane's list
    for (PLANE_LIST_PTR ptr = p_list; ptr != NULL; ptr = ptr->next){
        //printf("Km = %f\n",ptr->plane.Km);
        //printf("Ray plano.x = %f, plano.y = %f, plano.z = %f\n",ptr->plane.A,ptr->plane.B,ptr->plane.C);
        // Calculate the minor intersection
        t = - INFINITY;
        t = -1*( ( ptr->plane.A*eye.x + ptr->plane.B*eye.y + ptr->plane.C*eye.z + ptr->plane.D ) / ( ptr->plane.A*direction.x + ptr->plane.B*direction.y + ptr->plane.C*direction.z ) );
        if (t >= 0 && t < t_min) {
            t_min = t;
            intersection = ray_intersection(t_min,eye,direction);

            // Calculate the normal vector
            N = (VECTOR) {.x = ptr->plane.A, .y = ptr->plane.B, .z = ptr->plane.C};
            N = vector_norm(N);
            
            //Set the auxiliar variables
            color_aux = ptr->plane.p_color;
            aux_Kd = ptr->plane.Kd;
            aux_Ke = ptr->plane.Ke;
            aux_Km = ptr->plane.Km;

        }
    }

    //Intersection analysis for the cone's list
    for (CONE_LIST_PTR ptr = c_list; ptr != NULL; ptr = ptr->next){
        // Calculate the minor intersection
        t = - INFINITY;
        VECTOR c_dir = vector_norm(ptr->cone.c_direction);
        ray_cone_const (ptr, eye, direction, c_dir);
        t = ray_t_calculate(a,b,c);      
        if (t >= 0 && t < t_min) {
            t_min = t;
            intersection = ray_intersection(t_min,eye,direction);

            // Calculate the normal vector
            N = elements_cone_normal(ptr, intersection);
            
            //Set the auxiliar variables
            color_aux = ptr->cone.c_color;
            aux_Kd = ptr->cone.Kd;
            aux_Ke = ptr->cone.Ke;
            aux_Km = ptr->cone.Km;
        }
    }

    //Intersection analysis for the cylinder's list
    for (CYLINDER_LIST_PTR ptr = cy_list; ptr != NULL; ptr = ptr->next){
        // Calculate the minor intersection
        t = - INFINITY;
        VECTOR cy_dir = vector_norm(ptr->cylinder.cy_direction);
        ray_cylinder_const (ptr, eye, direction, cy_dir);
        t = ray_t_calculate(a,b,c);
        if (t >= 0 && t < t_min) {
            t_min = t;
            intersection = ray_intersection(t_min,eye,direction);

            // Calculate the normal vector
            N = elements_cylinder_normal(ptr, intersection, cy_dir);

            //Set the auxiliar variables
            color_aux = ptr->cylinder.cy_color;
            aux_Kd = ptr->cylinder.Kd;
            aux_Ke = ptr->cylinder.Ke;
            aux_Km = ptr->cylinder.Km;
        }
    }

    //Calculate the point intensity and the final color
    color_pixel = ray_point_color(t_min, color_aux, aux_Kd, aux_Ke, aux_Km, N);
    if (t_min < INFINITY) return color_pixel;
    return color_background;
}

///////////////////////////////////////////////////////////////////////////////
//
// ray_sphere_const -> Calculate the a,b,c constants for t's calculate of the sphere
//
///////////////////////////////////////////////////////////////////////////////
void ray_sphere_const (SPHERE_LIST_PTR ptr, VECTOR eye, VECTOR direction){
    a = (pow(direction.x,2) + pow(direction.y,2) + pow(direction.z,2));
    b = 2*(eye.x - ptr->sphere.s_center.x)*direction.x + 2*(eye.y - ptr->sphere.s_center.y)*direction.y + 2*(eye.z - ptr->sphere.s_center.z)*direction.z;
    c = pow(eye.x - ptr->sphere.s_center.x , 2) + pow(eye.y - ptr->sphere.s_center.y , 2) + pow(eye.z - ptr->sphere.s_center.z , 2) - pow(ptr->sphere.R , 2);
}

///////////////////////////////////////////////////////////////////////////////
//
// ray_cone_const -> Calculate the a,b,c constants for t's calculate of the cone
//
///////////////////////////////////////////////////////////////////////////////
void ray_cone_const (CONE_LIST_PTR ptr, VECTOR eye, VECTOR direction, VECTOR c_dir){
    VECTOR diff_eye_anchor = (VECTOR) {.x = eye.x - ptr->cone.c_anchor.x, .y = eye.y - ptr->cone.c_anchor.y, .z = eye.z - ptr->cone.c_anchor.z};
    a = vector_dot_prod(direction, direction) - (1 + pow(ptr->cone.c_open,2)) * pow(vector_dot_prod(direction, c_dir), 2);
    b = 2 * (vector_dot_prod(direction, diff_eye_anchor) - (1 + pow(ptr->cone.c_open, 2)) * vector_dot_prod(direction, c_dir) * vector_dot_prod(diff_eye_anchor, c_dir));
    c = vector_dot_prod(diff_eye_anchor, diff_eye_anchor) - (1 + pow(ptr->cone.c_open, 2)) * pow(vector_dot_prod(diff_eye_anchor, c_dir), 2);
}

///////////////////////////////////////////////////////////////////////////////
//
// ray_cylinder_const -> Calculate the a,b,c constants for t's calculate of the cylinder
//
///////////////////////////////////////////////////////////////////////////////
void ray_cylinder_const (CYLINDER_LIST_PTR ptr, VECTOR eye, VECTOR direction, VECTOR cy_dir){
    a = (direction.x * direction.x + direction.y * direction.y + direction.z * direction.z) - pow(direction.x * cy_dir.x + direction.y * cy_dir.y + direction.z * cy_dir.z, 2);
    b = 2 * ((eye.x - ptr->cylinder.cy_axis.x) * direction.x + (eye.y - ptr->cylinder.cy_axis.y) * direction.y + (eye.z - ptr->cylinder.cy_axis.z) * direction.z - ((eye.x - ptr->cylinder.cy_axis.x) * cy_dir.x + (eye.y - ptr->cylinder.cy_axis.y) * cy_dir.y + (eye.z - ptr->cylinder.cy_axis.z) * cy_dir.z) * (direction.x * cy_dir.x + direction.y * cy_dir.y + direction.z * cy_dir.z));
    c = pow(eye.x - ptr->cylinder.cy_axis.x, 2) + pow(eye.y - ptr->cylinder.cy_axis.y, 2) + pow(eye.z - ptr->cylinder.cy_axis.z, 2) - pow(((eye.x - ptr->cylinder.cy_axis.x) * cy_dir.x + (eye.y - ptr->cylinder.cy_axis.y) * cy_dir.y + (eye.z - ptr->cylinder.cy_axis.z) * cy_dir.z), 2) - pow(ptr->cylinder.R, 2);
}

///////////////////////////////////////////////////////////////////////////////
//
// ray_t_calculate -> Determinates the t constant for an intersection
//
///////////////////////////////////////////////////////////////////////////////
float ray_t_calculate (float a, float b, float c){
    disc = pow(b,2)-4*a*c;
    if (disc == 0 ) {
        t1 = (-b/(2*a));
        if (t1 >= 0) t = t1;
    } else if (disc > 0){
        disc = sqrt(disc);
        t1 = (-b + disc)/ (2*a);
        t2 = (-b - disc)/ (2*a);
        if (t1 < t2 && t1 >= 0) t = t1;
        else if (t2 < t1 && t2 >= 0) t = t2;
        else if (t2 > t1 && t1 < 0 && t2 >= 0) t = t2;
        else if (t1 > t2 && t2 < 0 && t1 >= 0) t = t1;
    }
    return t;
}

///////////////////////////////////////////////////////////////////////////////
//
// ray_intersection -> Evaluate the t intersection in a ray function
//
///////////////////////////////////////////////////////////////////////////////
VECTOR ray_intersection (float t, VECTOR origin, VECTOR direction){
    VECTOR ray_int = (VECTOR) {.x = origin.x + t*direction.x, .y = origin.y + t*direction.y, .z = origin.z + t*direction.z};
    return ray_int;
}

///////////////////////////////////////////////////////////////////////////////
//
// ray_point_intensity -> Determinates the intensity of light in one point
//
///////////////////////////////////////////////////////////////////////////////
COLOR ray_point_color(float t, COLOR figure_color, float Kd, float Ke, float Km, VECTOR N){
    //Initialize the attenuation factors
    I = 0;
    E = 0;
    shadow_intersection = 0;
    k0 = 1;
    k1 = 0.00125;
    k2 = 0;
    COLOR color_aux = (COLOR) {0,0,0};
    COLOR color_out = (COLOR) {0,0,0};

    //Illumination analysis for the lightspot's list
    for (SPOTLIGHT_LIST_PTR ptr = spot_list; ptr != NULL; ptr = ptr->next){
        
        //Determinates the attenuation factor for each light
        d = (VECTOR) {.x = intersection.x - ptr->spotlight.origin.x, .y = intersection.y - ptr->spotlight.origin.y, .z = intersection.z - ptr->spotlight.origin.z};
        mag = vector_magnitude(d);
        F_att = 1 / (k0 + k1*mag + k2*pow(mag,2));

        //Inverts the normal vector of the planes so that it faces the light source
        VECTOR light_direction = {ptr->spotlight.origin.x - intersection.x, ptr->spotlight.origin.y - intersection.y, ptr->spotlight.origin.z - intersection.z};
        dot_prod = vector_dot_prod(N,light_direction);
        //if (dot_prod < 0) { N.x = -N.x; N.y = -N.y; N.z = -N.z;}

        //Calculation of individual contributions of each light
        float I = ray_difuse_lighting(N,ptr)*Kd*ptr->spotlight.Is*F_att;
        float E = ray_specular_lighting (N, ptr) *Ke*ptr->spotlight.Is*F_att;

        //Evaluates if necessary a ray of light due to shadow
        shadow_intersection = ray_shadow_intersection(ptr);

        //Generates mirros colors
        if (reflex < 2 && Km > 0) color_mirror = ray_mirror_lighting(N,direction);

        //Determinate final color 
        color_aux = ray_calculate_color(figure_color,ptr->spotlight.color, I , E, Km, color_mirror);
        shadow_intersection = 0;
        color_out = (COLOR) {color_out.r + color_aux.r, color_out.g + color_aux.g, color_out.b + color_aux.b};
        
    }
    if (color_out.r > 1) color_out.r = 1;
    if (color_out.g > 1) color_out.g = 1;
    if (color_out.b > 1) color_out.b = 1;
    return color_out;
    
}

///////////////////////////////////////////////////////////////////////////////
//
//  ray_specular_lighting -> Calculates E (specular coeficient)
//
///////////////////////////////////////////////////////////////////////////////
float ray_difuse_lighting (VECTOR N, SPOTLIGHT_LIST_PTR ptr){
    //Evaluates the Lambert's Law
    L.x = ptr->spotlight.origin.x - intersection.x;
    L.y = ptr->spotlight.origin.y - intersection.y;
    L.z = ptr->spotlight.origin.z - intersection.z;
    L = vector_norm(L);

    dot_prod = vector_dot_prod(L,N);
    dot_prod = (dot_prod < 0) ? 0 : dot_prod;

    return dot_prod;
}

///////////////////////////////////////////////////////////////////////////////
//
//  ray_specular_lighting -> Calculates E (specular coeficient)
//
///////////////////////////////////////////////////////////////////////////////
float ray_specular_lighting (VECTOR N, SPOTLIGHT_LIST_PTR ptr){
    L = (VECTOR) {.x = ptr->spotlight.origin.x - intersection.x, .y = ptr->spotlight.origin.y - intersection.y, .z = ptr->spotlight.origin.z - intersection.z};
    L = vector_norm(L);
    R = (VECTOR) {.x = (2*N.x*vector_dot_prod(N,L)) - L.x, .y = (2*N.y*vector_dot_prod(N,L)) - L.y, .z = (2*N.z*vector_dot_prod(N,L)) - L.z};
    V = (VECTOR) {.x = -direction.x, .y = -direction.y, .z = -direction.z};

    dot_prod = vector_dot_prod(V,R);
    dot_prod = (dot_prod < 0) ? 0 : dot_prod;

    return dot_prod;
}

///////////////////////////////////////////////////////////////////////////////
//
//  ray_mirror_lighting -> Calculates the mirror color 
//
///////////////////////////////////////////////////////////////////////////////
COLOR ray_mirror_lighting (VECTOR N, VECTOR direction){
    reflex++;
    dot_prod = vector_dot_prod(direction,N);
    R = (VECTOR) {.x = direction.x - 2*dot_prod*N.x, .y = direction.y - 2*dot_prod*N.y, .z = direction.z - 2*dot_prod*N.z};
    R = vector_norm(R);
    VECTOR shift_intersection = (VECTOR) {intersection.x + epsilon*R.x, intersection.y + epsilon*R.y, intersection.z + epsilon*R.z};
    color_mirror = first_intersection(shift_intersection,R);
    return color_mirror; 
}

///////////////////////////////////////////////////////////////////////////////
//
//  ray_shadow_intersection -> Calculates if exist an intersection in the 
//  ligth ray path
//
///////////////////////////////////////////////////////////////////////////////
int ray_shadow_intersection (SPOTLIGHT_LIST_PTR ptr){
    //Set variables for shadow analysis
    float t_shadow = INFINITY;
    shadow_intersection = 0;

    //Determinate new eye and direction vector from the intersection and light incoming direction
    VECTOR dir = (VECTOR) {ptr->spotlight.origin.x - intersection.x, ptr->spotlight.origin.y - intersection.y, ptr->spotlight.origin.z - intersection.z};
    dir = vector_norm(dir);
    VECTOR shift_intersection = (VECTOR) {intersection.x + epsilon*dir.x, intersection.y + epsilon*dir.y, intersection.z + epsilon*dir.z};

    //Forms intersections search
    for (SPHERE_LIST_PTR ptr = s_list; ptr != NULL; ptr = ptr->next){
        t = - INFINITY;
        ray_sphere_const (ptr, shift_intersection, dir);
        t = ray_t_calculate(a,b,c);
        if (t > 0 && t < t_shadow) return shadow_intersection = 1;
    }
    
    for (PLANE_LIST_PTR p_ptr = p_list; p_ptr != NULL; p_ptr = p_ptr->next){
        t = - INFINITY;
        t = -1*( ( p_ptr->plane.A*shift_intersection.x + p_ptr->plane.B*shift_intersection.y + p_ptr->plane.C*shift_intersection.z + p_ptr->plane.D ) / ( p_ptr->plane.A*dir.x + p_ptr->plane.B*dir.y + p_ptr->plane.C*dir.z ) );
        if (t >= 0 && t < t_min) {
            VECTOR shadow_intersection_vector = ray_intersection(t,shift_intersection,dir);
            VECTOR aux1 = (VECTOR) {ptr->spotlight.origin.x - intersection.x, ptr->spotlight.origin.y - intersection.y, ptr->spotlight.origin.z - intersection.z};
            VECTOR aux2 = (VECTOR) {shadow_intersection_vector.x - intersection.x, shadow_intersection_vector.y - intersection.y, shadow_intersection_vector.z - intersection.z};
            float v1 = vector_magnitude(aux1);
            float v2 = vector_magnitude(aux2);
            float delta = fabs(v2-v1);
            if ((delta > 0.1) && (v2 < v1)){
                return shadow_intersection = 1;
            } 
        }
    }
    for (CONE_LIST_PTR c_ptr = c_list; c_ptr != NULL; c_ptr = c_ptr->next){
        t = - INFINITY;
        VECTOR c_dir = vector_norm(c_ptr->cone.c_direction);
        ray_cone_const (c_ptr, shift_intersection, dir, c_dir);
        t = ray_t_calculate(a,b,c);      
        if (t >= 0 && t < t_min) {
            VECTOR shadow_intersection_vector = ray_intersection(t,shift_intersection,dir);
            VECTOR aux1 = (VECTOR) {ptr->spotlight.origin.x - intersection.x, ptr->spotlight.origin.y - intersection.y, ptr->spotlight.origin.z - intersection.z};
            VECTOR aux2 = (VECTOR) {shadow_intersection_vector.x - intersection.x, shadow_intersection_vector.y - intersection.y, shadow_intersection_vector.z - intersection.z};
            float v1 = vector_magnitude(aux1);
            float v2 = vector_magnitude(aux2);
            float delta = fabs(v2-v1);
            if ((delta > 0.1) && (v2 < v1)){
                return shadow_intersection = 1;
            } 
        }
    }
    for (CYLINDER_LIST_PTR ptr = cy_list; ptr != NULL; ptr = ptr->next){
        t = - INFINITY;
        VECTOR cy_dir = vector_norm(ptr->cylinder.cy_direction);
        ray_cylinder_const (ptr, shift_intersection, dir, cy_dir);
        t = ray_t_calculate(a,b,c);
        if (t > 0 && t < t_shadow) return shadow_intersection = 1;
    }
    return shadow_intersection;
}

///////////////////////////////////////////////////////////////////////////////
//
//  ray_calculate_light -> Calculates the total illumination for a pixel
//
///////////////////////////////////////////////////////////////////////////////
COLOR ray_calculate_color (COLOR pixel, COLOR light, float I, float E, float Km, COLOR mirror_color){
    // Ambient light applied to the color of the object
    //printf("I = %f, E = %f\n", I,E);
    COLOR ambient_color;
    ambient_color.r = illumination->Ia * pixel.r;
    ambient_color.g = illumination->Ia * pixel.g;
    ambient_color.b = illumination->Ia * pixel.b;

    // Diffused light applied to the color of the object
    COLOR diffuse_color;
    diffuse_color.r = I * pixel.r;
    diffuse_color.g = I * pixel.g;
    diffuse_color.b = I * pixel.b;

    // Interpolation of specular light between the color of the object (diffuse) and the color of the light
    COLOR specular_color;
    specular_color.r = (1.0 - E) * diffuse_color.r + E * light.r;
    specular_color.g = (1.0 - E) * diffuse_color.g + E * light.g;
    specular_color.b = (1.0 - E) * diffuse_color.b + E * light.b;

    // Interpolation between the mirror color of the object (reflex) and the color of the object in terms 
    mirror_color.r = (1.0 - Km) * specular_color.r + Km * mirror_color.r;
    mirror_color.g = (1.0 - Km) * specular_color.g + Km * mirror_color.g;
    mirror_color.b = (1.0 - Km) * specular_color.b + Km * mirror_color.b;

    // Add the components to obtain the final color
    COLOR color_aux2;
    if (shadow_intersection) {
        color_aux2.r = ambient_color.r;
        color_aux2.g = ambient_color.g;
        color_aux2.b = ambient_color.b;
    } else {
        color_aux2.r = ambient_color.r + mirror_color.r;
        color_aux2.g = ambient_color.g + mirror_color.g;
        color_aux2.b = ambient_color.b + mirror_color.b;
    }
    return color_aux2;
}

///////////////////////////////////////////////////////////////////////////////