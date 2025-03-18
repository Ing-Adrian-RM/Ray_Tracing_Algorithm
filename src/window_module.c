///////////////////////////////////////////////////////////////////////////////
//
// Window Module
// Author: Adrián Rodríguez Murillo
//
///////////////////////////////////////////////////////////////////////////////

#include "window_module.h"

///////////////////////////////////////////////////////////////////////////////
//
// window_ini_buffer -> Initializes the buffer and adds background color
//
///////////////////////////////////////////////////////////////////////////////

void window_ini_buffer(){
    int i,j;

    //Reserve mem space 
    buffer = (COLOR **)malloc(x_res *sizeof(COLOR*));
    for (i=0; i < x_res; i++) {
        buffer[i] = (COLOR *)malloc(y_res * sizeof(COLOR));
    }

    //Initialize  buffer values
    for (i=0; i < x_res; i++){
        for(j=0; j < y_res; j++){
            buffer[i][j] = color_background;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
//
// window_plot_loop -> Paint a dot from buffer to frame
//
///////////////////////////////////////////////////////////////////////////////

void window_plot_loop(){
    int i,j;

    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //Color each pixel
    for (i=0; i < x_res; i++){
        for(j=0; j < y_res; j++){
            glColor3f(buffer[i][j].r, buffer[i][j].g, buffer[i][j].b);
            glBegin(GL_POINTS);
            glVertex2i(i,j);
            glEnd();
        }
    } 
}

///////////////////////////////////////////////////////////////////////////////
//
// window_plot_scene -> Paint the desire scene in the buffer
//
///////////////////////////////////////////////////////////////////////////////

void window_plot_scene(){
    
    //auto inicio = std::chrono::high_resolution_clock::now(); // Stopwatch start

    ray_tracing();

    //auto fin = std::chrono::high_resolution_clock::now(); // End of stopwatch
    //auto duracion = std::chrono::duration_cast<std::chrono::microseconds>(fin - inicio).count(); // Calculates elapsed time in microseconds
    //std::cout << "T = " << duracion << " microsegundos." << std::endl; // Print the execution time in microseconds 
}

///////////////////////////////////////////////////////////////////////////////
//
// window_gl -> Paint the scene save in the buffer
//
///////////////////////////////////////////////////////////////////////////////

int window_gl_manager(){
  
    GLFWwindow* window; //Define window

    if ( !glfwInit()) fprintf (stderr, "No se puede inicializar el GLFW"); //Initialize glfw
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE); //No resizing

    window = glfwCreateWindow(x_res, y_res, "Ventana Adrian", NULL, NULL); //Create window
    if ( !window) {
        fprintf (stderr, "No se abrir ventana");
        glfwTerminate();
        exit (EXIT_FAILURE);
    }

    glfwMakeContextCurrent(window); //Assign the context
    glfwSwapInterval(1);

    if ( !gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) fprintf (stderr, "No se puede incluir GLAD"); //Include Glad

    glViewport(0, 0, x_res, y_res); //Set point of view
    glad_glLoadIdentity();
    glOrtho(0.0, x_res, 0.0, y_res, 0.0, 1.0); //Set canvas

    window_ini_buffer(); //Initialize buffer
    window_plot_scene(); //Plot

    //Main loop
    while(!glfwWindowShouldClose(window)) {
        //Paint
        window_plot_loop();
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate(); //Close window
    exit( EXIT_SUCCESS ); //Exit
}

///////////////////////////////////////////////////////////////////////////////