#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;

    state.image_color= 0;
    state.image_depth=0;
    // std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;

    state.image_color = new pixel[width * height];
    for (unsigned i = 0; i < width * height; ++i){  //make all pixels black
        state.image_color[i] = make_pixel(0,0,0);
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    // std::cout<<"TODO: implement rendering."<<std::endl;
    data_geometry * data_geo = new data_geometry[state.num_vertices * state.floats_per_vertex];
    for (unsigned i = 0; i < state.num_vertices; i++){
        data_geo[i].data = new float[MAX_FLOATS_PER_VERTEX];
        data_vertex data_ver;
        data_ver.data = &state.vertex_data[i * state.floats_per_vertex];
        state.vertex_shader(data_ver, data_geo[i], state.uniform_data);
    }

    switch(type){
        case render_type::triangle:
            for(unsigned i = 0; i < state.num_vertices; i+= 3){
                rasterize_triangle(state, data_geo[i], data_geo[i+1], data_geo[i+2]);
            }
            break;
        case render_type::indexed:
            break;
        case render_type::fan:
            break;
        case render_type::strip:
            break;
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,v0,v1,v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    // std::cout<<"TODO: implement rasterization"<<std::endl;
    vec3 A = vec3(v0.gl_Position[0]/v0.gl_Position[3], v0.gl_Position[1]/v0.gl_Position[3], v0.gl_Position[2]/v0.gl_Position[3]);
    vec3 B = vec3(v1.gl_Position[0]/v1.gl_Position[3], v1.gl_Position[1]/v1.gl_Position[3], v1.gl_Position[2]/v1.gl_Position[3]);
    vec3 C = vec3(v2.gl_Position[0]/v2.gl_Position[3], v2.gl_Position[1]/v2.gl_Position[3], v2.gl_Position[2]/v2.gl_Position[3]);

    A[0] = (A[0] + 1 )/2.0 * (state.image_width -1);
    A[1] = (A[1] + 1 )/2.0 * (state.image_height -1);
    B[0] = (B[0] + 1 )/2.0 * (state.image_width -1);
    B[1] = (B[1] + 1 )/2.0 * (state.image_height -1);
    C[0] = (C[0] + 1 )/2.0 * (state.image_width -1);
    C[1] = (C[1] + 1 )/2.0 * (state.image_height -1);

    float areaABC = 0.5 * ( ( B[0] * C[1] - C[0] * B[1] ) + ( C[0] * A[1] - A[0] * C[1] ) + ( A[0] * B[1] - B[0] * A[1] ));
    
    float alpha, beta, gamma;
    for(unsigned i = 0; i < state.image_width; ++i){
        for(unsigned j = 0; j < state.image_height; ++j){
            vec3 P = vec3(i, j, 0);
            
            float areaPBC = 0.5 * ( ( B[0] * C[1] - C[0] * B[1] ) + ( C[0] * P[1] - P[0] * C[1] ) + ( P[0] * B[1] - B[0] * P[1] ) ); // might b fukt
            float areaAPC = 0.5 * ( ( P[0] * C[1] - C[0] * P[1] ) + ( C[0] * A[1] - A[0] * C[1] ) + ( A[0] * P[1] - P[0] * A[1] ) );
            float areaABP = 0.5 * ( ( B[0] * P[1] - P[0] * B[1] ) + ( P[0] * A[1] - A[0] * P[1] ) + ( A[0] * B[1] - B[0] * A[1] ) );


            // float alpha = (0.5 *((B[0]*C[1]-C[0]*B[1]) + (C[0]*P[1] - P[0]*C[1]) + (P[0]*B[1] - B[0]*P[1])) / areaABC);
            // float beta = (0.5 *((P[0]*C[1]-C[0]*P[1]) + (C[0]*A[1] - A[0]*C[1]) + (A[0]*P[1] - P[0]*A[1])) / areaABC);
            // float gamma = (0.5 *((B[0]*P[1]-P[0]*B[1]) + (P[0]*A[1] - A[0]*P[1]) + (A[0]*B[1] - B[0]*A[1])) / areaABC);

            alpha = areaPBC/areaABC;
            beta = areaAPC/areaABC;
            gamma = areaABP/areaABC;

            if(alpha >= 0 && beta >= 0 && gamma >= 0){
                state.image_color[i + j * state.image_width] = make_pixel(255,255,255);
            }
        }

    }
    
}

