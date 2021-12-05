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
    state.image_width = width;
    state.image_height = height;

    state.image_color = 0;
    state.image_depth = 0;
    // std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;

    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];
    for (unsigned i = 0; i < width * height; ++i){  //make all pixels black
        state.image_color[i] = make_pixel(0,0,0);
        state.image_depth[i] = std::numeric_limits<float>::max();
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
            for(int i = 0; i < state.num_vertices; i+= 3){
                clip_triangle(state, data_geo[i], data_geo[i+1], data_geo[i+2], 0);
            }
            break;
        case render_type::indexed:
            for(int i = 0; i < state.num_triangles; i++){
                clip_triangle(state, data_geo[state.index_data[i*state.floats_per_vertex]], data_geo[state.index_data[i*state.floats_per_vertex + 1]], data_geo[state.index_data[i*state.floats_per_vertex + 2]], 0);
            }
            break;
        case render_type::fan:
            for(int i = 0; i < state.num_vertices - 2; i++){
                clip_triangle(state, data_geo[0], data_geo[i+1], data_geo[i+2], 0);
            }
            break;
        case render_type::strip:
            for(int i = 0; i < state.num_vertices - 2; i++){
                if(i % 2 == 0)
                    clip_triangle(state, data_geo[i], data_geo[i+1], data_geo[i+2], 0);
                else
                    clip_triangle(state, data_geo[i], data_geo[i+2], data_geo[i+1], 0);
            }
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

    float lamAB,lamBC, lamCA;
    float k;

    int planedimension; //indicates x, y, or z 
    int planeval; // indicates positive or negative

    bool a_out; // indicates v0 lies outside, and therefore should be clipped
    bool b_out; // v1
    bool c_out; // v2
        
    data_geometry p;
    data_geometry q;
    p.data = new float[MAX_FLOATS_PER_VERTEX];
    q.data = new float[MAX_FLOATS_PER_VERTEX];

    if(face==0){ // -x direction
        planedimension = 0;
        planeval = -1;
    }
    else if(face==1){  //x
        planedimension = 0;
        planeval = 1;
    }
    else if(face==2){   //-y
        planedimension = 1; 
        planeval = -1;
    }
    else if(face==3){   //y
        planedimension = 1;
        planeval = 1;
    }
    else if(face==4){   //-z
        planedimension = 2;
        planeval = -1;
    }
    else if(face==5){   //z
        planedimension = 2;
        planeval = 1;
    }
    else if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    else{
        exit(0);
    }


    //compute which of the vertices lie in or out of the cube
    if (planeval == -1) {   //if negative direction
        a_out = (v0.gl_Position[planedimension] < -v0.gl_Position[3])? 1: 0;
        b_out = (v1.gl_Position[planedimension] < -v1.gl_Position[3])? 1: 0;
        c_out = (v2.gl_Position[planedimension] < -v2.gl_Position[3])? 1: 0;
    } 
    else {  //if positive face
        a_out = (v0.gl_Position[planedimension] > v0.gl_Position[3])? 1: 0;
        b_out = (v1.gl_Position[planedimension] > v1.gl_Position[3])? 1: 0;
        c_out = (v2.gl_Position[planedimension] > v2.gl_Position[3])? 1: 0;
    }

    // call clipper differently depending on the case
    if(a_out && b_out && !c_out){   // c in, a b out
        lamBC = (v2.gl_Position[3] * planeval - v2.gl_Position[planedimension]) / (v1.gl_Position[planedimension] - v2.gl_Position[planedimension]);
        lamCA = (v2.gl_Position[3] * planeval - v2.gl_Position[planedimension]) / (v0.gl_Position[planedimension] - v2.gl_Position[planedimension]);

        p.gl_Position = v2.gl_Position + lamBC * (v1.gl_Position - v2.gl_Position);
        q.gl_Position = v2.gl_Position + lamCA * (v0.gl_Position - v2.gl_Position);

        for (unsigned i = 0; i < MAX_FLOATS_PER_VERTEX; ++i){
            if (state.interp_rules[i] == interp_type::noperspective){
                p.data[i] = v2.data[i] + lamBC * (v1.data[i] - v2.data[i]);
                q.data[i] = v2.data[i] + lamCA * (v0.data[i] - v2.data[i]);
            }

            else if (state.interp_rules[i] == interp_type::smooth){
                k = lamBC / v2.gl_Position[3] + (1 - lamBC) * v1.gl_Position[3];
                p.data[i] = v2.data[i] + lamBC/(v2.gl_Position[3] * k) * (v1.data[i] - v2.data[i]);
                k = lamCA / v2.gl_Position[3] + (1 - lamCA) * v0.gl_Position[3];
                q.data[i] = v2.data[i] + lamCA/(v2.gl_Position[3] * k) * (v0.data[i] - v2.data[i]);
            }
        }

        clip_triangle(state, v2, q, p, face + 1);

    }

    else if(a_out && !b_out && c_out){  // b in, a c out
        lamAB = (v1.gl_Position[3] * planeval - v1.gl_Position[planedimension]) / (v0.gl_Position[planedimension] - v1.gl_Position[planedimension]);
        lamBC = (v1.gl_Position[3] * planeval - v1.gl_Position[planedimension]) / (v2.gl_Position[planedimension] - v1.gl_Position[planedimension]);

        p.gl_Position = v1.gl_Position + lamAB * (v0.gl_Position - v1.gl_Position);
        q.gl_Position = v1.gl_Position + lamBC * (v2.gl_Position - v1.gl_Position);

        for (unsigned i = 0; i < MAX_FLOATS_PER_VERTEX; ++i){
            if (state.interp_rules[i] == interp_type::noperspective){
                p.data[i] = v1.data[i] + lamAB * (v0.data[i] - v1.data[i]);
                q.data[i] = v1.data[i] + lamBC * (v2.data[i] - v1.data[i]);
            }

            else if (state.interp_rules[i] == interp_type::smooth){
                k = lamAB / v1.gl_Position[3] + (1 - lamAB) * v0.gl_Position[3];
                p.data[i] = v1.data[i] + lamAB/(v1.gl_Position[3] * k) * (v0.data[i] - v1.data[i]);
                k = lamBC / v1.gl_Position[3] + (1 - lamBC) * v2.gl_Position[3];
                q.data[i] = v1.data[i] + lamBC/(v1.gl_Position[3] * k) * (v2.data[i] - v1.data[i]);
            }
        }

        clip_triangle(state, v1, q, p, face + 1);
    }

    else if(!a_out && b_out && c_out){  // a in, b c out
        lamAB = (v0.gl_Position[3] * planeval - v0.gl_Position[planedimension]) / (v1.gl_Position[planedimension] - v0.gl_Position[planedimension]);
        lamCA = (v0.gl_Position[3] * planeval - v0.gl_Position[planedimension]) / (v2.gl_Position[planedimension] - v0.gl_Position[planedimension]);

        p.gl_Position = v0.gl_Position + lamCA * (v2.gl_Position - v0.gl_Position);
        q.gl_Position = v0.gl_Position + lamAB * (v1.gl_Position - v0.gl_Position);

        for (unsigned i = 0; i < MAX_FLOATS_PER_VERTEX; ++i){
            if (state.interp_rules[i] == interp_type::noperspective){
                p.data[i] = v0.data[i] + lamCA * (v2.data[i] - v0.data[i]);
                q.data[i] = v0.data[i] + lamAB * (v1.data[i] - v0.data[i]);
            }

            else if (state.interp_rules[i] == interp_type::smooth){
                k = lamCA / v0.gl_Position[3] + (1 - lamCA) * v2.gl_Position[3];
                p.data[i] = v0.data[i] + lamCA/(v0.gl_Position[3] * k) * (v2.data[i] - v0.data[i]);
                k = lamAB / v0.gl_Position[3] + (1 - lamAB) * v1.gl_Position[3];
                q.data[i] = v0.data[i] + lamAB/(v0.gl_Position[3] * k) * (v1.data[i] - v0.data[i]);
            }
        }

        clip_triangle(state, v0, q, p, face + 1);
    }

    else if(a_out && !b_out && !c_out){ // b c in, a out
        lamAB = (v1.gl_Position[3]*planeval - v1.gl_Position[planedimension]) / (v0.gl_Position[planedimension] - v1.gl_Position[planedimension]);
        lamCA = (v2.gl_Position[3]*planeval - v2.gl_Position[planedimension]) / (v0.gl_Position[planedimension] - v2.gl_Position[planedimension]);

        p.gl_Position = v1.gl_Position + lamAB * (v0.gl_Position - v1.gl_Position);
        q.gl_Position = v2.gl_Position + lamCA * (v0.gl_Position - v2.gl_Position);

        for (unsigned i = 0; i < MAX_FLOATS_PER_VERTEX; ++i){
            if (state.interp_rules[i] == interp_type::noperspective){
                p.data[i] = v1.data[i] + lamAB * (v0.data[i] - v1.data[i]);
                q.data[i] = v2.data[i] + lamCA * (v0.data[i] - v2.data[i]);
            }

            else if (state.interp_rules[i] == interp_type::smooth){
                k = lamAB / v1.gl_Position[3] + (1-lamAB) / v0.gl_Position[3];
                float lamABpersp = lamAB / (v1.gl_Position[3] * k);
                p.data[i] = v1.data[i] + lamABpersp * (v0.data[i] - v1.data[i]);
                k = lamCA / v2.data[3] + (1 - lamCA) / v0.data[3];
                float lamCApersp = lamCA / (v2.gl_Position[3] * k);
                q.data[i] = v2.data[i] + lamCApersp * (v0.data[i] - v2.data[i]);
            }
        }

        clip_triangle(state, v1, v2, p, face + 1);
        clip_triangle(state, p, v2, q, face + 1);
    }

    else if(!a_out && !b_out && c_out){ // a b in, c out
        lamCA = (v0.gl_Position[3]*planeval - v0.gl_Position[planedimension]) / (v2.gl_Position[planedimension] - v0.gl_Position[planedimension]);
        lamBC = (v1.gl_Position[3]*planeval - v1.gl_Position[planedimension]) / (v2.gl_Position[planedimension] - v1.gl_Position[planedimension]);

        p.gl_Position = v0.gl_Position + lamCA * (v2.gl_Position - v0.gl_Position);
        q.gl_Position = v1.gl_Position + lamBC * (v2.gl_Position - v1.gl_Position);

        for (unsigned i = 0; i < MAX_FLOATS_PER_VERTEX; ++i){
            if (state.interp_rules[i] == interp_type::noperspective){
                p.data[i] = v0.data[i] + lamCA * (v2.data[i] - v0.data[i]);
                q.data[i] = v1.data[i] + lamBC * (v2.data[i] - v1.data[i]);
            }

            else if (state.interp_rules[i] == interp_type::smooth){
                k = lamCA / v0.data[3] + (1-lamCA) / v2.data[3];
                float lamCApersp = lamCA / (v0.gl_Position[3] * k);
                p.data[i] = v0.data[i] + lamCApersp * (v2.data[i] - v0.data[i]);
                k = lamBC / v1.data[3] + (1 - lamBC) / v2.data[3];
                float lamBCpersp = lamBC / (v1.gl_Position[3] * k);
                q.data[i] = v1.data[i] + lamBCpersp * (v2.data[i] - v1.data[i]);
            }

        }
        clip_triangle(state, v0, v1, p, face+1);
        clip_triangle(state, p, v1, q, face+1);
    }

    else if(!a_out && b_out && !c_out){ // a c in, b out
        lamBC = (v2.gl_Position[3]*planeval - v2.gl_Position[planedimension]) / (v1.gl_Position[planedimension] - v2.gl_Position[planedimension]);
        lamAB = (v0.gl_Position[3]*planeval - v0.gl_Position[planedimension]) / (v1.gl_Position[planedimension] - v0.gl_Position[planedimension]);

        p.gl_Position = v2.gl_Position + lamBC * (v1.gl_Position - v2.gl_Position);
        q.gl_Position = v0.gl_Position + lamAB * (v1.gl_Position - v0.gl_Position);

        for (unsigned i = 0; i < MAX_FLOATS_PER_VERTEX; ++i){
            if (state.interp_rules[i] == interp_type::noperspective){
                p.data[i] = v2.data[i] + lamBC * (v1.data[i] - v2.data[i]);
                q.data[i] = v0.data[i] + lamAB * (v1.data[i] - v0.data[i]);
            }

            else if (state.interp_rules[i] == interp_type::smooth){
                k = lamBC / v2.data[3] + (1 - lamBC) / v1.data[3];
                float lamBCpersp = lamBC / (v2.gl_Position[3] * k);
                p.data[i] = v2.data[i] + lamBCpersp * (v1.data[i] - v2.data[i]);
                k = lamAB / v0.data[3] + (1 - lamAB) / v1.data[3];
                float lamABpersp = lamAB / (v0.gl_Position[3] * k);
                q.data[i] = v0.data[i] + lamABpersp * (v1.data[i] - v0.data[i]);
            }

        }

        clip_triangle(state, v2, v0, p, face+1);
        clip_triangle(state, p, v0, q, face+1);
    }

    else if(!a_out && !b_out && !c_out){ // all points lie inside
        clip_triangle(state, v0, v1, v2, face + 1);
    }


    if (p.data != nullptr)
        delete[] p.data;
    if (q.data != nullptr)
        delete[] q.data;
    // std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    // clip_triangle(state,v0,v1,v2,face+1);
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
                // if((MAX_FLOATS_PER_VERTEX) == 3){
                //     if(interp_rules[2] == flat){
                //         state.image_color[i + j * state.image_width] = make_pixel(v0.gl_Position[3]);
                //     }
                //     if(interp_rules[2] == smooth){  //check one vertex interp rules
                //         state.image_color[i + j * state.image_width] = make_pixel();
                //     }
                //     if(interp_rules[2] == noperspective){
                //         state.image_color[i + j * state.image_width] = make_pixel();
                //     }
                // }
                P[2] = alpha * A[2] + beta * B[2] + gamma * C[2] + 1;
                data_fragment frag_color;
                frag_color.data = new float[MAX_FLOATS_PER_VERTEX];
                data_output p_color;
                state.fragment_shader(frag_color, p_color, state.uniform_data);

                int r_val = p_color.output_color[0] * 255;
                int g_val = p_color.output_color[1] * 255;
                int b_val = p_color.output_color[2] * 255;

                

                if(state.floats_per_vertex > 3){
                    float k_val = (alpha / v0.gl_Position[3]) + (beta / v1.gl_Position[3]) + (gamma / v2.gl_Position[3]);

                    switch(state.interp_rules[3]){
                        case interp_type::flat:
                            r_val = v0.data[3]*255;
                            g_val = v0.data[4]*255;
                            b_val = v0.data[5]*255;
                        break;
                        case interp_type::smooth:  //check one vertex interp rules
                            alpha /= k_val * v0.gl_Position[3];
                            beta  /= k_val * v1.gl_Position[3];
                            gamma /= k_val * v2.gl_Position[3];


                            r_val = (v0.data[3]*alpha + v1.data[3]*beta + v2.data[3]*gamma)*255;
                            g_val = (v0.data[4]*alpha + v1.data[4]*beta + v2.data[4]*gamma)*255;
                            b_val = (v0.data[5]*alpha + v1.data[5]*beta + v2.data[5]*gamma)*255;
                        break;
                        case interp_type::noperspective:
                            r_val = (v0.data[3]*alpha + v1.data[3]*beta + v2.data[3]*gamma)*255;
                            g_val = (v0.data[4]*alpha + v1.data[4]*beta + v2.data[4]*gamma)*255; 
                            b_val = (v0.data[5]*alpha + v1.data[5]*beta + v2.data[5]*gamma)*255;
                        break;
                      
                    }
                }


                if(P[2] < state.image_depth[state.image_width *j + i]){
                    state.image_color[i + j * state.image_width] = make_pixel(r_val, g_val, b_val);
                    state.image_depth[state.image_width *j + i]= P[2];
                }

                delete[] frag_color.data;

            }

        }
    
    }

}