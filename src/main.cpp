// C++ include
#include <iostream>
#include <string>
#include <vector>

// Utilities for the Assignment
#include "utils.h"

// Image writing library
#define STB_IMAGE_WRITE_IMPLEMENTATION // Do not include this line twice in your project!
#include "stb_image_write.h"

// Shortcut to avoid Eigen:: everywhere, DO NOT USE IN .h
using namespace Eigen;

void raytrace_sphere()
{
    std::cout << "Simple ray tracer, one sphere with orthographic projection" << std::endl;

    const std::string filename("sphere_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the
    // unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // Intersect with the sphere
            // NOTE: this is a special case of a sphere centered in the origin and for orthographic rays aligned with the z axis
            Vector2d ray_on_xy(ray_origin(0), ray_origin(1));
            const double sphere_radius = 0.9;

            if (ray_on_xy.norm() < sphere_radius)
            {
                // The ray hit the sphere, compute the exact intersection point
                Vector3d ray_intersection(
                    ray_on_xy(0), ray_on_xy(1),
                    sqrt(sphere_radius * sphere_radius - ray_on_xy.squaredNorm()));

                // Compute normal at the intersection point
                Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_parallelogram()
{
    std::cout << "Simple ray tracer, one parallelogram with orthographic projection" << std::endl;

    const std::string filename("plane_orthographic.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is orthographic, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    


    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // Prepare the ray
            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = camera_view_direction;

            // TODO: Check if the ray intersects with the parallelogram
            double t = (pgram_origin - ray_origin).dot(pgram_v.cross(pgram_u)); 
            
            t = t / ray_direction.dot(pgram_v.cross(pgram_u));

            if (t>=0)
            {
                // TODO: The ray hit the parallelogram, compute the exact intersection
                // point
                //Vector3d ray_intersection = ray_origin + t * ray_direction;
                Vector3d ray_intersection = ray_origin + t * ray_direction;
                
                //COMMENTED OUT-> Vector3d ray_intersection(0, 0, 0);

                // TODO: Compute normal at the intersection point
                
                Vector3d ray_normal = (pgram_v.cross(pgram_u)).normalized();
                //COMMENTED OUT-> Vector3d ray_normal = ray_intersection.normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}

void raytrace_perspective()
{
    std::cout << "Simple ray tracer, one parallelogram with perspective projection" << std::endl;

    const std::string filename("plane_perspective.png");
    MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / C.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / C.rows(), 0);

    // TODO: Parameters of the parallelogram (position of the lower-left corner + two sides)
    const Vector3d pgram_origin(-0.5, -0.5, 0);
    const Vector3d pgram_u(0, 0.7, -10);
    const Vector3d pgram_v(1, 0.4, 0);

    // Single light source
    const Vector3d light_position(-1, 1, 1);

    for (unsigned i = 0; i < C.cols(); ++i)
    {
        for (unsigned j = 0; j < C.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;
 

            // TODO: Prepare the ray (origin point and direction)
            //const Vector3d ray_origin = pixel_center;
            const Vector3d ray_origin = camera_origin;

            //const Vector3d ray_direction = camera_view_direction;
            const Vector3d ray_direction = (pixel_center - camera_origin).normalized();

            // TODO: Check if the ray intersects with the parallelogram
            double t = (pgram_origin - ray_origin).dot(pgram_v.cross(pgram_u));
            
            t = t / ray_direction.dot(pgram_v.cross(pgram_u));
            if (t>=0)
            {
                // TODO: The ray hit the parallelogram, compute the exact intersection point
                //Vector3d ray_intersection(0, 0, 0);
                Vector3d ray_intersection = ray_origin + t * ray_direction;

                // TODO: Compute normal at the intersection point
                //Vector3d ray_normal = ray_intersection.normalized();
                Vector3d ray_normal = (pgram_v.cross(pgram_u)).normalized();

                // Simple diffuse model
                C(i, j) = (light_position - ray_intersection).normalized().transpose() * ray_normal;

                // Clamp to zero
                C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C, C, C, A, filename);
}


bool ray_sphere_intersection(const Vector3d &ray_origin, const Vector3d &ray_direction, const Vector3d &sphere_center, double sphere_radius, double &t)
{
    Vector3d distance = ray_origin - sphere_center;

    double x = ray_direction.dot(ray_direction);
    double y = 2.0 * distance.dot(ray_direction);
    double z = distance.dot(distance) - sphere_radius * sphere_radius;

    double discriminant = y * y - 4 * x * z;

    if (discriminant < 0) {
        return false;

    } else {
        double sqrt_discriminant = sqrt(discriminant);

        double root1 = (-y - sqrt_discriminant) / (2.0 * x);
        
        if (root1 < 0) {
            root1 = (-y + sqrt_discriminant) / (2.0 * x);
            if (root1 < 0) {
                return false;
            }
        }

        t = root1;
        return true;
    }
}

void raytrace_shading()
{
    std::cout << "Simple ray tracer, one sphere with different shading" << std::endl;

    const std::string filename("shading.png");
    //MatrixXd C = MatrixXd::Zero(800, 800); // Store the color
    MatrixXd C_R = MatrixXd::Zero(800, 800); //R, G, then B
    MatrixXd C_G = MatrixXd::Zero(800, 800); 
    MatrixXd C_B = MatrixXd::Zero(800, 800); 
   
    MatrixXd A = MatrixXd::Zero(800, 800); // Store the alpha mask

    const Vector3d camera_origin(0, 0, 3);
    const Vector3d camera_view_direction(0, 0, -1);

    // The camera is perspective, pointing in the direction -z and covering the unit square (-1,1) in x and y
    const Vector3d image_origin(-1, 1, 1);
    const Vector3d x_displacement(2.0 / A.cols(), 0, 0);
    const Vector3d y_displacement(0, -2.0 / A.rows(), 0);

    //Sphere setup
    const Vector3d sphere_center(0, 0, 0);
    const double sphere_radius = 0.9;

    //material params
    const Vector3d diffuse_color(1, 0, 1);
    const double specular_exponent = 100;
    const Vector3d specular_color(0., 0, 1);

    // Single light source
    const Vector3d light_position(-1, 1, 1);
    double ambient = 0.1;

    for (unsigned i = 0; i < C_R.cols(); ++i)
    {
        for (unsigned j = 0; j < C_R.rows(); ++j)
        {
            const Vector3d pixel_center = image_origin + double(i) * x_displacement + double(j) * y_displacement;

            // TODO: Prepare the ray (origin point and direction)
            //const Vector3d ray_origin = pixel_center;
            //const Vector3d ray_direction = camera_view_direction;

            const Vector3d ray_origin = pixel_center;
            const Vector3d ray_direction = (pixel_center - camera_origin).normalized();

            // Intersect with the sphere
            // TODO: implement the generic ray sphere intersection
            double t;
            const bool intersects = ray_sphere_intersection(ray_origin, ray_direction, sphere_center, sphere_radius, t);
            //if (true)
            if (intersects)
            {
                // TODO: The ray hit the sphere, compute the exact intersection point
                //Vector3d ray_intersection(0, 0, 0);
                Vector3d ray_intersection = ray_origin + t * ray_direction;

                // TODO: Compute normal at the intersection point
                //Vector3d ray_normal = ray_intersection.normalized();
                Vector3d ray_normal = (ray_intersection - sphere_center).normalized();

                // TODO: Add shading parameter here
                //const double diffuse = (light_position - ray_intersection).normalized().dot(ray_normal);
                //const double specular = (light_position - ray_intersection).normalized().dot(ray_normal);
                const double diffuse = std::max(0.0, ray_normal.dot((light_position - ray_intersection).normalized()));
                const double specular = pow(std::max(0.0, ray_normal.dot((ray_intersection - light_position).normalized())), specular_exponent);


                // Simple diffuse model
                C_R(i, j) = ambient * diffuse_color[0] + diffuse * diffuse_color[0] + specular * specular_color[0];
                C_G(i, j) = ambient * diffuse_color[1] + diffuse * diffuse_color[1] + specular * specular_color[1];
                C_B(i, j) = ambient * diffuse_color[2] + diffuse * diffuse_color[2] + specular * specular_color[2];

                // Clamp to [0, 1]
                C_R(i, j) = std::min(std::max(C_R(i, j), 0.0), 1.0);
                C_G(i, j) = std::min(std::max(C_G(i, j), 0.0), 1.0);
                C_B(i, j) = std::min(std::max(C_B(i, j), 0.0), 1.0);

                //C(i, j) = ambient + diffuse + specular;

                // Clamp to zero
                //C(i, j) = std::max(C(i, j), 0.);

                // Disable the alpha mask for this pixel
                A(i, j) = 1;
            }
        }
    }

    // Save to png
    write_matrix_to_png(C_R, C_G, C_B, A, filename);
}



int main()
{
    raytrace_sphere();
    raytrace_parallelogram();
    raytrace_perspective();
    raytrace_shading();

    return 0;
}
