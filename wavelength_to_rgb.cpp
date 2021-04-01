/*
 * The CIE, or Commission Internationale de l’Eclairage, 
 * has quantative data for the subjective realm of color
 * perception. This data can be used to map a wavelength 
 * of light to a sRGB color space, among other uses. In
 * this program there are two functions for approximating 
 * RGB colors from a wavelength.
 *
 *  Resources: http://www.fourmilab.ch/documents/specrend/
 *             https://stackoverflow.com/questions/1472514/convert-light-frequency-to-rgb
 *             http://cvrl.ioo.ucl.ac.uk/cmfs.htm
 *             https://www.oceanopticsbook.info/view/photometry-and-visibility/from-xyz-to-rgb
 *             http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html#WSMatrices
 */

#include <iostream>
#include <cmath>
#include <fstream>           //for writing output to a .ppm file
#include "cie_color_match.h" //for the CIE color matching functions
using namespace std;

double r, g, b;

/* utility functions */

// clamp for color values, default is the range [0, 1]
inline double clamp(double x, double min = 0.0, double max = 0.99)
{
    return (x < min) ? min : (x > max) ? max : x;
}

// ppm files only support RGB values in range [0, 255], so we clamp to those
inline void rgb_clamp(double &r, double &g, double &b)
{
    r = static_cast<int>(255.999 * clamp(r));
    g = static_cast<int>(255.999 * clamp(g));
    b = static_cast<int>(255.999 * clamp(b));
}

// The CIE color matching functions used are from the
// 10-deg XYZ CMFs transformed from the CIE (2006) 2-deg LMS cone fundamentals
// located at http://cvrl.ioo.ucl.ac.uk/cmfs.htm
// and stored in cie_color_match.h as cie_color_match[351][3]
void wavelength_to_rgb(double lambda, double &r, double &g, double &b)
{
    // convert wavelength lambda to X, Y, Z coordinates using CIE color matching functions
    double X, Y, Z;

    // shift range to match array, then get the respective val for X, Y, Z
    int idx = (lambda - 400);

    X = cie_colour_match[idx][0];
    Y = cie_colour_match[idx][1];
    Z = cie_colour_match[idx][2];

    // convert CIE XYZ coords. into x, y, z chromaticity coords.
    double magnitude = X+Y+Z;
    double x, y, z;
    x = X / magnitude;
    y = Y / magnitude;
    z = 1 - (x+y);
    x = X, y = Y, z = Z;

    /* convert x, y, z coordinates to linear sRGB
     * using transformation matrix
     * from: http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html#WSMatrices
     * 
     *  [ 3.2404542 -1.5371385 -0.4985314
     *   -0.9692660  1.8760108  0.0415560
     *   0.0556434 -0.2040259  1.0572252 ]
     * 
     * sRGB uses the IlluminantD65, where (x, y) = (0.3127, 0.3291) for the color white
     */
    r =  3.2404542*x - 1.5371385*y - 0.4985314*z;
    g = -0.9692660*x + 1.8760108*y + 0.0415560*z;
    b =  0.0556434*x - 0.2040259*y + 1.0572252*z;
    /* the difference in colors isn't linear, however
     * (e.g. the ends of the spectrum look darker)
     * so we apply gamma correction of 2.4
     * https://www.oceanopticsbook.info/view/photometry-and-visibility/from-xyz-to-rgb
     */
    r = (r <= 0.0031308) ? 12.92*r : 1.055*pow(r, (1/2.4)) - 0.055;
    g = (g <= 0.0031308) ? 12.92*g : 1.055*pow(g, (1/2.4)) - 0.055;
    b = (b <= 0.0031308) ? 12.92*b : 1.055*pow(b, (1/2.4)) - 0.055;
}

// method two is an approximation -
// a c++ port of Tarc's code from 
// https://stackoverflow.com/questions/1472514/convert-light-frequency-to-rgb
void approx_wavelength_to_rgb(double lambda, double &r, double &g, double &b)
{
    double falloff;
    static double gamma = 0.8;

    // convert wavelength lambda to x, y, z coordinates using approximations
    if ((lambda >= 380) && (lambda < 440))
    {
        r = -(lambda - 440) / (440 - 380);
        g = 0.0;
        b = 1.0;
    } else if ((lambda >= 440) && (lambda < 490))
    {
        r = 0.0;
        g = (lambda - 440) / (490 - 440);
        b = 1.0;
    
    } else if ((lambda >= 490) && (lambda < 510))
    {
        r = 0.0;
        g = 1.0;
        b = -(lambda - 510) / (510 - 490);
    
    } else if ((lambda >= 510) && (lambda < 580))
    {
        r = (lambda - 510) / (580-510);
        g = 1.0;
        b = 0.0;
    
    } else if ((lambda >= 580) && (lambda < 645))
    {
        r = 1.0;
        g = -(lambda - 645) / (645-580);
        b = 0.0;
    } else if ((lambda >= 645) && (lambda < 781))
    {
        r = 1.0;
        g = 0.0;
        b = 0.0;
    } else {
        r = 0.0;
        g = 0.0;
        b = 0.0;
    }
    // intensity falloff at the ends of the visible spectrum
    if ((lambda >= 380) && (lambda < 420))
    {
        falloff = 0.3 + 0.7*(lambda - 380) / (420 - 380);
    } else if ((lambda >= 420) && (lambda < 701))
    {
        falloff = 1.0;
    } else if ((lambda >= 701) && (lambda < 781))
    {
        falloff = 0.3 + 0.7*(780-lambda) / (780-700);
    } else {
        falloff = 0.0;
    }
    // gamma correction
    r = (r == 0.0) ? 0 : pow(r*falloff, gamma);
    g = (g == 0.0) ? 0 : pow(g*falloff, gamma);
    b = (b == 0.0) ? 0 : pow(b*falloff, gamma);
}

int main()
{
    int px_per_wavelength = 2; // change this for more px's 
    int vis_spectrum_range = 700 - 400;

    // file setup
    const char image_name[] = "images/spectrum.ppm";
    const int  image_width  = vis_spectrum_range * px_per_wavelength; 
    const int  image_height = 100;

    ofstream file;
    file.open(image_name);
    file << "P3\n" << image_width << " " << image_height << "\n255\n";
    
    
    for (int j=image_height; j>=0; --j)
    {
        // 380nm to 780 nm is the visible spectrum
        // 400 to 700 nm is a *good enough* approximation

        // sample wavelegth at every 1 nm
        for (int i=400; i<700; ++i)
        {
            wavelength_to_rgb(i, r, g, b);
            rgb_clamp(r, g, b);
            // generate multiple pixels per wavelength sample
            for (int s=0; s<px_per_wavelength; ++s)
            {
                file << r << " " << g << " " << b << "\n";
            }
        }
    }
    cerr << "\n" << image_name << " created successfully.\n";
    file.close(); // make sure to close the file!
    return 0;
}