/**
 * Copyright (c) 2016, David Stutz
 * Contact: david.stutz@rwth-aachen.de, davidstutz.de
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <fstream>
#include <opencv2/opencv.hpp>
#include <boost/filesystem.hpp>
#include <boost/timer.hpp>
#include <boost/program_options.hpp>
#include <boost/timer.hpp>
#include "ergc_opencv.h"
#include "io_util.h"
#include "superpixel_tools.h"
#include "visualization.h"

/** \brief Command line tool for running ERGC.
 * Usage:
 * \code{sh}
 *   $ ../bin/ergc_cli --help
 *   Allowed options:
 *     -h [ --help ]                   produce help message
 *     -i [ --input ] arg              the folder to process
 *     -s [ --superpixels ] arg (=400) number of superpixels
 *     -r [ --color-space ] arg (=1)   color space; 0 = RGB, >0 = Lab
 *     -p [ --perturb-seeds ] arg (=1) >0 for perturbing seeds
 *     -c [ --compacity ] arg (=0)     compacity
 *     -f [ --fair ]                   for a fair comparison with other algorithms, 
 *                                     quadratic blocks are used for initialization
 *     -o [ --csv ] arg                save segmentation as CSV file
 *     -v [ --vis ] arg                visualize contours
 *     -x [ --prefix ] arg             output file prefix
 *     -w [ --wordy ]                  verbose/wordy/debug
 * \endcode
 * \author David Stutz
 */

void WriteImageP5(int *labels, int num_cols, int num_rows, const char *filename)
{
    FILE *fp = NULL;
    int p, hi, lo;
    uchar *data8 = NULL;
    ushort *data16 = NULL;

    fp = fopen(filename, "wb");

    fprintf(fp, "P5\n");
    fprintf(fp, "%d %d\n", num_cols, num_rows);

    int maximum = 0, minimum = -1;
    for (int i = 0; i < num_cols * num_rows; i++)
    {
        if (labels[i] > maximum)
            maximum = labels[i];
        if (labels[i] < minimum || minimum == -1)
            minimum = labels[i];
    }

    if ((maximum < 256) && (minimum >= 0))
    {
        fprintf(fp, "%d\n", 255);
        data8 = (uchar *)calloc(num_cols * num_rows, sizeof(uchar));
        for (p = 0; p < num_cols * num_rows; p++)
            data8[p] = (uchar)labels[p];

        fwrite(data8, sizeof(uchar), num_cols * num_rows, fp);
        free(data8);
    }
    else if (maximum < 65536)
    {
        fprintf(fp, "%d\n", 65535);
        data16 = (ushort *)calloc(num_cols * num_rows, sizeof(ushort));
        for (p = 0; p < num_cols * num_rows; p++)
            data16[p] = (ushort)labels[p];

        {
#define HI(num) (((num)&0x0000FF00) >> 8)
#define LO(num) ((num)&0x000000FF)
            for (p = 0; p < num_cols * num_rows; p++)
            {
                hi = HI(data16[p]);
                lo = LO(data16[p]);
                fputc(hi, fp);
                fputc(lo, fp);
            }
        }

        free(data16);
    }
    else
    {
        printf("Cannot write image as P5 (%d/%d)", maximum, minimum);
    }
    fclose(fp);
}

void WriteImageP2(int *labels, int num_cols, int num_rows, const char *filename)
{
    FILE *fp = NULL;
    int p;

    int maximum = 0;
    for (int i = 0; i < num_cols * num_rows; i++)
        if (labels[i] > maximum)
            maximum = labels[i];

    fp = fopen(filename, "w");

    fprintf(fp, "P2\n");
    fprintf(fp, "%d %d\n", num_cols, num_rows);
    fprintf(fp, "%d\n", maximum);

    int k = 0;
    for (int i = 0; i < num_cols; i++)
    {
        for (int j = 0; j < num_rows; j++)
        {
            fprintf(fp, "%d ", labels[k]);
            k++;
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

void WriteLabels(int *labels, int num_cols, int num_rows, const char *filename)
{
    int maximum = 0;
    for (int i = 0; i < num_cols * num_rows; i++)
        if (labels[i] > maximum)
            maximum = labels[i];

    if (maximum > 255)
        WriteImageP2(labels, num_cols, num_rows, filename);
    else
        WriteImageP5(labels, num_cols, num_rows, filename);
}


int main(int argc, const char** argv) {
    
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help,h", "produce help message")
        ("input,i", boost::program_options::value<std::string>(), "the folder to process")
        ("superpixels,s", boost::program_options::value<int>()->default_value(400), "number of superpixels")
        ("color-space,r", boost::program_options::value<int>()->default_value(1), "color space; 0 = RGB, >0 = Lab")
        ("perturb-seeds,p", boost::program_options::value<int>()->default_value(1), ">0 for perturbing seeds")
        ("compacity,c", boost::program_options::value<int>()->default_value(0), "compacity")
        ("fair,f", "for a fair comparison with other algorithms, quadratic blocks are used for initialization")
        ("csv,o", boost::program_options::value<std::string>()->default_value(""), "save segmentation as CSV file")
        ("vis,v", boost::program_options::value<std::string>()->default_value(""), "visualize contours")
        ("prefix,x", boost::program_options::value<std::string>()->default_value(""), "output file prefix")
        ("time,time", boost::program_options::value<std::string>()->default_value(""), "Time file")
        ("dtime,dtime", boost::program_options::value<std::string>()->default_value(""), "Detailed time file")
        ("wordy,w", "verbose/wordy/debug");

    boost::program_options::positional_options_description positionals;
    positionals.add("input", 1);
    
    boost::program_options::variables_map parameters;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(desc).positional(positionals).run(), parameters);
    boost::program_options::notify(parameters);

    if (parameters.find("help") != parameters.end()) {
        std::cout << desc << std::endl;
        return 1;
    }
    
    boost::filesystem::path output_dir(parameters["csv"].as<std::string>());
    if (!output_dir.empty()) {
        if (!boost::filesystem::is_directory(output_dir)) {
            boost::filesystem::create_directories(output_dir);
        }
    }
    
    boost::filesystem::path vis_dir(parameters["vis"].as<std::string>());
    if (!vis_dir.empty()) {
        if (!boost::filesystem::is_directory(vis_dir)) {
            boost::filesystem::create_directories(vis_dir);
        }
    }
    
    boost::filesystem::path input_dir(parameters["input"].as<std::string>());
    if (!boost::filesystem::is_directory(input_dir)) {
        std::cout << "Image directory not found ..." << std::endl;
        return 1;
    }
    
    std::string prefix = parameters["prefix"].as<std::string>();
    
    bool wordy = false;
    if (parameters.find("wordy") != parameters.end()) {
        wordy = true;
    }
    
    int superpixels = parameters["superpixels"].as<int>();
    int color_space = parameters["color-space"].as<int>();
    char dtime_file[255];

    sprintf(dtime_file, "%s", parameters["dtime"].as<std::string>().c_str());

    bool lab = false;
    if (color_space > 0) {
        lab = true;
    }
    
    int perturb_seeds_int = parameters["perturb-seeds"].as<int>();
    bool perturb_seeds = false;
    if (perturb_seeds_int > 0) {
        perturb_seeds = true;
    }
    
    int compacity = parameters["compacity"].as<int>();
    
    std::multimap<std::string, boost::filesystem::path> images;
    std::vector<std::string> extensions;
    IOUtil::getImageExtensions(extensions);
    IOUtil::readDirectory(input_dir, extensions, images);
    
    clock_t start, end;
    double sum_time = 0;
    for (std::multimap<std::string, boost::filesystem::path>::iterator it = images.begin(); 
            it != images.end(); ++it) {
        
        cv::Mat image = cv::imread(it->first);
        
        int region_width;
        int region_height;
        SuperpixelTools::computeHeightWidthFromSuperpixels(image, superpixels, 
                region_height, region_width);
        
        // If a fair comparison is requested:
        if (parameters.find("fair") != parameters.end()) {
            region_width = SuperpixelTools::computeRegionSizeFromSuperpixels(image, 
                    superpixels);
            region_height = region_width;
        }
        
        double time = 0;
        cv::Mat labels;
        start = clock();
        ERGC_OpenCV::computeSuperpixels(image, region_height, region_width, 
                lab, perturb_seeds, compacity, labels);
        end = clock();
        time = ((double)(end - start)) / (double)(CLOCKS_PER_SEC);
        sum_time += time;

        SuperpixelTools::relabelSuperpixels(labels);

        if (strcmp(dtime_file, "") != 0)
        {
            FILE *fp = fopen(dtime_file, "a+");
            fprintf(fp, "%s %.5f\n", it->second.stem().string().c_str(), time);
            fclose(fp);
        }
        
        if (!output_dir.empty()) {
            boost::filesystem::path csv_file(output_dir 
                    / boost::filesystem::path(prefix + it->second.stem().string() + ".pgm"));
            WriteLabels((int *)labels.data, image.cols, image.rows, csv_file.string().c_str());
        }
        
        if (!vis_dir.empty()) {
            boost::filesystem::path contours_file(vis_dir 
                    / boost::filesystem::path(prefix + it->second.stem().string() + ".png"));
            cv::Mat image_contours;
            Visualization::drawContours(image, labels, image_contours);
            cv::imwrite(contours_file.string(), image_contours);
        }
    }
    
    sprintf(dtime_file, "%s", parameters["time"].as<std::string>().c_str());
    if (strcmp(dtime_file, "") != 0)
    {
        // get image name
        FILE *fp = fopen(dtime_file, "a+");
        fprintf(fp, "%d %.5f\n", superpixels, sum_time / images.size());
        fclose(fp);
    }
    
    return 0;
}
