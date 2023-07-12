#include "utils.h"

void imagesc(std::string name_s, cv::Mat* map)
{
        char* name = const_cast<char*>(name_s.c_str());
        double min;
        double max;
        cv::minMaxIdx(*map, &min, &max);
        cv::Mat adjMap;
        // Histogram Equalization
        float scale = 255 / (max-min);
        map->convertTo(adjMap,CV_8UC1, scale, -min*scale);

        cv::Mat resultMap;
        applyColorMap(adjMap, resultMap, cv::COLORMAP_PARULA);

        cv::namedWindow(name, cv::WINDOW_NORMAL);
        cv::imshow(name, resultMap);
        cv::waitKey(1);
}

template <typename T>
void imagesc(std::string name_s, T* tab, int width, int height) {
        cv::Mat map = cv::Mat(height, width, CV_64F);

        for(int i=0; i<height; i++) {
                for(int j=0; j<width; j++)
                        map.at<double>(i, j) = (double)tab[i*width + j];
        }

        imagesc(name_s, &map);
}

template
void imagesc(std::string name_s, double* tab, int width, int height);

template
void imagesc(std::string name_s, int* tab, int width, int height);

template
void imagesc(std::string name_s, unsigned char* tab, int width, int height);

