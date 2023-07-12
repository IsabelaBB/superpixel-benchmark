#ifndef UTILS_H
#define UTILS_H

#include <opencv2/opencv.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <string>

template <typename T>
void imagesc(std::string name_s, T* tab, int width, int height);

#endif
