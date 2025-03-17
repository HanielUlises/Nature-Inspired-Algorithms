#ifndef UTILS_H
#define UTILS_H

#include <opencv2/opencv.hpp>
#include <cmath>

inline cv::Mat sigmoid_transform(const cv::Mat &gray_img, float alpha, float delta) {
    cv::Mat new_image = gray_img.clone();
    for (int i = 0; i < gray_img.rows; i++) {
        for (int j = 0; j < gray_img.cols; j++) {
            new_image.at<float>(i, j) = 1.0f / (1.0f + std::exp(-alpha * (gray_img.at<float>(i, j) - delta)));
        }
    }
    return new_image;
}

inline double entropy(const cv::Mat& img) {
    cv::Mat hist;
    int histSize = 256;
    float range[] = {0, 256};
    const float* histRange = {range};
    cv::calcHist(&img, 1, 0, cv::Mat(), hist, 1, &histSize, &histRange);
    hist /= img.total();
    double entropy = 0;
    for (int i = 0; i < histSize; i++) {
        float p = hist.at<float>(i);
        if (p > 0) {
            entropy -= p * log2(p);
        }
    }
    return entropy;
}

inline double standard_deviation(const cv::Mat& img) {
    cv::Scalar mean, stddev;
    cv::meanStdDev(img, mean, stddev);
    return stddev[0];
}

#endif // UTILS_H