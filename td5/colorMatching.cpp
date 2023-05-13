#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../stb/stb_image.h"

#include "../gx_random.h"
#include "../gx_vector.h"

#include <vector>
#include <iostream>

#include <unistd.h>

#ifndef COLOR_TRANSFER_N_ITER
    #define COLOR_TRANSFER_N_ITER 10
#endif

typedef std::vector<Vector> VectorImage;
typedef unsigned char * ImageArray;

#define SIZE_OF_PIXEL_COMPONENT sizeof(unsigned char)

// Utils
static Vector randomVector(){
    Vector out;
    randomSpherePoint(out[0], out[1], out[2]);
    return out;
}

VectorImage imageA2V(int W, int H, const ImageArray image){
    VectorImage vectorImage(W * H);
    #pragma omp parallel for
    for (int i=0; i<H; i++){
        for(int j=0; j<W; j++){
            int index = i*W + j;
            vectorImage[index] = Vector(image[3*index], image[3*index+1], image[3*index+2]);
        }
    }

    return vectorImage;
}

ImageArray imageV2A(int W, int H, const VectorImage& VectorImage){
    ImageArray image = (ImageArray) malloc(3*W*H*SIZE_OF_PIXEL_COMPONENT);
    if (!image) {
        fprintf(stderr, "Malloc failed in imageV2A.\n");
        exit(1);
    }

    #pragma omp parallel for
    for (int i=0; i<H; i++){
        for(int j=0; j<W; j++){
            int index = i*W + j;
            for (int channel=0; channel<3; channel++){
                // Add +0.5 because compiler truncates double
                image[3*index + channel] = VectorImage[index][channel] + 0.5;
            }
        }
    }

    return image;
}

void writeVectorImage(const char *filename, int W, int H, const VectorImage& vecImage){
    ImageArray image = imageV2A(W, H, vecImage);
    stbi_write_png(filename, W, H, 3, &image[0], 0);
    free(image);
}

VectorImage slicedOptimalTransportColorTransfer(const VectorImage& sourceImage, const VectorImage& targetImage){
    if (sourceImage.size() != targetImage.size()) {
        fprintf(stderr, "Error in slicedOptimalTransportColorTransfer: ");
        fprintf(stderr, "The two images must have the same number of pixels.\n");
        exit(1);
    }

    int numPixels = sourceImage.size();
    std::vector<std::pair<double, int>> sourceProj(numPixels);
    std::vector<double> targetProj(numPixels);

    VectorImage newImage(sourceImage);

    for(int iteration=0; iteration < COLOR_TRANSFER_N_ITER; iteration++){
        Vector v = randomVector();

        #pragma omp parallel for
        for (int pixelId=0; pixelId < numPixels; pixelId++){
            sourceProj[pixelId] = std::make_pair(dot(newImage[pixelId], v), pixelId);
            targetProj[pixelId] = dot(targetImage[pixelId], v);
        }

        std::sort(sourceProj.begin(), sourceProj.end());
        std::sort(targetProj.begin(), targetProj.end());

        #pragma omp parallel for
        for (int pixelId=0; pixelId < numPixels; pixelId++){
            newImage[sourceProj[pixelId].second] = (
                newImage[sourceProj[pixelId].second]
                + (targetProj[pixelId] - sourceProj[pixelId].first) * v
                ).clip(0., 255.);
        }
    }

    return newImage;
}

int main(int argc, char **argv){

    if (argc != 4) {
        fprintf(stderr, "Usage: %s source_image target_image output_path.\n", argv[0]);
        exit(1);
    }

    const char *source_image_path = argv[1];
    const char *target_image_path = argv[2];
    const char *write_image_path = argv[3];

    if (access(source_image_path, F_OK)) {
        fprintf(stderr, "Error: no file %s \n", source_image_path);
        exit(1);
    }

    if (access(target_image_path, F_OK)) {
        fprintf(stderr, "Error: no file %s \n", target_image_path);
        exit(1);
    }

    int W, H;
    ImageArray a_SourceImage, a_TargetImage;

    {
        int sW, sH, sC;
        int tW, tH, tC;
        a_SourceImage = stbi_load(source_image_path, &sW, &sH, &sC, STBI_rgb);
        a_TargetImage = stbi_load(target_image_path, &tW, &tH, &tC, STBI_rgb);

        if(sW != tW || sH != tH) {
            fprintf(stderr, "Error: Dimensions of input images do not match.\n");
            exit(1);
        }

        W = sW;
        H = sH;
    }

    // VectorImage whiteImage(W * H);
    // for(auto &p : whiteImage){
    //     p = Vector(255, 255, 255);
    // }

    // writeVectorImage(write_image_path, W, H, whiteImage);
    // exit(1);

    VectorImage v_SourceImage = imageA2V(W, H, a_SourceImage);
    VectorImage v_TargetImage = imageA2V(W, H, a_TargetImage);

    VectorImage v_OutputImage = slicedOptimalTransportColorTransfer(v_SourceImage, v_TargetImage);

    writeVectorImage(write_image_path, W, H, v_OutputImage);

    exit(0);
}