#include <iostream>
#include <vector>
//
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <math.h>
#include <string>

class img{
    public:
        int h,w,c;
        std::vector<std::vector<std::vector<float>>> fi;

        img(){}

        img(int _h, int _w, int _c){
            w=_w;
            h=_h;
            c=_c;
            fi=std::vector<std::vector<std::vector<float>>>(
                h, std::vector<std::vector<float>>(
                    w, std::vector<float>(
                        c, 0.0)));
        }

};

void smooth_gaussian(img& original, img& smoothed, int r, float sigma) {
    int r2 = 2 * r + 1;
    // precompute spatial stencil
    float stencil[r2 * r2];
    for (int dy = -r; dy <= r; ++dy){
        for (int dx = -r; dx <= r; ++dx){
            float h = dx * dx + dy * dy;
            int idx = dx + r + r2 * (dy + r);
            stencil[idx] = std::exp(-h / (2 * sigma * sigma));

        }    
    }
    
    // apply filter
    for (int px = 0; px < original.h; px++){
        for (int py = 0; py < original.w;  py++){
            for(int c=0;c<original.c;c++){
                float w_sum = 0.0;
                for (int dx = -r; dx <= r; ++dx){
                    for (int dy = -r; dy <= r; ++dy)
                    {
                        int px1 = px + dx;
                        int py1 = py + dy;
                        if (0 <= px1 && 0 <= py1 && px1 < original.h && py1 < original.w) {
                            float w = stencil[dy + r + r2 * (dx + r)];

                            
                            smoothed.fi[px][py][c] += w*original.fi[px1][py1][c];
                            w_sum += w;

                        }
                    }
                }
                smoothed.fi[px][py][c] /= w_sum;
                smoothed.fi[px][py][c]=std::min(1.f,std::max(0.f,smoothed.fi[px][py][c]));}
        }    
    }
    
};

void smooth_bilateral(img &in1, img &out, int r, float sigma_space, float sigma_range) {
    int r2 = 2 * r + 1;
    // precompute spatial stencil
    float stencil[r2 * r2];
    for (int dy = -r; dy <= r; ++dy){
        for (int dx = -r; dx <= r; ++dx){
            float h = dx * dx + dy * dy;
            int idx = dx + r + r2 * (dy + r);
            stencil[idx] = std::exp(-h/ (2 * sigma_space * sigma_space));
        }
    }
    // apply filter
    for (int py = 0; py < in1.h; py++){
        for (int px = 0; px < in1.w;  px++)
        {
                // int idx0 = px + in1.w * py;
                float w_sum = 0.0;
                for (int dy = -r; dy <= r; ++dy){
                    for (int dx = -r; dx <= r; ++dx)
                    {
                        int px1 = px + dx;
                        int py1 = py + dy;
                        if (0 <= px1 && 0 <= py1 && px1 < in1.w && py1 < in1.h) {
                            float w_space = stencil[dx + r + r2 * (dy + r)];
                            float h1=0.0;
                            for(int c=0;c<in1.c;c++){
                                float dr=in1.fi[py1][px1][c]-in1.fi[py][px][c];
                                h1+=dr*dr;}

                            float w_range=std::exp(-h1/ (2 * sigma_range * sigma_range));                            
                            float w = w_space * w_range;
                            
                            for(int c=0;c<in1.c;c++){
                                out.fi[py][px][c] += w*in1.fi[py1][px1][c];}
                            w_sum += w;
                        }
                    }
                }
            
            for(int c=0;c<in1.c;c++){out.fi[py][px][c] /= w_sum;
            out.fi[py][px][c]=std::min(1.f,std::max(0.f,out.fi[py][px][c]));}
        } 
    }
};

void smooth_bilateral_joint(img &in1,img &in2, img &out, int r, float sigma_space, float sigma_range) {
    int r2 = 2 * r + 1;
    // precompute spatial stencil
    float stencil[r2 * r2];
    for (int dy = -r; dy <= r; ++dy){
        for (int dx = -r; dx <= r; ++dx){
            float h = dx * dx + dy * dy;
            int idx = dx + r + r2 * (dy + r);
            stencil[idx] = std::exp(-h/ (2 * sigma_space * sigma_space));
        }
    }
    // apply filter
    for (int py = 0; py < in1.h; py++){
        for (int px = 0; px < in1.w;  px++)
        {
                // int idx0 = px + in1.w * py;
                float w_sum = 0.0;
                for (int dy = -r; dy <= r; ++dy){
                    for (int dx = -r; dx <= r; ++dx)
                    {
                        int px1 = px + dx;
                        int py1 = py + dy;
                        if (0 <= px1 && 0 <= py1 && px1 < in2.w && py1 < in2.h) {
                            float w_space = stencil[dx + r + r2 * (dy + r)];
                            float h1=0.0;
                            for(int c=0;c<in2.c;c++){
                                float dr=in2.fi[py1][px1][c]-in2.fi[py][px][c];
                                h1+=dr*dr;}

                            float w_range=std::exp(-h1/ (2 * sigma_range * sigma_range));                            
                            float w = w_space * w_range;
                            
                            for(int c=0;c<in1.c;c++){
                                out.fi[py][px][c] += w*in1.fi[py1][px1][c];}
                            w_sum += w;
                        }
                    }
                }
            
            for(int c=0;c<in2.c;c++){out.fi[py][px][c] /= w_sum;
            out.fi[py][px][c]=std::min(1.f,std::max(0.f,out.fi[py][px][c]));}
        } 
    }
};

void bilinear_inter(img &inp, img &out,float a){
    for(int i=0; i<out.h;i++){
        for(int j=0; j<out.w;j++){
            for(int c=0;c<out.c;c++){
                float p_x=i/a;
                float p_y=j/a;
                int x1=floor(p_x);
                int x2=x1+1;
                int y1=floor(p_y);
                int y2=y1+1;
                float w11=(x2-p_x)*(y2-p_y);
                float w12=(x2-p_x)*(p_y-y1);
                float w21=(p_x-x1)*(y2-p_y);
                float w22=(p_x-x1)*(p_y-y1);
            
                if(x2<inp.h && y2<inp.w){
                    out.fi[i][j][c]+=w11*inp.fi[x1][y1][c] + w12*inp.fi[x1][y2][c] + w21*inp.fi[x2][y1][c] + w22*inp.fi[x2][y2][c];
                }
                else if(y2<inp.w){
                    out.fi[i][j][c]+=(y2-p_y)*inp.fi[x1][y1][c] + (p_y-y1)*inp.fi[x1][y2][c];
                }
                else if(x2<inp.h){
                    out.fi[i][j][c]+=(x2-p_x)*inp.fi[x1][y1][c] + (p_x-x1)*inp.fi[x2][y1][c];
                }
                else{
                    out.fi[i][j][c]=inp.fi[x1][y1][c];}
                out.fi[i][j][c]=std::min(1.f,std::max(0.f,out.fi[i][j][c]));
            }
        }
    }
};

void upsample(img &inp, img &out, float a) {
    float s = 1.0;
    img imge(inp.h*a,inp.w*a,inp.c);    
    
    bilinear_inter(inp, imge, a);
    int stencil = 2;
    smooth_gaussian(imge, out, stencil, s);
}

void downsample(img &inp, img &out, float a) {
    int stencil = 2;
    float s = 1.0;
    img imge(inp.h,inp.w,inp.c);
    
    smooth_gaussian(inp, imge, stencil, s);
    bilinear_inter(imge, out, a);
}

const float sigma_1 = 12.0;
const float sigma_2 = 0.05;

int main(){
    int width,height,ch;
    //read png
    unsigned char *image_ = stbi_load("../image.png", &height, &width, &ch, 0);    
    img datum=img(height,width,ch);
    for (int h = 0; h < height; h++) {
        for (int w = 0; w < width; w++) {
            int index = (h*width + w) * ch;
            for (int c = 0; c < ch; c++) 
                datum.fi[h][w][c] = static_cast<float>(image_[index + c]) / 255.0;
        }
    }
    stbi_image_free(image_);

    int n=floor(std::log2(std::max(height,width))) - 5;
    std::vector<img> P1(n + 1);
    std::vector<img> P2(n + 1);
    std::vector<img> P3(n + 1);
    std::vector<img> PSU(n + 1);
    std::vector<img> Gau(n + 1);
    std::vector<img> Lpl(n + 1);

    Gau[0] = datum;
    for (int l = 0; l < n; l++) {
        Gau[l + 1] = img(Gau[l].h*0.5, Gau[l].w*0.5, Gau[l].c);
        downsample(Gau[l], Gau[l + 1], 0.5);
        // std::cout << "Gau " << l << std::endl;        
    }

    Lpl[n] = Gau[n];
    for (int l = 0; l < n; l++) {
        img LG=img(Gau[l + 1].h*2.0,Gau[l + 1].w*2.0,Gau[l + 1].c);
        upsample(Gau[l + 1], LG, 2.0);
        // std::cout << LG.h<<" "<<Gau[l].h<<"  "<<Gau[l].w<<" "<<LG.w << std::endl;
        Lpl[l]=img(Gau[l].h,Gau[l].w,Gau[l].c);
        for (int i = 0; i < LG.h; ++i) {
            for (int j = 0; j < LG.w; ++j) {
                for (int c = 0; c < LG.c; ++c){
                    Lpl[l].fi[i][j][c]=std::max(0.f,Gau[l].fi[i][j][c]-LG.fi[i][j][c]);
                }
            }
        }
        std::cout << "Lap " << l << std::endl;
    }

    PSU[n] = Gau[n];
    float sigma_l = sigma_1/std::pow(2.0, n+1);
    for (int l = n; l > 0; l--) {
        sigma_l *= 2;
        int k1 = std::round(std::max(sigma_l, 1.f));
        int k2 = std::round(std::max(4.f*sigma_l, 1.f));
        
        PSU[l - 1]=img(PSU[l].h*2.0,PSU[l].w*2.0,PSU[l].c);
        P1[l - 1]=img(PSU[l].h*2.0,PSU[l].w*2.0,PSU[l].c);
        P2[l - 1]=img(PSU[l].h*2.0,PSU[l].w*2.0,PSU[l].c);
        P3[l - 1]=img(Lpl[l - 1].h,Lpl[l - 1].w,Lpl[l - 1].c);
        
        bilinear_inter(PSU[l], P1[l - 1], 2.0);
        smooth_bilateral_joint(P1[l - 1], Gau[l - 1], P2[l - 1], k1, sigma_l, sigma_2);

        for (int i = 0; i < P2[l - 1].h; i++) {
            for (int j = 0; j < P2[l - 1].w; j++) {
                for (int c = 0; c < P2[l - 1].c; c++){
                    P3[l - 1].fi[i][j][c]=std::min(1.f,P2[l - 1].fi[i][j][c]+Lpl[l - 1].fi[i][j][c]);
                }
            }
        }
        
        img P= img(PSU[l].h*2.0,PSU[l].w*2.0,PSU[l].c);
        // std::cout << P3[l-1].h<<" "<<P2[l-1].h<<"  "<<PSU[l-1].h<<" "<<P3[l-1].w << std::endl;

        smooth_bilateral_joint(P3[l - 1], P2[l - 1], PSU[l - 1], k2, sigma_l, sigma_2);       
        // std::cout << l << std::endl;
        smooth_bilateral(PSU[l - 1], P, k2, sigma_l, sigma_2);
        PSU[l - 1] = P;
        if(l==1){
            int hh=PSU[l-1].h;
            int ww=PSU[l-1].w;
            int cc=PSU[l-1].c;
            std::vector<unsigned char> image_(hh*ww*cc);
            for (int i = 0; i < hh; i++) {
                for (int j = 0; j < ww; j++) {
                    int index = (i*ww + j) * cc;
                    for (int z = 0; z < cc; z++) {
                        // if(j==37){std::cout << R[l-1].fi[i][j][z] <<std::endl;}
                        unsigned char pc = static_cast<unsigned char>(PSU[l-1].fi[i][j][z] * 255.0);
                        image_[index + z] = std::max((unsigned char)(0), std::min((unsigned char)255, pc));
                    }
                }
            }
            int op = stbi_write_png("output.png", height, width, ch, image_.data(), height*ch);
            if (!op){std::cerr << "Write Failure" << std::endl;}
        }
    }
}