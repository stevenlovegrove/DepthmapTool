#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <cvd/image_io.h>
#include <cvd/image.h>

using namespace std;
using namespace TooN;
using namespace CVD;

void readdoubleimage(istream& s, Image<double>& img)
{
    for(int i = 0 ; i < img.size().y ; i++)
    {
        for (int j = 0 ; j < img.size().x ; j++)
        {
            double val = 0;
            s >> val;
            img[i][j] = val;
        }
    }
}

void normalise(Image<double>& in, Image<byte>& out)
{
    assert(in.size() == out.size() );

    double mind = *min_element(in.data(),in.data()+in.size().area());
    double maxd = *max_element(in.data(),in.data()+in.size().area());
    const double range = maxd - mind;

    for(int i = 0 ; i < in.size().y ; i++ )
    {
        for(int j = 0 ; j < in.size().x ; j++ )
            out[i][j] = 255.0 * (in[i][j] - mind) / range;
    }
}

void PovrayMesh(Image<double>& depth, ostream& s )
{
    const int w = depth.size().x;
    const int h = depth.size().y;
    const float F = w*880/640;

    const int u0 = w/2;
    const int v0 = h/2;

    const float ms = 3;

    Vector<3>* vs = new Vector<3>[w*h];

    s << "#declare dm= mesh2{" << endl;

    { // Vertices
        s << "vertex_vectors{ " << w*h;
        for(int v = 0 ; v < h ; v++ )
            for(int u = 0 ; u < w ; u++ )
            {
                unsigned id = u+v*w;
                const double x = (u - u0) / F;
                const double y = -(v - v0) / F;
                const double z = -1;
                const double c = depth[v][u] / sqrt(x*x + y*y + 1);
                vs[id] = makeVector(c*x, c*y, c*z);
                s << "," << endl << "<" << c*x << "," << c*y << "," << c*z << ">";
            }
        s << endl << "}" << endl;
    }

    { // Indices
        s << "face_indices{ " << 2*(w-1)*(h-1);
        for(int v = 0 ; v < (h-1) ; v++ )
            for(int u = 0 ; u < (w-1) ; u++ )
            {
                const int id = u + v*w;

                const double s1 = abs(depth[v][u] - depth[v][u+1]); // top
                const double s2 = abs(depth[v+1][u] - depth[v+1][u+1]);   // bottom
                const double s3 = abs(depth[v][u] - depth[v+1][u]);   // left
                const double s4 = abs(depth[v][u+1] - depth[v+1][u+1]); // right
                const double s5 = abs(depth[v][u+1] - depth[v+1][u]);   // diag

                if( max(s1,max(s3,s5)) < 5E-2 )
                {
                    s << ",<" << id << "," << id+1 << "," << id+w << ">";
                }else{
                    s << ",<" << 0 << "," << 1 << "," << w << ">";
                }
                if( max(s2,max(s4,s5)) < 5E-2 )
                {
                    s << ",<" << id+1 << "," << id+w << "," << id+w+1 << ">";
                }else{
                    s << ",<" << 1 << "," << w << "," << w+1 << ">";
                }
            }
        s << endl << "}" << endl;
    }

    s << "}" << endl;
}

int main(int argc, char* argv[])
{
        if( argc != 4 )
        {
            cerr << "Usage: " << argv[0] << " width height rawdepthfile" << endl;
            return -1;
        }

        const int w   = atoi(argv[1]);
        const int h   = atoi(argv[2]);
        const string fin(argv[3]);
        const string fpng = fin + ".png";
        const string fpov = fin + ".inc";

        ifstream ifile(fin.c_str());
        if( ifile.is_open() )
        {
            CVD::Image<double>img(ImageRef(w,h));
            readdoubleimage(ifile, img);
            ifile.close();

            CVD::Image<byte>img_norm(ImageRef(w,h));
            normalise(img, img_norm);
            img_save(img_norm,fpng);

            ofstream ofile(fpov.c_str());
            PovrayMesh(img, ofile);
            ofile.close();

        }else{
            cerr << "Couldn't open file" << endl;
        }

        return 0;
}
