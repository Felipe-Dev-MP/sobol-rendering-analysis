// Adaptacao para o projeto de CG. Cena e transporte de luz de Kevin Beason.
// Original: https://www.kevinbeason.com/smallpt/ (MIT).
// Compilar no MSYS2 MSYS:
// g++ -std=gnu++17 -O2 smallpt_pesquisa.cpp -o smallpt_pesquisa.exe
// Exemplo: ./smallpt_pesquisa.exe sobol 16 sobol_16.ppm 320 240 1
// Baseline deste experimento: modo sobol. random e apenas um controle auxiliar.
// J/P: suplemento do artigo de Abdalla G. M. Ahmed, EGSR 2024,
// DOI 10.2312/sr.20241147. Mantidos em amostragem.h.
// LP/MMD e exportacao de radiancia linear ainda pendentes.
// Os dois modos usam sementes por amostra e nao reproduzem bit a bit o
// fluxo pseudoaleatorio por linha do smallpt original.
#include "amostragem.h"
#include <cstring>
#include <cerrno>
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
  double x, y, z;                  // position, also color (r,g,b)
  Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
  Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
  Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
  Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
  Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
  Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
  double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross:
  Vec operator%(Vec&b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
};
struct Ray { Vec o, d; Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
struct Sphere {
  double rad;       // radius
  Vec p, e, c;      // position, emission, color
  Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
  Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
    rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
  double intersect(const Ray &r) const { // returns distance, 0 if nohit
    Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
    double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
    if (det<0) return 0; else det=sqrt(det);
    return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
  }
};
Sphere spheres[] = {//Scene: radius, position, emission, color, material
  Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
  Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
  Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
  Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF),//Frnt
  Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
  Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
  Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
  Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
  Sphere(600, Vec(50,681.6-.27,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite
};
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }
inline bool intersect(const Ray &r, double &t, int &id){
  double n=sizeof(spheres)/sizeof(Sphere), d, inf=t=1e20;
  for(int i=int(n);i--;) if((d=spheres[i].intersect(r))&&d<t){t=d;id=i;}
  return t<inf;
}
Vec radiance(const Ray &r, int depth, unsigned short *Xi){
  double t;                               // distance to intersection
  int id=0;                               // id of intersected object
  if (!intersect(r, t, id)) return Vec(); // if miss, return black
  const Sphere &obj = spheres[id];        // the hit object
  Vec x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c;
  double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
  if (++depth>5) if (erand48(Xi)<p) f=f*(1/p); else return obj.e; //R.R.
  if (obj.refl == DIFF){                  // Ideal DIFFUSE reflection
    double r1=2*M_PI*erand48(Xi), r2=erand48(Xi), r2s=sqrt(r2);
    Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
    Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
    return obj.e + f.mult(radiance(Ray(x,d),depth,Xi));
  } else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
    return obj.e + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth,Xi));
  Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
  bool into = n.dot(nl)>0;                // Ray from outside going in?
  double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
  if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
    return obj.e + f.mult(radiance(reflRay,depth,Xi));
  Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
  double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
  double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
  return obj.e + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
    radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
    radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
}

// Argumentos: metodo spp arquivo largura altura semente
int inteiro(const char *texto, int limite) {
  char *fim = 0; errno = 0;
  long valor = strtol(texto, &fim, 10);
  if (errno || !*texto || *fim || valor < 1 || valor > limite) return -1;
  return int(valor);
}
int main(int argc, char *argv[]) {
  if (argc != 7 || (strcmp(argv[1], "sobol") && strcmp(argv[1], "random"))) {
    fprintf(stderr, "Uso: %s sobol|random spp saida.ppm largura altura semente\n", argv[0]);
    return 1;
  }
  int spp = inteiro(argv[2], 1048576), w = inteiro(argv[4], 4096);
  int h = inteiro(argv[5], 4096), seed = inteiro(argv[6], 2147483647);
  if (spp < 4 || (spp & (spp-1)) || w < 2 || h < 2 || seed < 1) {
    fprintf(stderr, "Use spp potencia de 2 >= 4, dimensoes entre 2 e 4096 e semente positiva.\n");
    return 1;
  }
  FILE *f = fopen(argv[3], "w");
  if (!f) { perror(argv[3]); return 1; }
  const bool usarSobol = strcmp(argv[1], "sobol") == 0;
  const int samps = spp / 4; // Original: quatro subpixels, com filtro tent.
  Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm());
  Vec cx=Vec(w*.5135/h), cy=(cx%cam.d).norm()*.5135;
  Vec *c = new Vec[w*h];
  fprintf(stderr, "Metodo=%s, %dx%d, %d spp (%d por subpixel), semente=%d\n",
          argv[1], w, h, spp, samps, seed);
  for (int y=0; y<h; ++y) {
    fprintf(stderr, "\rRenderizando %5.1f%%", 100.*y/(h-1));
    for (int x=0; x<w; ++x) {
      const int i=(h-y-1)*w+x;
      for (int sy=0; sy<2; ++sy) for (int sx=0; sx<2; ++sx) {
        Vec r;
        const std::uint32_t dominio = std::uint32_t((y*w+x)*4+sy*2+sx);
        const std::uint32_t chave = amostragem::mistura(dominio ^ amostragem::mistura(seed));
        for (int s=0; s<samps; ++s) {
          double u, v;
          if (usarSobol) {
            amostragem::sobol(s, chave, u, v);
          } else {
            unsigned short cameraXi[3];
            amostragem::estado(amostragem::mistura(chave ^ 0x12345678u) ^ amostragem::mistura(s), cameraXi);
            u = erand48(cameraXi); v = erand48(cameraXi);
          }
          // Mesma transformacao tent do smallpt original nos dois modos.
          double r1=2*u, dx=r1<1 ? sqrt(r1)-1 : 1-sqrt(2-r1);
          double r2=2*v, dy=r2<1 ? sqrt(r2)-1 : 1-sqrt(2-r2);
          Vec d = cx*(((sx+.5+dx)/2+x)/w-.5) + cy*(((sy+.5+dy)/2+y)/h-.5) + cam.d;
          // Transporte usa a mesma regra de sementes nos dois modos.
          // Separar este estado evita que trocar a camera consuma seus numeros.
          unsigned short Xi[3];
          amostragem::estado(amostragem::mistura(chave ^ 0x87654321u) ^ amostragem::mistura(s), Xi);
          r = r + radiance(Ray(cam.o+d*140,d.norm()),0,Xi)*(1./samps);
        }
        // Mantido do original. PPM tem clamp e gamma: nao usar para RMSE linear.
        c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
      }
    }
  }
  fprintf(f, "P3\n%d %d\n255\n", w, h);
  for (int i=0; i<w*h; ++i)
    fprintf(f, "%d %d %d\n", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
  const bool erro = ferror(f) != 0;
  const int fechou = fclose(f);
  delete[] c;
  if (erro || fechou) { fprintf(stderr, "Falha ao gravar imagem.\n"); return 1; }
  fprintf(stderr, "\nImagem salva em %s\n", argv[3]);
  return 0;
}
