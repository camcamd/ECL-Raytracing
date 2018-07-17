#include <iostream>
#include "CImg.h"
#include <vector>
#include <algorithm>
#include <random>



using namespace std;
std::default_random_engine engine((unsigned int)time(0));
std::uniform_real_distribution<double> distrib(0,1);


// Définition de la classe vecteur
class Vector {

public:
    Vector(){};
    Vector(double x, double y, double z) {
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
    };
    Vector(const Vector &b) {
        memcpy(xyz, b.xyz, 3*sizeof(double));
    }


    double operator[](int i) const {
        return xyz[i];
    }
    double squaredNorm() {
        return xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2];
    }
    void Normalize(){
        double n = sqrt(squaredNorm());
        xyz[0]/=n;
        xyz[1]/=n;
        xyz[2]/=n;
    }
    double xyz[3];

};

// Définition des opérations sur les vecteurs
Vector operator+(const Vector &a, const Vector &b) {
    return Vector(a[0]+b[0], a[1]+b[1], a[2]+b[2]) ;
};

Vector operator-(const Vector &a, const Vector &b) {
    return Vector(a[0]-b[0], a[1]-b[1], a[2]-b[2]) ;
};

Vector operator*(const Vector &a, const Vector &b) {
    return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]) ;
};

Vector operator*(double a, const Vector &b) {
    return Vector(a*b[0], a*b[1], a*b[2]) ;
};

Vector operator/(const Vector &b, double a) {
    return Vector(b[0]/a, b[1]/a, b[2]/a) ;
};


double dot(const Vector &a, const Vector &b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
};

// Définition de la classe Ray, défini par une origine et un vecteur directeur
class Ray {

public:
    Vector O; // L'origine
    Vector u; // La direction
    Ray(Vector O, Vector u){
        this->O = O;
        this->u = u;
    };

};


// Définition de la classe Material
class Material {

public:
    int specularite; // 0 = opaque, 1=spéculaire, 2=transparent
    double indice;
    Vector couleur;
    int emissivite;
    Material(int specul = 0, double ind = 1., Vector coul = Vector(1.,1.,1.), int emissiv = 0){
        this->specularite = specul;
        this->indice = ind;
        this->couleur = coul;
        this->emissivite = emissiv;
    };

};



// Définition de la classe Sphère, une origine, un rayon et une matière
class Sphere {

public:
    Vector O; // Le centre
    int R; // Le rayon
    Material mat; // la matière

    Sphere(){};
    Sphere(Vector O, int R, const Material& mat){
        this->O = O;
        this->R = R;
        this->mat = mat;

    };

    bool IntersectRay(Ray ray, double &t){
        Vector u = ray.u;
        Vector C = ray.O;
        Vector O = this->O;
        int R = this->R;

        double a = 1.;
        double b = 2.*dot(u, C-O);
        double c = (C-O).squaredNorm() - R*R;
        double delta = b*b - 4*a*c;

        if (delta>=0.) {
            t = (-b - sqrt(delta)) / (2. * a);

            if (t >= 0.) {
                return true;
            }
            else {
                t = (-b + sqrt(delta)) / (2. * a);
                if (t > 0.) {
                    return true;
                }
                else {
                    return false;
                }
            }
        }
        else {
            return false;
        }
    };

};

// Définition de la classe Scène, vecteur de Sphères
class Scene {

public:
    std::vector<Sphere> vsphere;

    Scene(std::vector<Sphere> vsphere){
        this->vsphere = vsphere;
    };


    bool IntersectScene(Ray ray, double &tmin, int &indice){
        int i;
        double test;
        indice = 0;
        tmin = 100000.;
        for (i=0;i<vsphere.size();i++){
            if (vsphere[i].IntersectRay(ray, test)==true){
                if (tmin>test){
                    tmin = test;
                    indice = i ;
                }
            }
        }
        if (tmin==100000.){
            return false;
        }
        else {return true;}
    };


};

class Lumiere{
public:
    Vector L; // Position
    double I; // Intensité
    double gamma; // coeff gamma

    Lumiere(Vector L, double I, double gamma){
        this->I=I;
        this->L=L;
        this->gamma = gamma;
    }
};


Ray Reflexion(Ray ray, Vector P, Vector n){
    Vector dirReflec = ray.u - 2.* dot(ray.u,n)*n;
    dirReflec.Normalize();
    Ray reflecP(P,dirReflec);
    return reflecP;
};

Ray Refraction(Ray ray, Vector P, Vector n, Sphere sph){

    // Check du sens de la normale de la sphère et adaptation sinon
    double sens = (dot(ray.u,n) < 0.) - (dot(ray.u,n) > 0.); // renvoie 1 si bon sens, -1 si mauvais
    n = sens*n;

    // Si sens normale = 1 on calcule nair/nsph si -1 on calcule nsph/nair
    double ratioInd = pow(1./sph.mat.indice, sens);
    double racine = 1 - ratioInd*ratioInd*(1 - dot(ray.u,n)*dot(ray.u,n));

    if (racine<0.){
        Vector dirReflexTot = ray.u - 2.*dot(ray.u,n)*n;
        dirReflexTot.Normalize();
        Ray absorbP(P+0.001*n,dirReflexTot);
        return absorbP;
    }

    else {
        Vector dirRefract = ratioInd * (ray.u - dot(ray.u,n)*n) - sqrt(racine)*n;
        dirRefract.Normalize();
        Ray refractP(P-0.1*n,dirRefract);
        return refractP;
    }


};


Vector RandomVector(Vector n){


    double approxPi = 3.1415926535;
    double r1 = distrib(engine);
    double r2 = distrib(engine);

    Vector I1 = {cos(2.*approxPi*r1)*sqrt(1.-r2),sin(2.*approxPi*r2)*sqrt(1.-r2),sqrt(r2)};
    Vector t1;
    if (n[0]!=0 || n[1]!=0){ // calcul première tangente
        t1 = {n[1],-n[0],0.};
    }
    else {
        t1 = {0.,n[2],-n[1]};
    }

    Vector t2 = {n[1]*t1[2]-n[2]*t1[1],n[2]*t1[0]-n[0]-t1[2],n[0]*t1[1]-n[1]*t1[0]}; // produit vectoriel pour deuxième tangente

    Vector I2 = {t1[0]*I1[0]+t2[0]*I1[1]+n[0]*I1[2],t1[1]*I1[0]+t2[1]*I1[1]+n[1]*I1[2],t1[2]*I1[0]+t2[2]*I1[1]+n[2]*I1[2]};
    I2.Normalize();

    return I2;

};


Vector getColor(Ray ray, Scene scene, Lumiere light, int n_rebond){

    // Initialisation
    Vector couleur(0.,0.,0.);
    int indice;
    double t = 100000.;

    if (scene.IntersectScene(ray, t, indice)==true){


        Vector P(t * ray.u[0] + ray.O[0], t * ray.u[1] + ray.O[1], t * ray.u[2] + ray.O[2]);
        Vector n(P - scene.vsphere[indice].O);
        n.Normalize();

        Vector alealum = RandomVector(n);

        Vector L = scene.vsphere[0].O + scene.vsphere[0].R * alealum;
        Vector l(L - P);
        double d = dot(l, l);
        l.Normalize();



        double epsilon = 0.001;
        Ray raybis(P+epsilon*n,l);
        double tbis;
        int indbis;
        double ombre = 1.;

        Vector alea = RandomVector(n);

        //Calcul de l'ombre, s'il y a une ombre on renvoie tout simplement du noir
        if (scene.IntersectScene(raybis,tbis,indbis)==true && tbis*tbis <d && indbis!=0){ // si indbis = 0 c'est la sphère lumière
            ombre = 0.;
        }


        // Ensuite, calcul de la couleur selon matériau opaque ou miroir
        if (n_rebond > 0) {
            if (scene.vsphere[indice].mat.specularite == 1){ // Si surface miroir
                couleur = getColor(Reflexion(ray,P+epsilon*n,n),scene,light,n_rebond-1);
            }
            else if (scene.vsphere[indice].mat.specularite == 0){ // Si surface diffuse
                double lum = (std::max(0., dot(l, n))) * light.I / d;

                Vector PL = P - scene.vsphere[0].O;
                PL.Normalize();

                couleur = lum * ombre * dot(PL,alealum)/dot(-1*l, alealum) * scene.vsphere[indice].mat.couleur
                          + lum * scene.vsphere[indice].mat.couleur * getColor(Ray(P+epsilon*n, alea),scene,light,n_rebond-1);

            }

            else if(scene.vsphere[indice].mat.specularite == 2){
                couleur = getColor(Refraction(ray,P,n,scene.vsphere[indice]),scene,light,n_rebond-1);
            }
        }
        else {
            couleur = {0.,0.,0.};
            return couleur;
        }
    }

    return couleur;

};



int main(int argc, const char* argv[]) {
    int W = 1024;
    int H = 1024;
    double I = 2000.;
    double gamma = 2.2;
    std::vector<unsigned char> pixels(W*H*3, 0);



    Vector C(0,0,55);
    Vector L(-10,20,40);

    Lumiere light(L,I,gamma);
    Material matlum(0,1.,Vector(1.,1.,1.),2000);
    Sphere lumiere(L,1,matlum);

    // On définit la première sphère (centre)
    int R = 10;
    Vector O(0,0,0);
    Material mat(0, 1., Vector(1.,1.,1.),0);
    Sphere S1(O,R,mat);

    // On définit la deuxième sphère (sol)
    Vector O2(0,-1000,0);
    int R2 = 990;
    Material mat2(0, 1.,Vector(0.,0.,1.),0);
    Sphere S2(O2, R2, mat2);

    //Les autres sphères
    Vector O3(0,0,1000); // Devant
    int R3 = 940;
    Material mat3(0, 1.,Vector(1.,0.,1.),0);
    Sphere S3(O3, R3, mat3);

    Vector O4(0,0,-1000); // Derrière
    int R4 = 940;
    Material mat4(0, 1., Vector(0.,1.,0.),0);
    Sphere S4(O4, R4, mat4);

    Vector O5(0,1000,0); // Plafond
    int R5 = 940;
    Material mat5(0, 1., Vector(1.,0.,0.),0);
    Sphere S5(O5, R5, mat5);

    Vector O6(1000,0,0);
    int R6 = 940;
    Material mat6(0, 1., Vector(1.,1.,0.),0);
    Sphere S6(O6, R6, mat6);

    Vector O7(-1000,0,0);
    int R7 = 940;
    Material mat7(0, 1., Vector(0.,1.,1.),0);
    Sphere S7(O7, R7, mat7);


    vector<Sphere> v(8);
    v[0] = lumiere;
    v[1] = S1;
    v[2] = S2;
    v[3] = S3;
    v[4] = S4;
    v[5] = S5;
    v[6] = S6;
    v[7] = S7;
    Scene scene(v);

    double fov = 60*3.14/180.;
    double approxPi = 3.1415926535;

    int n = 30; // nombre de rayons générés par pixel pour les BRDFs


        for (int i = 0; i < H; i++) {

            for (int j = 0; j < W; j++) {

                Vector couleur(0., 0., 0.);
                Vector u(j - W / 2., H / 2. - i, -W / (2. * tan(fov / 2.)));
                u.Normalize();

                Ray ray(C, u);
                for (int k = 0; k < n; k++) {

                    couleur = couleur + getColor(ray, scene, light, 3);
                }

                couleur = 1. / n * couleur;


                pixels[i * W + j] = (std::min(255., 255. * pow(couleur[0], 1. / light.gamma)));
                pixels[i * W + j + W * H] = (std::min(255., 255. * pow(couleur[1], 1. / light.gamma)));
                pixels[i * W + j + 2 * W * H] = (std::min(255., 255. * pow(couleur[2], 1. / light.gamma)));


            }


        }
        cimg_library::CImg<unsigned char> cimg(&pixels[0], W, H, 1, 3);
        cimg.save("DevauxCamille.bmp");


}