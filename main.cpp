//
//  main.cpp
//  DOHW1
//
//  Created by Jessica Hoffmann on 10/03/14.
//  Copyright (c) 2014 ENS. All rights reserved.
//

#include <iostream>
#include <limits.h>
#include <vector>
#include "graph.h"
#include "lodepng.h"

using namespace std;

const int lambdaEdge = 10;
const int lambdaNonEdge = 20;
const int deltaI = 8;
const int Inf = numeric_limits<int>::max();

const int n = 288;
const int m = 384;
const int L = 16;

typedef Graph<int,int,int> GraphType;

void readImage(const char* filename, vector<vector<int> > &imRead);

void writeImage(char* filename, const vector<vector<int> > &imToWrite);

int vp(const vector<vector<int> > &imRight, const vector<vector<int> > &imLeft, int i,int j,int d);

int vpq(const vector<vector<int> > &imRight, int i1,int j1, int i2, int j2, int d1, int d2);

int computeEnergy(const vector<vector<int> > &imRight, const vector<vector<int> > &imLeft, const vector<vector<int> > &disparity);

int computeWpq(int nbToCompare, int treshold, int lambda1, int lambda2);

void disparityFromFlow(GraphType* g, vector<vector<int> > &disparity);

void disparityFromFlowAlpha(GraphType* g, vector<vector<int> > &disparity, int label);

void exactOptim(const vector<vector<int> > &imRight, const vector<vector<int> > &imLeft, vector<vector<int> > &disparity);

void alphaExpansion(const vector<vector<int> > &imRight, const vector<vector<int> > &imLeft, vector<vector<int> > &disparity);





void readImage(const char* filename, vector<vector<int> > &imRead)
{
    std::vector<unsigned char> image; //the raw pixels
    unsigned width, height;
    
    //decode
    unsigned error = lodepng::decode(image, width, height, filename);
    
    //if there's an error, display it
    if(error) std::cout << "decoder error " << error << ": " << lodepng_error_text(error) << std::endl;
    
    
    //convert to gray scale
    for (int j=0; j<m ; j++) {
        for (int i=0; i<n ; i++) {
            imRead[i][j] = 0.2989*image[4*(i*m + j)] + 0.5870*image[4*(i*m + j)+1] + 0.1140*image[4*(i*m +j) + 2];
        }
    }
    cout << "Image loaded. \n";
}

void writeImage(char* filename, const vector<vector<int> > &imToWrite) {
    vector<unsigned char> raw;
    
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            for(int color = 0; color < 3; color++) {
                raw.push_back((char)(imToWrite[i][j]));
            }
            raw.push_back(255);
        }
    }
    
    unsigned error = lodepng::encode(filename, raw, m, n);
    if(error) {
        cout << "encoder error " << error << ": "<< lodepng_error_text(error) << endl;
        exit(0);
    }
    cout << "Image written. \n";
}


int vp(const vector<vector<int> > &imRight, const vector<vector<int> > &imLeft, int i,int j,int d) {
    return abs(imLeft[i][j+d] - imRight[i][j]);
}


int vpq(const vector<vector<int> > &imLeft, int i1,int j1, int i2, int j2, int d1, int d2) {
    return computeWpq(abs(imLeft[i1][j1]-imLeft[i2][j2]), deltaI, lambdaEdge, lambdaNonEdge)*min(abs(d1-d2),2);
}


int computeEnergy(const vector<vector<int> > &imRight, const vector<vector<int> > &imLeft, const vector<vector<int> > &disparity) {
    int energy = 0;
    for (int i=0; i<n ; i++) {
        for (int j=0; j<m ; j++) {
            energy += vp(imRight, imLeft, i, j, disparity[i][j]/16);
            if (i<n-1) {
                energy += vpq(imLeft, i, j, i+1, j, disparity[i][j], disparity[i+1][j]/16);
            }
            if (j<m-1) {
                energy += vpq(imLeft, i, j, i, j+1, disparity[i][j], disparity[i][j+1]/16);
            }
        }
    }
    return energy;
}


int computeWpq(int nbToCompare, int treshold, int lambda1, int lambda2)
{
    return nbToCompare > treshold ? lambda1 : lambda2;
}


void disparityFromFlow(GraphType* g, vector<vector<int> > &disparity) {
    int label(0);
    for (int i=0; i<n ; i++) {
        for (int j=0; j<m ; j++) {
            label=0;
            while ((g->what_segment(label*n*m + i*m + j) == GraphType::SOURCE) && (label<L-1)) {
                label++;
            }
            disparity[i][j] = (label-1)*16;
        }
    }
    return;
}


void disparityFromFlowAlpha(GraphType* g, vector<vector<int> > &disparity, int label) {
     for (int i=0; i<n ; i++) {
        for (int j=0; j<m ; j++) {
            if (g->what_segment(i*m + j) == GraphType::SINK) {
                disparity[i][j] = (label-1)*16;
            }
        }
    }
    return;
}


void exactOptim(const vector<vector<int> > &imRight, const vector<vector<int> > &imLeft, vector<vector<int> > &disparity) {
    int wpq(0);
    int diffI(0);
    
    GraphType *g = new GraphType(n*m*L+2, 2*(n-1)*(m-1)*L);
    
    for (int i=0; i<n*m*L; i++) {
        g -> add_node();
    }
    
    for (int label=0;label < L; label++) {
        for (int i=0;i < n; i++) {
            for (int j=0;j < m; j++) {
                diffI = abs(imLeft[i][j+label] - imRight[i][j]) + 2*m*n*L*lambdaEdge;
                assert(diffI >= 0);
                
                // add vertical edges
                if (label==0) { 
                    // connected to SOURCE
                    g -> add_tweights(i*m + j, Inf, 0);
                }
                if (label==L-1) {
                    // connected to SINK
                    g -> add_tweights(label*m*n + i*m + j, 0, diffI);
                }
                else {
                    // normal edge
                    g -> add_edge(label*m*n + i*m + j, (label+1)*m*n + i*m + j, diffI, diffI);
                }
                
                // add horizontal edges
                if (i < n-1) {
                    wpq = computeWpq(abs(imLeft[i][j] - imLeft[i+1][j]), deltaI, lambdaEdge, lambdaNonEdge);
                    assert(wpq >= 0);
                    g -> add_edge(label*m*n + i*m + j, label*m*n + (i+1)*m + j, wpq, wpq);
                }
                if (j < m-1) {
                    wpq = computeWpq(abs(imLeft[i][j] - imLeft[i][j+1]), deltaI, lambdaEdge, lambdaNonEdge);
                    assert(wpq >= 0);
                    g -> add_edge(label*m*n + i*m + j, label*m*n + i*m + j+1, wpq, wpq);
                }
                
            }
        }
    }
    
    
    int flow = g -> maxflow();
    cout << "Flow : " << flow << endl ;
    cout << "Max-flow computed. \n";
    
    disparityFromFlow(g, disparity);
}

void alphaExpansion(const vector<vector<int> > &imRight, const vector<vector<int> > &imLeft, vector<vector<int> > &disparity) {
    bool change = true;
    int energyAct = computeEnergy(imRight, imLeft, disparity);
    int newEnergy(0), compteur(0);
    vector<vector<int> > disparityPot = disparity;
    
    
    // alpha-expansion
    while (change) {
        compteur++;
        change = false;
        for (int label=0; label<L; label++) {
            // create the graph with the right number of nodes
            GraphType *g = new GraphType(n*m + 2*(n-1)*(m-1), n*m + 2*(n-1)*(m-1));
            
            for (int i=0; i<n*m + n*(m-1) + (n-1)*m; i++) {
                /* 
                 node between 0 and n*m - 1, denoted by i*m + j : node "p", corresponding to the pixel (i,j)
                 node between n*m and n*m + n*(m-1) - 1 : node "alpha_pq", between the pixel (i,j) and the pixel (i+1,j)
                 node between n*m + n*(m-1) and n*m + n*(m-1) + m*(n-1) - 1 : node "alpha_pq", between the pixel (i,j) and the pixel (i,j+1)
                 */
                g -> add_node();
            }
            
            for (int i=0; i<n ; i++) {
                for (int j=0; j<m ; j++) {
                    g -> add_tweights(i*m + j, vp(imRight, imLeft, i, j, label), vp(imRight, imLeft, i, j, disparity[i][j]/16));
                    if (i<n-1) {
                        g -> add_edge(i*m + j, n*m + i*(m-1) + j, vpq(imRight, i, j, i+1, j,disparity[i][j]/16, label), vpq(imRight, i, j, i+1, j,disparity[i][j]/16, label));
                        g -> add_edge(n*m + i*(m-1) + j, (i+1)*m + j, vpq(imRight, i, j, i+1, j, disparity[i+1][j]/16, label), vpq(imRight, i, j, i+1, j, disparity[i+1][j]/16, label));
                        g -> add_tweights(n*m + i*(m-1) + j, 0, vpq(imRight,i, j, i+1, j, disparity[i][j]/16, disparity[i+1][j]/16));
                    }
                    if (j<m-1) {
                        g -> add_edge(i*m + j, n*m + i*(m-1) + j, vpq(imRight, i, j, i, j+1,disparity[i][j]/16, label), vpq(imRight, i, j, i, j+1,disparity[i][j]/16, label));
                        g -> add_edge(n*m + i*(m-1) + j, i*m + j+1, vpq(imRight, i, j, i, j+1, disparity[i][j+1]/16, label), vpq(imRight, i, j, i, j+1, disparity[i][j+1]/16, label));
                        g -> add_tweights(n*m + n*(m-1) + j*(n-1) + i, 0, vpq(imRight,i, j, i, j+1, disparity[i][j]/16, disparity[i][j+1]/16));
                    }
                }
            }
            int flow = g -> maxflow();
            
            disparityFromFlowAlpha(g, disparityPot, label);
            
            newEnergy = computeEnergy(imRight, imLeft, disparityPot);
            if (energyAct > newEnergy) {
                change = true;
                energyAct = newEnergy;
                disparity = disparityPot;
            }
            
            delete g;
        }
        cout << "Iteration " << compteur << ".\n";
    }
    cout << "Alpha-expansion done.\n";
    cout << "Energy : " << energyAct << ".\n";
}


int main (int argc, const char * argv[])
{
    vector<vector<int> > imRight(n,0), imLeft(n,0), disparity(n,0), disparity2(n,0);
    for (int i=0; i<n ; i++) {
        imRight[i].resize(m);
        imLeft[i].resize(m+L);
        disparity[i] .resize(m);
        disparity2[i] .resize(m);
    }
    
    
    cout << "Start. \n";
    
    
    readImage("/Users/jessicahoffmann/Documents/Cours/L3info/Codage/Cpp/DOHW1/DOHW1/imRight.png", imRight);
    
    readImage("/Users/jessicahoffmann/Documents/Cours/L3info/Codage/Cpp/DOHW1/DOHW1/imLeft.png", imLeft);
    
    
    for (int label=0; label<L; label++) {
        for (int i=0; i<n ; i++) {
            imLeft[i][m+label] = imLeft[i][m-1];
        }
    }
    
    exactOptim(imRight, imLeft, disparity);
    
    cout << "Energy : " << computeEnergy(imRight, imLeft, disparity) << ".\n";
    
    writeImage("/Users/jessicahoffmann/Documents/Cours/L3info/Codage/Cpp/DOHW1/DOHW1/res1.png", disparity);
    
    cout << "End of exo 1. \nStart exo 2.\n";
    
    
    alphaExpansion(imRight, imLeft, disparity);
    
    writeImage("/Users/jessicahoffmann/Documents/Cours/L3info/Codage/Cpp/DOHW1/DOHW1/alphaExpension.png", disparity);
    
    alphaExpansion(imRight, imLeft, disparity2);
    
    writeImage("/Users/jessicahoffmann/Documents/Cours/L3info/Codage/Cpp/DOHW1/DOHW1/alphaExpension2.png", disparity2);
    
    
        
    return 0;
}

