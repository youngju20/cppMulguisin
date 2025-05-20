#include <string.h>
#include <math.h>
#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <time.h>
#include <numeric>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <map>
#include <set>

#include <H5Cpp.h>
#include "mulguisin_kd.h"

// Qhull
#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullPoint.h"
#include "libqhullcpp/QhullUser.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"

using std::cout;
using std::endl;
using std::cerr;

using orgQhull::Qhull;
using orgQhull::QhullError;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullFacetListIterator;
using orgQhull::QhullFacetSet;
using orgQhull::QhullFacetSetIterator;
using orgQhull::QhullPoint;
using orgQhull::QhullPoints;
using orgQhull::QhullPointsIterator;
using orgQhull::QhullQh;
using orgQhull::QhullUser;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexList;
using orgQhull::QhullVertexListIterator;
using orgQhull::QhullVertexSet;
using orgQhull::QhullVertexSetIterator;
using orgQhull::RboxPoints;

using namespace H5;

// Compute tetrahedron
double tetravol(const std::vector<std::vector<double>>& pos) {
    std::vector<double> a(3), b(3), c(3);
    for (int i = 0; i < 3; ++i) {
        a[i] = pos[1][i] - pos[0][i];
        b[i] = pos[2][i] - pos[0][i];
        c[i] = pos[3][i] - pos[0][i];
    }
    std::vector<double> cross(3);
    cross[0] = b[1]*c[2] - b[2]*c[1];
    cross[1] = b[2]*c[0] - b[0]*c[2];
    cross[2] = b[0]*c[1] - b[1]*c[0];
    double dot = a[0]*cross[0] + a[1]*cross[1] + a[2]*cross[2];
    return std::abs(dot) / 6.0;
}
// ---------------------------------------------------

// Compute volume
double computeVoronoiVolume_fromVerticesOnly(
    const std::vector<int>& region_indices,
    const std::vector<std::vector<double>>& vertices_pos)
{
    // If it has infinite vertice, then it can not calculate
    if (std::find(region_indices.begin(), region_indices.end(), 0) != region_indices.end()) {
        return 0.0;
    }

    // Extract vertices of region
    std::vector<std::vector<double>> region_vertices;
    for (int idx : region_indices) {
        region_vertices.push_back(vertices_pos[idx]);
    }

    int n = region_vertices.size();
    if (n < 4) return 0.0;  // can not make tetrahedron

    // calculate centroid in region
    std::vector<double> centroid(3, 0.0);
    for (const auto& v : region_vertices) {
        for (int d = 0; d < 3; ++d) {
            centroid[d] += v[d];
        }
    }
    for (int d = 0; d < 3; ++d) {
        centroid[d] /= n;
    }

    // Simple way : make tetraheron using centroid and three vertices
    double volume = 0.0;
    for (int i = 1; i < n - 1; ++i) {
        std::vector<std::vector<double>> tetra = {
            centroid,
            region_vertices[0],
            region_vertices[i],
            region_vertices[i+1]
        };
        try {
            volume += tetravol(tetra);
        } catch (...) {
            continue;
        }
    }

    return volume;
}
// ---------------------------------------------------

// Compute voronoi
void qvoronoi_o(const Qhull &qhull,
                std::vector<std::vector<double>> &voronoiVertices,
                std::vector<std::vector<int>> &voronoiRegions)
{
    int voronoiDimension = qhull.hullDimension() - 1;
    int numfacets = qhull.facetCount();
    size_t numpoints = qhull.points().size();

    // 1. Gather voronoi vertices
    std::map<int, int> facetIdToVertexIndex;
    voronoiVertices.clear();

    // Infinite vertex is in index 0
    std::vector<double> vertexAtInfinity(voronoiDimension, qh_INFINITE);
    voronoiVertices.push_back(vertexAtInfinity);

    int vertexIndex = 1;
    QhullFacetListIterator facetIt(qhull.facetList());
    while (facetIt.hasNext()) {
        QhullFacet facet = facetIt.next();
        if (!facet.isGood()) continue;

        voronoiVertices.push_back(facet.getCenter().toStdVector());
        facetIdToVertexIndex[facet.id()] = vertexIndex;
        ++vertexIndex;
    }

    // 2. Gather voronoi region
    voronoiRegions.clear();
    voronoiRegions.resize(numpoints);

    QhullVertexListIterator vertexIt(qhull.vertexList());
    while (vertexIt.hasNext()) {
        QhullVertex vertex = vertexIt.next();
        std::vector<int> voronoiRegion;
        bool hasInfinite = false;

        QhullFacetSetIterator neighborIt(vertex.neighborFacets());
        while (neighborIt.hasNext()) {
            QhullFacet neighbor = neighborIt.next();
            if (!neighbor.isGood()) continue;

            if (neighbor.isUpperDelaunay()) {
                // handle infinite region
                if (!hasInfinite) {
                    hasInfinite = true;
                    voronoiRegion.push_back(0);  // infinite vertex index
                }
            } else {
                int fid = neighbor.id();
                if (facetIdToVertexIndex.find(fid) != facetIdToVertexIndex.end()) {
                    voronoiRegion.push_back(facetIdToVertexIndex[fid]);
                }
            }
        }

        if (!voronoiRegion.empty()) {
            int siteId = vertex.point().id();
            if (siteId >= 0 && siteId < int(numpoints)) {
                voronoiRegions[siteId] = voronoiRegion;
            }
        }
    }
}
// ---------------------------------------------------


//==============================================================================
// Run MulGuisin without ROOT
//==============================================================================

using namespace std;

const double pi = 3.141592;

// -------------------------------------------------------------------------------
// Required one to READ HDF5 file
// -------------------------------------------------------------------------------

double* ReadDoubleData(char* strFile,char* strDataset); // Used for 8-byte data
float* ReadFloatData(char* strFile,char* strDataset);   // Used for 4-byte dat
int* ReadIntData(char* strFile,char* strDataset);   // Used for 4-byte dat
int ReadDatasize(char* strFile,char* strDataset);   // Used for 4-byte dat

// -------------------------------------------------------------------------------

// -------------------------------------------------------------------------------
// Prepare data array
// -------------------------------------------------------------------------------

typedef struct GalaxyArray{double x,y,z;} GalaxyArray; 
typedef struct MGSArray{double sx,sy,sz,w;} MGSArray;

std::vector<GalaxyArray>* ReadGalaxy(char filename[128], int &nGalaxies)
{
   std::vector<GalaxyArray>* galaxies=new vector<GalaxyArray>;

   //char filename[128] = "Mr20.0_1.hdf5";
   printf("filename : %s\n", filename);

   double* sx = ReadDoubleData(filename,(char*)"x");
   double* sy = ReadDoubleData(filename,(char*)"y");
   double* sz = ReadDoubleData(filename,(char*)"z");
   int Dsize = ReadDatasize(filename,(char*)"x");
   //int count = 0;

   int count = Dsize;

   nGalaxies=0;
   for(int i=0;i<count;i++)
   {
      if (i%10000==0) printf("Galaxy %8d : rsd=(% 06.2f,% 06.2f,% 06.2f)\n",\
            i,sx[i],sy[i],sz[i]);
      GalaxyArray galaxy;
      galaxy.x=sx[i];
      galaxy.y=sy[i];
      galaxy.z=sz[i];
      galaxies->push_back(galaxy);
      nGalaxies++;
   }
   delete[] sx;
   delete[] sy;
   delete[] sz;
   return galaxies;
}

  
// -------------------------------------------------------------------------------

int main(int argc, char *argv[])
{

  // ----------------------------------------------------------------------
  // Set initial parameters
  // ----------------------------------------------------------------------

  float fResolution1 = std::stof(argv[1]); // fResolution for 1st MGS
  //bool SortGrp = (std::stoi(argv[3])==1)?1:0; // sort groups by population  

  // ----------------------------------------------------------------------
  // READ data file (extension : txt)
  // ----------------------------------------------------------------------

  //int nElements;
  //char filename[128] = "/home/young/test_mgs/sim_sumi_200Mpc_2_230813.hdf5";
  //std::vector<GalaxyArray>* iarray = ReadGalaxy(filename, nElements);
  //std::vector<GalaxyArray>* oarray = new std::vector<GalaxyArray>;
  //printf("A total of %d galaxies were read in.\n",nElements);

  std::vector<GalaxyArray>* iarray = new std::vector<GalaxyArray>;

  double minx=1000.0;
  double miny=1000.0;
  double minz=1000.0;
  double maxx=0.0;
  double maxy=0.0;
  double maxz=0.0;
  int nGalaxies=0;
  cout << "Data READ" << endl;
  //std::ifstream file("N128L100/snapshot_006.txt");
  std::ifstream file(argv[2]);
  if (file.is_open()) {
     //for (int i=0; i<8; i++) { file.ignore(200, '\n'); }
     string line;
     int ic=0;
     while (getline(file, line)) {
        //cout << line << endl;
        std::istringstream iss(line);
        //string junk;
        double col1;
        double col2;
        double col3;
        double col4;
        //iss >> junk >> junk >> junk >> junk >> junk >> col1 >> col2 >> col3 >> junk >> junk >> junk >> junk;
        iss >> col1 >> col2 >> col3 >> col4;
        if (ic%1000000==0) {cout << ic << " " << col2 << " " << col3 << " " << col4 << endl;}
        if (col2 < minx) minx=col2;
        if (col2 > maxx) maxx=col2;
        if (col3 < miny) miny=col3;
        if (col3 > maxy) maxy=col3;
        if (col4 < minz) minz=col4;
        if (col4 > maxz) maxz=col4;
        GalaxyArray galaxy;
        galaxy.x = col2;
        galaxy.y = col3;
        galaxy.z = col4;
        iarray->push_back(galaxy);
        ic++;
     }
     nGalaxies = ic;
     file.close();
  }
  cout << "Total # of galaxies : " << nGalaxies << endl;

  cout << "min/max data : " << endl;
  cout << "x : " << minx << " " << maxx << endl;
  cout << "y : " << miny << " " << maxy << endl;
  cout << "z : " << minz << " " << maxz << endl;


  // =====================================================================
  // Process 1 : Estimate Density
  // =====================================================================
  // calculate voronoi volume
  // This process is for making density. In this code, the density is
  // inverse of voronoi volume. With this process, we can count a weight
  // of each galaxy.
  // ! Define 'container' to use 'voro++' library
  // 'voro++' is library for making voronoi tessellation and also can
  // ---------------------------------------------------------------------

  // -----------------------------------------------------------------------
  std::cout << "\nStart estimating Density\n" << std::endl;

  //float_t rx,ry,rz,sx,sy,sz,vx,vy,vz,m,u,w;
  float_t sx,sy,sz, w;
  std::vector<MGSArray>* arr=new std::vector<MGSArray>;

  //sub-boxes
  std::vector<std::vector<double>> galvor;
  //int newic = 0;
  for (int i=0; i<nGalaxies; i++){
      GalaxyArray gali=iarray->at(i);
      double x1=gali.x;
      double y1=gali.y;
      double z1=gali.z;
      galvor.push_back({x1, y1, z1});
      //if (x1 <= 10.0 && y1 <= 10.0 && z1 <= 10.0){
      //   galvor.push_back({x1, y1, z1});
      //   newic++;
      //}
  }
  delete iarray;

  cout << "TEST : size input : " << galvor.size() << endl;

  nGalaxies = galvor.size();

  // Put data into array about qhull
  clock_t start, end;
  double time_result;
  start = clock();

  const int numrow = galvor.size();
  const int numcol = galvor[0].size();
  std::vector<double> pointarray(numrow * numcol);
  for (int i = 0; i < numrow; ++i) {
      for (int j = 0; j < numcol; ++j) {
          pointarray[i * numcol + j] = galvor[i][j];
          //cout << pointarray[i * numcol + j] << " ";
      }
      //cout << endl;
  }

  cout << "array size " << pointarray.size() << endl;

  Qhull qhull;
  //qhull.runQhull("", 3, nGalaxies, pointarray.data(), "v Qbb FN");
  qhull.runQhull("voronoi", 3, nGalaxies, pointarray.data(), "d Qbb Qx");

  qhull.defineVertexNeighborFacets();

  // Get voronoi diagram information
  std::vector<std::vector<double>> vertices_pos;
  std::vector<std::vector<int>> vertices_id;
  qvoronoi_o(qhull, vertices_pos, vertices_id);

  for (int i = 0; i < nGalaxies; ++i) {
      if (i%10000==0) { cout << i << "  th galaxies is under density estimation" << endl; }

      const auto& region_indices = vertices_id[i];

    double volume = 0.0;
    // If it contain infinite vertex, then skip (unbounded)
    if (std::find(region_indices.begin(), region_indices.end(), 0) != region_indices.end()) {
        std::cout << "Region " << i << " is unbounded, skipping." << std::endl;
        volume = 0.0;
    } else {
        volume = computeVoronoiVolume_fromVerticesOnly(region_indices, vertices_pos);
        //std::cout << "Voronoi volume for point " << i << ": " << volume << std::endl;
    }

      MGSArray gal;
      gal.sx = galvor[i][0];
      gal.sy = galvor[i][1];
      gal.sz = galvor[i][2];
      gal.w = (volume > 0.0) ? 1.0 / volume : 1e-10;

      //cout << gal.w << endl;

      arr->push_back(gal);

  }

  end = clock();
  time_result = (double)(end-start);
  cout << "Time for estimating Voronoi volume : " << ((time_result)/CLOCKS_PER_SEC) << " Sec" << endl;
      
  // =====================================================================
  // Process 2 : Calculate MGS
  // =====================================================================

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // create Mulguising finder and initialize the MGS arrays 
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  std::cout << "\nStart finding MGS groups/clusters\n" << std::endl;
  cout << "Linking Length = " << fResolution1 << endl;

  MGS3DClusterFinder *mgs=new MGS3DClusterFinder();
  mgs->Init(nGalaxies);
  mgs->SetMode(2);	// grouping method (2: group points when d<Resolution)
  mgs->SetResolution(fResolution1); // distance cut to be grouped
  //mgs->SetMergingDistance(Resolution); // distance parameter for merging
  mgs->SetVerbose(0);	// print level (0: turn off, 1: turn on)

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // initialize input arrays
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  float *Gm=mgs->Gm;
  float *Gx=mgs->Gx;
  float *Gy=mgs->Gy;
  float *Gz=mgs->Gz;

  //double *Gm=mgs->Gm;
  //double *Gx=mgs->Gx;
  //double *Gy=mgs->Gy;
  //double *Gz=mgs->Gz;

  for (int i=0;i<nGalaxies;i++)
  {
    MGSArray gali=arr->at(i);
    double xi=gali.sx;
    double yi=gali.sy;
    double zi=gali.sz;
    double wi=gali.w;
    Gm[i]=wi;
    Gx[i]=xi;
    Gy[i]=yi;
    Gz[i]=zi;
  }

  delete arr;

  // ---------------------------------------------------------------------
  //clock_t start, end;
  //double time_result;
  start = clock();

  mgs->Run();

  end = clock();
  time_result = (double)(end-start);
  cout << "Time for calculating MGS : " << ((time_result)/CLOCKS_PER_SEC) << " Sec" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Print out each MGS cluster with the number of particles attached  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  int NbOfGrp=mgs->GetNbOfGrp();
  int *NbPntInGrp=mgs->NbOfPntInGrp;	// Number of points in a group
  std::vector<int> *IdOfPntInGrp=mgs->IdOfPntInGrp; //id of a point in a group

  printf("- - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf("Mulguisin found a total of %d groups.\n",NbOfGrp);
  printf("- - - - - - - - - - - - - - - - - - - - - - - - -\n");

  
  //for (int igrp=0; igrp<NbOfGrp;igrp++)
  for (int igrp=0; igrp<20;igrp++)
  {
    printf("The %d'th group has %d points\n",igrp,mgs->GetNbOfPntInGrp(igrp));
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Sort MGS clusters according to the number of particles associated
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  printf("- - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf("Groups are sorted by its population.\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - -\n");

  int *copyNbPntInGrp = new int[nGalaxies];
  for (int ii = 0; ii < nGalaxies; ii++){
    copyNbPntInGrp[ii] = NbPntInGrp[ii];
  }
  std::sort(copyNbPntInGrp, copyNbPntInGrp + nGalaxies, greater<int>());

  for (int isg=0; isg<20; isg++)
  {
       printf("The %d'th biggest group has %d points\n",isg,copyNbPntInGrp[isg]);
  }
  delete[] copyNbPntInGrp;

  // - - - - - - - - - - - - - - - - - - - - - -
  // Get label and link informations
  // - - - - - - - - - - - - - - - - - - - - - -

  // Get label for MGS, it is ordered as same as data order
  cout << "Estimate MGS label" << endl;
  mgs->GetLabel();
  int *label = mgs->MGSlabel;

  cout << "Estimate MGS link" << endl;
  mgs->GetLink();
  std::vector<GLink> *link = mgs->MGSLink;

  // -------------------------------------------------------------------------
  // Read MGS topological information
  // -------------------------------------------------------------------------

  cout << "Read MGS topological information" << endl;

  std::vector<float> Totlen;
  std::vector<float> Getlen;
  std::vector<float> GetOA;
  std::vector<float> GetPA;
  std::vector<int> GetGen;
  std::vector<int> Getchild;
  std::vector<int> GetBranch;
  std::vector<int> GetDegree;
  for (int i=0; i<NbOfGrp; i++){
      Totlen.push_back(mgs->GetTotalLength(i));
      GetBranch.push_back(mgs->GetNBranches(i));
      for (int j=0; j<NbPntInGrp[i]; j++){
          GetGen.push_back(mgs->GetGeneration(i,j));
          Getlen.push_back(mgs->GetLength(i,j));
          GetOA.push_back(mgs->GetOpeningAngle(i,j));
          GetPA.push_back(mgs->GetLinkPolarAngle(i,j));
          Getchild.push_back(mgs->GetNChildren(i,j));
          GetDegree.push_back(mgs->GetDegree(i,j));
      }
  }
  cout << "Done" << endl;


  // Define save file name
  H5File* wfile = new H5File(argv[3], H5F_ACC_TRUNC);

  // Save data for each dataspace

  // label
  hsize_t dims[1];
  dims[0] = nGalaxies;
  DataSpace dataspace(1,dims);

  const std::string dataname = "label";
  DataSet* dataset = new DataSet(wfile->createDataSet(
    dataname, PredType::NATIVE_INT, dataspace));
  dataset->write(label, PredType::NATIVE_INT);
  delete dataset;

  // TotLength
  hsize_t dims2[1];
  dims2[0] = NbOfGrp;
  DataSpace dataspace2(1,dims2);

  const std::string dataname2 = "Totlen";
  DataSet* dataset2 = new DataSet(wfile->createDataSet(
     dataname2, PredType::NATIVE_FLOAT, dataspace2));
  dataset2->write(Totlen.data(), PredType::NATIVE_FLOAT);
  delete dataset2;

  // Generation
  hsize_t dims3[1];
  dims3[0] = nGalaxies;
  DataSpace dataspace3(1,dims3);

  const std::string dataname3 = "Gen";
  DataSet* dataset3 = new DataSet(wfile->createDataSet(
     dataname3, PredType::NATIVE_INT, dataspace3));
  dataset3->write(GetGen.data(), PredType::NATIVE_INT);
  delete dataset3;

  // Branch
  hsize_t dim5[1];
  dim5[0] = NbOfGrp;
  DataSpace dataspace5(1,dim5);

  const std::string dataname5 = "nBranch";
  DataSet* dataset5 = new DataSet(wfile->createDataSet(
     dataname5, PredType::NATIVE_INT, dataspace5));
  dataset5->write(GetBranch.data(),PredType::NATIVE_INT);
  delete dataset5;

  // Child
  hsize_t dim6[1];
  dim6[0] = nGalaxies;
  DataSpace dataspace6(1,dim6);

  const std::string dataname6 = "nChild";
  DataSet* dataset6 = new DataSet(wfile->createDataSet(
     dataname6, PredType::NATIVE_INT, dataspace6));
  dataset6->write(Getchild.data(),PredType::NATIVE_INT);
  delete dataset6;

  // Degree
  hsize_t dim7[1];
  dim7[0] = nGalaxies;
  DataSpace dataspace7(1,dim7);

  const std::string dataname7 = "Deg";
  DataSet* dataset7 = new DataSet(wfile->createDataSet(
     dataname7, PredType::NATIVE_INT, dataspace7));
  dataset7->write(GetDegree.data(),PredType::NATIVE_INT);
  delete dataset7;

  // Opening Angle
  hsize_t dim8[1];
  dim8[0] = nGalaxies;
  DataSpace dataspace8(1,dim8);

  const std::string dataname8 = "OpenAngle";
  DataSet* dataset8 = new DataSet(wfile->createDataSet(
     dataname8, PredType::NATIVE_FLOAT, dataspace8));
  dataset8->write(GetOA.data(),PredType::NATIVE_FLOAT);
  delete dataset8;

  // Polar Angle
  hsize_t dim9[1];
  dim9[0] = nGalaxies;
  DataSpace dataspace9(1,dim9);

  const std::string dataname9 = "PolAngle";
  DataSet* dataset9 = new DataSet(wfile->createDataSet(
     dataname9, PredType::NATIVE_FLOAT, dataspace9));
  dataset9->write(GetPA.data(),PredType::NATIVE_FLOAT);
  delete dataset9;


  // link
  CompType mtype(sizeof(GLink));
  mtype.insertMember("idbra", HOFFSET(GLink, idbra), PredType::NATIVE_INT);
  mtype.insertMember("idket", HOFFSET(GLink, idket), PredType::NATIVE_INT);
  mtype.insertMember("imgs", HOFFSET(GLink, imgs), PredType::NATIVE_INT);

  hsize_t dim4[1];
  dim4[0] = link->size();
  int rank = sizeof(dim4) / sizeof(hsize_t);
  DataSpace dataspace4(rank, dim4);

  DataSet* dataset4 = new DataSet(wfile->createDataSet(
     "links", mtype, dataspace4));
  dataset4->write(link->data(), mtype);
  delete dataset4;

  delete wfile;

  
  /*
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Sort MGS clusters according to the number of particles associated
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if (SortGrp) mgs->SortGrp();

  printf("- - - - - - - - - - - - - - - - - - - - - - - - -\n");
  printf("Groups are sorted by its population.\n");
  printf("- - - - - - - - - - - - - - - - - - - - - - - - -\n");

  for (int isg=0; isg<NbOfGrp;isg++)
  {
    int igrp=mgs->GetIdOfGrpFromSortedGrp(isg);
    int npnt=mgs->GetNbOfPntInGrp(igrp);
    if (isg<10 ||
       (isg%10==0&&isg<100) ||
       (isg%100==0&&isg<1000))
       printf("The %d'th biggest group was originally the %d'th group and it has %d points\n",isg,igrp,npnt);
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Merge MGS clusters to create Super clusters
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  mgs->MergeGrp();

  int NbOfCls=mgs->GetNbOfCls();
  printf("Mulguisin found a total of %d clusters.\n",NbOfCls);

  for (int icls=0; icls<NbOfCls;icls++)
  {
    printf("Cluster %d: ",icls);
    int NbOfGrpInCls=mgs->GetNbOfGrpInCls(icls);
    for (int i=0;i<NbOfGrpInCls;i++)
    {
      int isg=mgs->GetIdOfGrpInCls(icls,i);
      int igrp=mgs->GetIdOfGrpFromSortedGrp(isg);
      int npnt=mgs->GetNbOfPntInGrp(igrp);
      for (int ipnt=0;ipnt<npnt;ipnt++)
      {
        int idp=mgs->GetIdOfPntInGrp(igrp,ipnt);
        double x=Gx[idp];
        double y=Gy[idp];
        double z=Gz[idp];
        double m=Gm[idp];
      }
      printf("%d(%d)-",igrp,npnt);
    }
    printf("\n");
  }
  */

  delete mgs;
  return 0;
}

//==============================================================================
// Read data modules
//==============================================================================

double* ReadDoubleData(char* strFile,char* strDataset) {

        H5File file(strFile,H5F_ACC_RDONLY);

        DataSet dataset = file.openDataSet(strDataset);
        DataSpace filespace = dataset.getSpace();
        int rank = filespace.getSimpleExtentNdims();

        hsize_t* dims = new hsize_t[rank];
        rank = filespace.getSimpleExtentDims(dims);

        DataSpace mspace(rank,dims);
        int fullsize = 1;
        for(int irank=0;irank<rank;irank++) {
                fullsize *= dims[irank];
        }
        delete[] dims;

        double* data_out = new double[fullsize];
        dataset.read(data_out,PredType::NATIVE_DOUBLE,mspace,filespace);

        return data_out;
}
float* ReadFloatData(char* strFile,char* strDataset) {

        H5File file(strFile,H5F_ACC_RDONLY);

        DataSet dataset = file.openDataSet(strDataset);
        DataSpace filespace = dataset.getSpace();
        int rank = filespace.getSimpleExtentNdims();

        hsize_t* dims = new hsize_t[rank];
        rank = filespace.getSimpleExtentDims(dims);

        DataSpace mspace(rank,dims);
        int fullsize = 1;
        for(int irank=0;irank<rank;irank++) {
                fullsize *= dims[irank];
        }
        delete[] dims;

        float* data_out = new float[fullsize];
        dataset.read(data_out,PredType::NATIVE_FLOAT,mspace,filespace);

        return data_out;
}

int* ReadIntData(char* strFile,char* strDataset) {

        H5File file(strFile,H5F_ACC_RDONLY);

        DataSet dataset = file.openDataSet(strDataset);
        DataSpace filespace = dataset.getSpace();
        int rank = filespace.getSimpleExtentNdims();

        hsize_t* dims = new hsize_t[rank];
        rank = filespace.getSimpleExtentDims(dims);

        DataSpace mspace(rank,dims);
        int fullsize = 1;
        for(int irank=0;irank<rank;irank++) {
                fullsize *= dims[irank];
        }
        delete[] dims;

        int* data_out = new int[fullsize];
        dataset.read(data_out,PredType::NATIVE_INT,mspace,filespace);

        return data_out;
}
int ReadDatasize(char* strFile,char* strDataset) {

        H5File file(strFile,H5F_ACC_RDONLY);

        DataSet dataset = file.openDataSet(strDataset);
        DataSpace filespace = dataset.getSpace();
        int rank = filespace.getSimpleExtentNdims();

        hsize_t* dims = new hsize_t[rank];
        rank = filespace.getSimpleExtentDims(dims);

        DataSpace mspace(rank,dims);
        int fullsize = 1;
        for(int irank=0;irank<rank;irank++) {
                fullsize *= dims[irank];
        }
        delete[] dims;

        return fullsize;
}
