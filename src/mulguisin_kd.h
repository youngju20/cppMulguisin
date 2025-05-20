#ifndef _MGS3DClusterFinder_
#define _MGS3DClusterFinder_

#include <iostream>
using namespace std;

typedef struct GLink{
    int idbra;
    int idket;
    int imgs;
} GLink;

class MGS3DClusterFinder  
{
private:
  int *Isort ;			// id of a point from the sorted array 
  int *Imgs;
  int *SortedMgsId;
  int *IdOfClsFromSortedGrp;

protected:
  int	Mode;			// Mulguisin mode
  float	Resolution;		// Mulguisin resolution for grouping points
  float	MergingDistance;	// Distance parameter for merging groups
  int	Verbose;		// level of comment print-out
  int	NbOfPnt;		// Number of input points
  int	NbOfGrp;		// number of groups
  int	NbOfCls;		// number of clusters

public:
  float *Gm;			// galaxy mass
  float *Gx;			// galaxy position x
  float *Gy;			// galaxy position y
  float *Gz;			// galaxy position z

  int   *NbOfPntInGrp;		// Number of points in a group
  std::vector<int> *IdOfPntInGrp;	// id of a point in a group
  std::vector<int> *IdOfOrgPntInGrp;	// id of the mother point 
  //std::vector<GLink> *Cll;

  int   *NbOfGrpInCls;		// Number of groups in a cluster
  std::vector<int> *IdOfGrpInCls;	// id of a group in a cluster

  int *MGSlabel;			// MGS label ordered as same as data order
  std::vector<GLink> *MGSLink;

  //
  //
  //

  MGS3DClusterFinder();
  virtual ~MGS3DClusterFinder();

  int	GetMode()           const { return Mode;}
  float	GetResolution()     const { return Resolution;}

  int	GetNbOfPnt()           const { return NbOfPnt;}

  int	GetNbOfGrp()           const { return NbOfGrp;}
  int	GetNbOfPntInGrp(int igrp) const {return NbOfPntInGrp[igrp];}
  int	GetIdOfPntInGrp(int igrp,int ipnt) const {return IdOfPntInGrp[igrp].at(ipnt);}

  int	GetNbOfCls()		const { return NbOfCls;}
  int	GetNbOfGrpInCls(int icls) const {return NbOfGrpInCls[icls];}
  int	GetIdOfGrpInCls(int icls,int igrp) const {return IdOfGrpInCls[icls].at(igrp);}

  int	GetIdOfGrpFromSortedGrp(int idx) const { return SortedMgsId[idx];}
  int   GetIdOfClsFromSortedGrp(int isg) const { return IdOfClsFromSortedGrp[isg];}

  float    GetTotalLength(int imgs);
  int      GetNode(int imgs,int igal);
  int      GetGeneration(int imgs,int ig);
  float    GetLength(int imgs,int ig);
  float    GetOpeningAngle(int imgs,int ig);
  float    GetLinkPolarAngle(int imgs,int ig);
  int      GetNChildren(int imgs,int ig);
  int      GetNBranches(int imgs);

  int      GetDegree(int imgs, int ig);

  // MGS main functions

  bool     Init(int n);
  void     Run();
  void     SortGrp();
  void     MergeGrp();

  void     Terminate();
  void     Qsort(float* numbers,int* partID,int left,int right);

  void     SetNbOfPnt(int n){NbOfPnt=n;}
  void     SetMode(int m) { Mode = m;}
  void     SetResolution(float r) { Resolution = r;}
  void     SetMergingDistance(float r) {MergingDistance = r;}
  void     SetVerbose(int m) { Verbose = m;}

  void     GetLabel();
  void     GetLink();
};

//void Qsort(float* numbers,int* partID,int left,int right);
float Dist2(float x1,float y1,float z1,float x2,float y2,float z2);
float Angle321(float x3,float y3,float z3,float x2,float y2,float z2,float x1,float y1,float z1);
#endif 
