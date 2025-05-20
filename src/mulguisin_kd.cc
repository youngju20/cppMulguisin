#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <math.h>
#include "mulguisin_kd.h"

#include <array>
#include <cmath>
#include <queue>
///////////////////////////////////////////////////////////////////////////
//
// MGS3DClusterFinder
//
// 3D Mulguisin clustering algorithm
//
// Authors: Inkyu Park
//
//      3D Mulguisin algorithm.
//
//      New features of the algorithm are
//
//      1) when Mode = 0
//              Mulguisin behaves like the convetional cone-based algorithm.
//
//      2) when Mode = 1
//              Mulguisin algorithm finds jets using a fixed cone size
//              of *Resol*. But, Jet direction is reestimated everytime
//              a cell is attached. So this mode should give better jet
//              resolution than other cone-based algorithm because
//
//      3) when Mode = 2
//              Mulguisin performs cluster finding by searching nearest cell
//              Now *Resol* means the minimum distance between the cell and
//              any cell inside the cluster
//
// 
// 20210331: merging algorithm added
//
//
///////////////////////////////////////////////////////////////////////////

//________________________________________________________________________
/**
 * Class for representing a point. coordinate_type must be a numeric type.
 */
template<typename coordinate_type, size_t dimensions>
class point {
public:
    point(std::array<coordinate_type, dimensions> c, size_t id) : coords_(c), id_(id) {}
    point(std::initializer_list<coordinate_type> list, size_t id) : id_(id) {
        size_t n = std::min(dimensions, list.size());
        std::copy_n(list.begin(), n, coords_.begin());
    }

    coordinate_type get(size_t index) const {
        return coords_[index];
    }

    double distance(const point& pt) const {
        double dist = 0;
        for (size_t i = 0; i < dimensions; ++i) {
            double d = get(i) - pt.get(i);
            dist += d * d;
        }
        return dist;
    }

    size_t getID() const {
       return id_;
    }

private:
    std::array<coordinate_type, dimensions> coords_;
    size_t id_;
};

/**
 * C++ k-d tree implementation, based on the C version at rosettacode.org.
 */
template<typename coordinate_type, size_t dimensions>
class kdtree {
public:
    typedef point<coordinate_type, dimensions> point_type;
private:
    // test
    //size_t nextID_ = 0;
    struct node {
        node(const point_type& pt) : point_(pt), left_(nullptr), right_(nullptr) {}
        coordinate_type get(size_t index) const {
            return point_.get(index);
        }
        double distance(const point_type& pt) const {
            return point_.distance(pt);
        }
        point_type point_;
        node* left_;
        node* right_;
    };
    node* root_ = nullptr;
    node* best_ = nullptr;
    double best_dist_ = 0;
    size_t visited_ = 0;
    std::vector<node> nodes_;

    struct node_cmp {
        node_cmp(size_t index) : index_(index) {}
        bool operator()(const node& n1, const node& n2) const {
            return n1.point_.get(index_) < n2.point_.get(index_);
        }
        size_t index_;
    };

    node* make_tree(size_t begin, size_t end, size_t index) {
        if (end <= begin)
            return nullptr;
        size_t n = begin + (end - begin)/2;
        auto i = nodes_.begin();
        std::nth_element(i + begin, i + n, i + end, node_cmp(index));
        index = (index + 1) % dimensions;

        // test
        //nodes_[n].point_.setID(nextID_);
        //++nextID_;

        nodes_[n].left_ = make_tree(begin, n, index);
        nodes_[n].right_ = make_tree(n + 1, end, index);
        return &nodes_[n];
    }

    void nearest(node* root, const point_type& point, size_t index) {
        if (root == nullptr)
            return;
        ++visited_;
        double d = root->distance(point);
        if (best_ == nullptr || d < best_dist_) {
            best_dist_ = d;
            best_ = root;
        }
        if (best_dist_ == 0)
            return;
        double dx = root->get(index) - point.get(index);
        index = (index + 1) % dimensions;
        nearest(dx > 0 ? root->left_ : root->right_, point, index);
        if (dx * dx >= best_dist_)
            return;
        nearest(dx > 0 ? root->right_ : root->left_, point, index);
    }

    // test
    std::priority_queue<std::pair<float, const node*>> nearest_neighbors_;

public:
    kdtree(const kdtree&) = delete;
    kdtree& operator=(const kdtree&) = delete;
    /**
     * Constructor taking a pair of iterators. Adds each
     * point in the range [begin, end) to the tree.
     *
     * @param begin start of range
     * @param end end of range
     */
    template<typename iterator>
    kdtree(iterator begin, iterator end) : nodes_(begin, end) {
        root_ = make_tree(0, nodes_.size(), 0);
    }

    /**
     * Constructor taking a function object that generates
     * points. The function object will be called n times
     * to populate the tree.
     *
     * @param f function that returns a point
     * @param n number of points to add
     */
    template<typename func>
    kdtree(func&& f, size_t n) {
        nodes_.reserve(n);
        for (size_t i = 0; i < n; ++i)
            nodes_.push_back(f());
        root_ = make_tree(0, nodes_.size(), 0);
    }

    /**
     * Returns true if the tree is empty, false otherwise.
     */
    bool empty() const { return nodes_.empty(); }

    /**
     * Returns the number of nodes visited by the last call
     * to nearest().
     */
    size_t visited() const { return visited_; }

    /**
     * Returns the distance between the input point and return value
     * from the last call to nearest().
     */
    double distance() const { return std::sqrt(best_dist_); }

    /**
     * Finds the nearest point in the tree to the given point.
     * It is not valid to call this function if the tree is empty.
     *
     * @param pt a point
     * @return the nearest point in the tree to the given point
     */
    const point_type& nearest(const point_type& pt) {
        if (root_ == nullptr)
            throw std::logic_error("tree is empty");
        best_ = nullptr;
        visited_ = 0;
        best_dist_ = 0;
        nearest(root_, pt, 0);
        return best_->point_;
    }
// test
public:
    std::vector<point_type> k_nearest(const point_type& point, size_t k) {
        best_ = nullptr;
        best_dist_ = 0;
        visited_ = 0;
        nearest_neighbors_ = std::priority_queue<std::pair<float, const node*>>();

        k_nearest(root_, point, 0, k);

        std::vector<point_type> result;
        while (!nearest_neighbors_.empty()) {
            result.push_back(nearest_neighbors_.top().second->point_);
            nearest_neighbors_.pop();
        }

        std::reverse(result.begin(), result.end()); // Reverse to get the closest points first
        return result;
    }
// test
private:
    void k_nearest(node* root, const point_type& point, size_t index, size_t k) {
        if (root == nullptr)
            return;

        float d = root->distance(point);
        if (best_ == nullptr || d < best_dist_) {
            best_dist_ = d;
            best_ = root;
        }

        nearest_neighbors_.push(std::make_pair(d, root));
        if (nearest_neighbors_.size() > k) {
            nearest_neighbors_.pop();
            best_dist_ = nearest_neighbors_.top().first;
        }

        if (best_dist_ == 0)
            return;

        float dx = root->get(index) - point.get(index);
        index = (index + 1) % dimensions;
        k_nearest(dx > 0 ? root->left_ : root->right_, point, index, k);
        if (dx * dx >= best_dist_)
            return;
        k_nearest(dx > 0 ? root->right_ : root->left_, point, index, k);
    }
private:
    void find_in_sphere(node* root, const point_type& point, size_t index, float radius) {
        if (root == nullptr)
            return;

        double d = root->distance(point);
        if (d <= radius) {
			points_in_sphere_.push_back(root->point_);
        }

        //if (best_dist_ == 0)
        //    return;

        double dx = root->get(index) - point.get(index);
        index = (index + 1) % dimensions;
        find_in_sphere(dx > 0 ? root->left_ : root->right_, point, index, radius);
        //if (dx * dx >= best_dist_)
        //    return;
        // Early exit
        if (dx * dx >= radius)
            return;
        find_in_sphere(dx > 0 ? root->right_ : root->left_, point, index, radius);
    }

public:
    std::vector<point_type> points_in_sphere(const point_type& point, float radius) {
        //best_ = nullptr;
        //best_dist_ = 0;
        //visited_ = 0;
		points_in_sphere_.clear();

        find_in_sphere(root_, point, 0, radius);

        std::sort(points_in_sphere_.begin(), points_in_sphere_.end(),
           [&point](const point_type& a, const point_type& b) {
                return a.distance(point) < b.distance(point);
           });

        return points_in_sphere_;
    }
private:
	std::vector<point_type> points_in_sphere_;
};

//________________________________________________________________________

//________________________________________________________________________
MGS3DClusterFinder::MGS3DClusterFinder()
{
   // Default constructor. 
}


//________________________________________________________________________
MGS3DClusterFinder::~MGS3DClusterFinder()
{
   // Destructor.
}

//________________________________________________________________________
bool MGS3DClusterFinder::Init(int n)
{
  //
  // prepare input arrays
  //
  NbOfPnt = n;
  if (NbOfPnt==0) return -99;
  Isort = new int[NbOfPnt];
  Gm = new float[NbOfPnt];  // galaxy mass
  Gx = new float[NbOfPnt];  // galaxy position x
  Gy = new float[NbOfPnt];  // galaxy position y
  Gz = new float[NbOfPnt];  // galaxy position z

  NbOfPntInGrp = new int[NbOfPnt];
  IdOfPntInGrp = new vector<int>[NbOfPnt];   // index of the galaxy ig
  IdOfOrgPntInGrp = new vector<int>[NbOfPnt];   // index of the mother galaxy
  //Cll = new vector<GLink>[NbOfPnt]; // link
  Imgs = new int[NbOfPnt];
  for (int i=0;i<NbOfPnt;i++)
  {
    Isort[i] = i;
    //Imgs[i]=0;
    Imgs[i]=-1;
    NbOfPntInGrp[i]=0;
  }
  printf("Mulguisin input arrays are ready (%d cells)\n",NbOfPnt);
  return 1;
}

//________________________________________________________________________
void MGS3DClusterFinder::Run()
{
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (1) Sort Cell array
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  vector <pair<float,int> > vi;
  for (int i=0;i<NbOfPnt;i++)
  {
    vi.push_back(make_pair(Gm[i],i));
  }
  sort(vi.rbegin(),vi.rend());
  for (int i=0;i<NbOfPnt;i++)
  {
    Isort[i]=vi[i].second;
  }

  // Check
  //for (size_t i=0; i<21; i++){
  //  std::cout << "id : " << Isort[i] << std::endl;
  //}

  // --------------------------------------------------------------

  std::vector<point<float,3>> arr_point;
  for (int i=0; i<NbOfPnt; i++)
  {
     int id_sorted = Isort[i]; 
     arr_point.push_back(point<float,3>({Gx[id_sorted], Gy[id_sorted], Gz[id_sorted]}, id_sorted));
  }
  
  typedef kdtree<float,3> tree3d;
  tree3d tree(std::begin(arr_point), std::end(arr_point));

  // --------------------------------------------------------------
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (1) Make 1 dim. sorted array from input common /JF_CELL/
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  // test
  typedef point<float,3> point3d;
  const size_t dim = 3;
  //size_t k = 2;

  NbOfGrp=0;
  for (int i=0; i<NbOfPnt; i++)
  {
    //if (i%10000==0) {cout << i << " th galaxies is under making MGS" << endl;}
    if (i%100000==0) {cout << i << " th galaxies is under making MGS" << endl;}

    // this cell will be tested
    int ic=Isort[i];

    // tag which indicates whether target belongs to existed MGS (0) or not (1)
    int findMGS = 1;
    // First target bocomes MGS group unconditionally
    if (i == 0)
    {
      Imgs[ic] = 0;
      // save target id in 'list for id which belongs to MGS group'
      IdOfPntInGrp[NbOfGrp].push_back(ic);
      // save target id in 'list for mother id'
      IdOfOrgPntInGrp[NbOfGrp].push_back(ic);
      // increase the number of member for MGS group
      NbOfPntInGrp[NbOfGrp]++;
	  NbOfGrp++;
    }
    else
    {
      // Make target as point struct. last value,'-1' is temporal value
      point<float,dim> queryP({Gx[ic],Gy[ic],Gz[ic]}, -1);

      // In the kdtree class, distance = dist*dist
      float Rcut=Resolution*Resolution;

      std::vector<point<float, dim>> pointsInS = tree.points_in_sphere(queryP, Rcut);
      
      // Search nearest point around target
      for (auto itP = std::next(pointsInS.begin()); itP != pointsInS.end(); itP++){
         size_t Pid = itP->getID();
         // check nearest point already belongs to other MGS group
         if (Imgs[Pid] != -1 ){
            Imgs[ic] = Imgs[Pid];
            // save target id in 'list for id which belongs to MGS group'
            IdOfPntInGrp[Imgs[Pid]].push_back(ic);
            // save nearest point id in 'list for mother id'
            IdOfOrgPntInGrp[Imgs[Pid]].push_back(Pid);
            // increase the number of member for MGS group
            NbOfPntInGrp[Imgs[Pid]]++;
            findMGS = 0;
            break;
         }
      }

      // If it cannot find existed MGS group, target becomes new MGS group
      if (findMGS == 1){
         Imgs[ic]=NbOfGrp;
         // save target id in 'list for id which belongs to MGS group'
         IdOfPntInGrp[NbOfGrp].push_back(ic);
         // save target id in 'list for mother id'
         IdOfOrgPntInGrp[NbOfGrp].push_back(ic);
         // increase the number of member for MGS group
         NbOfPntInGrp[NbOfGrp]++;
         NbOfGrp++;
      }

    }

  }
}
//________________________________________________________________________
void MGS3DClusterFinder::SortGrp()
{
  vector <pair<float,int> > vmgs;

  for (int imgs=0;imgs<NbOfGrp;imgs++)
      vmgs.push_back(make_pair(NbOfPntInGrp[imgs],imgs));

  sort(vmgs.rbegin(),vmgs.rend());

  for (int ism=0; ism<NbOfGrp;ism++)
  { 
    int imgs=vmgs[ism].second;
    SortedMgsId[ism]=imgs;
  }
}
//________________________________________________________________________
void MGS3DClusterFinder::MergeGrp()
{
  // create an array of merged MGS (Super cluster) 
  IdOfClsFromSortedGrp=new int[NbOfGrp];

  NbOfGrpInCls=new int[NbOfGrp];
  for (int i=0;i<NbOfGrp;i++) NbOfGrpInCls[i]=0;

  IdOfGrpInCls = new vector<int>[NbOfGrp];

  //

  NbOfCls=0;
  double x1,y1,z1,x2,y2,z2;
  double Dcut2=MergingDistance*MergingDistance;

  for (int ism1=0;ism1<NbOfGrp;ism1++)
  {
    int imgs1=SortedMgsId[ism1];
    int ngals1=NbOfPntInGrp[imgs1];
    int Iwc=-1;
    bool merged=false;
    for (int isc=0;isc<NbOfCls;isc++)
    {
      for (int icl=0;icl<NbOfGrpInCls[isc];icl++)
      {
        int ism2=IdOfGrpInCls[isc].at(icl);

        int imgs2=SortedMgsId[ism2];
        int ngals2=NbOfPntInGrp[imgs2];

        for (int ig1=0;ig1<ngals1;ig1++)
        {
          int idgal1=IdOfPntInGrp[imgs1].at(ig1);
          x1=Gx[idgal1];
          y1=Gy[idgal1];
          z1=Gz[idgal1];
          for (int ig2=0;ig2<ngals2;ig2++)
          {
            int idgal2=IdOfPntInGrp[imgs2].at(ig2);
            x2=Gx[idgal2];
            y2=Gy[idgal2];
            z2=Gz[idgal2];
            if (Dist2(x1,y1,z1,x2,y2,z2)<Dcut2) { merged=true; break; }
          }
          if (merged) break;
        }
        if (merged)
        {
          if (Verbose) printf("groups %d and %d are to be merged\n",ism1,ism2);
          break;
        }
      }
      if (merged)
      {
        Iwc=isc;
        break;
      }
    }
    //
    //
    //
    if (merged)
    { // this cluster is to be attached to the nearest cluster
      IdOfClsFromSortedGrp[ism1]=Iwc;
      IdOfGrpInCls[Iwc].push_back(ism1);
      NbOfGrpInCls[Iwc]++;
    }
    else
    { // this mgs becomes a new super cluster
 //     Isc[ism1]=NbOfCls;
      IdOfClsFromSortedGrp[ism1]=NbOfCls;
      IdOfGrpInCls[NbOfCls].push_back(ism1);
      NbOfGrpInCls[NbOfCls]++;
      NbOfCls++;
    }
  }
}

//________________________________________________________________________
float MGS3DClusterFinder::GetTotalLength(int imgs)
{ // this function returns the total length of whole links
  float lsum=0;
  for (int ig=0;ig<NbOfPntInGrp[imgs];ig++)
  { 
    int idbra=IdOfPntInGrp[imgs].at(ig);
    int idket=IdOfOrgPntInGrp[imgs].at(ig);
    float d2=Dist2(Gx[idbra],Gy[idbra],Gz[idbra],Gx[idket],Gy[idket],Gz[idket]);
    lsum+=sqrt(d2);
  }
  return lsum;
}
//________________________________________________________________________
int MGS3DClusterFinder::GetNode(int imgs,int idgal)
{ // this function returns the node galaxy id
  for (int ig=0;ig<NbOfPntInGrp[imgs];ig++)
  { 
    if (IdOfPntInGrp[imgs].at(ig)==idgal) return IdOfOrgPntInGrp[imgs].at(ig);
  }
  return -1;
}

//________________________________________________________________________
int MGS3DClusterFinder::GetGeneration(int imgs,int ig)
{ // this function returns the generation number of the galaxy
  int ngen = 0;
  int ileaf=IdOfPntInGrp[imgs].at(ig);
  int iseed=IdOfPntInGrp[imgs].at(0);
  while (ileaf != iseed)
  {
    int inode=GetNode(imgs,ileaf);
    ngen++;
    ileaf=inode;
  }
  return ngen;
}

//________________________________________________________________________
float MGS3DClusterFinder::GetLength(int imgs,int ig)
{ // this function returns the total length of whole links
  float lsum = 0;
  int ileaf=IdOfPntInGrp[imgs].at(ig);
  int iseed=IdOfPntInGrp[imgs].at(0);
  while (ileaf != iseed)
  {
    int inode=GetNode(imgs,ileaf);
    float d2=Dist2(Gx[ileaf],Gy[ileaf],Gz[ileaf],Gx[inode],Gy[inode],Gz[inode]);
    lsum+=sqrt(d2);
    ileaf=inode;
  }
  return lsum;
}

//________________________________________________________________________
float MGS3DClusterFinder::GetOpeningAngle(int imgs,int ig)
{ // this function returns the opening angle 
  int ileaf=IdOfPntInGrp[imgs].at(ig);
  int iseed=IdOfPntInGrp[imgs].at(0);
  //if (ileaf == iseed) return 0;
  if (ileaf == iseed) return -1000;
  int inode=GetNode(imgs,ileaf);
  //if (inode == iseed) return 0;
  if (inode == iseed) return -1000;
  int iaxle=GetNode(imgs,inode);
  return Angle321(Gx[ileaf],Gy[ileaf],Gz[ileaf],Gx[inode],Gy[inode],Gz[inode],Gx[iaxle],Gy[iaxle],Gz[iaxle]);
}

//________________________________________________________________________
float MGS3DClusterFinder::GetLinkPolarAngle(int imgs,int ig)
{ // this function returns the polar angle of a link 
  int i1=IdOfPntInGrp[imgs].at(ig);
  int i0=IdOfOrgPntInGrp[imgs].at(ig);
  if (i1 == i0) return -1000;
  return Angle321(Gx[i1],Gy[i1],Gz[i1],Gx[i0],Gy[i0],Gz[i0],0,0,0);
}

//________________________________________________________________________
int MGS3DClusterFinder::GetNChildren(int imgs,int ig)
{ // this function returns the number of children 
  int nchildren=0;
  for (int i=ig+1;i<NbOfPntInGrp[imgs];i++)
  { 
    if (IdOfOrgPntInGrp[imgs].at(i)==IdOfPntInGrp[imgs].at(ig)) nchildren++;
  }
  return nchildren;
}

//________________________________________________________________________
int MGS3DClusterFinder::GetNBranches(int imgs)
{ // this function returns the total number of branches
  int nbranches=0;
  for (int ig=0;ig<NbOfPntInGrp[imgs];ig++)
  { 
    if (GetNChildren(imgs,ig)==0) nbranches++;
  }
  return nbranches;
}


//________________________________________________________________________
void MGS3DClusterFinder::Terminate()
{

  delete [] Gm;
  delete [] Gx;
  delete [] Gy;
  delete [] Gz;

}

//________________________________________________________________________
// Aditional function
void MGS3DClusterFinder::GetLabel(){
	MGSlabel = new int[NbOfPnt];
	//for (auto &ele: MGSlabel){
	//	ele = -1;
	//}
	for (int ele=0; ele<NbOfPnt; ele++){
		MGSlabel[ele] = -1;
	}
	for (int i=0; i<NbOfGrp; i++){
		for (int j=0; j<NbOfPntInGrp[i]; j++){
			//MGSlabel[IdOfPntInGrp[i][j]] = i;
			//cout << IdOfPntInGrp[i][j] << endl;
			MGSlabel[IdOfPntInGrp[i].at(j)] = i;
		}
	}
}

//________________________________________________________________________
int MGS3DClusterFinder::GetDegree(int imgs, int ig)
{ // this function returns the Degree of each node
  int degree=0;
  int iseed=IdOfPntInGrp[imgs].at(0);
  for (int i=ig+1;i<NbOfPntInGrp[imgs];i++)
  {
    if (IdOfOrgPntInGrp[imgs].at(i)==IdOfPntInGrp[imgs].at(ig)) degree++;
  }
  if (IdOfPntInGrp[imgs].at(ig)!=iseed) degree++;
  return degree;
}

//________________________________________________________________________
void MGS3DClusterFinder::GetLink(){
  MGSLink = new vector<GLink>[NbOfGrp];
  for (int imgs=0; imgs<NbOfGrp; imgs++)
  {
    for (int ig=0; ig<NbOfPntInGrp[imgs]; ig++)
    {
      int idbra=IdOfOrgPntInGrp[imgs].at(ig);
      int idket=IdOfPntInGrp[imgs].at(ig);
      MGSLink->push_back({idbra,idket,imgs});
    }
  }
}
			

//________________________________________________________________________
//
//
// External utilities
//
//

void MGS3DClusterFinder::Qsort(float* numbers,int* partID,int left,int right)
//void Qsort(float* numbers,int* partID,int left,int right)
{
  //
  // Qsort algorithm. 
  //
  int pivotInd, pivotID, l_hold, r_hold;
  float pivot;

  l_hold = left;  // save values until the end
  r_hold = right;
  pivot = numbers[left];  // this will not change
  pivotID = partID[left];
  while (left < right) {
     while ((numbers[right] >= pivot) && (left < right))
         right--;
      if (left != right) {
         numbers[left] = numbers[right];
         partID[left] = partID[right];
         left++;
      }
      while ((numbers[left] <= pivot) && (left < right))
         left++;
      if (left != right) {
         numbers[right] = numbers[left];
         partID[right] = partID[left];
         right--;
      }
  }
  numbers[left] = pivot;
  partID[left] = pivotID;
  pivotInd = left;
  left = l_hold;
  right = r_hold;
  if (left < pivotInd) Qsort(numbers, partID, left, pivotInd-1);
  if (right > pivotInd) Qsort(numbers, partID, pivotInd+1, right);
}

float Dist2(float x1,float y1,float z1,float x2,float y2,float z2)
{
   return (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
}

float Angle321(float x3,float y3,float z3,float x2,float y2,float z2,float x1,float y1,float z1)
{
  float dx32=x3-x2;
  float dy32=y3-y2;
  float dz32=z3-z2;
  float dx21=x2-x1;
  float dy21=y2-y1;
  float dz21=z2-z1;
  float sprod=dx32*dx21+dy32*dy21+dz32*dz21;
  float rdx21=sqrt(dx21*dx21+dy21*dy21+dz21*dz21);
  float rdx32=sqrt(dx32*dx32+dy32*dy32+dz32*dz32);
  float cosalpha=sprod/(rdx21*rdx32);
  return acos(cosalpha);
}
