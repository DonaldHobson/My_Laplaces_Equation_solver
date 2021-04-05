using namespace std;
#ifndef P2_HEAT
#define P2_HEAT
#define Polysz 1
#define NearField 2
#define int_s long


template <class T>
class ArrayView;
template <class T>
class Array;
template <class T>
class Link;
template <class T>
class Llist;
class Complex;
class Poly;
class Node;
class Leaf;
class Calc;
class PCalc;
class Transformer;
class SVector;

template <class T>
class ArrayView {
public:
  /* A view of an array object
  Stores its own length, and a pointer to
 a piece of memory. It doesn't own the memory, something else does.
 Make sure not to use this after whatever it was pointing to is deleted.

  */

  int_s len;
  T* array;
  ArrayView(Array<T>& a);
  ArrayView();
  ArrayView(T* p,int_s l );
  T& operator[]( int_s value);
  std::ostream& operator<< (std::ostream &out);//, const Array &self);
};
template <class T>
class Array{
public:
  /*Basic array type. Stores its own length,
   clears up its own memory on deletion, passed by move semantics*/
  int_s len{0};
  T* array{0};
  Array(int_s l);
  Array();
  ~Array();
  Array(Array& a);
  Array<T>& operator=(Array<T>&& a);
  //void fromlist(Llist<T> l);
  Array(Llist<T> list);
  T& operator[]( int_s value);
  std::ostream& operator<< (std::ostream &out);//, const Array &self);
};

template <class T>
class Link {
public:
  /*A link in a singly linked list*/
  T data;
  class Link *next;
  Link(T& d);
  ~Link();
  std::ostream& operator<< (std::ostream &out);
};



template <class T>
class Llist {
public:
  /*the corresponding linked list. Stores pointers to both ends.
  Can be filter_strip into a pair of lists*/
  int_s len;
  class Link<T> *first;
  class Link<T> *last;
  Llist();
  T get_uniq();
  ~Llist();
  Llist<T>(Llist<T>& a);
  Llist<T>& operator=(Llist<T>& a);
  void Push(T data);
  void Push(Link<T>* n);
  Array<T> toArray();
  void append(Llist other);
  void filter_strip(function<bool (T)> filter,Llist<T> &other);
  template<class TT> friend std::ostream& operator<< (std::ostream &out, const Llist<TT> &self);
};

template <class T1, class T2>
class Either{public:
  /*Uses unions to make an either type. Only ever used to store pointers to either Node or Leaf in the quadtree*/
  bool isLeft;union{T1 left;T2 right;    };
  Either<T1, T2> Left(T1 x);
  Either<T1, T2> Right(T2 x);
  bool matchLeft(T1& x);
  bool matchRight(T2& x);
};


class Complex{
public:
  double re{},im{};
  double magsq();
  double mag();
  friend std::ostream& operator<< (std::ostream &out, const Complex &self);
  friend Complex operator+(const Complex &self, const Complex &other);
  friend Complex operator-(const Complex &self, const Complex &other);
  friend Complex operator*(const Complex &self, const Complex &other);
  friend Complex operator/(const Complex &self, const Complex &other);
};

class Poly{
public:
  /*represents polynomial. Not currently in use*/
  array<Complex,Polysz> coef;
  Poly shift(Complex i);
  Poly( );
  Poly(array<Complex,Polysz> c );
  Poly(array<Leaf*,4> leaves);
  Poly(array<Node*,4> nodes);
};


class NodePos{
public:
/*represents the position of a node. no particular reason this isn't part of the node class.*/
  int_s x;
  int_s y;
  int height;
  std::ostream& operator<< (std::ostream &out);
};


class Leaf{
public:
  /* the end leaves of a quadtree. Leaves are only put on the lowest level.
  The x and y positions are odd. Leaves are 2x2 in size. So every midpoint and edgepoint of any node or leaf is an integer.*/
  int_s x;
  int_s y;
  int_s index;
  double val;
  double get_sum(int x,int y)const;
  Leaf(int_s x,int_s y, double val);
  Leaf(const Leaf& l);
  Leaf();
};

//magic_number_gen1.py
const array<Complex,3>const_integrals{Complex{-1.0611754268825244,0.0},Complex{0.00400595283888705,0.0},Complex{0.34552084289931106,0.7853981633974483}};

class Node{
public:
  /*main node object in octree*/
  NodePos pos;
  int_s count;
  int_s index;
  //Poly val;//terms actually log(x),1/x, 1/x^2 of ect.
  Either<array<Node*,4>, array<Leaf*,4>> branches{};
  ~Node();
  Node(Llist<Leaf> luu,NodePos po);
  Node(Node& n);
  void println(int h);
  void print();
  void set_val();
  double get_sum(int x,int y);
  template <class V>
  V recurse_level(V start);
  void get_llist_level(Llist<Node*>& l,int h);
  void get_llist_leaf(Llist<Leaf*>& l);
  void add_pcalc(Llist<PCalc>& l,int_s x,int_s y);
  void add_calc(Llist<Calc>& l,int_s x,int_s y,int_s to);
  void print_debug();
  std::ostream& operator<< (std::ostream &out);
};
class PCalc{
public:
  //not in use
  //represents a calc that might not need doing.
double val;
int_s loc;
PCalc* prev;
bool used;
PCalc(double v);
void r_use();
};
class Calc{
public:
  /*rather than repeatedly recursively traversing a tree, calculations are flattened into a list of calcs.*/
  double val;
  int_s from;
  int_s to;
  //bool advance;
  std::ostream& operator<< (std::ostream &out);

};
class SVector{
public:
  /*an MPI parrallel vector. values are split between different processes. */
  int rank{};
  int size{};
  ArrayView<double> vals{};
  SVector(int r, int s,ArrayView<double> vals);
  SVector(ArrayView<double> vals);
  SVector();
  double sum(SVector other);
  void share_to(ArrayView<double> target);
  double& operator[]( int_s value);
};
class Transformer{
public:
  /*holds list of calcs, and uses them to transform start to my_end*/
  Array <Calc> coeffs{};//
  Transformer(Node& n,Llist<Leaf*>& lf, int mpi_size, int mpi_rank,int_s vlen);
  void run();
  Array<double> val;
  ArrayView<double> start;
  SVector my_start;
  SVector my_end;
};

//Stores the values, and the coefficients of the polynomials, all in one array.
/*class TVector{
public:
  Array<double> coeffs;
  ArrayView<double> start;
  SVector end;
  void set_views(int l);
};*/

#endif
