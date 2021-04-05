using namespace std;
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <mpi.h>
#include "p3.h"

#define FALSE bool(0);
#define TRUE bool(1);
//#define cerr //cerr

/*reasons for segmentation faults −

    Accessing an array out of bounds
    Dereferencing NULL pointers
    Dereferencing freed memory
    Dereferencing uninitialized pointers
    Incorrect use of the "&" (address of) and "*" (dereferencing) operators
    Improper formatting specifiers in printf and scanf statements
    Stack overflow
    Writing to read-only memory
    */

  template <class T>
  ArrayView<T>::ArrayView(){
      array=0;
      len=0;
    }
template <class T>
ArrayView<T>::ArrayView(Array<T> &a){
  array=a.array;
  len=a.len;
}
template <class T>
ArrayView<T>::ArrayView(T* p , int_s l){
  array=p;
  len=l;
}



template <class T>
T& ArrayView<T>::operator[] (int_s index){
  return array[index];
}


template <class T>
std::ostream& operator<< (std::ostream &out, const ArrayView<T> &self){
  out<<"(";
  for (int_s i=0;i<self.len;i++){
    out << self.array[i]<<",";
  }
  out<<")";
  return out;
}




template <class T>
Array<T>::Array(int_s l){
  len=l;
  array=new T[len]{};
}
template <class T>
Array<T>::Array(){
  len=0;
  array=0;
}
template <class T>
Array<T>::Array(Llist<T> l){

  array=new T[l.len]{};
  len=l.len;

  Link<T> *st=l.first;
  for (int_s i=0;i<len;i++){

    array[i]=st->data;


    if (st->next==0){
      st=l.last;
      cerr<<"st="<<(*st)<<" at i="<<i<<"of "<<len<<endl;
      cerr<<"st "<<st<<" last "<<l.last<<endl;
    }
    st=st->next;
  }
  cerr<<"Array(Llist) succeeded"<<endl;
}
template <class T>
Array<T>::~Array(){
  delete [] array;
  array=0;
}

template <class T>
Array<T>::Array(Array& a){
  array = a.array;
  a.array = 0;
  len=a.len;
  a.len=0;
}
template <class T>
Array<T>& Array<T>::operator=(Array<T>&& a)
	{
    cerr<<"op="<<endl;
		// Self-assignment detection
		if (&a == this){
      cerr<<"same"<<endl;
      return *this;
    }

		// Release any resource we're holding
    cerr<<"about to delete"<<endl;
		delete[] array;
    cerr<<"deleted"<<endl;
		// Copy the resource
		array=a.array;
		len=a.len;
    a.array=0;
    a.len=0;
    cerr<<this<<endl;
    cerr<<(*this)<<endl;
		//return *this;
	}
template <class T>
T& Array<T>::operator[] (int_s index){
  return array[index];
}


template <class T>
std::ostream& operator<< (std::ostream &out, const Array<T> &self){
  out<<"<"<<self.len<<">(";
  for (int_s i=0;i<self.len&&i<10;i++){
    //out << self.array[i]<<",";
  }
  out<<")";
  return out;
}




template <class T>
Link<T>::Link(T& d){
  data=d;
  next=0;
}
template <class T>
std::ostream& operator<< (std::ostream &out, const Link<T> &self){
  out<<"Data<"<<self.data<<">["<<self.next<<"]";
  return out;
}
template <class T>
Link<T>::~Link(){

  //delete next;
  //next=0;
}


template <class T>
Llist<T>::Llist(){
  len=0;
  first=0;
  last=0;
}
template <class T>
T Llist<T>::get_uniq(){
  if (len!=1){
    //cerr<<this<<endl;
    cerr<<this<<endl;
    throw runtime_error("get_uniq called on Llist containing the wrong number of variables");
  }
  T ret=first->data;
  delete first;
  len=0;
  first=0;
  last=0;
  return ret;
}
template <class T>
Llist<T>::~Llist(){

  Link<T> * f{first};
  //if (f==0) return;
  Link<T>  * g{first};
  while (f!=0){
    g=f->next;
    delete f;
    f=g;
  }
  //delete first;
  first=0;
}

template <class T>
Llist<T>::Llist(Llist<T>& a){

  first = a.first;
  a.first = 0;
  last = a.last;
  a.last = 0;
  len=a.len;
  a.len=0;
}
template <class T>
Llist<T>& Llist<T>::operator=(Llist<T>& a) // note: not const
{
if (&a == this)
return *this;

delete first;
first = a.first;
a.first = 0;
last = a.last;
a.last = 0;
len=a.len;
a.len=0;
return *this;
}
template <class T>
void Llist<T>::Push(T data){
  Link<T>* n=new Link<T>{data};
  if (last==0){
    first=n;
    last=n;

  }else{
    last->next=n;
    last=n;

  }
  len++;
}
template <class T>
void Llist<T>::Push(Link<T>* n){
  n->next=0;
  if (last==0){
    first=n;
    last=n;
  }else{
    last->next=n;
    last=n;
  }
  len++;
}

template <class T>
void Llist<T>::append(Llist other){
  len+=other.len;
  last->next=other.first;
  last=other.last;
}
template <class T>
void Llist<T>::filter_strip(function<bool (T)> filter,Llist<T> &other){
  Link<T>** nx=&first;
  Link<T>* t;
  Link<T>** ln=&other.first;
  while (*nx!=0){
    if(filter((*nx)->data)){
      t=*nx;
      *nx=t->next;
      other.Push(t);

    }else{
      nx=&((*nx)->next);
    }
  }
  len-=other.len;

}

template <class T>
std::ostream& operator<< (std::ostream &out, const Llist<T> &self){
  Link<T> *n=self.first;
  out<<"[";
  while(n!=0){
    out << (n->data)<<",";
    n=n->next;
  }
  out<<"]";
  return out;
}
template <class T>
std::ostream& operator<< (std::ostream &out, const Llist<T*> &self){
  Link<T*> *n=self.first;
  out<<"[";
  while(n!=0){
    if (n->data!=0){
      out <<"*" <<(*(n->data))<<",";
    }else{
      out<<"Null,";
    }

    n=n->next;
  }
  out<<"]";
  return out;
}

std::ostream& operator<< (std::ostream &out, const tuple<int_s,int_s> &self){
  out<<"Tuple("<<get<0>(self)<<","<<get<1>(self)<<")";
  return out;
}

//coppied Either from https://gist.github.com/jvranish/784615

//class Either{public:
//  bool isLeft;union{T1 left;T2 right;    };
//template<class T1_, class T2_> friend Either<T1_, T2_> Left(T1_ x);
//template<class T1_, class T2_> friend Either<T1_, T2_> Right(T1_ x);
template <class T1, class T2>
bool Either<T1,T2>::matchLeft(T1& x){
  if (isLeft)x = left;
  return isLeft;
}
template <class T1, class T2>
bool Either<T1,T2>::matchRight(T2& x){
  if (!isLeft)x = right;
  return !isLeft;
}

template <class T1, class T2>
Either<T1, T2> Either<T1,T2>::Left(T1 x){
  Either<T1, T2> e;
  e.isLeft = true;
  e.left = x;
  return e;
}
template <class T1, class T2>
Either<T1, T2> Either<T1,T2>::Right(T2 x){
  Either<T1, T2> e;
  e.isLeft = false;
  e.right = x;
  return e;
}


double Complex::magsq(){
  return re*re+im*im;
}
double Complex::mag(){
  return sqrt(magsq());
}

std::ostream& operator<< (std::ostream &out,const Complex &self){//, const Complex &self){
  out<<self.re<<(self.im>=0.?"+":"")<<self.im<<"i";
  return out;
}
Complex operator+ (const Complex &self,const Complex &other){
  return Complex{self.re+other.re,self.im+other.im};
}

Complex operator- (const Complex &self, const Complex &other){
  return Complex{self.re-other.re,self.im-other.im};
}
Complex operator* (const Complex &self, const Complex &other){
  return Complex{self.re*other.re-self.im*other.im,self.re*other.im+self.im*other.re};
}
Complex operator* (const Complex &self, const double &other){
  return Complex{self.re*other,self.im*other};
}
Complex operator* (const double &other, const Complex &self){
  return Complex{self.re*other,self.im*other};
}
Complex operator/ (const Complex &self, const Complex &other){
  double m=other.re*other.re+other.im*other.im;
  return Complex{(self.re*other.re+self.im*other.im)/m,(self.im*other.re-self.re*other.im)/m};
}

Poly Poly::shift(Complex y1){//see magic_number_gen2.py
  array<Complex,Polysz> newcoef{};
  Complex y2{y1*y1};
  newcoef[0]= coef[0];
  newcoef[1]= y1*coef[0]+coef[1];
  newcoef[2]= -0.5*y2*coef[0]-y1*coef[1]+coef[2];

}
Poly::Poly( ){
  //array<Complex,Polysz> coef{};
  //coef[0]=Complex{6,6};
}
Poly::Poly(array<Complex,Polysz> c ){
  coef=c;
}
Poly::Poly(array<Leaf*,4> leaves){
  //array<Complex,Polysz> coef{};
  for (int i =0;i<4;i++){
    if (leaves[i]!=0){
      coef[0]=coef[0]+Complex{leaves[i]->val,0};
      //cerr<<coef[0]<<", ";
    }
  }
  //return Poly(coef);
}
Poly::Poly(array<Node*,4> nodes){
  //array<Complex,Polysz> coef{};
  //coef[0]=Complex{3,7};
  for (int i =0;i<4;i++){
    if (nodes[i]!=0){
      //coef[0]=coef[0]+(nodes[i]->val.coef[0]);
    }
  }
  //cerr<<"PN "<<coef[0];
  //return Poly(coef);
}

Poly operator+ (const Poly &self, const Poly &other){
  array<Complex,Polysz> c{};
  for (int i=0;i<Polysz;i++)c[i]=self.coef[i]+other.coef[i];
  return Poly(c);
}
Poly operator- (const Poly &self, const Poly &other){
  array<Complex,Polysz> c{};
  for (int i=0;i<Polysz;i++)c[i]=self.coef[i]-other.coef[i];
  return Poly(c);
}
Poly operator* (const Poly &self, const Poly &other){
  array<Complex,Polysz> c{};
  for (int i=0;i<Polysz;i++)c[i]=self.coef[i]*other.coef[i];
  return Poly(c);
}
Poly operator/ (const Poly &self, const Poly &other){
  array<Complex,Polysz> c{};
  for (int i=0;i<Polysz;i++)c[i]=self.coef[i]/other.coef[i];
  return Poly(c);
}

std::ostream& operator<< (std::ostream &out,const NodePos &self){//, const Complex &self){
  out<<self.x<<","<<self.y<<","<<self.height;
  return out;
}


double Leaf::get_sum(int xv,int yv)const {

  int gx=abs(x-xv);
  int gy=abs(y-yv);
  int mx=max(gx,gy);
  int mn=min(gx,gy);
  if (mx<NearField){
    int v{};



    return const_integrals[mx*(mx+1)/2+mn].re-const_integrals[0].re;


  }else{
    return 0.5*log(mx*mx+mn*mn)-const_integrals[0].re;
  }
}

std::ostream& operator<< (std::ostream &out, const Leaf &self){
  out<<"Leaf("<<self.x<<","<<self.y<<","<<self.val<<")";
  return out;
}
Leaf::Leaf(int_s xx,int_s yy,double valv){
  x=(xx<<1)+1;
  y=(yy<<1)+1;
  val=valv;
}
Leaf::Leaf(const Leaf& l){
  x=l.x;
  y=l.y;
  val=l.val;
}
Leaf::Leaf(){
  x=0;
  y=0;
  val=0;

}
//const array<double,3>const_integrals{1.,0.5,0.4};
/*class Node{
public:
NodePos pos;
int_s count;
Poly val;//terms actually log(x),1/x, 1/x^2 of ect.
Either<array<Node*,4>, array<Leaf*,4>> branches{};
*/
Node::~Node(){
  if (branches.isLeft){
    for (int i = 0; i < 4; i++) {
      delete branches.left[i];
    }
  }else{
    for (int i = 0; i < 4; i++) {
      delete branches.right[i];
    }
  }
}
Node::Node(Llist<Leaf> luu,NodePos po){
  pos=po;
  count=luu.len;
  int_s x=pos.x<<1;
  int_s y=pos.y<<1;
  int h=pos.height-1;

  int_s filx=((x+1)<<h)-1;
  int_s fily=((y+1)<<h)-1;


  Llist<Leaf> lud{};
  luu.filter_strip([fily](Leaf l)->bool{return l.y>fily;},lud);
  //cerr<<luu<<endl;
  //cerr<<lud<<endl;

  Llist<Leaf> ldd{};
  lud.filter_strip([filx](Leaf l)->bool{return l.x>filx;},ldd);
  Llist<Leaf> ldu{};
  luu.filter_strip([filx](Leaf l)->bool{return l.x>filx;},ldu);
  Llist<Leaf>* all[4]={&luu,&ldu,&lud,&ldd};
  Llist<Leaf>* p;
  if (h==0){
    branches.isLeft=FALSE;
    for (int i=0;i<4;i++){
      p=all[i];
      if (p->len==1){
        branches.right[i]=new Leaf{p->get_uniq()};//new Leaf{NodePos{x|(i&1), y|(i>>1),h},0.};
      }else if(p->len>0){
        cerr<<*p<<endl;
        throw runtime_error("Error Multiple points in same position");
      }
    }
  }else{
    branches.isLeft=TRUE;
    for (int i=0;i<4;i++){
      p=all[i];
      if (p->len>0){
        branches.left[i]=new Node{*p,NodePos{x|(i&1), y|(i>>1),h}};
      }
    }
  }

}
Node::Node(Node& n){//cerr<<"EEEEEEEEEEE";
pos=n.pos;
count=n.count;
branches=n.branches;
n.branches.left={0,0,0,0};

}
void Node::println(int h){


  if (!branches.isLeft){
    h*=2;
    for (int j=h;j<h+2;j++){
      if (branches.right[j]!=0){
        cerr<<"x";//branches.right[j]->pos.x;//
      }else{
        cerr<<".";
      }
    }
  }else{
    int q=(3<<(pos.height-2))-1;
    //cerr<<q<<endl;
    //q=16;
    int f;
    if (h==q){
      for (int i=0;i<q;i++)cerr<<"-";
      cerr<<"+";
      for (int i=0;i<q;i++)cerr<<"-";
      return;
    }else if(h<q){
      f=0;
    }else{
      f=2;
      h-=q+1;
    }



    for (int j=f;j<f+2;j++){
      if (branches.left[j]!=0){
        branches.left[j]->println(h);

      }else{
        for (int i=0;i<q;i++)cerr<<" ";
      }
      if (j==f){cerr<<"|";}
    }


  }
}
void Node::print(){
  for (int i=0;i<(3<<(pos.height-1))-1;i++){
    println(i);cerr<<endl;
  }
}
void Node::set_val(){
  if (branches.isLeft){
    for (int i=0;i<4;i++){
      if (branches.left[i]!=0){
        branches.left[i]->set_val();
      }
    }
    //val=Poly{branches.left};
    //val.coef[0].im=9;
    //cerr<<"SD "<<val.coef[0]<<",";
  }else{
    //val=Poly{branches.right};

  }
  //cerr<<val.coef[0]<<",";
  //cerr<<"A val<<"
}
//template<class T>


template<class T, class V >
V recurse(array<T*,4> pts,V start,function<V (T,V)> func ){
  for (int i=0;i<4;i++){
    if(pts[i]!=0){
      start=func(*pts[i],start);
    }
  }
  return start;
}
template<class T, class V >
void recurse__void(array<T*,4> pts,V start,function<void (V,T)> func ){
  for (int i=0;i<4;i++){
    if(pts[i]!=0){
      func(*pts[i],start);
    }
  }
}
template<class T>
double recurse_sum(array<T*,4> pts,int x,int y){
  double sum=0.;
  for (int i=0;i<4;i++){
    if(pts[i]!=0){
      sum+=pts[i]->get_sum(x,y);
    }
  }
  return sum;
}
template< class V >
V Node::recurse_level(V start){
 return start;
}
void Node::get_llist_level(Llist<Node*>& l,int h){
  if (pos.height==h){
    l.Push(this);
  }else{
    if (!branches.isLeft){
      throw runtime_error("Bad h value in get_llist_level");
    }
    for (int i=0;i<4;i++){
      if(branches.left[i]!=0){
        branches.left[i]->get_llist_level(l,h);
      }
    }
    //recurse_void(this->branches.left,[h,l](Node n)->void{n.get_llist_level(l,h)})
  }
}
void Node::get_llist_leaf(Llist<Leaf*>& l){
  if (branches.isLeft){
    for (int i=0;i<4;i++){
      if(branches.left[i]!=0){
        branches.left[i]->get_llist_leaf(l);
      }
    }
  }else{
    for (int i=0;i<4;i++){
      if(branches.right[i]!=0){
        l.Push(branches.right[i]);
      }
    }
  }
}
void Node::add_calc(Llist<Calc>& l,int_s x,int_s y,int_s to){
  if (branches.isLeft){
    for (int i=0;i<4;i++){
      if(branches.left[i]!=0){
        branches.left[i]->add_calc(l,x,y,to);
      }
    }
  }else{
    for (int i=0;i<4;i++){
      if(branches.right[i]!=0){
        l.Push(Calc{branches.right[i]->get_sum(x,y),branches.right[i]->index,to});
      }
    }
  }
}
/*void Node::add_pcalc(Llist<PCalc>& l,int_s x,int_s y){
  if (branches.isLeft){
    for (int i=0;i<4;i++){
      if(branches.left[i]!=0){
        branches.left[i]->add_pcalc(l,x,y);
      }
    }
  }else{
    for (int i=0;i<4;i++){
      if(branches.right[i]!=0){
        l.Push(PCalc{branches.right[i]->get_sum(x,y),branches.right[i]->index,0});
      }
    }
  }
}*/
double Node::get_sum(int x,int y){
  double sum=0.;
  if (branches.isLeft){
    sum=recurse_sum(branches.left,x,y);
  }else{
    sum=recurse_sum(branches.right,x,y);
  }
  return sum;
}


std::ostream& operator<< (std::ostream &out,const Node& self){
  out<<"Node("<<self.pos<<")<";
  if (self.branches.isLeft){

    for (int i=0;i<4;i++){
      if (self.branches.left[i]!=0){
        out<<*self.branches.left[i];
      }
    }
  }else{
    for (int i=0;i<4;i++){
      if (self.branches.right[i]!=0){
        out<<(*(self.branches.right[i]))<<",";
      }
    }
    out<<">";
  }
  return out;
}

std::ostream& operator<< (std::ostream &out,const Calc& self){
  out<<"Calc{"<<self.val<<","<<self.from<<","<<self.to<<"}";
}



double SVector::sum(SVector other){
  double local_sum{},global_sum{};
  for (int_s i=0;i<vals.len;i++){
    local_sum+=vals[i]*other.vals[i];
  }
  MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return global_sum;
}

  void SVector::share_to(ArrayView<double> target){
    Array<double> copy{vals.len};
    for (int_s i=0;i<vals.len;i++){
      copy[i]=vals[i];
    }
    MPI_Allgather(copy.array,copy.len,MPI_DOUBLE,target.array,target.len,MPI_DOUBLE,MPI_COMM_WORLD);

  }

SVector::SVector(int r, int s,ArrayView<double> val){
  vals=val;
  rank=r;
  size=s;
}
SVector::SVector(ArrayView<double> val){
  vals=val;
}
SVector::SVector(){
}
double& SVector::operator[]( int_s value){
  //cerr<<vals.len<<"££"<<value<<endl;
  return vals[value];
}




















Transformer::Transformer(Node& n,Llist<Leaf*>& lf, int mpi_size, int mpi_rank,int_s vlen){
  //Llist<Leaf*> lf{};
  n.get_llist_leaf(lf);
  int_s total_len=mpi_size*vlen;
  cerr<<"total_len ="<<total_len<<",mpi_rank ="<<mpi_rank<<", vlen="<<vlen<<endl;
  val=Array<double>{total_len*2};
  //b.set_views(lf.len);
  start=ArrayView<double>{val.array,total_len};
  my_start=SVector{ArrayView<double>{&val.array[mpi_rank*vlen],vlen}};
  my_end=SVector{ArrayView<double>{&val.array[total_len],vlen}};
  //cerr<<"lf"<<lf<<endl;
  Llist<Node*> nd{};
  for (int i=1;i<n.pos.height+1;i++){
    n.get_llist_level(nd,i);
  }
  int i=0;
  cerr<<"got node list"<<endl;
  Link<Leaf*>* st{lf.first};
  while (st!=0){
    st->data->index=i;
    val[i]=st->data->val;
    i++;
    st=st->next;
  }
  cerr<<"leaf loop"<<endl;
  Link<Node*>* st2{nd.first};
  while (st2!=0){
    st2->data->index=i;
    i++;
    st2=st2->next;
  }
  cerr<<"Node loop"<<endl;
  Llist<Calc> calc{};
  //for (int_s i=0;i<total_len;i++){
    //calc.Push(PCalc{0,0,0,1});
  //}
  st=lf.first;
  int_s lfc=0;
  int_s my_n=vlen*mpi_size;
  while (st!=0){
    if (vlen*mpi_rank<=lfc&& lfc<vlen*(mpi_rank+1)){
      /*if (st->data==0){
        cerr<<"error data ptr is 0"<<endl;
      }*/
      n.add_calc(calc,st->data->x,st->data->y,my_n);
      my_n++;
      //calc.last->data.advance=TRUE;
    }
    st=st->next;
    lfc++;
  }
  //cerr<<endl<<"CALCS"<<calc<<endl<<endl;
  cerr<<"Calc Llist made"<<endl;//*(calc.first)<<","<<*(calc.last)<<","<<calc.len<<endl;
  coeffs=Array<Calc> {calc};
  cerr<<"Transformer made"<<endl;
  for (int i=0;i<1000000;i++){}
  cerr<<"Transformer made"<<endl;
  //cerr<<coeffs.len<<","<<coeffs<<endl<<endl;
}
void Transformer::run(){
for (int_s i=start.len;i<val.len;i++){
  val[i]=0;
}
//int j=start.len;
//cerr<<"transforming ..."<<coeffs.len<<endl;
for(int i=0;i<coeffs.len;i++){

  val[coeffs[i].to]+=coeffs[i].val*val[coeffs[i].from];
  //cerr<<"TTT<"<<j<<", "<<coeffs[i].from<<">";
  //j+=coeffs[i].advance;
}

}

/*void TVector::set_views(int l,int mpi_size, int mpi_rank){
  start.array=coeffs.array;
  start.len=l;
  end=ArrayView<double>{&coeffs.array[coeffs.len-l]};

}*/
bool test(Llist<int> l){
  cerr<<"test len"<<l.len<<endl;
  return 0;
}
void rapply(const Leaf& l,Array<Array<Poly>>& polys, Array<double>& base,int_s x,int_s y,int d, int h){
int_s ofs=1<<(h);
int_s gx=x+ofs-l.x;
int_s gy=y+ofs-l.y;
int_s agx=abs(gx);
int_s agy=abs(gy);
int_s dst=agx*agx+agy*agy;
//if ()
if (h>0){
  rapply(l,polys, base,x,y,d+1,h-1);
  rapply(l,polys, base,x+ofs,y,d+1,h-1);
  rapply(l,polys, base,x,y+ofs,d+1,h-1);
  rapply(l,polys, base,x+ofs,y+ofs,d+1,h-1);
}else{
  base[(x>>1)*(1<<d)+(y>>1)]+=l.get_sum(x,y)*l.val;

}
}
int main(int argc, char** argv){

  int mpi_rank, mpi_size;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);//cerr<<"kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk";

  mpi_size=1;
  cerr<<"rank"<<mpi_rank<<endl;
  for (int i=0;i<argc;i++){


    cerr<<"argv"<<i<<" :"<<argv[i]<<endl;
  }
  int logs=stoi(argv[1]);
  int s=1<<logs;
  int intsz=4;
  int doublesz=8;
  double rval{};
  int xl{},yl{};
  Llist<Leaf> ll{};
  //Array<char> temp{intsz*2+doublesz};

  int_s lla_len{};
  if (mpi_rank==0){
    if (*argv[2]=='y'){
      for (int_s i=0;i<s;i++){
        ll.Push(Leaf{0,i,1});
        ll.Push(Leaf{s-1,i,-1});
      }
    }else{


      for (int i=0;i<s*s&&!isnan(rval);i++){
        fread(&xl,sizeof(int),1,stdin);
        fread(&yl,sizeof(int),1,stdin);
        fread(&rval,sizeof(double),1,stdin);
        //rval=*(double*)(&temp[intsz*2]);
        if (!isnan(rval)){
          ll.Push(Leaf{xl,yl,rval});
        }
        //cerr<<xl<<","<<yl<<","<<rval<<endl;


      }
    }
    lla_len=ll.len;
  }
  MPI_Bcast(&lla_len,1,MPI_LONG_INT,0,MPI_COMM_WORLD);
  cerr<<"Parsed input of "<<lla_len<<" points"<<endl;
Array<Leaf> lla{lla_len};
if (mpi_rank==0){
  lla=Array<Leaf>(ll);
}
  MPI_Bcast((void*)lla.array,
    lla_len*sizeof(Leaf),
    MPI_CHAR,
    0,
    MPI_COMM_WORLD);

  for (int_s i=0;i<lla.len;i++){
    ll.Push(lla[i]);
  }
cerr<<"Leaves shared "<<lla[0]<<","<<lla[1]<<" ... "<<endl;

//cerr<<mpi_rank<<endl;



  /*if(mpi_rank>0){
    cerr<<"qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq";
    MPI_Finalize();
    return 0;
  }*/


  //for (int i=0;i<s;i++){
  //  ll.Push(Leaf{i,15-i,-1.});
  //}

  //ll.Push(Leaf{3,4,-5.});

  int_s tlen=ll.len;
  int_s vlen= (tlen+mpi_size-1)/mpi_size;

  //cerr<<"here-1 ";
  Node n{ll,NodePos{0,0,logs+1}};
  cerr<<"Node created "<<endl;
  //cerr<<"here0 ";
  //Llist<Node*> q{};
  //cerr<<"here1 ";
  //n.get_llist_level(q,2);
  //cerr<<"here2 ";
  //cerr<<n<<endl;
  //cerr<<q<<endl;
  //cerr<<"here3 ";
  //TVector p{};//Array<double>{}};
  Llist<Leaf*> lf{};

  Transformer tr{n,lf,mpi_size,mpi_rank,vlen};
  cerr<<"Transformer created "<<endl;
  //cerr<<"wwwwwwwwwww"<<tr.coeffs<<endl;
  //int_s vlen=p.start.len;
  Array<double> r_{vlen};
  Array<double> x_{vlen};

  double rsold{},rsnew{},rsrat{};
  SVector x{ArrayView<double>{x_}};
  SVector r{ArrayView<double>{r_}};
  //cerr<<"xxxxxxx"<<endl;
  //SVector pstart{ArrayView<double>{&p.start[mpi_rank*vlen],vlen}}
  for (int_s i=0;i<vlen;i++){
    r_[i]=tr.my_start[i];
    //rsold+=p.start[i]*p.start[i];
  }
  //cerr<<"yyyyyyy"<<endl;
  rsold=tr.my_start.sum(tr.my_start);
  //cerr<<"r"<<r.vals<<endl;
  double dotp{},alpha{};
  int loopcount;
  if (mpi_rank==0){
    loopcount=stoi(argv[3]);
  }
  MPI_Bcast(&loopcount,1,MPI_INT,0,MPI_COMM_WORLD);

  for (int loop=0;loop<loopcount;loop++){
    //cerr<<"x"<<loop<<" ,"<<x_<<endl;
    //cerr<<"r"<<loop<<" ,"<<r<<endl;
    //cerr<<"p.start"<<loop<<" ,"<<p.start<<endl;
    //cerr<<"p.end"<<loop<<" ,"<<p.end<<endl;

    //cerr<<rsold<<","<<rsnew<<","<<alpha<<","<<rsrat<<","<<dotp<<endl;
    cerr<<"start_run"<<endl;fflush(stderr);
    tr.run();
    cerr<<"end_run"<<endl;fflush(stderr);
    dotp=0;
    for (int_s i=0;i<vlen;i++){
      //dotp+=p.start[i]*p.end[i];
    }
    //cerr<<"st"<<tr.my_start.vals<<endl;
    //cerr<<"end"<<tr.my_end.vals<<endl;
    dotp=tr.my_start.sum(tr.my_end);

    alpha=rsold/dotp;

    //cerr<<"a"<<loop<<" ,"<<alpha<<"="<<rsold<<"/"<<dotp<<endl;
    rsnew=0;
    for (int_s i=0;i<vlen;i++){
      x[i]+=alpha*tr.my_start[i];
      r[i]-=alpha*tr.my_end[i];
      //rsnew+=r[i]*r[i];
    }
    rsnew=r.sum(r);

    rsrat=rsnew/rsold;
    cerr<<"rsnew "<<rsnew<<","<<mpi_rank<<endl;fflush(stderr);
    if (rsnew<1e-5){
      break;
    }
    //cerr<<"rsrat"<<loop<<" ,"<<rsrat<<"="<<rsnew<<"/"<<rsold<<endl;
    for (int_s i=0;i<vlen;i++){
      tr.my_start[i]=r[i]+rsrat*tr.my_start[i];
    }
    cerr<<"start_share"<<endl;fflush(stderr);
    tr.my_start.share_to(tr.start);
    cerr<<"end_share"<<endl;fflush(stderr);
    rsold=rsnew;

  }
//cerr<<"tr.end"<<" ,"<<tr.my_end.vals<<endl;
//cerr<<"x"<<x.vals<<endl;
cerr<<"finished loop"<<endl;fflush(stderr);

Link<Leaf*>* lfptr=lf.first;
  for (int_s i=0;i<vlen;i++){
    tr.my_start[i]=x[i];
    lfptr->data->val=x[i];
    lfptr=lfptr->next;

  }
  cerr<<"coppied data"<<endl;fflush(stderr);
  tr.run();
  //cerr<<"tr.end"<<tr.my_end.vals<<endl;

  //cerr<<"here4 ";


  /*
  function x = conjgrad(A, b, x)
      r = b - A * x;
      p = r;
      rsold = r' * r;

      for i = 1:length(b)
          Ap = A * p;
          alpha = rsold / (p' * Ap);
          x = x + alpha * p;
          r = r - alpha * Ap;
          rsnew = r' * r;
          if sqrt(rsnew) < 1e-10
                break
          end
          p = r + (rsnew / rsold) * p;
          rsold = rsnew;
      end
  end*/




  //int s=16;
  //int_s this_up=((s*s-1)/mpi_size+1);
  //int_s squaresup{mpi_size*this_up};

  Array<double> base{s*s};
  //cerr<<base<<endl;
  Array<Array<Poly>> polys{logs};
  for (int i=0;i<logs;i++){
    polys[i]=Array<Poly>(1<<i);
    //cerr<<polys[i].len<<"len"<<endl;
  }
  lfptr=lf.first;
  int_s count=0;
  cerr<<"about to rapply"<<endl;fflush(stderr);
  while (lfptr!=0){
    //cerr<<lfptr->data->x<<",";
    if (lf.len*mpi_rank<=count*mpi_size & count*mpi_size<lf.len*(mpi_rank+1)){
      rapply(*(lfptr->data),polys,base, 1,1,0,logs);
    }

    lfptr=lfptr->next;
    count++;
  }
  cerr<<"Done"<<endl;fflush(stderr);
  Array<double> baseall{};
  if (mpi_rank==0){
    baseall=Array<double>{s*s};
  }

  MPI_Reduce(base.array,baseall.array,s*s,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  //rapply(const Leaf& l,Array<Array<Poly>>& polys, Array<double>& base,int_s x,int_s y,int d, int h){
  /*double rerr=0;
  for (int i=1;i<s-1;i++){
    for (int j=1;j<s-1;j++){
      if ((j!=15-i)&&((j!=5) || (i!=5))){
        rerr+=abs(base[i*s+j]*4-base[i*s+j-1]-base[i*s+j+1]-base[i*s+j-s]-base[i*s+j+s]);
      }

    }
  }*/
  //cerr<<endl<<"rerr "<<rerr<<endl;
  //cerr<<base;
  if (mpi_rank==0){
    fwrite(baseall.array,sizeof(double),s*s,stdout);
  }

  MPI_Finalize();
  return 0;
}
