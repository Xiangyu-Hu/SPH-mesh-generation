#ifndef MY_MEMORY_POOL_H
#define MY_MEMORY_POOL_H

#define TBB_PREVIEW_MEMORY_POOL 1
#include "tbb/memory_pool.h"  //HAN
#include "tbb/mutex.h"
#include <list>

using namespace std;
using namespace tbb;

//------------------------------
//my memory pool
//------------------------------
template<class T>
class Mypool{

  T sample;
  tbb::memory_pool< std::allocator<T> > my_pool;           //memory pool
  typedef tbb::memory_pool_allocator<T> pool_allocator_t;  //memory allocator
  std::list<T, pool_allocator_t> data_list;                //list of all nodes allocated
  std::list<T*> free_list;                                 //list of all free nodes
  tbb::mutex gMutex;

public:

  //constructor
  Mypool():data_list( (pool_allocator_t( my_pool )) ) {};
  //deconstructor
  ~Mypool(){
    //my_pool.recycle();
  };
  //prepare an avaliable node
  T* malloc()
  {
    tbb::mutex::scoped_lock lock(gMutex);
    if(free_list.empty()){
      data_list.push_back(sample);
      return (&data_list.back());
    }
    else
    {
      T* result = free_list.front();
      free_list.pop_front();
      return result;
    }
  };
  //relinquish an unused node
  void free(T* ptr)
  {
    tbb::mutex::scoped_lock lock(gMutex);
    free_list.push_back(ptr);
  };
  //return the total number of nodes allocated
  int capacity()
  {
    tbb::mutex::scoped_lock lock(gMutex);
    return data_list.size();
  };
  //return the number of current available nodes
  int available_node()
  {
    tbb::mutex::scoped_lock lock(gMutex);
    return free_list.size();
  };
};

#endif
