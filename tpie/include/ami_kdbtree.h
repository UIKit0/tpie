// Copyright (C) 2001,2002 Octavian Procopiuc
//
// File:    ami_kdbtree.h
// Author:  Octavian Procopiuc <tavi@cs.duke.edu>
//
// K-D-B-tree definition and implementation. 
//
// $Id: ami_kdbtree.h,v 1.5 2003-04-21 04:07:26 tavi Exp $
//

#ifndef _KDBTREE_H
#define _KDBTREE_H

#include <iostream>
#include <ami_point.h>
#include <ami_kdtree.h>
#include <ami_kdbtree_base.h>

#define KDBTREE_HEADER_MAGIC_NUMBER 0xA9542F

// Forward declarations.
template<class coord_t, size_t dim, class BTECOLL> class Kdbtree_node;
template<class coord_t, size_t dim, class BTECOLL> class Kdbtree_leaf;

// The Kdbtree class. 
template<class coord_t, size_t dim, class Bin_node=Bin_node_default<coord_t, dim>, class BTECOLL = BTE_COLLECTION>
class Kdbtree {
public:

  //  typedef IdPoint<coord_t, dim> point_t;
  typedef Record<coord_t, size_t, dim> point_t;
  typedef Record<coord_t, size_t, dim> record_t;
  typedef Point<coord_t, dim> key_t;
  typedef AMI_STREAM<point_t> stream_t;
  typedef AMI_collection_single<BTECOLL> collection_t;
  typedef Kdbtree_node<coord_t, dim, BTECOLL> node_t;
  typedef Kdbtree_leaf<coord_t, dim, BTECOLL> leaf_t;
  typedef kdb_item_t<coord_t, dim> item_t;

  Kdbtree(const char *base_file_name, AMI_collection_type type, 
	  const Kdbtree_params& params);

  // Transform a kdtree into a kdbtree, in place. Returns true if
  // succesful, false otherwise (i.e., the status is not
  // KDBTREE_STATUS_KDTREE or the kdtree nodes contain too many keys).
  bool kd2kdb();

  size_t window_query(const point_t& p1, const point_t& p2, 
		      stream_t* stream);

  // Find a point.
  bool find(const point_t& p);

  // Insert p into the kdbtree. Return true if successful.
  bool insert(const point_t& p);

  // Traverse the tree in dfs preorder. Return next node and its level
  // (root is on level 0). Start the process by calling this with
  // level=-1.
  item_t dfs_preorder(int& level);

  // Set persistence. It passes per along to the two collections.
  void persist(persistence per);

  // Inquire the (real) parameters.
  const Kdbtree_params& params() const { return params_; }

  // Inquire the status.
  Kdbtree_status status() const { return status_; };

  // Inquire the size (number of points stored).
  size_t size() const { return header_.size; }

  // Inquire the mbr_lo point.
  point_t mbr_lo() const { return header_.mbr_lo; }

  // Inquire the mbr_hi point.
  point_t mbr_hi() const { return header_.mbr_hi; }

  // Inquire the statistics.
  const tpie_stats_tree &stats();

  // Inquire the leaf block size (in bytes)
  size_t leaf_block_size() const { return pcoll_leaves_->block_size(); }

  // Inquire the node block size (in bytes)
  size_t node_block_size() const { return pcoll_nodes_->block_size(); }

  // Destructor.
  ~Kdbtree();

  node_t* fetch_node(AMI_bid bid = 0);
  leaf_t* fetch_leaf(AMI_bid bid = 0);
  void release_node(node_t* q);
  void release_leaf(leaf_t* q);

  class header_type {
  public:
    unsigned int magic_number;
    point_t mbr_lo;
    point_t mbr_hi;
    AMI_bid root_bid;
    Node_type root_type;
    size_t size;
    
    header_type():
      magic_number(KDBTREE_HEADER_MAGIC_NUMBER), mbr_lo(0), mbr_hi(0), 
      root_bid(0), root_type(BLOCK_LEAF), size(0) {}
  };

  typedef header_type header_t;

protected:

  // Function object for the node cache write out.
  class remove_node {
  public:
    void operator()(node_t* p) { delete p; }
  };
  // Function object for the leaf cache write out.
  class remove_leaf { 
  public:
    void operator()(leaf_t* p) { delete p; }
  };

  // The node cache.
  AMI_CACHE_MANAGER<node_t*, remove_node>* node_cache_;
  // The leaf cache.
  AMI_CACHE_MANAGER<leaf_t*, remove_leaf>* leaf_cache_;

  // The collection storing the leaves.
  collection_t * pcoll_leaves_;

  // The collection storing the internal nodes (could be the same).
  collection_t * pcoll_nodes_;

  // Critical information: root bid and type, mbr, size (will be
  // stored into the header of the nodes collection).
  header_type header_;

  // The status.
  Kdbtree_status status_;

  // Run-time parameters.
  Kdbtree_params params_;

  // Stack to store the path to a leaf.
  stack<path_stack_item_t<coord_t, dim> > path_stack_;

  // Stack for dfs_preorder
  stack<path_stack_item_t<coord_t, dim> > dfs_stack_;

  // Statistics object.
  tpie_stats_tree stats_;

  bool insert_empty(const point_t& p);

  size_t window_query(const item_t& ki, const region_t<coord_t, dim>& r,
		      stream_t* stream);

  // Various initialization common to all constructors.
  void shared_init(const char* base_file_name, AMI_collection_type type);
  void kd2kdb(const item_t& ki, b_vector<item_t >& bv);
  void kd2kdb_node( Kdtree_node<coord_t, dim, Bin_node, BTECOLL>* bn, 
		    size_t i, region_t<coord_t, dim> r, 
		    b_vector<item_t >& bv, size_t& bv_pos);
  inline size_t split_dim_longest_span(const path_stack_item_t<coord_t, dim>& top);
  bool split_leaf_and_insert(const path_stack_item_t<coord_t, dim>& top, 
			     leaf_t* bl, item_t& ki1, 
			     item_t& ki2, const point_t& p);
  void split_leaf(coord_t sp, size_t d, const item_t& kis, 
		  item_t& ki1, item_t& ki2);
  void split_node_and_insert(const path_stack_item_t<coord_t, dim>& top,
			     item_t& ki1, item_t& ki2, const item_t& ki);
  void split_node(coord_t sp, size_t d, const item_t& kis, 
		  item_t& ki1, item_t& ki2);
  void find_split_position(const path_stack_item_t<coord_t, dim>& top, 
			   coord_t& sp, size_t& d);
  // Empty the path stack.
  inline void empty_stack(bool update_weight = false);
};

struct Kdbtree_leaf_info {
  size_t size;
  AMI_bid next;
  size_t split_dim;
};

template<class coord_t, size_t dim, class BTECOLL>
class Kdbtree_leaf: public AMI_block<Record<coord_t, size_t, dim>, Kdbtree_leaf_info, BTECOLL> {
public:
  typedef Record<coord_t, size_t, dim> point_t;
  typedef Record<coord_t, size_t, dim> record_t;
  typedef AMI_STREAM<point_t> stream_t;
  typedef AMI_collection_single<BTECOLL> collection_t;
  typedef Kdbtree_leaf_info info_t;

  static size_t el_capacity(size_t block_size);

  Kdbtree_leaf(collection_t* pcoll, AMI_bid bid = 0): 
    AMI_block<Record<coord_t, size_t, dim>, Kdbtree_leaf_info, BTECOLL>(pcoll, 0, bid) {
    if (bid == 0) {
      size() = 0;
      next() = 0;
      split_dim() = 0;
    }      
  }

  // Number of points stored in this leaf.
  size_t& size() { return info()->size; }
  const size_t& size() const { return info()->size; }

  // The weight of a leaf is the size. Just for symmetry with the
  // nodes.
  const size_t& weight() const { return info()->size; }
  
  // Next leaf. All leaves of a tree are chained togther for easy
  // retrieval.
  const AMI_bid& next() const { return info()->next; }
  AMI_bid& next() { return info()->next; }

  size_t& split_dim() { return info()->split_dim; }
  const size_t& split_dim() const { return info()->split_dim; }
 
  // Maximum number of points that can be stored in this leaf.
  size_t capacity() const { return el.capacity(); }

  // Find a point. Return the index of the point found in the el
  // vector (if not found, return size()).
  size_t find(const point_t &p) const {
    size_t i = 0;
    while (i < size()) {
      if (p == el[i])
	break;
      i++;
    }
    return i; 
  }

  size_t window_query(const point_t &lop, const point_t &hip,
		      stream_t* stream) const {
    size_t i, di, result = 0;
    for (i = 0; i < size(); i++) {
      // Test on all dimensions.
      if (lop < el[i] && el[i] < hip) {
	result++;
	if (stream != NULL)
	  stream->write_item(el[i]);
      }
    }
    return result;
  }

  // Insert a point, assuming the leaf is not full.
  bool insert(const point_t &p) {
    assert(size() < el.capacity());
    if (size() > 0 && find(p) < size())
      return false;

    el[size()] = p;
    size()++;
    dirty() = 1;
    return true;
  }

  bool erase(const point_t &p) {
    bool ans = false;
    size_t idx;
    if ((idx = find(p)) < size()) {
      if (idx < size() - 1) {
	// Copy the last item indo pos idx. We could use el.erase() as
	// well, but that's slower. Here order is not important.
	el[idx] = el[size()-1];
      }
      size()--;
      ans = true;
      dirty() = 1;
    }
    return ans;
  }

  // Sort points on the given dimension.
  void sort(size_t d) {
    typename POINT::cmp cmpd(d);
    std::sort(&el[0], &el[0] + size(), cmpd);
  }

  // Find median point on the given dimension. Return the index of the
  // median in the el vector.
  size_t find_median(size_t d) {
    sort(d);
    size_t ans = (size() - 1) / 2; // preliminary median.
    ///    while ((ans + 1 < size()) && (cmpd.compare(el[ans], el[ans+1]) == 0))
    while ((ans + 1 < size()) && (el[ans][d] == el[ans+1][d]))
      ans++;
    return ans;
  }
};


struct Kdbtree_node_info {
  size_t size;
  size_t weight;
  size_t split_dim;
};

// The Kdbtree_node class. Inherits from Kdtree_node and adds new
// find() and insert() procedures.
template<class coord_t, size_t dim, class BTECOLL>
class Kdbtree_node: public AMI_block<kdb_item_t<coord_t, dim>, Kdbtree_node_info, BTECOLL> {
public:
  
  typedef Record<coord_t, size_t, dim> point_t;
  typedef AMI_STREAM<point_t> stream_t;
  typedef AMI_collection_single<BTECOLL> collection_t;
  typedef kdb_item_t<coord_t, dim> item_t;
  typedef Kdbtree_node_info info_t;

  static size_t el_capacity(size_t block_size);

  // A node is an AMI_block containing kdb_item_t's as elements and no links.
  Kdbtree_node(collection_t* pcoll, AMI_bid bid = 0):
    AMI_block<kdb_item_t<coord_t, dim>, Kdbtree_node_info, BTECOLL>(pcoll, 0, bid) {
    if (bid == 0) {
      size() = 0;
      weight() = 0;
      //      split_dim() = 0;
    }
  }

  // Number of kdb_item_t's stored in this node.
  size_t& size() { return info()->size; }
  const size_t& size() const { return info()->size; }

  // Weight (ie, number of points stored in the subtree rooted at this
  // node).
  size_t& weight() { return info()->weight; }
  const size_t& weight() const { return info()->weight; }

  // Splitting dimension.
  size_t& split_dim() { return info()->split_dim; }
  const size_t& split_dim() const { return info()->split_dim; }

  // Maximum number of kdb_item_t's that can be stored in this node.
  size_t capacity() const { return el.capacity(); }

  // Find the index of the kdb_item_t containing the given point. If
  // no item contains the point, return size().
  size_t find(const point_t& p) {
    size_t i;
    for (i = 0; i < size(); i++) {
      if (el[i].region.contains(p.key))
	break;
    }
    return i;
  }

  // Insert a kdb_item_t in this node.
  bool insert(const item_t& ki) {
    //    assert(size() < el.capacity());
    el[size()] = ki;
    size()++;
    dirty() = 1;
    return true;
  }
};


////////////////////////////////////////////////////////////
///////////////     ***Implementation***    ////////////////
////////////////////////////////////////////////////////////

#define KDBTREE       Kdbtree<coord_t, dim, Bin_node, BTECOLL>
#define KDBTREE_NODE  Kdbtree_node<coord_t, dim, BTECOLL>
#define KDBTREE_LEAF  Kdbtree_leaf<coord_t, dim, BTECOLL>
#define REGION        region_t<coord_t, dim>
#define KDB_ITEM      kdb_item_t<coord_t, dim>
#define STACK_ITEM    path_stack_item_t<coord_t, dim>
#undef TPLOG
#define TPLOG(msg) 
//   (LOG_APP_DEBUG(msg),LOG_FLUSH_LOG)

//////////////////////////////////////
////////// **Kdbtree_leaf** //////////
//////////////////////////////////////

template<class coord_t, size_t dim, class BTECOLL>
size_t KDBTREE_LEAF::el_capacity(size_t block_size) {
  return AMI_block<Record<coord_t, size_t, dim>, Kdbtree_leaf_info, BTECOLL>::el_capacity(block_size, 0);
}

//////////////////////////////////////
////////// **Kdbtree_node** //////////
//////////////////////////////////////

template<class coord_t, size_t dim, class BTECOLL>
size_t KDBTREE_NODE::el_capacity(size_t block_size) {
  return AMI_block<KDB_ITEM, Kdbtree_node_info, BTECOLL>::el_capacity(block_size, 0);
}


//////////////////////////////////////
///////////// **Kdbtree** ////////////
//////////////////////////////////////

//// *Kdbtree::Kdbtree* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
KDBTREE::Kdbtree(const char *base_file_name, AMI_collection_type type, 
		 const Kdbtree_params& params): header_(), params_(params) {

  shared_init(base_file_name, type);
}

//// *Kdbtree::shared_init* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::shared_init(const char* base_file_name, AMI_collection_type type) {

  assert(base_file_name != NULL);
  char collname[124];

  // Open the two block collections.
  strncpy(collname, base_file_name, 124 - 2);
  strcat(collname, ".l");
  pcoll_leaves_ = new collection_t(collname, type, params_.leaf_block_factor);

  strncpy(collname, base_file_name, 124 - 2);
  strcat(collname, ".n");
  pcoll_nodes_ = new collection_t(collname, type, params_.node_block_factor);

  if (pcoll_nodes_->status() != AMI_COLLECTION_STATUS_VALID ||
      pcoll_leaves_->status() != AMI_COLLECTION_STATUS_VALID) {
    status_ = KDBTREE_STATUS_INVALID;
    delete pcoll_leaves_;
    delete pcoll_nodes_;
    return;
  }

  // Read the header info, if relevant.
  if (pcoll_leaves_->size() != 0) {
    unsigned int magic = *((unsigned int *) pcoll_nodes_->user_data());
    if (magic == KDTREE_HEADER_MAGIC_NUMBER) {
      status_ = KDBTREE_STATUS_KDTREE;
    } else if (magic == KDBTREE_HEADER_MAGIC_NUMBER) {
      status_ = KDBTREE_STATUS_VALID;
      //      header_ = *((header_type *) pcoll_nodes_->user_data());
      memcpy((void *)(&header_), pcoll_nodes_->user_data(), sizeof(header_));
      // TODO: sanity checks on the header.
    } else {
      status_ = KDBTREE_STATUS_INVALID;
      LOG_WARNING_ID("Invalid kdbtree magic number:"<<magic);
      delete pcoll_leaves_;
      delete pcoll_nodes_;
      return;
    }
  }

  // Initialize the caches.
  leaf_cache_ = new AMI_CACHE_MANAGER<KDBTREE_LEAF*, remove_leaf>(params_.leaf_cache_size, 4);
  node_cache_ = new AMI_CACHE_MANAGER<KDBTREE_NODE*, remove_node>(params_.node_cache_size, 1);
  
  // Give meaningful values to parameters, if necessary.
  size_t leaf_capacity = KDBTREE_LEAF::el_capacity(pcoll_leaves_->block_size());
  if (params_.leaf_size_max == 0 || params_.leaf_size_max > leaf_capacity)
    params_.leaf_size_max = leaf_capacity;
  TPLOG("  Kdbtree::shared_init leaf_size_max="<<params_.leaf_size_max<<"\n");
  
  size_t node_capacity = KDBTREE_NODE::el_capacity(pcoll_nodes_->block_size());
  if (params_.node_size_max == 0 || params_.node_size_max > node_capacity)
    params_.node_size_max = node_capacity;
  TPLOG("  Kdbtree::shared_init node_size_max="<<params_.node_size_max<<"\n");
  TPLOG("  sizeof(kdb_item_t)="<<sizeof(KDB_ITEM)<<"\n");

  // Set the right block factor parameters for the case of an existing tree.
  params_.leaf_block_factor = pcoll_leaves_->block_factor();
  params_.node_block_factor = pcoll_nodes_->block_factor();
}


//// *Kdbtree::kd2kdb* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
bool KDBTREE::kd2kdb() {
  TPLOG("Kdbtree::kd2kdb() Entering " << "\n");
  if (status_ != KDBTREE_STATUS_KDTREE) {
    LOG_WARNING_ID("  kd2kdb: status is not KDBTREE_STATUS_KDTREE. operation aborted.");
    return false;
  }

  bool ans = true;
  typename KDTREE::header_t kdheader;
  // = *((KDTREE::header_t *) pcoll_nodes_->user_data());
  memcpy((void *)(&kdheader), pcoll_nodes_->user_data(), sizeof(kdheader));
  header_.root_bid = kdheader.root_bid;
  header_.root_type = kdheader.root_type;
  header_.size = kdheader.size;
  header_.mbr_lo = kdheader.mbr_lo;
  header_.mbr_hi = kdheader.mbr_hi;

  REGION r; // Unbounded region, corresponding to the root.
  KDB_ITEM ki(r, header_.root_bid, header_.root_type);
  if (ki.type == BLOCK_NODE) {

    // Create temporary buffer.
    KDBTREE_NODE* buffer = new KDBTREE_NODE(pcoll_nodes_);
    // Do the job.
    kd2kdb(ki, buffer->el);
    // Dispose of the temporary buffer.
    buffer->persist(PERSIST_DELETE);
    delete buffer;
    
  } else {
    // Just a leaf. Nothing to do.
  }

  if (status_ == KDBTREE_STATUS_KDTREE)
    status_ = KDBTREE_STATUS_VALID;
  else
    ans = false;

  TPLOG("Kdbtree::kd2kdb() Exiting ans=" << ans << "\n");
  return ans;
}


//// *Kdbtree::window_query* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
size_t KDBTREE::window_query(const POINT &p1, const POINT& p2, 
			     POINT_STREAM* stream) {
  TPLOG("  query window: "<<p1[0]<<" "<<p1[1]<<" "<<p2[0]<<" "<<p2[1]<<"\n");
  // The number of points found.
  size_t result = 0;

  // We can afford to do some error checking, since this is usually a
  // lengthy operation.
  if (status_ != KDBTREE_STATUS_VALID) {
    LOG_WARNING_ID("  window_query: tree is invalid or not loaded. query aborted.");
    return result;
  }

  // TODO...
  REGION r(p1.key, p2.key);
  REGION rr;
  KDB_ITEM ki(rr, header_.root_bid, header_.root_type);
  result = window_query(ki, r, stream);

  return result;
}

//// *Kdbtree::window_query* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
size_t KDBTREE::window_query(const KDB_ITEM& ki, const REGION& r, POINT_STREAM* stream) {
  TPLOG("  window query recusion: " << ki << "\n");
  size_t result = 0;
  if (ki.type == BLOCK_NODE) {
    KDBTREE_NODE* bn = fetch_node(ki.bid);
    size_t i;
    for (i = 0; i < bn->size(); i++) {
      if (bn->el[i].region.intersects(r))
	result += window_query(bn->el[i], r, stream);
    }
    release_node(bn);
  } else {
    assert(ki.type == BLOCK_LEAF);
    KDBTREE_LEAF* bl = fetch_leaf(ki.bid);
    result += bl->window_query(r.point_lo(), r.point_hi(), stream);
    release_leaf(bl);
  }
  return result;
}

//// *Kdbtree::kd2kdb* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::kd2kdb(const KDB_ITEM& ki, b_vector<KDB_ITEM >& bv) {
  TPLOG("Kdbtree::kd2kdb Entering bid=" << ki.bid << "\n");

  KDTREE_NODE *bno;
  KDBTREE_NODE *bn;
  Node_type ni_type;
  bno = new KDTREE_NODE(pcoll_nodes_, ki.bid);
  if (bno->size() + 1 > params_.node_size_max) {
    LOG_FATAL_ID("  kd2kdb: wrong kdtree node size;");
    LOG_FATAL_ID("  kd2kdb: kdbtree node capacity: " << params_.node_size_max);
    LOG_FATAL_ID("  kd2kdb: max kdtree node size allowed: " << params_.node_size_max-1);
    LOG_FATAL_ID("  kd2kdb: found kdtree node with size: " << bno->size());
    LOG_FATAL_ID("  kd2kdb: operation aborted.");
    delete bno;
    status_ = KDBTREE_STATUS_INVALID;
  } else {
    size_t free_pos = 0;
    kd2kdb_node(bno, 0, ki.region, bv, free_pos);
    delete bno;
    TPLOG("  Kdbtree_node size="<<free_pos<<"\n");

    // Open the same block, but as a KDBTREE_NODE.
    bn = new KDBTREE_NODE(pcoll_nodes_, ki.bid);
    bn->el.copy(0, free_pos, bv, 0);
    bn->size() = free_pos;
    bn->split_dim() = 0;
    
    for (size_t i = 0; i < free_pos; i++) {
      if (bn->el[i].type == BLOCK_NODE && status_ != KDBTREE_STATUS_INVALID)
	kd2kdb(bn->el[i], bv);
    }
    delete bn;
  }

  TPLOG("Kdbtree::kd2kdb Exiting bid=" << ki.bid << "\n");
}

//// *Kdbtree::kd2kdb_node* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::kd2kdb_node(KDTREE_NODE* bn, size_t i, REGION r, 
			  b_vector<KDB_ITEM >& bv, size_t& bv_pos) {
  TPLOG("Kdbtree::kd2kdb_node Entering " << "\n");

  REGION ni_r;
  size_t ni;
  Node_type ni_type;

  ni_r = r;
  bn->el[i].get_low_child(ni, ni_type);
  ni_r.cutout_hi(bn->el[i].get_discriminator_val(), bn->el[i].get_discriminator_dim());
  if (ni_type == BIN_NODE) {
    kd2kdb_node(bn, ni, ni_r, bv, bv_pos);
  } else {
    assert(ni_type == BLOCK_LEAF || ni_type == BLOCK_NODE);
    bv[bv_pos].region = ni_r;
    bv[bv_pos].bid = bn->lk[ni];
    bv[bv_pos].type = ni_type;
    bv_pos++;
  }

  ni_r = r;
  bn->el[i].get_high_child(ni, ni_type);
  ni_r.cutout_lo(bn->el[i].get_discriminator_val(), bn->el[i].get_discriminator_dim());
  if (ni_type == BIN_NODE) {
    kd2kdb_node(bn, ni, ni_r, bv, bv_pos);
  } else {
    assert(ni_type == BLOCK_LEAF || ni_type == BLOCK_NODE);
    bv[bv_pos].region = ni_r;
    bv[bv_pos].bid = bn->lk[ni];
    bv[bv_pos].type = ni_type;
    bv_pos++;
  }
  TPLOG("Kdbtree::kd2kdb_node Exiting " << "\n");  
}

//// *Kdbtree::insert_empty* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
bool KDBTREE::insert_empty(const Record<coord_t, size_t, dim>& p) {
  bool ans;
  KDBTREE_LEAF* bl = fetch_leaf();
  ans = bl->insert(p);
  assert(ans);
  bl->split_dim() = 0;
  bl->next() = 0;

  header_.size = 1;
  header_.root_bid = bl->bid();
  header_.root_type = BLOCK_LEAF;
  header_.mbr_lo = p;
  header_.mbr_lo.id() = 1;
  header_.mbr_hi = p;
  header_.mbr_hi.id() = 1;

  status_ = KDBTREE_STATUS_VALID;
  release_leaf(bl);
  return ans;
}

//// *Kdbtree::find* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
bool KDBTREE::find(const Record<coord_t, size_t, dim>& p) {

  TPLOG("Kdbtree::find Entering "<<"\n");

  if (header_.size == 0)
    return false;

  bool ans;
  size_t i;

  KDBTREE_NODE* bn;
  REGION r;
  KDB_ITEM ki(r, header_.root_bid, header_.root_type);
  STACK_ITEM si(ki, 0);

  while (si.item.type == BLOCK_NODE) {
    bn = fetch_node(si.item.bid);

    i = bn->find(p);
    assert(i < bn->size());
    si.item = bn->el[i];
    //    si.el_idx = i;
    release_node(bn);
  }
  
  assert(si.item.type == BLOCK_LEAF);

  // Fetch the leaf.
  KDBTREE_LEAF* bl = fetch_leaf(si.item.bid);
  // Check whether item is in the leaf.
  ans = (bl->find(p) < bl->size());
  // Release the leaf.
  release_leaf(bl);

  TPLOG("Kdbtree::find Exiting ans="<<ans<<"\n");
  return ans;
}


//// *Kdbtree::insert* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
bool KDBTREE::insert(const Record<coord_t, size_t, dim>& p) {

  TPLOG("Kdbtree::insert Entering "<<"\n");

  // The first insertion is treated separately.
  if (header_.size == 0)
    return insert_empty(p);

  // Update the MBR.
  for (size_t i = 0; i < dim; i++) {
    header_.mbr_lo[i] = std::min(header_.mbr_lo[i], p[i]);
    header_.mbr_hi[i] = std::max(header_.mbr_hi[i], p[i]);
  }

  bool ans;
  KDBTREE_NODE* bn;
  REGION r; // Infinite region.
  KDB_ITEM ki(r, header_.root_bid, header_.root_type);
  size_t i;

  // Stack item; initially unbounded, corresponding to the root node.
  STACK_ITEM si(ki, 0);
  path_stack_.push(si);

  // Go down the tree until the appropriate leaf is found.
  while (si.item.type == BLOCK_NODE) {
    bn = fetch_node(si.item.bid);

    i = bn->find(p);
    assert(i < bn->size());
    si.item = bn->el[i];
    si.d = (si.d + 1) % dim;
    si.el_idx = i;

    release_node(bn);
    path_stack_.push(si);
  }

  // Make sure we reached a leaf.
  assert(si.item.type == BLOCK_LEAF);

  // Fetch the leaf.
  KDBTREE_LEAF* bl = fetch_leaf(si.item.bid);

  // Check for duplicate key. For now, just exit with false answer if found.
  if (bl->find(p) < bl->size()) { // Found.

    ans = false;
    release_leaf(bl);
    empty_stack(ans);

  } else if (bl->size() < params_.leaf_size_max) {
    // The very easy case. Just insert into the leaf.

    ans = bl->insert(p);
    release_leaf(bl);
    // TODO [****] this is wrong! Should be empty_stack(ans);
    //    empty_stack(false);
    empty_stack(ans);

  } else {
    // Need to split the leaf. Maybe some nodes, too.

    KDB_ITEM ki1, ki2;
    STACK_ITEM top;

    // Pop the leaf from the stack. 
    top = path_stack_.top();
    path_stack_.pop();

    // Split bl into two leaves, one of which is bl, insert p into the
    // appropriate leaf, store the kdb_item_t's pointing to the two
    // leaves in ki1 and ki2, and return true if insertion was
    // successful.
    ans = split_leaf_and_insert(top, bl, ki1, ki2, p);

    release_leaf(bl);

    bool done = false;
    size_t el_idx;

    while (!path_stack_.empty() && !done) {
      // Save top.el_idx.
      el_idx = top.el_idx;

      // Pop the next node on the path to the root.
      top = path_stack_.top();
      path_stack_.pop();

      // Fetch the node.
      bn = fetch_node(top.item.bid);
      // Update bn.
      bn->el[el_idx] = ki1;
      if (ans) bn->weight()++;
  
      if (bn->size() == params_.node_size_max) {

	// Split bn into two nodes, one of which is bn, insert ki2
	// into the appropriate node, store the resulting regions into
	// ki1 and ki2 (for next iteration).
	ki = ki2;
	release_node(bn);
	split_node_and_insert(top, ki1, ki2, ki);

      } else {
	// Insert ki2 into bn and exit the while loop.
	bn->insert(ki2);
	done = true;
	release_node(bn);
      }
    } // end of while.

    // Check if root was split.
    if (path_stack_.empty() && !done) {
      // Create new root node.
      bn = fetch_node();
      // Insert the two kdb_item's into this new root node.
      bn->insert(ki1);
      bn->insert(ki2);
      // Update the header information.
      header_.root_bid = bn->bid();
      header_.root_type = BLOCK_NODE;
      bn->split_dim() = 0;
      release_node(bn);
    } else {
      empty_stack(ans);
    }
  }

  if (ans) header_.size++;

  TPLOG("Kdbtree::insert Exiting ans="<<ans<<"\n");
  return ans;
}

//// *Kdbtree::dfs_preorder* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
KDB_ITEM KDBTREE::dfs_preorder(int& level) {

  // level signals the start/end of the traversal. All necessary state
  // information is kept on dfs_stack_.

  if (level == -1) {

    // Empty the stack. This allows restarts in the middle of a
    // traversal. All previous state information is lost.
    while (!dfs_stack_.empty())
      dfs_stack_.pop();
    REGION r; // Infinite region.
    KDB_ITEM ki(r, header_.root_bid, header_.root_type);
    dfs_stack_.push(STACK_ITEM(ki, 0, 0));

    level = dfs_stack_.size() - 1;
    return ki;

  } else {

    KDBTREE_NODE* bn;
    KDB_ITEM ki;
    if (dfs_stack_.top().item.type == BLOCK_NODE) {
      // Fetch the node ...
      bn = fetch_node(dfs_stack_.top().item.bid);
      // ... and get the appropriate child.
      ki = bn->el[dfs_stack_.top().el_idx];
      dfs_stack_.push(STACK_ITEM(ki, 0, 0));
      release_node(bn);
    } else { // i.e., dfs_stack_.top().item.type == BLOCK_LEAF
      // Remove the leaf from the stack.
      dfs_stack_.pop();
      bool done = false;
      while (!dfs_stack_.empty() && !done) {
	// Fetch the node ...
	bn = fetch_node(dfs_stack_.top().item.bid);

	(dfs_stack_.top().el_idx)++;
	if (dfs_stack_.top().el_idx < bn->size()) {
	  ki = bn->el[dfs_stack_.top().el_idx];
	  dfs_stack_.push(STACK_ITEM(ki, 0, 0));
	  done = true;
	} else
	  dfs_stack_.pop();

	release_node(bn);
      }
    }

    // Note: if stack is empty, level will be -1, signaling the end of
    // the traversal.
    level = dfs_stack_.size() - 1;
    return ki;
  }
}

//// *Kdbtree::split_leaf_and_insert* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
size_t KDBTREE::split_dim_longest_span(const STACK_ITEM& top) {
  size_t d;

  coord_t longest_span =
    (top.item.region.is_bounded_hi(0) ? top.item.region.hi(0): header_.mbr_hi[0]) - 
    (top.item.region.is_bounded_lo(0) ? top.item.region.lo(0): header_.mbr_lo[0]);
  d = 0;
  for (size_t i = 1; i < dim; i++)
    if ((top.item.region.is_bounded_hi(i) ? top.item.region.hi(i): header_.mbr_hi[i]) -
	(top.item.region.is_bounded_lo(i) ? top.item.region.lo(i): header_.mbr_lo[i]) > longest_span) {
      longest_span = 
	(top.item.region.is_bounded_hi(i) ? top.item.region.hi(i): header_.mbr_hi[i]) -
	(top.item.region.is_bounded_lo(i) ? top.item.region.lo(i): header_.mbr_lo[i]);
      d = i;
    }
  
  return d;
}

//// *Kdbtree::split_leaf_and_insert* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
bool KDBTREE::split_leaf_and_insert(const STACK_ITEM& top, KDBTREE_LEAF* bl,
				    KDB_ITEM& ki1, KDB_ITEM& ki2, const POINT& p) {
  // Get the median point.
  // Move points higher than the median in bl_hi.
  KDBTREE_LEAF* bl_hi = fetch_leaf();
  bool ans;
  size_t d;
  if (params_.split_heuristic == CYCLICAL)
    d = bl->split_dim();
  else if (params_.split_heuristic == LONGEST_SPAN)
    d = split_dim_longest_span(top);   
  else if (params_.split_heuristic == RANDOM)
    d = random() % dim;
  size_t med = bl->find_median(d); // the index of the median point in bl->el.
  POINT sp = bl->el[med];

  if (med + 1 >= bl->size()) {
    std::cerr << "\nbl->bid()=" << bl->bid() << ", bl->size()=" <<bl->size() << ", med=" << med << "\n";
    std::cerr << "bl: ";
    for (size_t i = 0; i < bl->size(); i++) {
      std::cerr << "[" << bl->el[i][0] << "," << bl->el[i][1] << "] ";
    }
    std::cerr << "\n";
  }
  assert(med + 1 < bl->size());
  bl_hi->size() = bl->size() - (med + 1); // the size of bl_hi
  bl->size() = med + 1; // the new size of bl

  // Cycle through dimensions.
  bl->split_dim() = (d + 1) % dim;
  bl_hi->split_dim() = (d + 1) % dim;
  // Update next pointers.
  bl_hi->next() = bl->next();
  bl->next() = bl_hi->bid();

  bl_hi->el.copy(0, bl_hi->size(), bl->el, med + 1); // copy points from bl to bl_hi.

  assert(top.item.type == BLOCK_LEAF);
  ki1 = top.item;
  ki1.region.cutout_hi(sp[d], d);
  assert(ki1.bid == bl->bid()); 

  ki2 = top.item;
  ki2.region.cutout_lo(sp[d], d);
  ki2.bid = bl_hi->bid();
  
  ans = (ki1.region.contains(p.key) ? bl: bl_hi)->insert(p);
  release_leaf(bl_hi);
  return ans;
}

//// *Kdbtree::split_leaf* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::split_leaf(coord_t sp, size_t d,  const KDB_ITEM& kis, 
			 KDB_ITEM& ki1, KDB_ITEM& ki2) {

  assert(kis.type == BLOCK_LEAF);
  KDBTREE_LEAF* bl = fetch_leaf(kis.bid);
  KDBTREE_LEAF* bl_hi = fetch_leaf();
  size_t bl_size = bl->size();
  POINT p;

  bl->size() = 0;
  assert(bl_hi->size() == 0);
  for (size_t i = 0; i < bl_size; i++) {
    p = bl->el[i];
    (sp < p[d] ? bl_hi: bl)->insert(p);
  }
  
  ki1 = kis;
  ki1.region.cutout_hi(sp, d);
  assert(ki1.bid == bl->bid());

  ki2 = kis;
  ki2.region.cutout_lo(sp, d);
  ki2.bid = bl_hi->bid();

  release_leaf(bl_hi);
  release_leaf(bl);
}

//// *Kdbtree::split_node_and_insert* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::split_node_and_insert(const STACK_ITEM& top,
				    KDB_ITEM& ki1, KDB_ITEM& ki2, const KDB_ITEM& ki) {
  // The split position.
  coord_t sp;
  // The split dimension.
  size_t d;

  // Find a split position and store it in sp.
  find_split_position(top, sp, d);

  // Split.
  split_node(sp, d, top.item, ki1, ki2);

  // Insert ki.
  KDBTREE_NODE *bn, *bn_hi;
  bn = fetch_node(ki1.bid);
  bn_hi = fetch_node(ki2.bid);
  int pos = ki.region.relative_to_plane(sp, d);
  if (pos == -1) {
    assert(bn->size() < params_.node_size_max);
    bn->insert(ki);
  } else if (pos == 1) {
    assert(bn_hi->size() < params_.node_size_max);
    bn_hi->insert(ki);
  } else {
    assert(bn->size() < params_.node_size_max && bn_hi->size() < params_.node_size_max);
    KDB_ITEM lki1, lki2;

    if (ki.type == BLOCK_LEAF)
      split_leaf(sp, d, ki, lki1, lki2);
    else
      split_node(sp, d, ki, lki1, lki2);

    bn->insert(lki1);
    bn_hi->insert(lki2);
  }

  bn->split_dim() = (d + 1) % dim;
  bn_hi->split_dim() = (d + 1) % dim;
  release_node(bn_hi);
  release_node(bn);
}

//// *Kdbtree::find_split_position* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::find_split_position(const STACK_ITEM& top, coord_t& sp, size_t& d) {
  KDBTREE_NODE *bn;
  bn = fetch_node(top.item.bid);

  // Find split dimension.
#if 0
  if (params_.split_heuristic == CYCLICAL)
    d = bn->split_dim();
  else if (params_.split_heuristic == LONGEST_SPAN)
    d = split_dim_longest_span(top);   
  else if (params_.split_heuristic == RANDOM)
    d = random() % dim;
#else
  d = bn->split_dim();
#endif

  vector<coord_t > cv(0);
  size_t unbounded = 0;
  // Collect all low boundaries from bn and, if they are bounded,
  // store them in cv. The unbounded ones are counted only.
  for (size_t i = 0; i < bn->size(); i++) {
    if (bn->el[i].region.is_bounded_lo(d))
      cv.push_back(bn->el[i].region.lo(d));
    else
      unbounded++;
  }
  assert(cv.size() > 0);
  // Sort.
  sort(cv.begin(), cv.end());
  // Get median value.
  size_t median = (bn->size() / 2 > unbounded ? bn->size() / 2 - unbounded: 0);
  // Make sure we don't return the leftmost boundary.
  if (unbounded == 0)
    while (cv[median] == cv[0])
      median++;
  assert(median < cv.size());
  sp = cv[median];
  release_node(bn);
}

//// *Kdbtree::split_node* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::split_node(coord_t sp, size_t d, const KDB_ITEM& kis, 
			 KDB_ITEM& ki1, KDB_ITEM& ki2) {
  assert(kis.type == BLOCK_NODE);
  KDBTREE_NODE *bn, *bn_hi;

  bn = fetch_node(kis.bid);
  bn_hi = fetch_node();
  KDB_ITEM ki;
  ///  size_t d = bn->split_dim();

  ki1 = ki2 = kis;
  ki1.region.cutout_hi(sp, d);
  assert(ki1.bid == bn->bid());
  ki2.region.cutout_lo(sp, d);
  ki2.bid = bn_hi->bid();

  size_t next_free_lo = 0, next_free_hi = 0, i;
  size_t bn_size = bn->size();
  int pos;
  bn->size() = 0;
  assert(bn_hi->size() == 0);
  // Cycle through dimensions.
  ///  bn->split_dim() = bn_hi->split_dim() = (d + 1) % dim;

  for (i = 0; i < bn_size; i++) {
    pos = bn->el[i].region.relative_to_plane(sp, d);
    ki = bn->el[i];
    if (pos == -1) {
      // Move to bn.
      bn->insert(ki);
    } else if (pos == 1) {
      // Move to bn_hi;
      bn_hi->insert(ki);
    } else {
      // The hard case: intersection.
      if (ki.type == BLOCK_LEAF)
	split_leaf(sp, d, ki, bn->el[bn->size()], bn_hi->el[bn_hi->size()]);
      else // BLOCK_NODE
	split_node(sp, d, ki, bn->el[bn->size()], bn_hi->el[bn_hi->size()]);
      bn->size()++;
      bn_hi->size()++;
    }
  }
  
  release_node(bn_hi);
  release_node(bn);
}

//// *Kdbtree::empty_stack* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::empty_stack(bool update_weight) {
  node_t* bn;
  while (!path_stack_.empty()) {
    if (update_weight && path_stack_.top().item.type == BLOCK_NODE) {
      bn = fetch_node(path_stack_.top().item.bid);
      bn->weight()++;
      release_node(bn);
    }
    path_stack_.pop(); 
  }
}

//// *Kdbtree::persist* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::persist(persistence per) {
  pcoll_leaves_->persist(per);
  pcoll_nodes_->persist(per);
}

//// *Kdbtree::fetch_node* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
KDBTREE_NODE* KDBTREE::fetch_node(AMI_bid bid) {
  KDBTREE_NODE* q;
  stats_.record(NODE_FETCH);
  // Warning: using short-circuit evaluation. Order is important.
  if ((bid == 0) || !node_cache_->read(bid, q)) {
    q = new KDBTREE_NODE(pcoll_nodes_, bid);
  }
  return q;
}

//// *Kdbtree::fetch_leaf* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
KDBTREE_LEAF* KDBTREE::fetch_leaf(AMI_bid bid) {
  KDBTREE_LEAF* q;
  stats_.record(LEAF_FETCH);
  // Warning: using short-circuit evaluation. Order is important.
  if ((bid == 0) || !leaf_cache_->read(bid, q)) {
    q = new KDBTREE_LEAF(pcoll_leaves_, bid);
  }
  return q;
}

//// *Kdbtree::release_node* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::release_node(KDBTREE_NODE* q) {
  stats_.record(NODE_RELEASE);
  if (q->persist() == PERSIST_DELETE)
    delete q;
  else
    node_cache_->write(q->bid(), q);
}

//// *Kdbtree::release_leaf* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
void KDBTREE::release_leaf(KDBTREE_LEAF* q) {
  stats_.record(LEAF_RELEASE);
  if (q->persist() == PERSIST_DELETE)
    delete q;
  else
    leaf_cache_->write(q->bid(), q);
}

//// *Kdbtree::stats* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
const tpie_stats_tree &KDBTREE::stats() {
  node_cache_->flush();
  leaf_cache_->flush();
  stats_.set(LEAF_READ, pcoll_leaves_->stats().get(BLOCK_GET));
  stats_.set(LEAF_WRITE, pcoll_leaves_->stats().get(BLOCK_PUT));
  stats_.set(LEAF_CREATE, pcoll_leaves_->stats().get(BLOCK_NEW));
  stats_.set(LEAF_DELETE, pcoll_leaves_->stats().get(BLOCK_DELETE));
  stats_.set(LEAF_COUNT, pcoll_leaves_->size());
  stats_.set(NODE_READ, pcoll_nodes_->stats().get(BLOCK_GET));
  stats_.set(NODE_WRITE, pcoll_nodes_->stats().get(BLOCK_PUT));
  stats_.set(NODE_CREATE, pcoll_nodes_->stats().get(BLOCK_NEW));
  stats_.set(NODE_DELETE, pcoll_nodes_->stats().get(BLOCK_DELETE));
  stats_.set(NODE_COUNT, pcoll_nodes_->size());
  return stats_;
}

//// *Kdbtree::~Kdbtree* ////
template<class coord_t, size_t dim, class Bin_node, class BTECOLL>
KDBTREE::~Kdbtree() {
  
  if (status_ == KDBTREE_STATUS_VALID) {
    // Write initialization info into the pcoll_nodes_ header.
    //    *((header_type *) pcoll_nodes_->user_data()) = header_;
    memcpy(pcoll_nodes_->user_data(), (void *)(&header_), sizeof(header_));
  }

  delete node_cache_;
  delete leaf_cache_;

  // Delete the two collections.
  delete pcoll_leaves_;
  delete pcoll_nodes_;
 
}

#endif // _KDBTREE_H
