// Copyright (c) 1994 Darren Erik Vengroff
//
// File: mm_register.cpp
// Author: Darren Erik Vengroff <dev@cs.duke.edu>
// Created: 5/31/94
//

// A simple registration based memory manager.

static char mm_register_id[] = "$Id: mm_register.cpp,v 1.3 1994-09-26 19:28:08 darrenv Exp $";

#include <tpie_log.h>
#include <tpie_assert.h>

#define MM_IMP_REGISTER
#include <mm.h>


MM_register::MM_register()
{
    instances++;

    //*tpl << setpriority(TP_LOG_ASSERT);

    tp_assert(instances == 1,
              "Only 1 instance of MM_register_base should exist.");

    // Why does this cause seg faults?
    // LOG_INFO("Created MM_register object.\n");
}


MM_register::~MM_register(void)
{
    tp_assert(instances == 1,
              "Only 1 instance of MM_register_base should exist.");

    instances--;
}


MM_err MM_register::register_allocation(size_t sz)
{
    if (sz > remaining) {
        return MM_ERROR_INSUFFICIENT_SPACE;
    }

    remaining -= sz;

    LOG_INFO("Allocated ");
    LOG_INFO((unsigned int)sz);
    LOG_INFO("; ");
    LOG_INFO((unsigned int)remaining);
    LOG_INFO(" remaining.\n");
    
    return MM_ERROR_NO_ERROR;
}


MM_err MM_register::register_deallocation(size_t sz)
{
    if (sz + remaining > max_sz) {
        return MM_ERROR_UNDERFLOW;
    }

    remaining += sz;

    LOG_INFO("De-allocated ");
    LOG_INFO((unsigned int)sz);
    LOG_INFO("; ");
    LOG_INFO((unsigned int)remaining);
    LOG_INFO(" now available.\n");
    
    return MM_ERROR_NO_ERROR;
}


MM_err MM_register::available(size_t *sz)
{
    *sz = remaining;
    return MM_ERROR_NO_ERROR;    
}


MM_err MM_register::resize_heap(size_t sz)
{
    if (max_sz - remaining > sz) {
        return MM_ERROR_EXCESSIVE_ALLOCATION;
    } else {
        // These are unsigned, so be careful.
        if (sz < max_sz) {
            remaining -= max_sz - sz;
        } else {
            remaining += sz - max_sz;
        }
        max_sz = sz;
        return MM_ERROR_NO_ERROR;
    }
}



// The number of instances.  Implicitly set to zero.
int MM_register::instances;

// The actual memory manager.
MM_register MM_manager;

// The counter of mm_register_init instances.  It is implicity set to 0.
unsigned int mm_register_init::count;

// The constructor and destructor that ensure that the log files are
// created exactly once, and destroyed when appropriate.
mm_register_init::mm_register_init(void)
{
    if (count++ == 0) {
        MM_manager.remaining = MM_manager.max_sz = MM_DEFAULT_MM_SIZE;
    }
}


mm_register_init::~mm_register_init(void)
{
    --count;
}

