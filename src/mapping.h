#ifndef MAPPING_H
#define MAPPING_H

#include <stddef.h>

struct mapped_region {
    void *addr;
    size_t length;
};

int map_file(const char *pathname, struct mapped_region *mapped_region);

int unmap_file(struct mapped_region *mapped_region);

#endif
