#include "mapping.h"

#include <fcntl.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/mman.h>
#endif

int map_file(const char *pathname, struct mapped_region *mapped_region) {
    int retval = 0;
 // Get file status
    struct stat sb;
    if (stat(pathname, &sb) == -1) {
        return -1;
    }
 // Test if file is a regular file
    if (!S_ISREG(sb.st_mode)) {
        return -1;
    }
 // Get file length
    mapped_region->length = sb.st_size;
 // Map file
#ifdef _WIN32
    HANDLE hFile = CreateFileA(
        pathname,
        GENERIC_READ,
        FILE_SHARE_READ,
        NULL,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL,
        NULL
    );
    if (hFile == INVALID_HANDLE_VALUE) {
        return -1;
    }
    HANDLE hMem = CreateFileMappingA(
        hFile,
        NULL,
        PAGE_READONLY,
        0,
        0,
        NULL
    );
    if (hMem == NULL) {
        retval = -1;
        goto close_file;
    }
    mapped_region->addr = MapViewOfFile(
        hMem,
        FILE_MAP_READ,
        0,
        0,
        0
    );
    if (mapped_region->addr == NULL) {
        retval = -1;
    }
#else
    int fd = open(
        pathname,
        O_RDONLY
    );
    if (fd == -1) {
        return -1;
    }
    mapped_region->addr = mmap(
        NULL,
        sb.st_size,
        PROT_READ,
        MAP_SHARED,
        fd,
        0
    );
    if (mapped_region->addr == MAP_FAILED) {
        retval = -1;
    }
#endif
#ifdef _WIN32
    if (CloseHandle(hMem) == 0) {
        retval = -1;
    }
close_file:
    if (CloseHandle(hFile) == 0) {
        retval = -1;
    };
#else
    if (close(fd) == -1) {
        retval = -1;
    }
#endif
    return retval;
}

int unmap_file(struct mapped_region *mapped_region) {
    int retval = 0;
 // Check if region is already unmapped
    if (!mapped_region->addr) {
        return -1;
    }
#ifdef _WIN32
    if (UnmapViewOfFile(mapped_region->addr) == 0) {
        return -1;
    }
#else
    if (munmap(mapped_region->addr, mapped_region->length) == -1) {
        return -1;
    }
#endif
 // Reset address and length
    mapped_region->addr = NULL;
    mapped_region->length = 0;
    return retval;
}
