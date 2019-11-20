#include "mapping.h"

#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <sys/mman.h>
#endif

int map_file(const char *pathname, struct mapped_region *mapped_region, char mode) {
    int retval = 0;
 // Check mode
    if (mode != 'r' && mode != 'w') {
        errno = 8;
        return -1;
    }
 // Get file status
    struct stat sb;
    if (stat(pathname, &sb) == -1) {
        errno = 1;
        return -1;
    }
 // Test if file is a regular file
    if (!S_ISREG(sb.st_mode)) {
        errno = 7;
        return -1;
    }
 // Get file length
    mapped_region->length = sb.st_size;
 // Map file
#ifdef _WIN32
    HANDLE hFile = CreateFileA(
        pathname,
        mode == 'w' ? GENERIC_READ | GENERIC_WRITE : GENERIC_READ,
        mode == 'w' ? 0 : FILE_SHARE_READ,
        NULL,
        OPEN_EXISTING,
        FILE_ATTRIBUTE_NORMAL,
        NULL
    );
    if (hFile == INVALID_HANDLE_VALUE) {
        errno = 2;
        return -1;
    }
    HANDLE hMem = CreateFileMappingA(
        hFile,
        NULL,
        mode == 'w' ? PAGE_READWRITE : PAGE_READONLY,
        0,
        0,
        NULL
    );
    if (hMem == NULL) {
        errno = 3;
        retval = -1;
        goto close_file;
    }
    mapped_region->addr = MapViewOfFile(
        hMem,
        mode == 'w' ? FILE_MAP_ALL_ACCESS : FILE_MAP_READ,
        0,
        0,
        0
    );
    if (mapped_region->addr == NULL) {
        errno = 4;
        retval = -1;
    }
#else
    int fd = open(
        pathname,
        mode == 'w' ? O_RDWR : O_RDONLY
    );
    if (fd == -1) {
        errno = 2;
        return -1;
    }
    mapped_region->addr = mmap(
        NULL,
        sb.st_size,
        mode == 'w' ? PROT_READ | PROT_WRITE : PROT_READ,
        MAP_SHARED,
        fd,
        0
    );
    if (mapped_region->addr == MAP_FAILED) {
        errno = 4;
        retval = -1;
    }
#endif
#ifdef _WIN32
    if (CloseHandle(hMem) == 0) {
        errno = 5;
        retval = -1;
    }
close_file:
    if (CloseHandle(hFile) == 0) {
        errno = 6;
        retval = -1;
    };
#else
    if (close(fd) == -1) {
        errno = 6;
        retval = -1;
    }
#endif
    return retval;
}

int unmap_file(struct mapped_region *mapped_region) {
    int retval = 0;
 // Check if region is already unmapped
    if (!mapped_region->addr) {
        errno = 1;
        return -1;
    }
#ifdef _WIN32
    if (UnmapViewOfFile(mapped_region->addr) == 0) {
        errno = 2;
        return -1;
    }
#else
    if (munmap(mapped_region->addr, mapped_region->length) == -1) {
        errno = 2;
        return -1;
    }
#endif
 // Reset address and length
    mapped_region->addr = NULL;
    mapped_region->length = 0;
    return retval;
}
