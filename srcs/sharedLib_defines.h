#ifndef SHAREDLIB_DEFINES_H
#define SHAREDLIB_DEFINES_H

#if defined(_WIN32)
	#define PFEM_EXPORT __declspec(dllexport)
	#define PFEM_IMPORT __declspec(dllimport)
	#define PFEM_NOEXPORT
	#if defined(_MSC_VER)
		#define PFEM_DEPRECATED __declspec(deprecated)
	#else
		#define PFEM_DEPRECATED __attribute__ ((__deprecated__))
	#endif
	#define PFEM_DEPRECATED_EXPORT PFEM_EXPORT PFEM_DEPRECATED
	#define PFEM_DEPRECATED_NOEXPORT PFEM_NOEXPORT PFEM_DEPRECATED
#elif defined(__linux__) or defined(__unix__)
	#define PFEM_EXPORT __attribute__((visibility("default")))
	#define PFEM_IMPORT __attribute__((visibility("default")))
	#define PFEM_NOEXPORT __attribute__((visibility("hidden")))
	#define PFEM_DEPRECATED __attribute__ ((__deprecated__))
	#define PFEM_DEPRECATED_EXPORT PFEM_EXPORT PFEM_DEPRECATED
	#define PFEM_DEPRECATED_NOEXPORT PFEM_NOEXPORT PFEM_DEPRECATED
#elif defined(__APPLE__) && defined(__MACH__)
	#define PFEM_EXPORT __attribute__((visibility("default")))
	#define PFEM_IMPORT __attribute__((visibility("default")))
	#define PFEM_NOEXPORT __attribute__((visibility("hidden")))
	#define PFEM_DEPRECATED __attribute__ ((__deprecated__))
	#define PFEM_DEPRECATED_EXPORT PFEM_EXPORT PFEM_DEPRECATED
	#define PFEM_DEPRECATED_NOEXPORT PFEM_NOEXPORT PFEM_DEPRECATED
#else
	#error("Unsupported compiler!")
#endif

#endif /* SHAREDLIB_DEFINES_H */
