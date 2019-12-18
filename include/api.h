/*
 * Copyright 2016-2019 CNRS-UM LIRMM, CNRS-AIST JRL
 */

#pragma once

// Handle portable symbol export.
// Defining manually which symbol should be exported is required
// under Windows whether MinGW or MSVC is used.
//
// The headers then have to be able to work in two different modes:
// - dllexport when one is building the library,
// - dllimport for clients using the library.
//
// On Linux, set the visibility accordingly. If C++ symbol visibility
// is handled by the compiler, see: http://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
// On Microsoft Windows, use dllimport and dllexport to tag symbols.
#define COPRA_DLLIMPORT __declspec(dllimport)
#define COPRA_DLLEXPORT __declspec(dllexport)
#define COPRA_DLLLOCAL
#else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#if __GNUC__ >= 4
#define COPRA_DLLIMPORT __attribute__((visibility("default")))
#define COPRA_DLLEXPORT __attribute__((visibility("default")))
#define COPRA_DLLLOCAL __attribute__((visibility("hidden")))
#else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#define COPRA_DLLIMPORT
#define COPRA_DLLEXPORT
#define COPRA_DLLLOCAL
#endif // __GNUC__ >= 4
#endif // defined _WIN32 || defined __CYGWIN__

#ifdef COPRA_STATIC
// If one is using the library statically, get rid of
// extra information.
#define COPRA_DLLAPI
#define COPRA_LOCAL
#else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#ifdef copra_EXPORTS
#define COPRA_DLLAPI COPRA_DLLEXPORT
#else
#define COPRA_DLLAPI COPRA_DLLIMPORT
#endif // COPRA_EXPORTS
#define COPRA_LOCAL COPRA_DLLLOCAL
#endif // COPRA_STATIC
